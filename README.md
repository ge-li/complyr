
<!-- README.md is generated from README.Rmd. Please edit that file -->

# complyr

<!-- badges: start -->
<!-- badges: end -->

The goal of complyr is to create a one-stop package for analyzing
compiler causal effects.

Unignorable confounding is no stranger even for randomized controlled
trials (RCTs) in the presence of treatment non-compliance. One fallback
is the intention-to-treat (ITT) analysis, which unfortunately only
reflects the assignment-induced causal effect. The celebrated (Angrist,
Imbens & Rubin 1996) framework provides an excellent platform to tackle
unmeasured confounding for estimating the complier (local) average
treatment effect (ATE). This package aims to provide a suite of tools
for analyzing different complier causal effects including ATE, quantile
treatment effect (QTE), and probabilistic index.

## Installation

You can install the development version of complyr from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ge-li/complyr")
```

## Example

This is a basic example which shows you how to estimate a complier
probabilistic index model, with comparison to intention-to-treat and
per-protocol methods.

``` r
library(complyr)
# Simulate some RCT data with non-compliance, see details in function docs. 
set.seed(123)
df <- dgp_rct(
  n = 1000,
  p_c = 0.8,
  alpha = 1,
  beta_1 = 0.5,
  beta_2 = -0.7,
  error_dist = "gumbel"
)
# Estimate the complier probabilistic index models
# We use three methods to analyze this data set: ITT, PP, IV
# ITT: intention-to-treat
itt_fit <- upim::pim_fit(y = df$y, X = df[, c("z", "x1", "x2")], link = "logit")
# PP: per-protocol
df_pp <- df[df$z == df$a, ]
pp_fit <- upim::pim_fit(y = df_pp$y, X = df_pp[, c("a", "x1", "x2")], link = "logit")
# IV: instrumental variable
ps_model <- glm(df$z ~ 1, family = binomial(link = "logit"), x = TRUE)
iv_fit <- complyr::ivpim(y = df$y, z = df$z, a = df$a, X = df[, c("x1", "x2")], 
                         ps_model =  ps_model, link = "logit")
sum_stat <- function(model_fit) {
  # get summary stats for downstream analysis
  ss <- c(model_fit$coef, sqrt(diag(model_fit$vcov)))
  names(ss) <- c("alpha_hat", "beta_1_hat", "beta_2_hat",
                 "alpha_se", "beta_1_se", "beta_2_se")
  round(ss, 3)
}
results <- as.data.frame(rbind(sum_stat(itt_fit),
                               sum_stat(pp_fit),
                               sum_stat(iv_fit)))
results$methods <- c("itt", "pp", "iv")
results
#>   alpha_hat beta_1_hat beta_2_hat alpha_se beta_1_se beta_2_se methods
#> 1     0.807      0.375     -0.469    0.082     0.042     0.046     itt
#> 2     1.229      0.464     -0.584    0.089     0.043     0.046      pp
#> 3     1.186      0.545     -0.674    0.111     0.065     0.072      iv
```
