#' Generate Data for Simulation Studies
#'
#' @param n_obs numeric Number of observations.
#' @param ps_type character Type of IV propensity score DGP.
#' @param p_c numeric Compliance rate.
#' @param alpha numeric Coefficient associated with the treatment difference.
#' @param beta_1 numeric Coefficients associated with `x1`.
#' @param beta_2 numeric Coefficients associated with `x2`.
#' @param beta_3 numeric Coefficients associated with `x1 * x2`.
#' @param error_dist "gumbel" or "normal".
#'
#' @return A data frame
#' @export
#'
#' @examples
#' set.seed(42)
#' generate_data(n_obs = 6, ps_type = "la")
generate_data <- function(n_obs = 200,
                          ps_type = c("la", "nla", "nlna", "omit", "complex", "rct"),
                          p_c = 0.7,
                          alpha = 1,
                          beta_1 = 0.5,
                          beta_2 = -0.7,
                          beta_3 = 0,
                          error_dist = c("gumbel", "normal"))
{
  ps_type <- match.arg(ps_type)
  error_dist <- match.arg(error_dist)
  # Check if p_c lies in the range [0,1]
  if (p_c < 0 | p_c > 1) {
    stop("Compliance rate (p_c) must be between 0 and 1.")
  }

  # Generate covariates and propensity scores
  x1 <- stats::runif(n_obs)
  x2 <- stats::runif(n_obs)
  ps_type = match.arg(ps_type)
  if (ps_type == "la") {
    eta <- 1 + 2 * x1 - 3 * x2
  } else if (ps_type == "nla") {
    eta <- -2 + sin(2 * pi * x1) + exp(x2)
  } else if (ps_type == "nlna") {
    eta <- sin(2 * pi * x1) * exp(x2)
  } else if (ps_type == "omit") {
    eta <- 1 + 2 * x1 - 3 * x2 + 0.3 * stats::rnorm(n_obs)
  } else if (ps_type == "complex") {
    eta <- sin(1 / x1) * cos(1 / x2)
  } else if (ps_type == "rct") {
    eta <- rep(0, n_obs)
  }
  p <- stats::plogis(eta)
  z <- stats::rbinom(n_obs, 1, p)
  df <- data.frame(x1, x2, eta, p, z)

  # Assign each subject to a group: always-taker, complier, never-taker
  principle_strata_prob <- c((1 - p_c) / 2, p_c, (1 - p_c) / 2)
  c <- sample(c(0, 1, 2), n_obs, replace = T, prob = principle_strata_prob)

  a <- (c == 0) + df$z * (c == 1) # actual treatment

  # Calculate expected outcome
  ey <- alpha * a + beta_1 * df$x1 + beta_2 * df$x2 + beta_3 * df$x1 * df$x2
  ey[c == 0] <- (ey[c == 0] - alpha) * 0.7 + alpha + 1
  ey[c == 2] <- (ey[c == 2] ) * (- 0.7) - 1

  # Add error term
  if (error_dist == "normal") {
    err <- stats::rnorm(n_obs, sd = 1 / sqrt(2))
  } else if (error_dist == "gumbel") {
    err <- evd::rgumbel(n_obs, scale = 1)
  } else {
    stop("Error term must be either normal or gumbel!")
  }

  y <- ey + err

  # add new columns to the data frame
  df$c <- c
  df$a <- a
  df$ey <- ey
  df$y <- y
  # return the data frame
  df
}

