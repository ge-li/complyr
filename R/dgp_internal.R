#' Data generating process: linear transformation model with non-compliance
#'
#' TODO: ADD MORE DESCRIPTION HERE.
#'
#' @param n numeric Number of samples (observations).
#' @param x1_dist character {"binary", "uniform", "normal"}.
#' @param x2_dist character {"binary", "uniform", "normal"}.
#' @param sigma_x numeric Standard deviation of x (if not binary).
#' @param nonlinearity character Apply what nonlinear transformation to `x`?
#' @param ps_link character link function for propensity score DGP, see ?make.link.
#' @param ps_trim numeric This enforces positivity assumption.
#'                  (ps_trim, 1 - ps_trim) will be the bound of generated
#'                  propensity score, implemented by rejection sampling.
#' @param eta_0 numeric Coefficients for the propensity score DGP.
#' @param eta_1 numeric Coefficients for the propensity score DGP.
#' @param eta_2 numeric Coefficients for the propensity score DGP.
#' @param eta_3 numeric Coefficients for the propensity score DGP.
#' @param eta_4 numeric Coefficients for the propensity score DGP.
#' @param p_c numeric Compliance rate.
#' @param rho_a numeric Proportion of always-takers among non-compliers.
#' @param nc_scale numeric Scale non-complier coefficients.
#' @param nc_shift numeric Shift non-complier outcomes.
#' @param alpha numeric Coefficient associated with the treatment difference
#' @param beta_1 numeric Coefficients associated with `x1`.
#' @param beta_2 numeric Coefficients associated with `x2`.
#' @param beta_3 numeric Coefficients associated with `x1 * x2`.
#' @param sigma_e standard deviation of the error term, default is 1.
#' @param error_dist "gumbel" or "normal".
#'
#' @return
#' df - generated data frame with columns {x1, x2, x2_t, lp, ps, z, c, a, ey, y} \cr
#' df$... - see "dgp_prop_score()" \cr
#' df$z - randomized treatment assignment based on the propensity scores \cr
#' df$c - complier status {0, 1, 2} = {always-taker, complier, never-taker} \cr
#' df$a - actual treatment received \cr
#' df$ey - expected value of y \cr
#' df$y - outcome \cr
#' @export
dgp_internal <- function(n = 200, x1_dist = "normal", x2_dist = "normal", sigma_x = 1,
                         nonlinearity = "sin", ps_link = "logit", ps_trim = 0.05,
                         eta_0 = 0, eta_1 = 0, eta_2 = 0, eta_3 = 0, eta_4 = 0,
                         p_c = 0.7, rho_a = 0.5, nc_scale = 0.5, nc_shift = 1,
                         alpha = 1, beta_1 = 0.5, beta_2 = -0.7, beta_3 = 0,
                         sigma_e = 1, error_dist = "gumbel") {

  # generate covariates and propensity scores
  df <- dgp_prop_score(n, x1_dist, x2_dist, sigma_x,
                       nonlinearity, ps_link, ps_trim,
                       eta_0, eta_1, eta_2, eta_3, eta_4)

  # c = 0, 1, 2: always-taker, complier, never-taker
  principle_strata_prob <- c((1 - p_c) * rho_a, p_c, (1 - p_c) * (1 - rho_a))
  c <- sample(c(0, 1, 2), n, replace = T, prob = principle_strata_prob)

  a <- (c == 0) + df$z * (c == 1) # actual treatment

  # compliers: E[Y|A,X]
  ey <- alpha * a + beta_1 * df$x1 + beta_2 * df$x2 + beta_3 * df$x1 * df$x2
  # always-takers: scale and shift randomly
  ey[c == 0] <- (ey[c == 0] - alpha) * nc_scale + nc_shift
  # never takers: scale and shift randomly
  ey[c == 2] <- (ey[c == 2] ) / nc_scale - nc_shift

  # add error term
  if (error_dist == "normal") {
    # for normal error, divide sigma_e by sqrt(2) so that
    err <- stats::rnorm(n, sd = sigma_e / sqrt(2))
  } else if (error_dist == "gumbel") {
    err <- evd::rgumbel(n, scale = sigma_e)
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
