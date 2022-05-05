#' Data generating process: randomized controlled trials
#'
#' This function inherit dgp_internal(). TODO: ADD MORE DESCRIPTION.
#'
#' @param n numeric Number of samples (observations).
#' @param p_c numeric Compliance rate.
#' @param alpha numeric Coefficient associated with the treatment difference
#' @param beta_1 numeric Coefficients associated with `x1`.
#' @param beta_2 numeric Coefficients associated with `x2`.
#' @param beta_3 numeric Coefficients associated with `x1 * x2`.
#' @param error_dist "gumbel" or "normal".
#'
#' @return TODO
#' @export
#'
#' @examples
#' set.seed(42)
#' dgp_rct(n = 6)
dgp_rct <- function(n = 200,
                    p_c = 0.7,
                    alpha = 1,
                    beta_1 = 0.5,
                    beta_2 = -0.7,
                    beta_3 = 0,
                    error_dist = "gumbel") {
  dgp_internal(
    n = n,
    p_c = p_c,
    alpha = alpha,
    beta_1 = beta_1,
    beta_2 = beta_2,
    beta_3 = beta_3,
    error_dist = error_dist
  )
}
