#' Abadie's kappa
#'
#' This function computes the \eqn{\kappa = \kappa (Z, A, X) = 1 -
#'   \frac{A(1-Z)}{1 - \Pr(Z = 1 \mid X)} - \frac{(1-A)Z}{\Pr(Z = 1 \mid X)}}
#'   defined in Abadie, Angrist & Imbens (2002), and Abadie (2003). The expected
#'   value of \eqn{\kappa} is the probability of complier.
#'
#' @param z numeric A binary instrument, e.g., randomized treatment assignment.
#' @param a numeric The actual received treatment, must be the same length as `z`.
#' @param ps numeric The propensity scores \eqn{\Pr(Z = 1 \mid X)}.
#'
#' @return numeric The formula defined above.
#' @export
#'
#' @examples
#' abadie_k(1, 1, 0.5)
#' abadie_k(1, 0, 0.5)
#' abadie_k(0, 1, 0.5)
#' abadie_k(0, 0, 0.5)
#' abadie_k(c(1, 1, 0, 0), c(1, 0, 1, 0), 0.5)
#' abadie_k(c(1, 1), c(1, 1), c(0.4, 0.6))
#' abadie_k(c(1, 1), c(0, 0), c(0.4, 0.6))
abadie_k <- function(z, a, ps) {
  1 - a * (1 - z) / (1 - ps) - (1 - a) * z / ps
}
