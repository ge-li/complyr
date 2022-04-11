#' Abadie's nu
#'
#' This function computes the derivative of Abadie's kappa,
#' \deqn{\nu(\pi)=\dot{\kappa}(\pi)=\frac{Z(1-A)}{\pi^{2}}-\frac{A(1-Z)}{(1-\pi)^{2}}}.
#'
#' @param z numeric A binary instrument, e.g., randomized treatment assignment.
#' @param a numeric The actual received treatment, must be the same length as `z`.
#' @param ps numeric The propensity scores \eqn{\Pr(Z = 1 \mid X)}.
#'
#' @return numeric The formula defined above.
#' @export
#'
#' @examples
#' abadie_v(1, 1, 0.5)
#' abadie_v(1, 0, 0.5)
#' abadie_v(0, 1, 0.5)
#' abadie_v(0, 0, 0.5)
#' abadie_v(c(1, 1, 0, 0), c(1, 0, 1, 0), 0.5)
#' abadie_v(c(1, 1), c(1, 1), c(0.4, 0.6))
#' abadie_v(c(1, 1), c(0, 0), c(0.4, 0.6))
abadie_v <- function(z, a, ps) {
  (1 - a) * z / (ps^2) - a * (1 - z) / ((1 - ps)^2)
}
