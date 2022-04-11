#' Empirical influence function of glm object
#'
#' @param glm_obj The object returned by glm().
#'
#' @return matrix A `n` by `p` matrix.
#' @export
#'
#' @examples
#' x1 <- rnorm(1e3)
#' x2 <- rnorm(1e3)
#' y <- rbinom(n = 1e3, size = 1, plogis(-1 + x1 + 2 * x2))
#' glm_obj <- glm(y ~ x1 + x2, family = binomial(link = "logit"))
#' all.equal(crossprod(glm_IF(glm_obj)) / 1e6, sandwich::vcovHC(glm_obj, type = "HC"))
glm_IF <- function(glm_obj) {
  # https://github.com/ohines/plmed/issues/1
  # Thanks to Oliver Hines.
  # This is slightly slower but won't affect performance much
  # since it's still O(n).
  sandwich::estfun(glm_obj)%*%sandwich::bread(glm_obj)
}
