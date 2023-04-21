#' Instrumental variable probabilistic index model
#'
#' This model targets the compatible complier probabilistic index model.
#'
#' @param y numeric Outcome vector.
#' @param z numeric Binary instrument, e.g., randomized treatment assignment.
#' @param a numeric Actual received treatment, must be the same length as `z`.
#' @param X numeric Covariates.
#' @param ps_model glm First step propensity score `glm()` object.
#' @param link character Link function: "logit", or "probit".
#' @param init numeric Initial guess of Newton's method.
#' @param tol numeric Numeric tolerance of `nleqslv()`.
#' @param max.iter numeric Maximum iteration of Newton's method.
#' @param nleqslv.global character The global strategy for Newton's method. See ?nleqslv::nleqslv.
#' @param trace logical Show Newton's method iteration report if TRUE.
#' @param test.nleqslv logical Test different global strategies for Newton's method if TRUE. See ?nleqslv::testnslv.
#' @param keep.data logical Should the returned object keep original data?
#'
#' @return A list containing the estimated coefficients and their covaraince
#'           matrix. It also contains the diagnostics of `nleqlsv()` procedure.
#'           If `keep.data` is `TRUE`, then the inputs `y`, `z`, `a`, `X`, `link`, `w`
#'           will also be returned.
#' @export
ivpim <- function(y,
                  z,
                  a,
                  X = NULL,
                  ps_model,
                  link = "logit",
                  init = NULL,
                  tol = sqrt(.Machine$double.eps),
                  max.iter = 100,
                  nleqslv.global = "none",
                  trace = FALSE,
                  test.nleqslv = FALSE,
                  keep.data = FALSE
) {
  n <- NROW(y)
  if (is.null(X)) {
    p <- 1
  } else {
    p <- NCOL(X) + 1
  }
  k <- abadie_k(z, a, ps_model$fitted.values)
  v <- abadie_v(z, a, ps_model$fitted.values)
  w <- as.vector(t(outer(k, k))) # iv weight
  W <- as.matrix(cbind(a, X)) # covariates
  pim_obj <- upim::pim_fit(
    y = y,
    X = W,
    link = link,
    w = w,
    init = init,
    tol = tol,
    max.iter = max.iter,
    nleqslv.global = nleqslv.global,
    trace = trace,
    test.nleqslv = test.nleqslv,
    keep.data = keep.data
  )

  # calculate drift matrix
  sigma_dot <- ps_model$family$mu.eta(ps_model$linear.predictors)
  kvs <- as.vector(t(outer(k, v * sigma_dot))) # kappa_i * v_j * sigma_dot_j
  vsk <- as.vector(t(outer(v * sigma_dot, k))) # v_i * sigma_dot_i * kappa_j
  ind <- seq_len(n) # c(1, 2, ..., n)
  drift_right <- kvs * model.matrix(ps_model)[rep(ind, n), , drop = F] + # (j: 12...n, 12...n, ...)
                 vsk * model.matrix(ps_model)[rep(ind, rep(n, n)), , drop = F] # (i: 1...1, 2...2, ...)
  U_raw <- upim::est_fun(upim::create_pseudo_obs(y),
                         upim::new_design_matrix(W),
                         pim_obj$coef,
                         link = link) # n^2 x p
  drift_matrix <- crossprod(U_raw, drift_right) / (n^2) # p x p
  drifts <- tcrossprod(glm_IF(ps_model), drift_matrix)

  # adjusted ivpim sandwich meat
  if (p == 1) {
    U_i <- 2 * sapply(split.data.frame(U_raw * w, rep(seq(n), each = n)), colSums) / (n - 1) # n x p
  } else {
    U_i <- 2 * t(sapply(split.data.frame(U_raw * w, rep(seq(n), each = n)), colSums) / (n - 1)) # n x p
  }

  modified_U_i <- U_i + drifts
  ivpim_sandwich_meat <- crossprod(modified_U_i) / n

  # correct asymptotic variance
  ivpim_vcov <- solve(pim_obj$jac) %*% ivpim_sandwich_meat %*% solve(pim_obj$jac) / n
  colnames(ivpim_vcov) <- colnames(pim_obj$vcov)
  rownames(ivpim_vcov) <- rownames(pim_obj$vcov)
  pim_obj$vcov_naive <- pim_obj$vcov
  pim_obj$vcov <- ivpim_vcov

  # return ivpim object
  if (keep.data) {
    pim_obj$z <- z
    pim_obj$a <- a
  }
  pim_obj
}
