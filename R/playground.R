#
# df <- dgp_obs(
#   n = 400,
#   p_c = 0.7,
#   alpha = 1,
#   beta_1 = 0.5,
#   beta_2 = -0.7,
#   beta_3 = 0,
#   error_dist = "gumbel"
# )
#
# ps_model_correct <- glm(z ~ x1_t + x2_t, family = binomial(link = "logit"), data = df, x = T)
#
# iv_correct_fit <- ivpim_obs(
#   y = df$y,
#   z = df$z,
#   a = df$a,
#   X = model.matrix(~ x1 + x2 - 1, df),
#   ps_model = ps_model_correct,
#   link = "logit"
# )
#
# pim_obj <- upim::pim_fit(y = df$y, X = df[, c("a", "x1", "x2")], link = "logit",
#               w = w,
#               tol = 1e-8, keep.data = TRUE)
#
# pim_obj$vcov
# iv_correct_fit$vcov_naive
#
