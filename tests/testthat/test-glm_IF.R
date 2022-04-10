test_that("glm_IF works", {
  x1 <- rnorm(1e3)
  x2 <- rnorm(1e3)
  y <- rbinom(n = 1e3, size = 1, plogis(-1 + x1 + 2 * x2))
  glm_obj <- glm(y ~ x1 + x2, family = binomial(link = "logit"))
  expect_equal(crossprod(glm_IF(glm_obj)) / 1e6,
               sandwich::vcovHC(glm_obj, type = "HC"))
})
