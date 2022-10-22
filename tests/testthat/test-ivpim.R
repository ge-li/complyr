test_that("marginal ivpim works", {
  set.seed(42)
  n <- 1e2
  p_c <- 0.8

  c <- sample(c(0, 1, 2), n, replace = T, prob = c(0.1, 0.8, 0.1))
  z <- sample(c(0, 1), n, replace = T)
  a <- (c == 0) + z * (c == 1) # actual treatment

  x <- rnorm(n)
  y <- rnorm(n)

  ps_model <- glm(z ~ 1, family = binomial(link = "logit"), x = TRUE)

  expect_error(
    iv_pim_marginal <- ivpim(
      y = y,
      z = z,
      a = a,
      X = NULL,
      link = "logit",
      ps_model = ps_model
    ), NA # means no error is expected!!
  )
})

test_that("ivpim works", {
  set.seed(42)
  n <- 1e2
  p_c <- 0.8

  c <- sample(c(0, 1, 2), n, replace = T, prob = c(0.1, 0.8, 0.1))
  z <- sample(c(0, 1), n, replace = T)
  a <- (c == 0) + z * (c == 1) # actual treatment

  x <- rnorm(n)
  y <- rnorm(n)

  ps_model <- glm(z ~ 1, family = binomial(link = "logit"), x = TRUE)

  expect_error(
    iv_pim_marginal <- ivpim(
      y = y,
      z = z,
      a = a,
      X = x,
      link = "logit",
      ps_model = ps_model
    ), NA # means no error is expected!!
  )
})
