test_that("abadie_k works", {
  expect_equal(abadie_k(1, 1, 0.5), 1)
  expect_equal(abadie_k(1, 0, 0.5), -1)
  expect_equal(abadie_k(0, 1, 0.5), -1)
  expect_equal(abadie_k(0, 0, 0.5), 1)
  expect_equal(abadie_k(c(1, 1, 0, 0), c(1, 0, 1, 0), 0.5), c(1, -1, -1, 1))
  expect_equal(abadie_k(c(1, 1), c(1, 1), c(0.4, 0.6)), c(1, 1))
  expect_equal(abadie_k(c(1, 1), c(0, 0), c(0.4, 0.6)), c(1 - 1 / 0.4, 1 - 1 / 0.6))
})
