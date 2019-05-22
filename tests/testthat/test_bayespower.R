context("Errors work")
library(BayesianPower)

test_that("No invalid entries are allowed", {
  n <- 12.3
  h1 <- matrix(c(1,2,3,3,2,1), nrow = 2, byrow = TRUE)
  h2 <- "u"
  m1 <- c(1.3, 23, 1)
  m2 <- c(1, 0, 21)
  expect_error(bayes_power(n, h1, h2, m1, m2))
})
