omega <- 1.4
y <- 2
cum.h <- 0.8
x <- c(1, 2, 1)
theta <- c(0.7, 0.8, 0.1)
gamma <- 1.2

test_that("sample_omega returns a double", {
  delta.omega <- 0.8
  s1 <- sample_omega(omega, y, cum.h, x, theta, gamma, delta.omega)
  s2 <- sample_omega(omega, y, cum.h, x, theta, gamma)

  expect_type(s1, "double")
  expect_type(s2, "double")
})

test_that("sample_omega correctly assigns delta", {
  delta.omega <- y + 1
  set.seed(42)
  s1 <- sample_omega(omega, y, cum.h, x, theta, gamma, delta.omega)
  set.seed(42)
  s2 <- sample_omega(omega, y, cum.h, x, theta, gamma)

  expect_equal(s1, s2)
})
