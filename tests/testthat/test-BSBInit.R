test_that("initial setup only accepts one of df or t1 and t2", {
  df <- data.frame(a = 1:3)
  t1 <- survival::Surv(c(1, 2, 3))
  t2 <- survival::Surv(c(5, 2, 1))
  
  expect_error(
    BSBInit(df, t1, t2),
    "Only one of df or t1 and t2 must be supplied"
  )
  expect_error(
    BSBInit(df, t1), "Only one of df or t1 and t2 must be supplied"
  )
  expect_error(
    BSBInit(df, t2 = t2),
    "Only one of df or t1 and t2 must be supplied"
  )
})

test_that("initial setup throws an error if t1 or t2 are not of class 'Surv'", {
  expect_error(BSBInit(t1 = c(1, 2, 3), t2 = c(1, 2, 3)))
})

test_that("initial setup throws an error if t1, t2 dont have the same length", {
  t1 <- survival::Surv(c(1, 1, 1))
  t2 <- survival::Surv(c(1, 1))
  expect_error(
    BSBInit(t1 = t1, t2 = t2), "t1 and t2 must have the same length"
  )
})

test_that("initial setup fails if df is not a data frame", {
  expect_error(
    BSBInit(df = c(1, 2, 3)), "df must be of class 'data.frame'"
  )
})

test_that("t1 and t2 must be right-censored", {
  t1 <- survival::Surv(c(1, 1, 1), c(1, 1, 1), type = "left")
  t2 <- survival::Surv(c(1, 1, 1))
  expect_error(
    BSBInit(t1 = t1, t2 = t2), "t1 has censoring of type: left"
  )
  expect_error(
    BSBInit(t1 = t2, t2 = t1), "t2 has censoring of type: left"
  )
})

test_that("initial setup returns an object of class 'BSBinit'", {
  t1 <- survival::Surv(c(1, 2, 3))
  t2 <- survival::Surv(c(1, 2, 3))
  df1 <- data.frame(t1 = c(1, 2, 3), t2 = c(1, 2, 2))
  df2 <- data.frame(t1 = c(1, 2), t2 = c(1, 2),
                    delta1 = c(1, 1), delta2 = c(0, 0))
  
  expect_s3_class(BSBInit(t1 = t1, t2 = t2), class = "BSBinit")
  expect_s3_class(BSBInit(df1), class = "BSBinit")
  expect_s3_class(BSBInit(df2), class = "BSBinit")
})

