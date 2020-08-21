test_that("partitions are created correctly", {
  t1 <- c(0.1, 0.2, 1, 3, 5)
  p1 <- c(0, 2, 4, 6)
  p2 <- c(0, 5)
  p3 <- c(0, 1, 2, 3, 4, 5)
  p4 <- c(0, 10)
  
  expect_equal(partition(t1, 2), p1)
  expect_equal(partition(t1, 5), p2)
  expect_equal(partition(t1, 1), p3)
  expect_equal(partition(t1, 10), p4)
})

test_that("times are located correctly", {
  t <- c(0.1, 1.1, 2.1, 3, 4, 5)
  p1 <- c(0, 1, 2, 3, 4, 5)
  p2 <- c(0, 1, 2, 3, 4, 5, 6)
  p3 <- c(0, 0.5, 1, 1.5, 2)
  
  expect_equal(partition_location(t, p1), c(1, 2, 3, 3, 4, 5))
  expect_equal(partition_location(t, p2), c(1, 2, 3, 3, 4, 5))
  expect_equal(partition_location(t, p3), c(1, 3, 5, 5, 5, 5))
})

test_that("times by interval are counted correctly", {
  t <- c(0.1, 1.1, 1.2, 2)
  p1 <- c(0, 1, 2)
  p2 <- c(0, 0.5, 1, 1.5, 2)
  
  expect_equal(partition_count(t, p1), c(1, 3))
  expect_equal(partition_count(t, p2), c(1, 0, 2, 1))
})

test_that("cumulative hazards are computed correctly", {
  t1 <- 3
  p1 <- c(0, 1, 2, 3)
  lambda <- c(0.1, 0.2, 0.3)
  t2 <- c(0.2, 0.7)
  p2 <- c(0, 0.5, 1, 1.5)
  out2 <- c(0.2 * 0.1, 0.5 * 0.1 + 0.2 * 0.2)
  
  expect_equal(cum_h(t1, p1, lambda), 0.1 + 0.2 + 0.3)
  expect_equal(cum_h(t2, p2, lambda), out2)
})
