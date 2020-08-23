test_that("sample_y returns a double", {
  expect_type(sample_y(0.1, 0.1, 0.1), "double")
})
