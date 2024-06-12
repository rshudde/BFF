library(testthat)
library(ggplot2)

###Testing plot.BFF function

example_BFF_generic <- list(
  BFF = list(tau2 = seq(0, 1, length.out = 100), log_bf = rnorm(100)),
  generic_test = TRUE,
  test_type = "t_test"
)

example_BFF_specific <- list(
  BFF = list(omega = seq(0, 1, length.out = 100), log_bf = rnorm(100)),
  generic_test = FALSE,
  test_type = "t_test"
)

test_that("plot.BFF stops if BFF is NULL", {
  expect_error(plot.BFF(list(BFF = NULL)),
               "Bayes factor function can be plotted only if a specific omega/tau2 is not user set")
})

test_that("plot.BFF returns a data.frame when plot = FALSE", {
  df <- plot.BFF(example_BFF_generic, plot = FALSE)
  expect_s3_class(df, "data.frame")
  expect_named(df, c("x", "log_BF"))
})

test_that("plot.BFF returns a ggplot object when plot = TRUE", {
  p <- plot.BFF(example_BFF_generic, plot = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("plot.BFF handles additional arguments correctly", {
  p <- plot.BFF(example_BFF_generic, plot = TRUE, title = "Title", ylab = "Ylab")
  expect_equal(p$labels$title, "Title")
  expect_equal(p$labels$y, "Ylab")
})

test_that("plot.BFF sets effect size cutpoints and colors correctly", {
  cutpoints <- .get_effect_size_cutpoints("t_test")
  expect_equal(cutpoints, c(0.2, 0.5, 0.8))

  range <- .get_effect_size_range("t_test")
  expect_equal(range, c(0.0, 1.0))
})
