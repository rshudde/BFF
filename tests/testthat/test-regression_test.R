test_that("two-sample: basic functionality", {

  ### test fixed omega -- non-local prior regression (two-sided)
  fit <- regression_test_BFF(
    t_stat = 1.5,
    alternative = "two.sided",
    n = 50,
    k = 3,
    omega = 0.5)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf_h1, -0.8319082, tolerance = 1e-5)
  testthat::expect_equal(fit$omega_h1,  0.5)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local regression test"  ,
      ""                                        ,
    "log Bayes factor = -0.83"                 ,
      "omega = 0.50 (signed Cohen's f)"                ,
      "alternative = two.sided"
    )
  )
  testthat::expect_error(plot(fit), "Bayes factor function can be plotted only if a specific omega/tau2 is not user set")

  # # vdiffr::expect_doppelganger("t_test-two_sample-two_sided-posterior",           posterior_plot(fit))
  # # vdiffr::expect_doppelganger("t_test-two_sample-two_sided-posterior_and_prior", posterior_plot(fit, prior = TRUE))


  ### test fixed omega -- non-local prior regression test (one-sided, also set r)
  fit <- regression_test_BFF(
    t_stat = 1.5,
    alternative = "greater",
    n = 50,
    k = 2,
    r = 3,
    omega = 0.5)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf_h1, -0.27897, tolerance = 1e-4)
  testthat::expect_equal(fit$omega_h1,  0.50)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local regression test"  ,
      ""                                        ,
      "log Bayes factor = -0.28"                 ,
      "omega = 0.50 (signed Cohen's f)"                ,
      "alternative = greater"
    )
  )
  # vdiffr::expect_doppelganger("t_test-two_sample-one_sided-posterior",           posterior_plot(fit))
  # vdiffr::expect_doppelganger("t_test-two_sample-one_sided-posterior_and_prior", posterior_plot(fit, prior = TRUE))


  ### test unspecified omega -- BFF (one-sided; also change n and k)
  fit <- regression_test_BFF(
    t_stat = 0.5,
    alternative = "two.sided",
    n = 25,
    k = 1)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf_h1, 0.00, tolerance = 1e-5)
  testthat::expect_equal(fit$omega_h1,  0.00)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local regression test"  ,
      ""                                        ,
      "maximized (in favor of alternative) log Bayes factor = 0.00",
      "maximized (in favor of alternative) omega = 0.00 (signed Cohen's f)",
      "minimized (in favor of null for medium/large effect sizes) log Bayes factor = -3.45",
      "minimized (in favor of null for medium/large effect sizes) omega = 1.00 (signed Cohen's f)",
      "alternative = two.sided"
    )
  )
  # vdiffr::expect_doppelganger("t_test_BFF-two_sample-two_sided-BFF",                 plot(fit))
  # vdiffr::expect_doppelganger("t_test_BFF-two_sample-two_sided-posterior",           posterior_plot(fit))
  # vdiffr::expect_doppelganger("t_test_BFF-two_sample-two_sided-posterior_and_prior", posterior_plot(fit, prior = TRUE, color = c("red", "blue"), linetype = c(3,5),
  #                                                                                                   linewidth = c(2, 1), x_limit = c(-2, 2)))


  # check that the data.frame plot output also works
  # no_plot_plot <- posterior_plot(fit, plot = FALSE, prior = TRUE)
  # testthat::expect_true(is.data.frame(no_plot_plot))
  # testthat::expect_equal(colnames(no_plot_plot), c("x", "prior", "posterior"))

  ### test unspecified omega -- BFF (one-sided, also set r)
  fit <- regression_test_BFF(
    t_stat = 0.5,
    alternative = "less",
    r = 3,
    n = 50,
    k = 4)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf_h1, 0.0, tolerance = 1e-5)
  testthat::expect_equal(fit$omega_h1,  0.0)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(summary(fit), print = TRUE, width = 100),
    c(
      "\tBayesian non-local regression test",
      ""                                        ,
      "maximized (in favor of alternative) log Bayes factor = 0.00",
      "maximized (in favor of alternative) omega = 0.00 (signed Cohen's f)",
      "minimized (in favor of null for medium/large effect sizes) log Bayes factor = -8.63",
      "minimized (in favor of null for medium/large effect sizes) omega = 1.00 (signed Cohen's f)",
      "alternative = less"
    )
  )
  # vdiffr::expect_doppelganger("regression_test_BFF-two_sample-one_sided-BFF", plot(fit))
  # testthat::expect_error(posterior_plot(fit), "There is no non-local prior distribution")
})

