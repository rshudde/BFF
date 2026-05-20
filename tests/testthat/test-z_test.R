test_that("two-sample: basic functionality", {

  ### test fixed omega -- non-local prior z-test (two-sided)
  fit <- z_test_BFF(
    z_stat = 1.5,
    alternative = "two.sided",
    n1 = 50,
    n2 = 50,
    omega = 0.5)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf_h1, -0.27839, tolerance = 1e-5)
  testthat::expect_equal(fit$omega_h1,  0.50)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local two-sample z test"  ,
      ""                                        ,
      "log Bayes factor = -0.28"                 ,
      "omega = 0.50 (Cohen's d)"                ,
      "alternative = two.sided"
    )
  )
  testthat::expect_error(plot(fit), "Bayes factor function can be plotted only if a specific omega/tau2 is not user set")

  # vdiffr::expect_doppelganger("t_test-two_sample-two_sided-posterior",           posterior_plot(fit))
  # vdiffr::expect_doppelganger("t_test-two_sample-two_sided-posterior_and_prior", posterior_plot(fit, prior = TRUE))


  ### test fixed omega -- non-local prior z-test (one-sided)
  fit <- z_test_BFF(
    z_stat = 1.5,
    alternative = "greater",
    n1 = 50,
    n2 = 50,
    omega = 0.5)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf_h1, 0.4009388, tolerance = 1e-4)
  testthat::expect_equal(fit$omega_h1,  0.50)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local two-sample z test"  ,
      ""                                        ,
      "log Bayes factor = 0.40"                 ,
      "omega = 0.50 (Cohen's d)"                ,
      "alternative = greater"
    )
  )
  # vdiffr::expect_doppelganger("z_test-two_sample-one_sided-posterior",           posterior_plot(fit))
  # vdiffr::expect_doppelganger("z_test-two_sample-one_sided-posterior_and_prior", posterior_plot(fit, prior = TRUE))


  ### test unspecified omega -- BFF (one-sided; also change n1/n2)
  fit <- z_test_BFF(
    z_stat = 2.5,
    alternative = "two.sided",
    n1 = 25,
    n2 = 75)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf_h1, 2.07857, tolerance = 1e-5)
  testthat::expect_equal(fit$omega_h1,  0.45)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local two-sample z test"  ,
      ""                                        ,
      "maximized (in favor of alternative) log Bayes factor = 2.08",
      "maximized (in favor of alternative) omega = 0.45 (Cohen's d)",
      "minimized (in favor of null for medium/large effect sizes) log Bayes factor = 1.21",
      "minimized (in favor of null for medium/large effect sizes) omega = 1.00 (Cohen's d)",
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

  ### test unspecified omega -- BFF (one-sided)
  fit <- z_test_BFF(
    z_stat = 1.5,
    alternative = "less",
    n1 = 50,
    n2 = 50)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf_h1, 0.00, tolerance = 1e-5)
  testthat::expect_equal(fit$omega_h1,  0.00)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(summary(fit), print = TRUE, width = 100),
    c(
      "\tBayesian non-local two-sample z test",
      ""                                        ,
      "maximized (in favor of alternative) log Bayes factor = 0.00",
      "maximized (in favor of alternative) omega = 0.00 (Cohen's d)",
      "minimized (in favor of null for medium/large effect sizes) log Bayes factor = -5.80",
      "minimized (in favor of null for medium/large effect sizes) omega = 1.00 (Cohen's d)",
      "alternative = less"
    )
  )
  # vdiffr::expect_doppelganger("z_test_BFF-two_sample-one_sided-BFF", plot(fit))
  # testthat::expect_error(posterior_plot(fit), "There is no non-local prior distribution")
})

