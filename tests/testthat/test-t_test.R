test_that("two-sample: basic functionality", {

  ### test fixed omega -- non-local prior t-test (two-sided)
  fit <- t_test_BFF(
    t_stat = 1.5,
    alternative = "two.sided",
    n1 = 50,
    n2 = 50,
    omega = 0.5)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf, 0.1172115, tolerance = 1e-5)
  testthat::expect_equal(fit$omega,  0.5)

  # test S3 methods
  testthat::expect_equal(
    capture_output_lines(fit, print = TRUE, width = 100),
    c(
    "\tBayesian non-local two-sample t test"  ,
    ""                                        ,
    "log Bayes factor = 0.12"                 ,
    "omega = 0.50 (Cohen's d)"                ,
    "alternative = two.sided"
    )
  )
  testthat::expect_error(plot(fit), "Bayes factor function can be plotted only if a specific omega/tau2 is not user set")
  vdiffr::expect_doppelganger("t_test-two_sample-two_sided-posterior",           posterior_plot(fit))
  vdiffr::expect_doppelganger("t_test-two_sample-two_sided-posterior_and_prior", posterior_plot(fit, prior = TRUE))


  ### test fixed omega -- non-local prior t-test (one-sided, also set r)
  fit <- t_test_BFF(
    t_stat = 1.5,
    alternative = "greater",
    n1 = 50,
    n2 = 50,
    r = 2,
    omega = 0.5)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf, 0.4952852, tolerance = 1e-5)
  testthat::expect_equal(fit$omega,  0.50)

  # test S3 methods
  testthat::expect_equal(
    capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local two-sample t test"  ,
      ""                                        ,
      "log Bayes factor = 0.50"                 ,
      "omega = 0.50 (Cohen's d)"                ,
      "alternative = greater"
    )
  )
  vdiffr::expect_doppelganger("t_test-two_sample-one_sided-posterior",           posterior_plot(fit))
  vdiffr::expect_doppelganger("t_test-two_sample-one_sided-posterior_and_prior", posterior_plot(fit, prior = TRUE))


  ### test unspecified omega -- BFF (one-sided; also change n1/n2)
  fit <- t_test_BFF(
    t_stat = 0.5,
    alternative = "two.sided",
    n1 = 25,
    n2 = 75)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf, 0.03230874, tolerance = 1e-5)
  testthat::expect_equal(fit$omega,  0.05)

  # test S3 methods
  testthat::expect_equal(
    capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local two-sample t test"  ,
      ""                                        ,
      "maximized log Bayes factor = 0.03"       ,
      "maximized omega = 0.05 (Cohen's d)"      ,
      "alternative = two.sided"
    )
  )
  vdiffr::expect_doppelganger("t_test_BFF-two_sample-two_sided-BFF",                 plot(fit))
  vdiffr::expect_doppelganger("t_test_BFF-two_sample-two_sided-posterior",           posterior_plot(fit))
  vdiffr::expect_doppelganger("t_test_BFF-two_sample-two_sided-posterior_and_prior", posterior_plot(fit, prior = TRUE, color = c("red", "blue"), linetype = c(3,5),
                                                                                                    linewidth = c(2, 1), x_limit = c(-2, 2)))


  # check that the data.frame plot output also works
  no_plot_plot <- posterior_plot(fit, plot = FALSE, prior = TRUE)
  testthat::expect_true(is.data.frame(no_plot_plot))
  testthat::expect_equal(colnames(no_plot_plot), c("x", "prior", "posterior"))

  ### test unspecified omega -- BFF (one-sided, also set r)
  fit <- t_test_BFF(
    t_stat = 0.5,
    alternative = "less",
    r = 3,
    n1 = 50,
    n2 = 50)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf, 0, tolerance = 1e-5)
  testthat::expect_equal(fit$omega,  0)

  # test S3 methods
  testthat::expect_equal(
    capture_output_lines(summary(fit), print = TRUE, width = 100),
    c(
      "\tBayesian non-local two-sample t test",
      ""                                        ,
      "maximized log Bayes factor = 0.00"       ,
      "maximized omega = 0.00 (Cohen's d)"      ,
      "alternative = less"
    )
  )
  vdiffr::expect_doppelganger("t_test_BFF-two_sample-one_sided-BFF", plot(fit))
  testthat::expect_error(posterior_plot(fit), "There is no non-local prior distribution")
})
