test_that("chi2: basic functionality", {

  ### test fixed omega -- non-local prior chi2-test (LRT = TRUE)
  fit <- chi2_test_BFF(
    chi2_stat = 1.5,
    n = 25,
    df = 24,
    LRT = TRUE,
    omega = 0.5)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf_h1, -24.60179, tolerance = 1e-5)
  testthat::expect_equal(fit$omega_h1,  0.5)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local chi2 test"  ,
      ""                                        ,
      "log Bayes factor = -24.60"                 ,
      "omega = 0.50 (RMSES)"
    )
  )
  testthat::expect_error(plot(fit), "Bayes factor function can be plotted only if a specific omega/tau2 is not user set")

  # # vdiffr::expect_doppelganger("t_test-two_sample-two_sided-posterior",           posterior_plot(fit))
  # # vdiffr::expect_doppelganger("t_test-two_sample-two_sided-posterior_and_prior", posterior_plot(fit, prior = TRUE))
  #

  ### test fixed omega -- non-local prior chi2-test (LRT = FALSE)
  fit <- chi2_test_BFF(
    chi2_stat = 5.5,
    n = 25,
    df = 24,
    LRT = FALSE,
    omega = 0.5)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf_h1, -22.76035, tolerance = 1e-5)
  testthat::expect_equal(fit$omega_h1,  0.50)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local chi2 test"  ,
      ""                                        ,
      "log Bayes factor = -22.76"                 ,
      "omega = 0.50 (RMSES)"
    )
  )
  #MODIFY FOR CHI2
  # vdiffr::expect_doppelganger("z_test-two_sample-one_sided-posterior",           posterior_plot(fit))
  # vdiffr::expect_doppelganger("z_test-two_sample-one_sided-posterior_and_prior", posterior_plot(fit, prior = TRUE))


  ### test unspecified omega -- chi2 - test (LRT = FALSE)
  fit <- chi2_test_BFF(
    chi2_stat = 9.5,
    n = 45,
    df = 44,
    LRT = FALSE)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf_h1, 0.00, tolerance = 1e-5)
  testthat::expect_equal(fit$omega_h1,  0.00)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local chi2 test"  ,
      ""                                        ,
      "maximized (in favor of alternative) log Bayes factor = 0.00"       ,
      "maximized (in favor of alternative) omega = 0.00 (RMSES)"       ,
      "minimized (in favor of null for medium/large effect sizes) log Bayes factor = -82.72"       ,
      "minimized (in favor of null for medium/large effect sizes) omega = 1.00 (RMSES)"
    )
  )
  #MODIFY FOR CHI2
  # vdiffr::expect_doppelganger("t_test_BFF-two_sample-two_sided-BFF",                 plot(fit))
  # vdiffr::expect_doppelganger("t_test_BFF-two_sample-two_sided-posterior",           posterior_plot(fit))
  # vdiffr::expect_doppelganger("t_test_BFF-two_sample-two_sided-posterior_and_prior", posterior_plot(fit, prior = TRUE, color = c("red", "blue"), linetype = c(3,5),
  #                                                                                                   linewidth = c(2, 1), x_limit = c(-2, 2)))


  # check that the data.frame plot output also works
  # no_plot_plot <- posterior_plot(fit, plot = FALSE, prior = TRUE)
  # testthat::expect_true(is.data.frame(no_plot_plot))
  # testthat::expect_equal(colnames(no_plot_plot), c("x", "prior", "posterior"))

  ### test unspecified omega -- chi2 test (LRT = TRUE)
  fit <- chi2_test_BFF(
    chi2_stat = 7.5,
    n = 45,
    df = 44,
    LRT = TRUE)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf_h1, 0.00, tolerance = 1e-5)
  testthat::expect_equal(fit$omega_h1,  0.00)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local chi2 test"  ,
      ""                                        ,
      "maximized (in favor of alternative) log Bayes factor = 0.00"       ,
      "maximized (in favor of alternative) omega = 0.00 (RMSES)" ,
      "minimized (in favor of null for medium/large effect sizes) log Bayes factor = -83.73",
      "minimized (in favor of null for medium/large effect sizes) omega = 1.00 (RMSES)"
    )
  )
  #MODIFY FOR CHI2
 # vdiffr::expect_doppelganger("z_test_BFF-two_sample-one_sided-BFF", plot(fit))
 # testthat::expect_error(posterior_plot(fit), "There is no non-local prior distribution")

  #Checking with different r
  fit <- chi2_test_BFF(
    chi2_stat = 7.5,
    n = 45,
    df = 44,
    LRT = FALSE,
    r = 5)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf_h1, 0.00, tolerance = 1e-5)
  testthat::expect_equal(fit$omega_h1,  0.00)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local chi2 test"  ,
      ""                                        ,

      "maximized (in favor of alternative) log Bayes factor = 0.00",
      "maximized (in favor of alternative) omega = 0.00 (RMSES)",
      "minimized (in favor of null for medium/large effect sizes) log Bayes factor = -94.06",
      "minimized (in favor of null for medium/large effect sizes) omega = 1.00 (RMSES)"
    )
  )
})

test_that("chi-square table dimensions must match degrees of freedom", {
  testthat::expect_error(
    chi2_test_BFF(
      chi2_stat = 12,
      n = 60,
      df = 5,
      omega = 0.2,
      effect_size = "cramers_v",
      table_dim = c(3, 4)
    ),
    "`df` must equal"
  )
})

test_that("statistic-level chi-square BFF rejects 2x2-only effect sizes", {
  testthat::expect_error(
    chi2_test_BFF(
      chi2_stat = 4.5,
      n = 100,
      df = 1,
      omega = 0.8,
      effect_size = "logOR"
    ),
    "Unsupported effect size for chi-square"
  )
})

test_that("vectorized chi-square test uses study-specific degrees of freedom", {
  fit <- chi2_test_BFF(
    chi2_stat = c(6, 12),
    n = c(80, 120),
    df = c(3, 5),
    omega = 0.2
  )

  expected <- sum(
    chi2_test_BFF(chi2_stat = 6, n = 80, df = 3, omega = 0.2)$log_bf_h1,
    chi2_test_BFF(chi2_stat = 12, n = 120, df = 5, omega = 0.2)$log_bf_h1
  )

  testthat::expect_equal(fit$log_bf_h1, expected, tolerance = 1e-10)
})

test_that("vectorized chi-square test handles conventional scales only with common df", {
  fit <- chi2_test_BFF(
    chi2_stat = c(6, 12),
    n = c(80, 120),
    df = c(3, 3),
    omega = 0.4,
    effect_size = "cohens_w"
  )

  expected <- sum(
    chi2_test_BFF(chi2_stat = 6, n = 80, df = 3, omega = 0.4, effect_size = "cohens_w")$log_bf_h1,
    chi2_test_BFF(chi2_stat = 12, n = 120, df = 3, omega = 0.4, effect_size = "cohens_w")$log_bf_h1
  )

  testthat::expect_length(fit$log_bf_h1, 1)
  testthat::expect_equal(fit$log_bf_h1, expected, tolerance = 1e-10)

  testthat::expect_error(
    chi2_test_BFF(
      chi2_stat = c(6, 12),
      n = c(80, 120),
      df = c(3, 5),
      omega = 0.4,
      effect_size = "cohens_w"
    ),
    "requires a common `df`"
  )
})

test_that("chi-square test rejects invalid inputs", {
  testthat::expect_error(
    chi2_test_BFF(chi2_stat = -1, n = 25, df = 5, omega = 0.2)
  )
  testthat::expect_error(
    chi2_test_BFF(chi2_stat = 5, n = 25, df = 0, omega = 0.2)
  )
  testthat::expect_error(
    chi2_test_BFF(
      chi2_stat = c(5, 8),
      n = c(25, 30),
      df = c(5, 6, 7),
      omega = 0.2
    )
  )
})
