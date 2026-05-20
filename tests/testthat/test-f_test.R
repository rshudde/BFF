test_that("two-sample: basic functionality", {

  ### test fixed omega -- non-local prior f-test (two-sample)
  fit <- f_test_BFF(
    f_stat = 1.75,
    n = 25,
    df1 = 50,
    df2 = 50,
    omega = 0.5)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf_h1, -2.73527, tolerance = 1e-5)
  testthat::expect_equal(fit$omega_h1,  0.5)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local f test"  ,
      ""                                        ,
      "log Bayes factor = -2.74"                 ,
      "omega = 0.50 (RMSES)"
    )
  )
  testthat::expect_error(plot(fit), "Bayes factor function can be plotted only if a specific omega/tau2 is not user set")

  # # vdiffr::expect_doppelganger("t_test-two_sample-two_sided-posterior",           posterior_plot(fit))
  # # vdiffr::expect_doppelganger("t_test-two_sample-two_sided-posterior_and_prior", posterior_plot(fit, prior = TRUE))
  #


  ### test unspecified omega -- BFF (two-sample; also change n1/n2)
  fit <- f_test_BFF(
    f_stat = 1.5,
    n = 50,
    df1 = 25,
    df2 = 75)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf_h1, 0.82374, tolerance = 1e-5)
  testthat::expect_equal(fit$omega_h1,  0.14)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local f test"  ,
      ""                                        ,
      "maximized (in favor of alternative) log Bayes factor = 0.82",
      "maximized (in favor of alternative) omega = 0.14 (RMSES)",
      "minimized (in favor of null for medium/large effect sizes) log Bayes factor = -23.13",
      "minimized (in favor of null for medium/large effect sizes) omega = 1.00 (RMSES)"
    )
  )
  #Modify for F test
  # vdiffr::expect_doppelganger("f_test_BFF-two_sample-two_sided-BFF",                 plot(fit))
  # vdiffr::expect_doppelganger("f_test_BFF-two_sample-two_sided-posterior",           posterior_plot(fit))
  # vdiffr::expect_doppelganger("f_test_BFF-two_sample-two_sided-posterior_and_prior", posterior_plot(fit, prior = TRUE, color = c("red", "blue"), linetype = c(3,5),
  #                                                                                                   linewidth = c(2, 1), x_limit = c(-2, 2)))


  # check that the data.frame plot output also works
  # no_plot_plot <- posterior_plot(fit, plot = FALSE, prior = TRUE)
  # testthat::expect_true(is.data.frame(no_plot_plot))
  # testthat::expect_equal(colnames(no_plot_plot), c("x", "prior", "posterior"))


  ### test different r -- non-local prior f-test (two-sample)
  fit <- f_test_BFF(
    f_stat = 1.75,
    n = 25,
    df1 = 50,
    df2 = 50,
    r = 3)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf_h1, 1.93649, tolerance = 1e-5)
  testthat::expect_equal(fit$omega_h1,  0.24)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local f test"  ,
      ""                                        ,
      "maximized (in favor of alternative) log Bayes factor = 1.94",
      "maximized (in favor of alternative) omega = 0.24 (RMSES)",
      "minimized (in favor of null for medium/large effect sizes) log Bayes factor = -22.46",
      "minimized (in favor of null for medium/large effect sizes) omega = 1.00 (RMSES)"
    )
  )

})

test_that("vectorized F-test uses study-specific degrees of freedom", {
  fit <- f_test_BFF(
    f_stat = c(1.4, 2.1),
    n = c(40, 65),
    df1 = c(2, 4),
    df2 = c(37, 60),
    omega = 0.25
  )

  expected <- sum(
    f_test_BFF(f_stat = 1.4, n = 40, df1 = 2, df2 = 37, omega = 0.25)$log_bf_h1,
    f_test_BFF(f_stat = 2.1, n = 65, df1 = 4, df2 = 60, omega = 0.25)$log_bf_h1
  )

  testthat::expect_equal(fit$log_bf_h1, expected, tolerance = 1e-10)
})

test_that("vectorized F-test handles conventional scales only with common numerator df", {
  fit <- f_test_BFF(
    f_stat = c(1.4, 2.1),
    n = c(40, 65),
    df1 = c(2, 2),
    df2 = c(37, 62),
    omega = 0.4,
    effect_size = "cohens_f"
  )

  expected <- sum(
    f_test_BFF(f_stat = 1.4, n = 40, df1 = 2, df2 = 37, omega = 0.4, effect_size = "cohens_f")$log_bf_h1,
    f_test_BFF(f_stat = 2.1, n = 65, df1 = 2, df2 = 62, omega = 0.4, effect_size = "cohens_f")$log_bf_h1
  )

  testthat::expect_length(fit$log_bf_h1, 1)
  testthat::expect_equal(fit$log_bf_h1, expected, tolerance = 1e-10)

  testthat::expect_error(
    f_test_BFF(
      f_stat = c(1.4, 2.1),
      n = c(40, 65),
      df1 = c(2, 4),
      df2 = c(37, 60),
      omega = 0.4,
      effect_size = "cohens_f"
    ),
    "requires a common `df1`"
  )
})

test_that("F-test rejects invalid inputs", {
  testthat::expect_error(
    f_test_BFF(f_stat = -1, n = 25, df1 = 5, df2 = 20, omega = 0.2)
  )
  testthat::expect_error(
    f_test_BFF(f_stat = 1.5, n = 25, df1 = 0, df2 = 20, omega = 0.2)
  )
  testthat::expect_error(
    f_test_BFF(f_stat = 1.5, n = 25, df1 = 5, df2 = 0, omega = 0.2)
  )
  testthat::expect_error(
    f_test_BFF(
      f_stat = c(1.5, 2.0),
      n = c(25, 30),
      df1 = c(5, 6),
      df2 = c(20, 21, 22),
      omega = 0.2
    )
  )
})
