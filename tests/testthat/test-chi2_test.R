test_that("chi2: basic functionality", {

  ### test fixed omega -- non-local prior chi2-test (LRT = TRUE)
  fit <- chi2_test_BFF(
    chi2_stat = 1.5,
    n = 25,
    LRT = TRUE,
    omega = 0.5)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf, -25.6263, tolerance = 1e-5)
  testthat::expect_equal(fit$omega,  0.5)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local chi2 test"  ,
      ""                                        ,
      "log Bayes factor = -25.63"                 ,
      "omega = 0.50 (Standardized effect size)"
    )
  )
  testthat::expect_error(plot(fit), "Bayes factor function can be plotted only if a specific omega/tau2 is not user set")

  # TODO: fix posterior plots
  # - I fixed the arguments not being properly passed
  # - however, the posterior distribution does not integrate to 1
  # (I remember that I raised this issue when I was in US, and Saptati was working on fixing it)

  # this is how the functions should work
  # posterior_plot(fit)
  # posterior_plot(fit, prior = TRUE)
  #
  # # this highlights the issue (run `devtools::load_all()` first)
  # tau2 <- get_two_sample_tau2(n1 = fit$input$n1, n2 = fit$input$n2, w = fit$omega, r = fit$r)
  #
  # # does not integrate to 1
  # integrate(
  #   f = function(x) .z_test.posterior(
  #     z_stat = fit$input$z_stat, tau2 = tau2, r = fit$r, effect_size = x,
  #     n = fit$input$n, n1 = fit$input$n1, n2 = fit$input$n2, one_sample = fit$one_sample, one_sided = fit$alternative != "two.sided"),
  #   lower = -Inf,
  #   upper = Inf
  # )
  #
  # # prior seems to work just fine (i.e., integrates to one)
  # integrate(
  #   f = function(x) .z_test.prior(
  #     tau2 = tau2, r = fit$r, effect_size = x,
  #     n = fit$input$n, n1 = fit$input$n1, n2 = fit$input$n2, one_sample = fit$one_sample, one_sided = fit$alternative != "two.sided"),
  #   lower = -Inf,
  #   upper = Inf
  # )
  # #MODIFY FOR CHI2
  # # <\TODO>
  #
  # # vdiffr::expect_doppelganger("t_test-two_sample-two_sided-posterior",           posterior_plot(fit))
  # # vdiffr::expect_doppelganger("t_test-two_sample-two_sided-posterior_and_prior", posterior_plot(fit, prior = TRUE))
  #

  ### test fixed omega -- non-local prior chi2-test (LRT = FALSE)
  fit <- chi2_test_BFF(
    chi2_stat = 5.5,
    n = 25,
    LRT = FALSE,
    omega = 0.5)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf, -25.29, tolerance = 1e-2)
  testthat::expect_equal(fit$omega,  0.50)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local chi2 test"  ,
      ""                                        ,
      "log Bayes factor = -25.29"                 ,
      "omega = 0.50 (Standardized effect size)"
    )
  )
  #MODIFY FOR CHI2
  # vdiffr::expect_doppelganger("z_test-two_sample-one_sided-posterior",           posterior_plot(fit))
  # vdiffr::expect_doppelganger("z_test-two_sample-one_sided-posterior_and_prior", posterior_plot(fit, prior = TRUE))


  ### test unspecified omega -- chi2 - test (LRT = FALSE)
  fit <- chi2_test_BFF(
    chi2_stat = 9.5,
    n = 45,
    LRT = FALSE)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf, 0.00, tolerance = 1e-5)
  testthat::expect_equal(fit$omega,  0.00)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local chi2 test"  ,
      ""                                        ,
      "maximized log Bayes factor = 0.00"       ,
      "maximized omega = 0.00 (Standardized effect size)"
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
    LRT = TRUE)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf, 0.00, tolerance = 1e-5)
  testthat::expect_equal(fit$omega,  0.00)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local chi2 test"  ,
      ""                                        ,
      "maximized log Bayes factor = 0.00"       ,
      "maximized omega = 0.00 (Standardized effect size)"
    )
  )
  #MODIFY FOR CHI2
 # vdiffr::expect_doppelganger("z_test_BFF-two_sample-one_sided-BFF", plot(fit))
 # testthat::expect_error(posterior_plot(fit), "There is no non-local prior distribution")

  #Checking with different r
  fit <- chi2_test_BFF(
    chi2_stat = 7.5,
    n = 45,
    LRT = FALSE,
    r = 5)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf, 0.00, tolerance = 1e-5)
  testthat::expect_equal(fit$omega,  0.00)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local chi2 test"  ,
      ""                                        ,
      "maximized log Bayes factor = 0.00"       ,
      "maximized omega = 0.00 (Standardized effect size)"
    )
  )
})
