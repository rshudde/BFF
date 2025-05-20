test_that("two-sample: basic functionality", {

  ### test fixed omega -- non-local prior t-test (two-sided)
  fit <- t_test_BFF(
    t_stat = 2.5,
    alternative = "two.sided",
    n1 = 50,
    n2 = 50,
    omega = 0.5)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf_h1, 1.917536, tolerance = 1e-5)
  testthat::expect_equal(fit$omega_h1,  0.5)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(fit, print = TRUE, width = 100),
    c(
    "\tBayesian non-local two-sample t test"  ,
    ""                                        ,
    "log Bayes factor = 1.92"                 ,
    "omega = 0.50 (Cohen's d)"                ,
    "alternative = two.sided"
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
  # tau2 <- get_two_sample_tau2(n1 = fit$input$n1, n2 = fit$input$n2, w = fit$omega_h1, r = fit$r)
  #
  # # does not integrate to 1
  # integrate(
  #   f = function(x) .t_test.posterior(
  #     t_stat = fit$input$t_stat, tau2 = tau2, r = fit$r, effect_size = x,
  #     n = fit$input$n, n1 = fit$input$n1, n2 = fit$input$n2, one_sample = fit$one_sample, one_sided = fit$alternative != "two.sided"),
  #   lower = -Inf,
  #   upper = Inf
  # )
  #
  # # prior seems to work just fine (i.e., integrates to one)
  # integrate(
  #   f = function(x) .t_test.prior(
  #     tau2 = tau2, r = fit$r, effect_size = x,
  #     n = fit$input$n, n1 = fit$input$n1, n2 = fit$input$n2, one_sample = fit$one_sample, one_sided = fit$alternative != "two.sided"),
  #   lower = -Inf,
  #   upper = Inf
  # )
  # # <\TODO>
  #
  # # vdiffr::expect_doppelganger("t_test-two_sample-two_sided-posterior",           posterior_plot(fit))
  # # vdiffr::expect_doppelganger("t_test-two_sample-two_sided-posterior_and_prior", posterior_plot(fit, prior = TRUE))
  #

  ### test fixed omega -- non-local prior t-test (one-sided, also set r)
  fit <- t_test_BFF(
    t_stat = 2.5,
    alternative = "greater",
    n1 = 50,
    n2 = 50,
    r = 2,
    omega = 0.5)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf_h1, 2.7725, tolerance = 1e-5)
  testthat::expect_equal(fit$omega_h1,  0.50)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local two-sample t test"  ,
      ""                                        ,
      "log Bayes factor = 2.77"                 ,
      "omega = 0.50 (Cohen's d)"                ,
      "alternative = greater"
    )
  )
  # vdiffr::expect_doppelganger("t_test-two_sample-one_sided-posterior",           posterior_plot(fit))
  # vdiffr::expect_doppelganger("t_test-two_sample-one_sided-posterior_and_prior", posterior_plot(fit, prior = TRUE))


  ### test unspecified omega -- BFF (one-sided; also change n1/n2)
  fit <- t_test_BFF(
    t_stat = 0.5,
    alternative = "two.sided",
    n1 = 25,
    n2 = 75)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf_h1, 0.0, tolerance = 1e-2)
  testthat::expect_equal(fit$omega_h1,  0.00)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(fit, print = TRUE, width = 100),
    c(
      "\tBayesian non-local two-sample t test"  ,
      ""                                        ,
      "maximized (in favor of alternative) log Bayes factor = 0.00",
      "maximized (in favor of alternative) omega = 0.00 (Cohen's d)",
      "minimized (in favor of null for medium/large effect sizes) log Bayes factor = -3.19",
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

  ### test unspecified omega -- BFF (one-sided, also set r)
  fit <- t_test_BFF(
    t_stat = 0.5,
    alternative = "less",
    r = 3,
    n1 = 50,
    n2 = 50)

  # check that the BF and omega is consistent
  testthat::expect_equal(fit$log_bf_h1, 0.00, tolerance = 1e-5)
  testthat::expect_equal(fit$omega_h1,  0.00)

  # test S3 methods
  testthat::expect_equal(
    testthat::capture_output_lines(summary(fit), print = TRUE, width = 100),
    c(
      "\tBayesian non-local two-sample t test",
      ""                                        ,
      "maximized (in favor of alternative) log Bayes factor = 0.00",
      "maximized (in favor of alternative) omega = 0.00 (Cohen's d)",
      "minimized (in favor of null for medium/large effect sizes) log Bayes factor = -6.85",
      "minimized (in favor of null for medium/large effect sizes) omega = 1.00 (Cohen's d)",
      "alternative = less"
    )
  )
  # vdiffr::expect_doppelganger("t_test_BFF-two_sample-one_sided-BFF", plot(fit))
  # testthat::expect_error(posterior_plot(fit), "There is no non-local prior distribution")
})

