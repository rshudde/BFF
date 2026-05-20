test_that("one-sided t and regression BFFs stay finite with large degrees of freedom", {
  testthat::expect_warning(
    t_fit <- t_test_BFF(
      t_stat = 2.5,
      n = 10000,
      one_sample = TRUE,
      alternative = "greater",
      omega = 0.2
    ),
    NA
  )
  testthat::expect_true(is.finite(t_fit$log_bf_h1))

  testthat::expect_warning(
    reg_fit <- regression_test_BFF(
      t_stat = 2.5,
      n = 10000,
      k = 3,
      alternative = "greater",
      omega = 0.2
    ),
    NA
  )
  testthat::expect_true(is.finite(reg_fit$log_bf_h1))
})

test_that("fractional-r chi-square BFF matches quadrature and normalizes the posterior", {
  cases <- list(
    list(chi2_stat = 6.5, n = 80, df = 3, omega = 0.2, r = 1.5),
    list(chi2_stat = 18, n = 120, df = 5, omega = 0.25, r = 2.25)
  )

  for(case in cases){
    tau2 <- BFF:::get_count_tau2(
      n = case$n,
      k = case$df,
      w = case$omega,
      r = case$r
    )
    prior <- function(effect_size) BFF:::.chi2_test.prior(
      tau2 = tau2,
      r = case$r,
      effect_size = effect_size,
      n = case$n,
      df = case$df
    )
    marginal <- integrate_density(
      f = function(effect_size) stats::dchisq(
        x = case$chi2_stat,
        df = case$df,
        ncp = BFF:::.chi2_test_ncp(
          effect_size = effect_size,
          n = case$n,
          df = case$df
        )
      ) * prior(effect_size),
      lower = 0
    )
    closed_form <- exp(BFF:::BFF_chi2_test(
      tau2 = tau2,
      chi2_stat = case$chi2_stat,
      k = case$df,
      r = case$r
    ))
    posterior_area <- integrate_density(
      f = function(effect_size) BFF:::.chi2_test.posterior(
        chi2_stat = case$chi2_stat,
        tau2 = tau2,
        r = case$r,
        effect_size = effect_size,
        n = case$n,
        df = case$df
      ),
      lower = 0
    )

    testthat::expect_equal(
      closed_form,
      marginal / stats::dchisq(case$chi2_stat, df = case$df),
      tolerance = 1e-6,
      info = paste("r =", case$r)
    )
    testthat::expect_equal(
      posterior_area,
      1,
      tolerance = 1e-5,
      info = paste("r =", case$r)
    )
  }
})

test_that("nonfinite omega values are rejected for F and regression scales", {
  testthat::expect_error(
    f_test_BFF(f_stat = 1.5, n = 50, df1 = 5, df2 = 40, omega = Inf),
    "omega"
  )
  testthat::expect_error(
    f_test_BFF(
      f_stat = 1.5,
      n = 50,
      df1 = 5,
      df2 = 40,
      omega_sequence = c(0.1, NA_real_),
      effect_size = "cohens_f2"
    ),
    "omega"
  )
  testthat::expect_error(
    regression_test_BFF(t_stat = 2, n = 50, k = 3, omega = NaN),
    "omega"
  )
  testthat::expect_error(
    regression_test_BFF(
      t_stat = 2,
      n = 50,
      k = 3,
      omega_sequence = c(0.1, Inf),
      effect_size = "partial_r"
    ),
    "omega"
  )
})

test_that("transformed effect-size scales reject nonfinite omega values", {
  testthat::expect_error(
    chi2_test_BFF(
      chi2_stat = 12,
      n = 60,
      df = 6,
      omega = Inf,
      effect_size = "cohens_w"
    ),
    "omega"
  )
  testthat::expect_error(
    prop_test_BFF(
      x = c(30, 20),
      n = c(100, 100),
      omega_sequence = c(0.2, NaN),
      effect_size = "logOR"
    ),
    "omega"
  )
})

test_that("F-test exposes tau2 using the fitted effect-size prior-mode convention", {
  fit <- f_test_BFF(
    f_stat = 1.75,
    n = 80,
    df1 = 5,
    df2 = 74,
    omega = 0.3,
    r = 2
  )

  expected_tau2 <- 80 * 5 * 0.3^2 / (2 * (5 + 2 * 2 - 1))
  testthat::expect_equal(unname(fit$tau2_h1[[1]]), expected_tau2, tolerance = 1e-12)
})

test_that("2x2 signed posterior plots preserve fitted direction across scales", {
  fit <- prop_test_BFF(
    x = c(20, 30),
    n = c(100, 100),
    omega = -0.4,
    effect_size = "logOR"
  )
  plot_data <- posterior_plot(
    fit,
    prior = TRUE,
    plot = FALSE,
    effect_size = "risk_difference",
    x_limit = c(-0.8, 0.8)
  )

  expect_posterior_plot_data(plot_data)
  testthat::expect_true(any(plot_data$prior[plot_data$x < 0] > 0))
  testthat::expect_true(any(plot_data$posterior[plot_data$x < 0] > 0))
  testthat::expect_true(all(plot_data$prior[plot_data$x > 0] == 0))
  testthat::expect_true(all(plot_data$posterior[plot_data$x > 0] == 0))
})
