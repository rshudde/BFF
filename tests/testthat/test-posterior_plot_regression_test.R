test_that("regression-test r=1 two-sided shortcut matches the closed-form implementation", {
  grid <- expand.grid(
    tau2 = c(0.1, 0.5, 2),
    t_stat = c(-2, 2),
    df = c(20, 50)
  )

  for(i in seq_len(nrow(grid))){
    testthat::expect_equal(
      BFF:::BFF_reg_test(
        tau2 = grid$tau2[i],
        t_stat = grid$t_stat[i],
        df = grid$df[i],
        r = 1,
        two_sided = TRUE
      ),
      BFF:::reg_t_val_r1(
        tau2 = grid$tau2[i],
        t_stat = grid$t_stat[i],
        df = grid$df[i]
      ),
      tolerance = 1e-10
    )
  }
})

test_that("regression-test prior and posterior densities integrate to one", {
  cases <- list(
    list(
      label = "two-sided r1",
      t_stat = 1.5, n = 50, k = 3,
      one_sided = FALSE, omega = 0.5, r = 1,
      lower = -Inf
    ),
    list(
      label = "greater r3",
      t_stat = 1.5, n = 50, k = 2,
      one_sided = TRUE, omega = 0.5, r = 3,
      lower = 0
    ),
    list(
      label = "greater r2",
      t_stat = 2.4, n = 80, k = 5,
      one_sided = TRUE, omega = 0.25, r = 2,
      lower = 0
    )
  )

  for(case in cases){
    tau2 <- BFF:::get_regression_tau2(n = case$n, k = case$k, w = case$omega, r = case$r)

    prior_integral <- integrate_density(
      f = function(delta) BFF:::.regression_test.prior(
        tau2 = tau2,
        r = case$r,
        effect_size = delta,
        n = case$n,
        k = case$k,
        one_sided = case$one_sided
      ),
      lower = case$lower
    )
    posterior_integral <- integrate_density(
      f = function(delta) BFF:::.regression_test.posterior(
        t_stat = case$t_stat,
        tau2 = tau2,
        r = case$r,
        effect_size = delta,
        n = case$n,
        k = case$k,
        one_sided = case$one_sided
      ),
      lower = case$lower
    )

    testthat::expect_equal(prior_integral, 1, tolerance = 1e-6, info = case$label)
    testthat::expect_equal(posterior_integral, 1, tolerance = 1e-5, info = case$label)
  }
})

test_that("posterior_plot works for a selected regression-test non-local prior", {
  fixed_fit <- regression_test_BFF(
    t_stat = 1.5,
    alternative = "two.sided",
    n = 50,
    k = 3,
    omega = 0.5
  )
  bff_fit <- regression_test_BFF(
    t_stat = 2.5,
    alternative = "two.sided",
    n = 50,
    k = 3
  )

  expect_posterior_plot_data(posterior_plot(fixed_fit, plot = FALSE, prior = TRUE))
  expect_posterior_plot_data(posterior_plot(bff_fit, plot = FALSE, prior = TRUE))
  testthat::expect_true(inherits(posterior_plot(fixed_fit), "ggplot"))
})

test_that("posterior_plot respects less-than regression-test support", {
  fit <- regression_test_BFF(
    t_stat = -2.5,
    n = 50,
    k = 3,
    alternative = "less",
    omega = 0.4
  )

  default_data <- posterior_plot(fit, plot = FALSE, prior = TRUE)
  signed_data <- posterior_plot(fit, plot = FALSE, prior = TRUE, x_limit = c(-1, 1))

  testthat::expect_equal(range(default_data$x), c(-3, 0))
  testthat::expect_true(any(signed_data$prior[signed_data$x < 0] > 0))
  testthat::expect_true(any(signed_data$posterior[signed_data$x < 0] > 0))
  testthat::expect_true(all(signed_data$prior[signed_data$x > 0] == 0))
  testthat::expect_true(all(signed_data$posterior[signed_data$x > 0] == 0))
})

test_that("regression-test two-sided BFF is sign symmetric", {
  fit_pos <- regression_test_BFF(
    t_stat = 2.5,
    alternative = "two.sided",
    n = 50,
    k = 3,
    omega = 0.4
  )
  fit_neg <- regression_test_BFF(
    t_stat = -2.5,
    alternative = "two.sided",
    n = 50,
    k = 3,
    omega = 0.4
  )

  testthat::expect_equal(fit_pos$log_bf_h1, fit_neg$log_bf_h1, tolerance = 1e-12)
})

test_that("posterior_plot rejects ambiguous or null selected regression-test priors", {
  null_fit <- regression_test_BFF(
    t_stat = 0.5,
    alternative = "less",
    r = 3,
    n = 50,
    k = 4
  )
  multi_fit <- regression_test_BFF(
    t_stat = 1.5,
    n = 50,
    k = 3,
    omega = c(0.3, 0.5)
  )
  vector_fit <- regression_test_BFF(
    t_stat = c(2.5, 2),
    n = c(50, 60),
    k = c(3, 4),
    omega = 0.3
  )

  testthat::expect_equal(null_fit$omega_h1, 0)
  testthat::expect_error(
    posterior_plot(null_fit),
    "There is no non-local prior distribution"
  )
  testthat::expect_error(
    posterior_plot(multi_fit),
    "requires a single selected omega"
  )
  testthat::expect_error(
    posterior_plot(vector_fit),
    "single t statistic"
  )
})
