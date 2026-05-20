test_that("z-test prior and posterior densities integrate to one", {
  cases <- list(
    list(
      label = "one-sample two-sided",
      z_stat = 2.1, n = 40, n1 = NULL, n2 = NULL,
      one_sample = TRUE, one_sided = FALSE, omega = 0.4, r = 1,
      lower = -Inf
    ),
    list(
      label = "one-sample one-sided",
      z_stat = 2.1, n = 60, n1 = NULL, n2 = NULL,
      one_sample = TRUE, one_sided = TRUE, omega = 0.3, r = 2,
      lower = 0
    ),
    list(
      label = "two-sample two-sided",
      z_stat = 1.5, n = NULL, n1 = 50, n2 = 50,
      one_sample = FALSE, one_sided = FALSE, omega = 0.5, r = 1,
      lower = -Inf
    ),
    list(
      label = "two-sample one-sided",
      z_stat = 2.5, n = NULL, n1 = 25, n2 = 75,
      one_sample = FALSE, one_sided = TRUE, omega = 0.35, r = 2,
      lower = 0
    )
  )

  for(case in cases){
    tau2 <- if(case$one_sample){
      BFF:::get_one_sample_tau2(n = case$n, w = case$omega, r = case$r)
    }else{
      BFF:::get_two_sample_tau2(n1 = case$n1, n2 = case$n2, w = case$omega, r = case$r)
    }

    prior_integral <- integrate_density(
      f = function(delta) BFF:::.z_test.prior(
        tau2 = tau2,
        r = case$r,
        effect_size = delta,
        n = case$n,
        n1 = case$n1,
        n2 = case$n2,
        one_sample = case$one_sample,
        one_sided = case$one_sided
      ),
      lower = case$lower
    )
    posterior_integral <- integrate_density(
      f = function(delta) BFF:::.z_test.posterior(
        z_stat = case$z_stat,
        tau2 = tau2,
        r = case$r,
        effect_size = delta,
        n = case$n,
        n1 = case$n1,
        n2 = case$n2,
        one_sample = case$one_sample,
        one_sided = case$one_sided
      ),
      lower = case$lower
    )

    testthat::expect_equal(prior_integral, 1, tolerance = 1e-6, info = case$label)
    testthat::expect_equal(posterior_integral, 1, tolerance = 1e-5, info = case$label)
  }
})

test_that("posterior_plot works for a selected z-test non-local prior", {
  fixed_fit <- z_test_BFF(
    z_stat = 1.5,
    alternative = "two.sided",
    n1 = 50,
    n2 = 50,
    omega = 0.5
  )
  bff_fit <- z_test_BFF(
    z_stat = 2.5,
    alternative = "two.sided",
    n1 = 25,
    n2 = 75
  )

  expect_posterior_plot_data(posterior_plot(fixed_fit, plot = FALSE, prior = TRUE))
  expect_posterior_plot_data(posterior_plot(bff_fit, plot = FALSE, prior = TRUE))
  testthat::expect_true(inherits(posterior_plot(fixed_fit), "ggplot"))
})

test_that("posterior_plot respects less-than z-test support", {
  fit <- z_test_BFF(
    z_stat = -2.5,
    n = 40,
    one_sample = TRUE,
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

test_that("posterior_plot rejects ambiguous or null selected z-test priors", {
  null_fit <- z_test_BFF(
    z_stat = 1.5,
    alternative = "less",
    n1 = 50,
    n2 = 50
  )
  multi_fit <- z_test_BFF(
    z_stat = 2.5,
    n = 50,
    one_sample = TRUE,
    omega = c(0.3, 0.5)
  )
  vector_fit <- z_test_BFF(
    z_stat = c(2.5, 2),
    n = c(50, 60),
    one_sample = TRUE,
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
    "single z statistic"
  )
})
