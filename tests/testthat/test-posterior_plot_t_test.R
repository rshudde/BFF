integrate_t_density <- function(f, lower, upper){
  suppressWarnings(stats::integrate(
    f = f,
    lower = lower,
    upper = upper,
    rel.tol = 1e-7,
    subdivisions = 1000
  )$value)
}

test_that("t-test prior and posterior densities integrate to one", {
  cases <- list(
    list(
      label = "one-sample two-sided",
      t_stat = 2.5, n = 50, n1 = NULL, n2 = NULL,
      one_sample = TRUE, one_sided = FALSE, omega = 0.28, r = 1,
      lower = -Inf
    ),
    list(
      label = "one-sample one-sided",
      t_stat = 2.5, n = 50, n1 = NULL, n2 = NULL,
      one_sample = TRUE, one_sided = TRUE, omega = 0.5, r = 2,
      lower = 0
    ),
    list(
      label = "two-sample two-sided",
      t_stat = 2.5, n = NULL, n1 = 50, n2 = 50,
      one_sample = FALSE, one_sided = FALSE, omega = 0.5, r = 1,
      lower = -Inf
    ),
    list(
      label = "two-sample one-sided",
      t_stat = 2.5, n = NULL, n1 = 20, n2 = 70,
      one_sample = FALSE, one_sided = TRUE, omega = 0.35, r = 3,
      lower = 0
    )
  )

  for(case in cases){
    tau2 <- if(case$one_sample){
      BFF:::get_one_sample_tau2(n = case$n, w = case$omega, r = case$r)
    }else{
      BFF:::get_two_sample_tau2(n1 = case$n1, n2 = case$n2, w = case$omega, r = case$r)
    }

    prior_integral <- integrate_t_density(
      f = function(delta) BFF:::.t_test.prior(
        tau2 = tau2,
        r = case$r,
        effect_size = delta,
        n = case$n,
        n1 = case$n1,
        n2 = case$n2,
        one_sample = case$one_sample,
        one_sided = case$one_sided
      ),
      lower = case$lower,
      upper = Inf
    )
    posterior_integral <- integrate_t_density(
      f = function(delta) BFF:::.t_test.posterior(
        t_stat = case$t_stat,
        tau2 = tau2,
        r = case$r,
        effect_size = delta,
        n = case$n,
        n1 = case$n1,
        n2 = case$n2,
        one_sample = case$one_sample,
        one_sided = case$one_sided
      ),
      lower = case$lower,
      upper = Inf
    )

    testthat::expect_equal(prior_integral, 1, tolerance = 1e-6, info = case$label)
    testthat::expect_equal(posterior_integral, 1, tolerance = 1e-5, info = case$label)
  }
})

test_that("posterior_plot works for a selected t-test non-local prior", {
  fit <- t_test_BFF(t_stat = 2.5, n = 50, one_sample = TRUE)

  plot_data <- posterior_plot(fit, plot = FALSE, prior = TRUE)

  testthat::expect_true(is.data.frame(plot_data))
  testthat::expect_equal(colnames(plot_data), c("x", "prior", "posterior"))
  testthat::expect_true(all(is.finite(plot_data$prior)))
  testthat::expect_true(all(is.finite(plot_data$posterior)))
  testthat::expect_true(all(plot_data$prior >= 0))
  testthat::expect_true(all(plot_data$posterior >= 0))
  testthat::expect_true(inherits(posterior_plot(fit), "ggplot"))
})

test_that("posterior_plot respects less-than one-sided support", {
  fit <- t_test_BFF(
    t_stat = -2.5,
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

test_that("posterior_plot rejects ambiguous or null selected t-test priors", {
  null_fit <- t_test_BFF(
    t_stat = 0.5,
    alternative = "less",
    r = 3,
    n1 = 50,
    n2 = 50
  )
  multi_fit <- t_test_BFF(
    t_stat = 2.5,
    n = 50,
    one_sample = TRUE,
    omega = c(0.3, 0.5)
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
})
