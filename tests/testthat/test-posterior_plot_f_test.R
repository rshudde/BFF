test_that("F-test r=1 shortcut matches the closed-form implementation", {
  grid <- expand.grid(
    tau2 = c(0.1, 0.5, 2),
    f_stat = c(1.5, 4.05),
    df1 = c(1, 2, 5),
    df2 = c(20, 80)
  )

  for(i in seq_len(nrow(grid))){
    testthat::expect_equal(
      BFF:::BFF_f_test(
        tau2 = grid$tau2[i],
        f_stat = grid$f_stat[i],
        k = grid$df1[i],
        m = grid$df2[i],
        r = 1
      ),
      BFF:::f_val_r1(
        tau2 = grid$tau2[i],
        f_stat = grid$f_stat[i],
        df1 = grid$df1[i],
        df2 = grid$df2[i]
      ),
      tolerance = 1e-10
    )
  }
})

test_that("F-test prior and posterior densities integrate to one", {
  cases <- list(
    list(
      label = "r1 small df1",
      f_stat = 4.05, n = 50, df1 = 2, df2 = 80,
      omega = 0.25, r = 1
    ),
    list(
      label = "r3 larger df1",
      f_stat = 1.75, n = 25, df1 = 5, df2 = 50,
      omega = 0.5, r = 3
    ),
    list(
      label = "unbalanced df",
      f_stat = 2.4, n = 120, df1 = 10, df2 = 90,
      omega = 0.18, r = 2
    )
  )

  for(case in cases){
    tau2 <- BFF:::get_linear_tau2(n = case$n, w = case$omega, k = case$df1, r = case$r)

    prior_integral <- integrate_density(
      f = function(effect_size) BFF:::.f_test.prior(
        tau2 = tau2,
        r = case$r,
        effect_size = effect_size,
        n = case$n,
        df1 = case$df1
      ),
      lower = 0
    )
    posterior_integral <- integrate_density(
      f = function(effect_size) BFF:::.f_test.posterior(
        f_stat = case$f_stat,
        tau2 = tau2,
        r = case$r,
        effect_size = effect_size,
        n = case$n,
        df1 = case$df1,
        df2 = case$df2
      ),
      lower = 0
    )

    testthat::expect_equal(prior_integral, 1, tolerance = 1e-6, info = case$label)
    testthat::expect_equal(posterior_integral, 1, tolerance = 1e-5, info = case$label)
  }
})

test_that("posterior_plot works for a selected F-test non-local prior", {
  fixed_fit <- f_test_BFF(
    f_stat = 1.75,
    n = 25,
    df1 = 50,
    df2 = 50,
    omega = 0.5
  )
  bff_fit <- f_test_BFF(
    f_stat = 1.5,
    n = 50,
    df1 = 25,
    df2 = 75
  )

  fixed_data <- posterior_plot(fixed_fit, plot = FALSE, prior = TRUE)
  bff_data <- posterior_plot(bff_fit, plot = FALSE, prior = TRUE)

  expect_posterior_plot_data(fixed_data)
  expect_posterior_plot_data(bff_data)
  testthat::expect_equal(range(fixed_data$x), c(0, 3))
  testthat::expect_true(inherits(posterior_plot(fixed_fit), "ggplot"))
})

test_that("posterior_plot rejects ambiguous or null selected F-test priors", {
  null_fit <- f_test_BFF(
    f_stat = 0.5,
    n = 50,
    df1 = 25,
    df2 = 75
  )
  multi_fit <- f_test_BFF(
    f_stat = 1.75,
    n = 25,
    df1 = 50,
    df2 = 50,
    omega = c(0.2, 0.5)
  )
  vector_fit <- f_test_BFF(
    f_stat = c(1.75, 1.5),
    n = c(25, 50),
    df1 = c(50, 25),
    df2 = c(50, 75),
    omega = 0.5
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
    "single F statistic"
  )
})
