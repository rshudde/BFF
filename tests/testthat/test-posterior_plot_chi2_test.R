test_that("chi-square BFF r=1 shortcut matches the closed-form implementation", {
  grid <- expand.grid(
    tau2 = c(0.1, 0.5, 2),
    chi2_stat = c(1.5, 6.5),
    df = c(1, 3, 6)
  )

  for(i in seq_len(nrow(grid))){
    testthat::expect_equal(
      BFF:::BFF_chi2_test(
        tau2 = grid$tau2[i],
        chi2_stat = grid$chi2_stat[i],
        k = grid$df[i],
        r = 1
      ),
      BFF:::G_val_r1(
        tau2 = grid$tau2[i],
        chi2_stat = grid$chi2_stat[i],
        df = grid$df[i]
      ),
      tolerance = 1e-10
    )
  }
})

test_that("chi-square prior and posterior densities integrate to one", {
  cases <- list(
    list(
      label = "count r1",
      chi2_stat = 60, n = 45, df = 44, LRT = FALSE,
      omega = 0.2, r = 1
    ),
    list(
      label = "count r5",
      chi2_stat = 30, n = 45, df = 10, LRT = FALSE,
      omega = 0.3, r = 5
    ),
    list(
      label = "LRT r1",
      chi2_stat = 12, n = 60, df = 6, LRT = TRUE,
      omega = 0.15, r = 1
    ),
    list(
      label = "LRT r3",
      chi2_stat = 20, n = 80, df = 8, LRT = TRUE,
      omega = 0.2, r = 3
    )
  )

  for(case in cases){
    tau2 <- if(case$LRT){
      BFF:::get_LRT_tau2(n = case$n, k = case$df, w = case$omega, r = case$r)
    }else{
      BFF:::get_count_tau2(n = case$n, k = case$df, w = case$omega, r = case$r)
    }

    prior_integral <- integrate_density(
      f = function(effect_size) BFF:::.chi2_test.prior(
        tau2 = tau2,
        r = case$r,
        effect_size = effect_size,
        n = case$n,
        df = case$df
      ),
      lower = 0
    )
    posterior_integral <- integrate_density(
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

    testthat::expect_equal(prior_integral, 1, tolerance = 1e-6, info = case$label)
    testthat::expect_equal(posterior_integral, 1, tolerance = 1e-5, info = case$label)
  }
})

test_that("posterior_plot works for a selected chi-square non-local prior", {
  fixed_fit <- chi2_test_BFF(
    chi2_stat = 12.65,
    n = 707,
    df = 6,
    LRT = FALSE,
    omega = 0.05
  )
  bff_fit <- chi2_test_BFF(
    chi2_stat = 35,
    n = 200,
    df = 6,
    LRT = TRUE
  )

  fixed_data <- posterior_plot(fixed_fit, plot = FALSE, prior = TRUE)
  bff_data <- posterior_plot(bff_fit, plot = FALSE, prior = TRUE)

  expect_posterior_plot_data(fixed_data)
  expect_posterior_plot_data(bff_data)
  testthat::expect_equal(range(fixed_data$x), c(0, 3))
  testthat::expect_true(inherits(posterior_plot(fixed_fit), "ggplot"))
})

test_that("posterior_plot rejects ambiguous or null selected chi-square priors", {
  null_fit <- chi2_test_BFF(
    chi2_stat = 7.5,
    n = 45,
    df = 44,
    LRT = TRUE
  )
  multi_fit <- chi2_test_BFF(
    chi2_stat = 12.65,
    n = 707,
    df = 6,
    omega = c(0.03, 0.05)
  )
  vector_fit <- chi2_test_BFF(
    chi2_stat = c(12.65, 10),
    n = c(707, 600),
    df = 6,
    omega = 0.05
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
    "single chi-square statistic"
  )
})
