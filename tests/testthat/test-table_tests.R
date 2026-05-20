test_that("prop_test_BFF matches uncorrected two-sample prop.test statistic", {
  x <- c(8, 22)
  n <- c(40, 40)
  ref_stat <- unname(stats::prop.test(x = x, n = n, correct = FALSE)$statistic)

  fit <- prop_test_BFF(x = x, n = n, omega = 0.2, effect_size = "cohens_w")
  ref <- chi2_test_BFF(chi2_stat = ref_stat, n = sum(n), df = 1, omega = 0.2)

  testthat::expect_equal(fit$input$chi2_stat, ref_stat, tolerance = 1e-12)
  testthat::expect_equal(fit$input$df, 1)
  testthat::expect_equal(fit$input$n, sum(n))
  testthat::expect_equal(fit$log_bf_h1, ref$log_bf_h1, tolerance = 1e-12)
  testthat::expect_equal(fit$omega_h1, ref$omega_h1, tolerance = 1e-12)
  testthat::expect_equal(fit$input$observed_effects$sign, -1)
  testthat::expect_lt(fit$input$observed_effects$logOR, 0)
  testthat::expect_lt(fit$input$observed_effects$risk_difference, 0)
})

test_that("contingency_table_BFF matches uncorrected chisq.test statistic", {
  tab <- matrix(c(10, 20, 30, 6, 14, 25, 18, 10, 7), nrow = 3, byrow = TRUE)
  ref_stat <- unname(stats::chisq.test(tab, correct = FALSE)$statistic)

  fit <- contingency_table_BFF(tab, omega = 0.2)
  ref <- chi2_test_BFF(chi2_stat = ref_stat, n = sum(tab), df = 4, omega = 0.2)

  testthat::expect_equal(fit$input$chi2_stat, ref_stat, tolerance = 1e-12)
  testthat::expect_equal(fit$input$df, 4)
  testthat::expect_equal(fit$input$table_dim, c(3, 3))
  testthat::expect_equal(fit$input$n, sum(tab))
  testthat::expect_equal(fit$log_bf_h1, ref$log_bf_h1, tolerance = 1e-12)
})

test_that("contingency_table_BFF computes likelihood-ratio G-squared statistic", {
  tab <- matrix(c(10, 20, 30, 6, 14, 25, 18, 10, 7), nrow = 3, byrow = TRUE)
  expected <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  positive <- tab > 0
  g2 <- 2 * sum(tab[positive] * log(tab[positive] / expected[positive]))

  fit <- contingency_table_BFF(tab, LRT = TRUE, omega = 0.2)
  ref <- chi2_test_BFF(chi2_stat = g2, n = sum(tab), df = 4, LRT = TRUE, omega = 0.2)

  testthat::expect_equal(fit$input$chi2_stat, g2, tolerance = 1e-12)
  testthat::expect_true(fit$input$LRT)
  testthat::expect_equal(fit$log_bf_h1, ref$log_bf_h1, tolerance = 1e-12)
})

test_that("2x2 contingency table stores margins and supports logOR plots", {
  tab <- matrix(c(8, 32, 22, 18), nrow = 2, byrow = TRUE)
  fit <- contingency_table_BFF(tab, omega = 0.8, effect_size = "logOR")

  testthat::expect_equal(fit$input$table_dim, c(2, 2))
  testthat::expect_equal(
    fit$input$table_margins,
    c(row1 = sum(tab[1, ]) / sum(tab), col1 = sum(tab[, 1]) / sum(tab)),
    tolerance = 1e-12
  )
  testthat::expect_equal(fit$input$observed_effects$sign, -1)
  expect_posterior_plot_data(posterior_plot(fit, prior = TRUE, plot = FALSE))
})

test_that("observed 2x2 effects are finite for zero interior cells", {
  tab <- matrix(c(0, 10, 5, 15), nrow = 2, byrow = TRUE)
  ref_stat <- suppressWarnings(unname(stats::chisq.test(tab, correct = FALSE)$statistic))
  haldane_anscombe <- tab + 0.5
  corrected_p <- haldane_anscombe[, 1] / rowSums(haldane_anscombe)
  corrected_or <- (haldane_anscombe[1, 1] * haldane_anscombe[2, 2]) /
    (haldane_anscombe[1, 2] * haldane_anscombe[2, 1])
  corrected_rr <- corrected_p[1] / corrected_p[2]

  fit <- contingency_table_BFF(tab, omega = 0.2, effect_size = "cohens_w")
  ref <- chi2_test_BFF(chi2_stat = ref_stat, n = sum(tab), df = 1, omega = 0.2)
  effects <- fit$input$observed_effects

  testthat::expect_true(all(rowSums(tab) > 0))
  testthat::expect_true(all(colSums(tab) > 0))
  testthat::expect_true(all(is.finite(unlist(effects, use.names = FALSE))))
  testthat::expect_equal(fit$input$chi2_stat, ref_stat, tolerance = 1e-12)
  testthat::expect_equal(fit$log_bf_h1, ref$log_bf_h1, tolerance = 1e-12)
  testthat::expect_equal(effects$sign, -1)
  testthat::expect_equal(effects$logOR, log(corrected_or), tolerance = 1e-12)
  testthat::expect_equal(effects$OR, corrected_or, tolerance = 1e-12)
  testthat::expect_equal(effects$logRR, log(corrected_rr), tolerance = 1e-12)
  testthat::expect_equal(effects$risk_ratio, corrected_rr, tolerance = 1e-12)
  testthat::expect_equal(effects$risk_difference, 0 / 10 - 5 / 20, tolerance = 1e-12)

  prop_fit <- prop_test_BFF(x = c(0, 5), n = c(10, 20), omega = 0.2, effect_size = "logOR")
  prop_ref_stat <- suppressWarnings(unname(stats::prop.test(
    x = c(0, 5), n = c(10, 20), correct = FALSE
  )$statistic))

  testthat::expect_true(all(is.finite(unlist(prop_fit$input$observed_effects, use.names = FALSE))))
  testthat::expect_equal(prop_fit$input$chi2_stat, prop_ref_stat, tolerance = 1e-12)
  testthat::expect_equal(unname(prop_fit$input$observed_effects$logOR), log(corrected_or), tolerance = 1e-12)
  testthat::expect_equal(unname(prop_fit$input$observed_effects$logRR), log(corrected_rr), tolerance = 1e-12)
})

test_that("2x2 wrapper metadata supports transformed density integration", {
  fit <- prop_test_BFF(x = c(30, 20), n = c(100, 100), omega = 0.5)
  tau2 <- BFF:::get_count_tau2(n = fit$input$n, k = fit$input$df, w = fit$omega_h1, r = fit$r)

  expect_transformed_density_integral(
    label = "prop_test prior logOR",
    lower = -Inf,
    upper = Inf,
    f = function(x) BFF:::.effect_size_density(
      test_type = fit$test_type,
      effect_size = "logOR",
      x = x,
      density = function(effect_size) BFF:::.chi2_test.prior(
        tau2 = tau2, r = fit$r, effect_size = effect_size, n = fit$input$n, df = fit$input$df
      ),
      input = fit$input
    )
  )
  expect_transformed_density_integral(
    label = "prop_test posterior logOR",
    lower = -Inf,
    upper = Inf,
    f = function(x) BFF:::.effect_size_density(
      test_type = fit$test_type,
      effect_size = "logOR",
      x = x,
      density = function(effect_size) BFF:::.chi2_test.posterior(
        chi2_stat = fit$input$chi2_stat, tau2 = tau2, r = fit$r,
        effect_size = effect_size, n = fit$input$n, df = fit$input$df
      ),
      input = fit$input
    )
  )
})

test_that("contingency-table wrapper metadata supports Cramer's V integration", {
  tab <- matrix(c(10, 20, 30, 5, 15, 25, 8, 12, 18, 20, 10, 7), nrow = 3, byrow = TRUE)
  fit <- contingency_table_BFF(tab, omega = 0.2)
  tau2 <- BFF:::get_count_tau2(n = fit$input$n, k = fit$input$df, w = fit$omega_h1, r = fit$r)

  expect_posterior_plot_data(posterior_plot(fit, prior = TRUE, plot = FALSE, effect_size = "cramers_v"))

  expect_transformed_density_integral(
    label = "contingency_table prior Cramer's V",
    lower = 0,
    upper = Inf,
    f = function(x) BFF:::.effect_size_density(
      test_type = fit$test_type,
      effect_size = "cramers_v",
      x = x,
      density = function(effect_size) BFF:::.chi2_test.prior(
        tau2 = tau2, r = fit$r, effect_size = effect_size, n = fit$input$n, df = fit$input$df
      ),
      input = fit$input
    )
  )
  expect_transformed_density_integral(
    label = "contingency_table posterior Cramer's V",
    lower = 0,
    upper = Inf,
    f = function(x) BFF:::.effect_size_density(
      test_type = fit$test_type,
      effect_size = "cramers_v",
      x = x,
      density = function(effect_size) BFF:::.chi2_test.posterior(
        chi2_stat = fit$input$chi2_stat, tau2 = tau2, r = fit$r,
        effect_size = effect_size, n = fit$input$n, df = fit$input$df
      ),
      input = fit$input
    )
  )
})

test_that("table-test wrappers validate invalid inputs", {
  testthat::expect_error(prop_test_BFF(x = 5, n = 4), "`x` and `n` must be numeric vectors of length two")
  testthat::expect_error(prop_test_BFF(x = c(1, 2), n = c(10, 20, 30)), "`x` and `n` must be numeric vectors of length two")
  testthat::expect_error(prop_test_BFF(x = c(5, 20), n = c(4, 20)), "`x` must contain success counts")
  testthat::expect_error(prop_test_BFF(x = c(1.5, 2), n = c(10, 20)), "`x` and `n` must contain integer counts")

  testthat::expect_error(contingency_table_BFF(matrix(c(1, 2, 3), nrow = 1)), "`table` must have at least two rows and two columns")
  testthat::expect_error(contingency_table_BFF(matrix(c(1, -1, 2, 3), nrow = 2)), "`table` counts must be nonnegative")
  testthat::expect_error(contingency_table_BFF(matrix(c(1.2, 2, 3, 4), nrow = 2)), "`table` counts must be integers")
  testthat::expect_error(contingency_table_BFF(matrix(c(1, 0, 2, 0), nrow = 2)), "`table` must not contain empty rows or columns")
})
