expect_transformed_density_integral <- function(label, lower, upper, f, tolerance = 1e-5){
  testthat::expect_equal(
    integrate_density(f = f, lower = lower, upper = upper),
    1,
    tolerance = tolerance,
    info = label
  )
}

test_that("minimum BFF helper respects the supplied effect-size cutoff", {
  testthat::expect_equal(
    BFF:::get_min_omega_bff(
      omega = c(0.05, 0.08, 0.20),
      bff = c(3, 1, 2),
      cutoff = 0.10
    ),
    c(2, 0.20)
  )
  testthat::expect_equal(
    BFF:::get_min_omega_bff(
      omega = c(0.05, 0.08, 0.20),
      bff = c(3, 1, 2),
      cutoff = 0.07
    ),
    c(1, 0.08)
  )
})

test_that("chi-square effect-size modes are converted to the internal omega scale", {
  internal_omega <- 0.2
  df <- 6
  table_dim <- c(3, 4)
  base <- chi2_test_BFF(chi2_stat = 12, n = 60, df = df, omega = internal_omega)

  cases <- list(
    list(effect_size = "cohens_w", value = sqrt(df) * internal_omega, table_dim = NULL),
    list(effect_size = "phi", value = sqrt(df) * internal_omega, table_dim = NULL),
    list(effect_size = "cramers_v", value = sqrt(df) * internal_omega / sqrt(min(table_dim - 1)), table_dim = table_dim),
    list(effect_size = "tschuprow_t", value = sqrt(df) * internal_omega / ((table_dim[1] - 1) * (table_dim[2] - 1))^(1/4), table_dim = table_dim),
    list(effect_size = "contingency_coefficient", value = sqrt(df) * internal_omega / sqrt(1 + df * internal_omega^2), table_dim = NULL)
  )

  for(case in cases){
    fit <- chi2_test_BFF(
      chi2_stat = 12,
      n = 60,
      df = df,
      omega = case$value,
      effect_size = case$effect_size,
      table_dim = case$table_dim
    )

    testthat::expect_equal(fit$omega_h1, internal_omega, tolerance = 1e-12, info = case$effect_size)
    testthat::expect_equal(fit$log_bf_h1, base$log_bf_h1, tolerance = 1e-12, info = case$effect_size)
  }
})

test_that("2x2 table effect-size modes are converted to the internal omega scale", {
  x <- c(30, 20)
  n <- c(100, 100)
  internal_omega <- 0.2
  base <- prop_test_BFF(x = x, n = n, omega = internal_omega, effect_size = "cohens_w")

  margins <- c(row1 = 0.5, col1 = 0.25)
  log_or_scale <- sqrt(prod(margins * (1 - margins)))
  log_rr_scale <- sqrt(margins[1] * (1 - margins[1]) * margins[2] / (1 - margins[2]))
  rd_scale <- sqrt(margins[1] * (1 - margins[1]) / (margins[2] * (1 - margins[2])))
  h_scale <- sqrt(margins[1] * (1 - margins[1]))

  cases <- list(
    list(effect_size = "logOR", value = internal_omega / log_or_scale),
    list(effect_size = "logOR", value = -internal_omega / log_or_scale),
    list(effect_size = "OR", value = exp(internal_omega / log_or_scale)),
    list(effect_size = "OR", value = exp(-internal_omega / log_or_scale)),
    list(effect_size = "logRR", value = internal_omega / log_rr_scale),
    list(effect_size = "risk_ratio", value = exp(internal_omega / log_rr_scale)),
    list(effect_size = "risk_difference", value = internal_omega / rd_scale),
    list(effect_size = "arcsine_h", value = internal_omega / h_scale)
  )

  for(case in cases){
    fit <- prop_test_BFF(
      x = x,
      n = n,
      omega = case$value,
      effect_size = case$effect_size
    )

    testthat::expect_equal(fit$omega_h1, internal_omega, tolerance = 1e-12, info = case$effect_size)
    testthat::expect_equal(fit$log_bf_h1, base$log_bf_h1, tolerance = 1e-12, info = case$effect_size)
  }
})

test_that("bounded 2x2 default effect-size sequences avoid excluded endpoints", {
  rd_sequence <- BFF:::.effect_size_default_sequence("prop_test", "risk_difference")
  h_sequence <- BFF:::.effect_size_default_sequence("prop_test", "arcsine_h")

  testthat::expect_lt(max(rd_sequence), 1)
  testthat::expect_lt(max(h_sequence), pi)

  testthat::expect_no_error(
    prop_test_BFF(x = c(30, 20), n = c(100, 100), effect_size = "risk_difference")
  )
  testthat::expect_no_error(
    prop_test_BFF(x = c(30, 20), n = c(100, 100), effect_size = "arcsine_h")
  )

  fit <- prop_test_BFF(x = c(30, 20), n = c(100, 100), omega = 0.4)
  plot_data <- posterior_plot(fit, prior = TRUE, plot = FALSE, effect_size = "arcsine_h")
  testthat::expect_true(all(abs(plot_data$x) < pi))
  expect_posterior_plot_data(plot_data)
})

test_that("2x2 bounded local density transforms are normalized on their legal scales", {
  n <- 100
  df <- 1
  r <- 1
  omega <- 0.8
  tau2 <- BFF:::get_count_tau2(n = n, k = df, w = omega, r = r)
  input <- list(
    n = n,
    df = df,
    table_dim = c(2, 2),
    table_margins = c(row1 = 0.4, col1 = 0.3)
  )
  cases <- list(
    list(effect_size = "risk_difference", lower = -1, upper = 1),
    list(effect_size = "arcsine_h", lower = -pi, upper = pi)
  )

  for(case in cases){
    area <- integrate_density(
      lower = case$lower,
      upper = case$upper,
      f = function(x) BFF:::.effect_size_density(
        test_type = "prop_test",
        effect_size = case$effect_size,
        x = x,
        density = function(effect_size) BFF:::.chi2_test.prior(tau2 = tau2, r = r, effect_size = effect_size, n = n, df = df),
        input = input
      )
    )

    testthat::expect_equal(area, 1, tolerance = 1e-5, info = case$effect_size)
  }
})

test_that("F-test and regression effect-size modes are converted to the internal omega scale", {
  f_base <- f_test_BFF(f_stat = 1.75, n = 25, df1 = 5, df2 = 50, omega = 0.5)
  conventional_f <- sqrt(5 / 2) * 0.5
  f_cases <- list(
    list(effect_size = "cohens_f", value = conventional_f, internal = 0.5),
    list(effect_size = "cohens_f2", value = conventional_f^2, internal = 0.5),
    list(effect_size = "partial_eta2", value = conventional_f^2 / (1 + conventional_f^2), internal = 0.5),
    list(effect_size = "partial_r2", value = conventional_f^2 / (1 + conventional_f^2), internal = 0.5)
  )

  for(case in f_cases){
    fit <- f_test_BFF(
      f_stat = 1.75,
      n = 25,
      df1 = 5,
      df2 = 50,
      omega = case$value,
      effect_size = case$effect_size
    )
    base <- f_test_BFF(f_stat = 1.75, n = 25, df1 = 5, df2 = 50, omega = case$internal)

    testthat::expect_equal(fit$omega_h1, case$internal, tolerance = 1e-12, info = case$effect_size)
    testthat::expect_equal(fit$log_bf_h1, base$log_bf_h1, tolerance = 1e-12, info = case$effect_size)
  }

  testthat::expect_equal(f_base$omega_h1, 0.5)

  reg_cases <- list(
    list(effect_size = "partial_r", value = -0.3, internal = 0.3 / sqrt(1 - 0.3^2)),
    list(effect_size = "partial_r2", value = 0.09, internal = sqrt(0.09 / 0.91)),
    list(effect_size = "cohens_f2", value = 0.16, internal = 0.4)
  )

  for(case in reg_cases){
    fit <- regression_test_BFF(
      t_stat = 2.5,
      n = 50,
      k = 3,
      omega = case$value,
      effect_size = case$effect_size
    )
    base <- regression_test_BFF(t_stat = 2.5, n = 50, k = 3, omega = case$internal)

    testthat::expect_equal(fit$omega_h1, case$internal, tolerance = 1e-12, info = case$effect_size)
    testthat::expect_equal(fit$log_bf_h1, base$log_bf_h1, tolerance = 1e-12, info = case$effect_size)
  }
})

test_that("plot.BFF transforms the Bayes factor function axis", {
  chi_fit <- chi2_test_BFF(
    chi2_stat = 12,
    n = 60,
    df = 6,
    omega_sequence = c(0.1, 0.2),
    effect_size = "cohens_w"
  )
  chi_data <- plot(chi_fit, plot = FALSE)
  testthat::expect_equal(chi_data$x, c(0, 0.1, 0.2), tolerance = 1e-12)

  log_or_fit <- prop_test_BFF(
    x = c(30, 20),
    n = c(100, 100),
    omega_sequence = c(0.4, 0.8),
    effect_size = "logOR"
  )
  log_or_data <- plot(log_or_fit, plot = FALSE)
  testthat::expect_equal(log_or_data$x, c(0, 0.4, 0.8), tolerance = 1e-12)

  negative_log_or_fit <- prop_test_BFF(
    x = c(30, 20),
    n = c(100, 100),
    omega_sequence = c(-0.4, -0.8),
    effect_size = "logOR"
  )
  negative_log_or_data <- plot(negative_log_or_fit, plot = FALSE)
  testthat::expect_equal(negative_log_or_data$x, c(0, -0.4, -0.8), tolerance = 1e-12)

  f_fit <- f_test_BFF(
    f_stat = 1.75,
    n = 25,
    df1 = 5,
    df2 = 50,
    omega_sequence = c(0.04, 0.16),
    effect_size = "cohens_f2"
  )
  f_data <- plot(f_fit, plot = FALSE)
  testthat::expect_equal(f_data$x, c(0, 0.04, 0.16), tolerance = 1e-12)

  reg_fit <- regression_test_BFF(
    t_stat = 2.5,
    n = 50,
    k = 3,
    omega_sequence = c(0.04, 0.16),
    effect_size = "partial_r2"
  )
  reg_data <- plot(reg_fit, plot = FALSE)
  testthat::expect_equal(reg_data$x, c(0, 0.04, 0.16), tolerance = 1e-12)

  signed_reg_fit <- regression_test_BFF(
    t_stat = 2.5,
    n = 50,
    k = 3,
    omega_sequence = c(-0.2, -0.4),
    effect_size = "cohens_f"
  )
  signed_reg_data <- plot(signed_reg_fit, plot = FALSE)
  testthat::expect_equal(signed_reg_data$x, c(0, -0.2, -0.4), tolerance = 1e-12)

  less_reg_fit <- regression_test_BFF(
    t_stat = 2.5,
    n = 50,
    k = 3,
    alternative = "less",
    omega_sequence = c(0.2, 0.4),
    effect_size = "partial_r"
  )
  less_reg_data <- plot(less_reg_fit, plot = FALSE)
  testthat::expect_equal(less_reg_data$x, c(0, -0.2, -0.4), tolerance = 1e-12)
})

test_that("print.BFF preserves signed transformed prior modes", {
  log_or_fit <- prop_test_BFF(
    x = c(30, 20),
    n = c(100, 100),
    omega = -0.4,
    effect_size = "logOR"
  )
  testthat::expect_match(
    paste(testthat::capture_output_lines(log_or_fit, print = TRUE, width = 100), collapse = "\n"),
    "prior mode = -0.40 \\(log odds ratio\\)"
  )

  reg_fit <- regression_test_BFF(
    t_stat = 2.5,
    n = 50,
    k = 3,
    omega = -0.3,
    effect_size = "cohens_f"
  )
  testthat::expect_match(
    paste(testthat::capture_output_lines(reg_fit, print = TRUE, width = 100), collapse = "\n"),
    "prior mode = -0.30 \\(signed Cohen's f\\)"
  )
})

test_that("chi-square transformed prior and posterior densities integrate to one", {
  n <- 60
  df <- 6
  chi2_stat <- 12
  r <- 1
  omega <- 0.2
  tau2 <- BFF:::get_count_tau2(n = n, k = df, w = omega, r = r)
  input <- list(n = n, df = df)
  table_dim <- c(3, 4)
  cases <- list(
    list(effect_size = "cohens_w", lower = 0, upper = Inf, table_dim = NULL),
    list(effect_size = "phi", lower = 0, upper = Inf, table_dim = NULL),
    list(effect_size = "cramers_v", lower = 0, upper = Inf, table_dim = table_dim),
    list(effect_size = "tschuprow_t", lower = 0, upper = Inf, table_dim = table_dim),
    list(effect_size = "contingency_coefficient", lower = 0, upper = 1, table_dim = NULL)
  )

  for(case in cases){
    expect_transformed_density_integral(
      label = paste("chi-square prior", case$effect_size),
      lower = case$lower,
      upper = case$upper,
      f = function(x) BFF:::.effect_size_density(
        test_type = "chi2_test",
        effect_size = case$effect_size,
        x = x,
        density = function(effect_size) BFF:::.chi2_test.prior(tau2 = tau2, r = r, effect_size = effect_size, n = n, df = df),
        input = input,
        table_dim = case$table_dim
      )
    )
    expect_transformed_density_integral(
      label = paste("chi-square posterior", case$effect_size),
      lower = case$lower,
      upper = case$upper,
      f = function(x) BFF:::.effect_size_density(
        test_type = "chi2_test",
        effect_size = case$effect_size,
        x = x,
        density = function(effect_size) BFF:::.chi2_test.posterior(chi2_stat = chi2_stat, tau2 = tau2, r = r, effect_size = effect_size, n = n, df = df),
        input = input,
        table_dim = case$table_dim
      )
    )
  }
})

test_that("2x2 table transformed prior and posterior densities integrate to one", {
  n <- 100
  df <- 1
  chi2_stat <- 4.5
  r <- 1
  omega <- 0.2
  tau2 <- BFF:::get_count_tau2(n = n, k = df, w = omega, r = r)
  input <- list(
    n = n,
    df = df,
    table_dim = c(2, 2),
    table_margins = c(row1 = 0.4, col1 = 0.3)
  )
  cases <- list(
    list(effect_size = "logOR", lower = -Inf, upper = Inf),
    list(effect_size = "OR", lower = 0, upper = Inf),
    list(effect_size = "logRR", lower = -Inf, upper = Inf),
    list(effect_size = "risk_ratio", lower = 0, upper = Inf),
    list(effect_size = "risk_difference", lower = -1, upper = 1),
    list(effect_size = "arcsine_h", lower = -pi, upper = pi)
  )

  for(case in cases){
    expect_transformed_density_integral(
      label = paste("2x2 prior", case$effect_size),
      lower = case$lower,
      upper = case$upper,
      f = function(x) BFF:::.effect_size_density(
        test_type = "prop_test",
        effect_size = case$effect_size,
        x = x,
        density = function(effect_size) BFF:::.chi2_test.prior(tau2 = tau2, r = r, effect_size = effect_size, n = n, df = df),
        input = input
      )
    )
    expect_transformed_density_integral(
      label = paste("2x2 posterior", case$effect_size),
      lower = case$lower,
      upper = case$upper,
      f = function(x) BFF:::.effect_size_density(
        test_type = "prop_test",
        effect_size = case$effect_size,
        x = x,
        density = function(effect_size) BFF:::.chi2_test.posterior(chi2_stat = chi2_stat, tau2 = tau2, r = r, effect_size = effect_size, n = n, df = df),
        input = input
      )
    )
  }
})

test_that("F-test transformed prior and posterior densities integrate to one", {
  n <- 25
  df1 <- 5
  df2 <- 50
  f_stat <- 1.75
  r <- 2
  omega <- 0.5
  tau2 <- BFF:::get_linear_tau2(n = n, k = df1, w = omega, r = r)
  input <- list(n = n, df1 = df1, df2 = df2)
  cases <- list(
    list(effect_size = "cohens_f", lower = 0, upper = Inf),
    list(effect_size = "cohens_f2", lower = 0, upper = Inf),
    list(effect_size = "partial_eta2", lower = 0, upper = 1),
    list(effect_size = "partial_r2", lower = 0, upper = 1)
  )

  for(case in cases){
    expect_transformed_density_integral(
      label = paste("F prior", case$effect_size),
      lower = case$lower,
      upper = case$upper,
      f = function(x) BFF:::.effect_size_density(
        test_type = "f_test",
        effect_size = case$effect_size,
        x = x,
        density = function(effect_size) BFF:::.f_test.prior(tau2 = tau2, r = r, effect_size = effect_size, n = n, df1 = df1),
        input = input
      )
    )
    expect_transformed_density_integral(
      label = paste("F posterior", case$effect_size),
      lower = case$lower,
      upper = case$upper,
      f = function(x) BFF:::.effect_size_density(
        test_type = "f_test",
        effect_size = case$effect_size,
        x = x,
        density = function(effect_size) BFF:::.f_test.posterior(f_stat = f_stat, tau2 = tau2, r = r, effect_size = effect_size, n = n, df1 = df1, df2 = df2),
        input = input
      )
    )
  }
})

test_that("regression transformed prior and posterior densities integrate to one", {
  n <- 50
  k <- 3
  t_stat <- 2.5
  r <- 2
  omega <- 0.4
  tau2 <- BFF:::get_regression_tau2(n = n, k = k, w = omega, r = r)
  input <- list(n = n, k = k)
  cases <- list(
    list(effect_size = "partial_r", lower = -1, upper = 1, alternative = "two.sided", one_sided = FALSE),
    list(effect_size = "partial_r2", lower = 0, upper = 1, alternative = "two.sided", one_sided = FALSE),
    list(effect_size = "cohens_f2", lower = 0, upper = Inf, alternative = "two.sided", one_sided = FALSE),
    list(effect_size = "partial_r", lower = -1, upper = 0, alternative = "less", one_sided = TRUE),
    list(effect_size = "partial_r2", lower = 0, upper = 1, alternative = "less", one_sided = TRUE),
    list(effect_size = "cohens_f2", lower = 0, upper = Inf, alternative = "less", one_sided = TRUE)
  )

  for(case in cases){
    expect_transformed_density_integral(
      label = paste("regression prior", case$effect_size, case$alternative),
      lower = case$lower,
      upper = case$upper,
      f = function(x) BFF:::.effect_size_density(
        test_type = "regression_test",
        effect_size = case$effect_size,
        x = x,
        density = function(effect_size) BFF:::.regression_test.prior(tau2 = tau2, r = r, effect_size = effect_size, n = n, k = k, one_sided = case$one_sided),
        input = input,
        alternative = case$alternative
      )
    )
    expect_transformed_density_integral(
      label = paste("regression posterior", case$effect_size, case$alternative),
      lower = case$lower,
      upper = case$upper,
      f = function(x) BFF:::.effect_size_density(
        test_type = "regression_test",
        effect_size = case$effect_size,
        x = x,
        density = function(effect_size) BFF:::.regression_test.posterior(t_stat = t_stat, tau2 = tau2, r = r, effect_size = effect_size, n = n, k = k, one_sided = case$one_sided),
        input = input,
        alternative = case$alternative
      )
    )
  }
})

test_that("posterior_plot returns finite data on transformed effect-size scales", {
  chi_fit <- chi2_test_BFF(chi2_stat = 12, n = 60, df = 6, omega = 0.2)
  expect_posterior_plot_data(posterior_plot(chi_fit, prior = TRUE, plot = FALSE, effect_size = "cohens_w"))
  expect_posterior_plot_data(posterior_plot(chi_fit, prior = TRUE, plot = FALSE, effect_size = "cramers_v", table_dim = c(3, 4)))
  testthat::expect_error(
    posterior_plot(chi_fit, prior = TRUE, plot = FALSE, effect_size = "cramers_v"),
    "`table_dim` must be supplied"
  )

  log_or_fit <- prop_test_BFF(
    x = c(30, 20),
    n = c(100, 100),
    omega = 0.8,
    effect_size = "logOR"
  )
  expect_posterior_plot_data(posterior_plot(log_or_fit, prior = TRUE, plot = FALSE))

  f_fit <- f_test_BFF(f_stat = 1.75, n = 25, df1 = 5, df2 = 50, omega = 0.5)
  expect_posterior_plot_data(posterior_plot(f_fit, prior = TRUE, plot = FALSE, effect_size = "cohens_f"))
  expect_posterior_plot_data(posterior_plot(f_fit, prior = TRUE, plot = FALSE, effect_size = "cohens_f2"))
  expect_posterior_plot_data(posterior_plot(f_fit, prior = TRUE, plot = FALSE, effect_size = "partial_eta2"))

  reg_fit <- regression_test_BFF(t_stat = 2.5, n = 50, k = 3, omega = 0.4)
  expect_posterior_plot_data(posterior_plot(reg_fit, prior = TRUE, plot = FALSE, effect_size = "partial_r"))
  expect_posterior_plot_data(posterior_plot(reg_fit, prior = TRUE, plot = FALSE, effect_size = "partial_r2"))
})
