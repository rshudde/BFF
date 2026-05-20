expect_prior_mode <- function(label, density, expected, interval, tolerance = 1e-5){
  mode <- stats::optimize(
    f = function(x) -density(x),
    interval = interval,
    tol = 1e-12
  )$minimum

  testthat::expect_equal(mode, expected, tolerance = tolerance, info = label)
}

selected_tau2 <- function(fit){
  tau2 <- fit$tau2_h1
  if(is.list(tau2)){
    tau2 <- tau2[[1]]
  }
  unname(tau2)
}

chi2_prior_density <- function(fit, effect_size, table_dim = NULL, branch_sign = NULL){
  tau2 <- selected_tau2(fit)
  function(x) BFF:::.effect_size_density(
    test_type   = fit$test_type,
    effect_size = effect_size,
    x           = x,
    density     = function(internal) BFF:::.chi2_test.prior(
      tau2       = tau2,
      r          = fit$r,
      effect_size = internal,
      n          = fit$input$n,
      df         = fit$input$df
    ),
    input       = fit$input,
    table_dim   = table_dim,
    table_margins = fit$input$table_margins,
    branch_sign = branch_sign
  )
}

f_prior_density <- function(fit, effect_size){
  tau2 <- selected_tau2(fit)
  function(x) BFF:::.effect_size_density(
    test_type   = "f_test",
    effect_size = effect_size,
    x           = x,
    density     = function(internal) BFF:::.f_test.prior(
      tau2       = tau2,
      r          = fit$r,
      effect_size = internal,
      n          = fit$input$n,
      df1        = fit$input$df1
    ),
    input       = fit$input
  )
}

regression_prior_density <- function(fit, effect_size){
  tau2 <- selected_tau2(fit)
  alternative <- fit$input$alternative.original
  one_sided <- fit$input$alternative != "two.sided"
  function(x) BFF:::.effect_size_density(
    test_type   = "regression_test",
    effect_size = effect_size,
    x           = x,
    density     = function(internal) BFF:::.regression_test.prior(
      tau2       = tau2,
      r          = fit$r,
      effect_size = internal,
      n          = fit$input$n,
      k          = fit$input$k,
      one_sided  = one_sided
    ),
    input       = fit$input,
    alternative = alternative
  )
}

test_that("z and t priors peak at requested Cohen's d values", {
  r <- 2
  omega <- 0.35
  tau2_z <- BFF:::get_one_sample_tau2(n = 40, w = omega, r = r)
  z_prior <- function(x) BFF:::.z_test.prior(
    tau2 = tau2_z, r = r, effect_size = x, n = 40,
    one_sample = TRUE, one_sided = FALSE
  )

  expect_prior_mode("z Cohen's d positive", z_prior,  omega, c(0.05, 0.8))
  expect_prior_mode("z Cohen's d negative", z_prior, -omega, c(-0.8, -0.05))

  tau2_t <- BFF:::get_two_sample_tau2(n1 = 35, n2 = 45, w = omega, r = r)
  t_prior <- function(x) BFF:::.t_test.prior(
    tau2 = tau2_t, r = r, effect_size = x, n1 = 35, n2 = 45,
    one_sample = FALSE, one_sided = FALSE
  )

  expect_prior_mode("t Cohen's d positive", t_prior,  omega, c(0.05, 0.8))
  expect_prior_mode("t Cohen's d negative", t_prior, -omega, c(-0.8, -0.05))
})

test_that("chi-square priors peak on every supported effect-size scale", {
  table_dim <- c(3, 4)
  cases <- list(
    list(effect_size = "omega", value = 0.20, interval = c(0.05, 0.50), table_dim = NULL),
    list(effect_size = "cohens_w", value = 0.40, interval = c(0.10, 0.80), table_dim = NULL),
    list(effect_size = "phi", value = 0.40, interval = c(0.10, 0.80), table_dim = NULL),
    list(effect_size = "cramers_v", value = 0.25, interval = c(0.05, 0.60), table_dim = table_dim),
    list(effect_size = "tschuprow_t", value = 0.25, interval = c(0.05, 0.60), table_dim = table_dim),
    list(effect_size = "contingency_coefficient", value = 0.30, interval = c(0.05, 0.80), table_dim = NULL)
  )

  for(case in cases){
    fit <- chi2_test_BFF(
      chi2_stat = 12,
      n = 60,
      df = 6,
      omega = case$value,
      effect_size = case$effect_size,
      table_dim = case$table_dim
    )
    expect_prior_mode(
      label = paste("chi-square", case$effect_size),
      density = chi2_prior_density(fit, case$effect_size, table_dim = case$table_dim),
      expected = case$value,
      interval = case$interval
    )
  }
})

test_that("contingency-table and two-proportion priors peak on their requested scales", {
  table <- matrix(
    c(16, 10, 8, 14, 17, 12, 9, 11, 13, 7, 15, 18),
    nrow = 3,
    byrow = TRUE
  )
  contingency_fit <- contingency_table_BFF(
    table = table,
    omega = 0.25,
    effect_size = "cramers_v"
  )
  expect_prior_mode(
    "contingency Cramer's V",
    chi2_prior_density(contingency_fit, "cramers_v"),
    0.25,
    c(0.05, 0.60)
  )

  cases <- list(
    list(effect_size = "logOR", value = 0.80, interval = c(0.05, 1.50)),
    list(effect_size = "logOR", value = -0.80, interval = c(-1.50, -0.05)),
    list(effect_size = "OR", value = 2.00, interval = c(1.05, 4.00)),
    list(effect_size = "OR", value = 0.50, interval = c(0.15, 0.95)),
    list(effect_size = "logRR", value = 0.40, interval = c(0.05, 1.20)),
    list(effect_size = "risk_ratio", value = 1.50, interval = c(1.05, 3.00)),
    list(effect_size = "risk_difference", value = 0.20, interval = c(0.02, 0.80)),
    list(effect_size = "arcsine_h", value = 0.40, interval = c(0.05, 1.20))
  )

  for(case in cases){
    fit <- prop_test_BFF(
      x = c(30, 20),
      n = c(100, 100),
      omega = case$value,
      effect_size = case$effect_size
    )
    expect_prior_mode(
      label = paste("two-proportion", case$effect_size),
      density = chi2_prior_density(
        fit,
        case$effect_size,
        branch_sign = fit$effect_size_sign_h1
      ),
      expected = case$value,
      interval = case$interval
    )
  }
})

test_that("F-test priors peak on every supported effect-size scale", {
  cases <- list(
    list(effect_size = "omega", value = 0.30, interval = c(0.05, 0.70)),
    list(effect_size = "cohens_f", value = 0.45, interval = c(0.10, 1.00)),
    list(effect_size = "cohens_f2", value = 0.16, interval = c(0.03, 0.50)),
    list(effect_size = "partial_eta2", value = 0.30, interval = c(0.05, 0.80)),
    list(effect_size = "partial_r2", value = 0.30, interval = c(0.05, 0.80))
  )

  for(case in cases){
    fit <- f_test_BFF(
      f_stat = 1.75,
      n = 80,
      df1 = 5,
      df2 = 74,
      omega = case$value,
      effect_size = case$effect_size,
      r = 2
    )
    expect_prior_mode(
      label = paste("F", case$effect_size),
      density = f_prior_density(fit, case$effect_size),
      expected = case$value,
      interval = case$interval
    )
  }
})

test_that("regression priors peak on every supported effect-size scale", {
  fit <- regression_test_BFF(
    t_stat = 2.5, n = 50, k = 3, omega = 0.35, effect_size = "cohens_f", r = 2
  )
  cohen_f_prior <- regression_prior_density(fit, "cohens_f")
  expect_prior_mode("regression Cohen's f positive", cohen_f_prior, 0.35, c(0.05, 0.80))
  expect_prior_mode("regression Cohen's f negative", cohen_f_prior, -0.35, c(-0.80, -0.05))

  fit <- regression_test_BFF(
    t_stat = 2.5, n = 50, k = 3, omega = 0.30, effect_size = "partial_r", r = 2
  )
  partial_r_prior <- regression_prior_density(fit, "partial_r")
  expect_prior_mode("regression partial r positive", partial_r_prior, 0.30, c(0.05, 0.70))
  expect_prior_mode("regression partial r negative", partial_r_prior, -0.30, c(-0.70, -0.05))

  fit <- regression_test_BFF(
    t_stat = 2.5, n = 50, k = 3, omega = 0.09, effect_size = "partial_r2", r = 2
  )
  expect_prior_mode(
    "regression partial R2",
    regression_prior_density(fit, "partial_r2"),
    0.09,
    c(0.01, 0.35)
  )

  fit <- regression_test_BFF(
    t_stat = 2.5, n = 50, k = 3, omega = 0.16, effect_size = "cohens_f2", r = 2
  )
  expect_prior_mode(
    "regression Cohen's f2",
    regression_prior_density(fit, "cohens_f2"),
    0.16,
    c(0.03, 0.50)
  )

  fit <- regression_test_BFF(
    t_stat = 2.5, n = 50, k = 3, alternative = "less",
    omega = 0.30, effect_size = "partial_r", r = 2
  )
  expect_prior_mode(
    "regression one-sided partial r",
    regression_prior_density(fit, "partial_r"),
    -0.30,
    c(-0.70, -0.05)
  )
})
