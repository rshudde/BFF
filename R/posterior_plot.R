#' @title Plot Prior and Posterior Distribution
#'
#' @description Creates a prior and posterior plot
#' of an BFF object. If no specific omega was set
#' when fitting the model, the omega resulting in
#' maximum BF against the null hypothesis is selected.
#'
#' @param x a BFF object
#' @param prior whether prior distribution should
#' be added to the figure
#' @param plot whether plot should be generated.
#' Default to \code{TRUE}. If \code{FALSE} a data
#' frame with the support, prior ordinates, and
#' posterior ordinates is returned instead.
#' @param ... additional arguments to the plotting
#' function. These include: \describe{
#'  \item{x_limit}{vector defining the plotting range,
#'  defaults to \code{c(-3, 3)}.}
#'  \item{color}{vector with color for the posterior and
#'  prior line. Defaults to \code{c("black", "grey")}}
#'  \item{linetype}{vector with linetype for the posterior and
#'  prior line. Defaults to \code{c(2, 1)}}
#'  \item{linewidth}{vector with linewidth for the posterior and
#'  prior line. Defaults to \code{c(1, 1)}}
#'  \item{effect_size}{effect-size scale for plotting the prior and
#'  posterior density. Defaults to the scale used to fit the BFF object,
#'  or to the package's internal omega scale.}
#'  \item{table_dim}{integer vector \code{c(rows, columns)} for
#'  chi-square transformations that require table dimensions.}
#'  \item{table_margins}{two marginal probabilities for 2x2
#'  table effect-size transformations.}
#' }
#'
#' @return either a ggplot2 object if \code{plot = TRUE} or a data.frame
#' with prior and posterior densities if \code{plot = FALSE}
#'
#' @seealso [z_test_BFF()], [t_test_BFF()], [chi2_test_BFF()], [f_test_BFF()]
#'
#' @export
posterior_plot <- function(x, prior = FALSE, plot = TRUE, ...){

  # this is the generic function for prior and posterior plot for all tests
  # it requires a data.frame with prior and posterior distributions from the
  # different tests that are provided by the .posterior_plot_data. + test_type functions

  if(!is.BFF(x))
    stop("`posterior_plot` is implemented only for `BFF` objects.")


  # catch and set additional arguments
  dots <- list(...)
  color     <- if(is.null(dots[["color"]]))     c("black", "grey") else dots[["color"]]
  linetype  <- if(is.null(dots[["linetype"]]))  c(1, 2)            else dots[["linetype"]]
  linewidth <- if(is.null(dots[["linewidth"]])) c(1, 1)            else dots[["linewidth"]]
  alternative <- if(x[["test_type"]] %in% c("z_test", "t_test", "regression_test")) .posterior_plot_alternative(x) else NULL
  effect_size <- .effect_size_for_object(x, dots[["effect_size"]])
  table_dim <- if(is.null(dots[["table_dim"]])) x$input$table_dim else dots[["table_dim"]]
  table_margins <- if(is.null(dots[["table_margins"]])) x$input$table_margins else dots[["table_margins"]]
  x_limit   <- if(is.null(dots[["x_limit"]])) .posterior_plot_default_x_limit(x, alternative, effect_size) else dots[["x_limit"]]
  # TODO: deal with positive/negative only plotting: maybe just implement everything for greater
  # (as in the t_test_BFF function) and then flip the support around y-axis


  # manual dispatching for different tests (could have been done with classes)
  if(x[["test_type"]] == "z_test"){
    plot_data <- .posterior_plot_data.z_test(x, prior, x_limit)
  }else if(x[["test_type"]] == "t_test"){
    plot_data <- .posterior_plot_data.t_test(x, prior, x_limit)
  }else if(x[["test_type"]] %in% c("chi2_test", "contingency_table", "prop_test")){
    plot_data <- .posterior_plot_data.chi2_test(x, prior, x_limit, effect_size, table_dim, table_margins)
  }else if(x[["test_type"]] == "f_test"){
    plot_data <- .posterior_plot_data.f_test(x, prior, x_limit, effect_size)
  }else if(x[["test_type"]] == "regression_test"){
    plot_data <- .posterior_plot_data.regression_test(x, prior, x_limit, effect_size)
  }else{
    stop("`posterior_plot` is not implemented for this BFF test type.")
  }


  # return data if no plot requested
  if(!plot){
    return(plot_data)
  }

  # create plot otherwise
  out <- ggplot2::ggplot(
    data      = plot_data,
    mapping   = ggplot2::aes(x = .data[["x"]])
  ) + ggplot2::geom_line(
    mapping   = ggplot2::aes(y = .data[["posterior"]]),
    color     = color[1],
    linetype  = linetype[1],
    linewidth = linewidth[1]
  ) + ggplot2::labs(
    x = if(x$generic_test) expression(tau^2) else .effect_size_plot_label(x, dots[["effect_size"]]),
    y = "Density")

  if(prior){
    out <- out + ggplot2::geom_line(
      mapping   = ggplot2::aes(y = .data[["prior"]]),
      color     = color[2],
      linetype  = linetype[2],
      linewidth = linewidth[2]
    )
  }

  return(out)
}

# helpers shared by test-specific plotting functions
.posterior_plot_alternative <- function(x){
  if(!is.null(x$input$alternative.original)){
    return(x$input$alternative.original)
  }
  if(!is.null(x$input$alternative)){
    return(x$input$alternative)
  }
  if(!is.null(x$alternative)){
    return(x$alternative)
  }
  stop("The BFF object does not contain an alternative specification.")
}

.posterior_plot_selected_omega <- function(x){
  omega <- x$omega_h1

  if(is.null(omega)){
    omega <- x$omega
  }

  if(is.null(omega) || length(omega) != 1 || is.na(omega)){
    stop("`posterior_plot` requires a single selected omega. Fit a single omega or a full BFF so the maximizing omega can be selected.")
  }

  return(omega)
}

.posterior_plot_selected_tau2 <- function(x){
  tau2 <- x$tau2_h1

  if(is.null(tau2)){
    return(NULL)
  }
  if(is.list(tau2)){
    if(length(tau2) != 1){
      stop("`posterior_plot` requires a single selected omega. Fit a single omega or a full BFF so the maximizing omega can be selected.")
    }
    tau2 <- tau2[[1]]
  }
  if(length(tau2) != 1 || is.na(tau2)){
    stop("`posterior_plot` requires a single selected omega. Fit a single omega or a full BFF so the maximizing omega can be selected.")
  }

  tau2
}

.posterior_plot_branch_sign <- function(x, effect_size){
  if(!.effect_size_chi2_family(x$test_type) ||
     !effect_size %in% c("logOR", "OR", "logRR", "risk_ratio", "risk_difference", "arcsine_h")){
    return(NULL)
  }
  if(is.null(x$input$effect_size) || .effect_size_normalize(x$test_type, x$input$effect_size) != effect_size){
    return(NULL)
  }
  if(is.null(x$effect_size_sign_h1) || length(x$effect_size_sign_h1) != 1){
    return(NULL)
  }

  x$effect_size_sign_h1
}

.posterior_plot_default_x_limit <- function(x, alternative, effect_size = NULL){
  transformed_x_limit <- .effect_size_default_x_limit(x[["test_type"]], effect_size, alternative)
  if(!is.null(transformed_x_limit)){
    return(transformed_x_limit)
  }

  switch(
    x[["test_type"]],
    "chi2_test" = c(0, 3),
    "f_test" = c(0, 3),
    switch(
      alternative,
      "two.sided" = c(-3, 3),
      "greater"   = c(0, 3),
      "less"      = c(-3, 0)
    )
  )
}

# test specific plotting functions
.posterior_plot_data.z_test <- function(x, prior, x_limit){

  x_seq <- seq(x_limit[1], x_limit[2], length.out = 301)

  z_stat <- x$input$z_stat
  r      <- x$r
  omega  <- .posterior_plot_selected_omega(x)

  if(length(z_stat) != 1)
    stop("`posterior_plot` for z-test BFF objects is currently implemented only for a single z statistic.")

  if(x$input$one_sample){
    if(length(x$input$n) != 1)
      stop("`posterior_plot` for z-test BFF objects is currently implemented only for a single z statistic.")
    tau2 <- get_one_sample_tau2(n = x$input$n, w = omega, r = r)
  }else{
    if(length(x$input$n1) != 1 || length(x$input$n2) != 1)
      stop("`posterior_plot` for z-test BFF objects is currently implemented only for a single z statistic.")
    tau2 <- get_two_sample_tau2(n1 = x$input$n1, n2 = x$input$n2, w = omega, r = r)
  }

  if(tau2 <= 0)
    stop("There is no non-local prior distribution that provides more evidence for the null hypothesis than the null prior distribution.")

  alternative <- .posterior_plot_alternative(x)
  one_sided  <- x$input$alternative != "two.sided"
  effect_size <- if(alternative == "less") -x_seq else x_seq

  lik.prior     <- .z_test.prior(tau2 = tau2, r = r, effect_size = effect_size, n = x$input$n, n1 = x$input$n1, n2 = x$input$n2, one_sample = x$input$one_sample, one_sided = one_sided)
  lik.posterior <- .z_test.posterior(z_stat = z_stat, tau2 = tau2, r = r, effect_size = effect_size, n = x$input$n, n1 = x$input$n1, n2 = x$input$n2, one_sample = x$input$one_sample, one_sided = one_sided)

  df <- data.frame(
    x         = x_seq,
    prior     = lik.prior,
    posterior = lik.posterior
  )

  return(df)
}

.posterior_plot_data.t_test <- function(x, prior, x_limit){

  x_seq <- seq(x_limit[1], x_limit[2], length.out = 301)

  # extract fitting information from the model object
  t_stat <- x$input$t_stat
  r      <- x$r
  omega  <- .posterior_plot_selected_omega(x)

  if(length(t_stat) != 1)
    stop("`posterior_plot` for t-test BFF objects is currently implemented only for a single t statistic.")

  if(x$input$one_sample){
    if(length(x$input$n) != 1)
      stop("`posterior_plot` for t-test BFF objects is currently implemented only for a single t statistic.")
    tau2 <- get_one_sample_tau2(n = x$input$n, w = omega, r = r)
  }else{
    if(length(x$input$n1) != 1 || length(x$input$n2) != 1)
      stop("`posterior_plot` for t-test BFF objects is currently implemented only for a single t statistic.")
    tau2 <- get_two_sample_tau2(n1 = x$input$n1, n2 = x$input$n2, w = omega, r = r)
  }

  # terminate if tau2 is equal to 0 -> BFF leads to omega 0 (e.g., opposite direction in one-sided test)
  if(tau2 <= 0)
    stop("There is no non-local prior distribution that provides more evidence for the null hypothesis than the null prior distribution.")

  alternative <- .posterior_plot_alternative(x)
  one_sided  <- x$input$alternative != "two.sided"
  effect_size <- if(alternative == "less") -x_seq else x_seq

  # compute prior and posterior
  lik.prior     <- .t_test.prior(tau2 = tau2, r = r, effect_size = effect_size, n = x$input$n, n1 = x$input$n1, n2 = x$input$n2, one_sample = x$input$one_sample, one_sided = one_sided)
  lik.posterior <- .t_test.posterior(t_stat = t_stat, tau2 = tau2, r = r, effect_size = effect_size, n = x$input$n, n1 = x$input$n1, n2 = x$input$n2, one_sample = x$input$one_sample, one_sided = one_sided)

  # create data.frame with values
  posterior = NULL
  df <- data.frame(
    x         = x_seq,
    prior     = lik.prior,
    posterior = lik.posterior
  )

  return(df)
}

.posterior_plot_data.chi2_test <- function(x, prior, x_limit, effect_size, table_dim = NULL, table_margins = NULL){

  x_seq <- seq(x_limit[1], x_limit[2], length.out = 301)

  chi2_stat <- x$input$chi2_stat
  r         <- x$r
  omega     <- .posterior_plot_selected_omega(x)

  if(length(chi2_stat) != 1 || length(x$input$n) != 1 || length(x$input$df) != 1)
    stop("`posterior_plot` for chi-square BFF objects is currently implemented only for a single chi-square statistic.")

  tau2 <- .posterior_plot_selected_tau2(x)
  if(is.null(tau2)){
    tau2 <- if(x$input$LRT){
      get_LRT_tau2(n = x$input$n, k = x$input$df, w = omega, r = r)
    }else{
      get_count_tau2(n = x$input$n, k = x$input$df, w = omega, r = r)
    }
  }

  if(tau2 <= 0)
    stop("There is no non-local prior distribution that provides more evidence for the null hypothesis than the null prior distribution.")

  branch_sign <- .posterior_plot_branch_sign(x, effect_size)
  lik.prior <- .effect_size_density(
    test_type   = x$test_type,
    effect_size = effect_size,
    x           = x_seq,
    density     = function(effect_size) .chi2_test.prior(tau2 = tau2, r = r, effect_size = effect_size, n = x$input$n, df = x$input$df),
    input       = x$input,
    table_dim   = table_dim,
    table_margins = table_margins,
    branch_sign = branch_sign
  )
  lik.posterior <- .effect_size_density(
    test_type   = x$test_type,
    effect_size = effect_size,
    x           = x_seq,
    density     = function(effect_size) .chi2_test.posterior(chi2_stat = chi2_stat, tau2 = tau2, r = r, effect_size = effect_size, n = x$input$n, df = x$input$df),
    input       = x$input,
    table_dim   = table_dim,
    table_margins = table_margins,
    branch_sign = branch_sign
  )

  df <- data.frame(
    x         = x_seq,
    prior     = lik.prior,
    posterior = lik.posterior
  )

  return(df)
}

.posterior_plot_data.f_test <- function(x, prior, x_limit, effect_size){

  x_seq <- seq(x_limit[1], x_limit[2], length.out = 301)

  f_stat <- x$input$f_stat
  r      <- x$r
  omega  <- .posterior_plot_selected_omega(x)

  if(length(f_stat) != 1 || length(x$input$n) != 1 || length(x$input$df1) != 1 || length(x$input$df2) != 1)
    stop("`posterior_plot` for F-test BFF objects is currently implemented only for a single F statistic.")

  tau2 <- .posterior_plot_selected_tau2(x)
  if(is.null(tau2)){
    tau2 <- get_linear_tau2(n = x$input$n, w = omega, k = x$input$df1, r = r)
  }

  if(tau2 <= 0)
    stop("There is no non-local prior distribution that provides more evidence for the null hypothesis than the null prior distribution.")

  lik.prior <- .effect_size_density(
    test_type   = "f_test",
    effect_size = effect_size,
    x           = x_seq,
    density     = function(effect_size) .f_test.prior(tau2 = tau2, r = r, effect_size = effect_size, n = x$input$n, df1 = x$input$df1),
    input       = x$input
  )
  lik.posterior <- .effect_size_density(
    test_type   = "f_test",
    effect_size = effect_size,
    x           = x_seq,
    density     = function(effect_size) .f_test.posterior(f_stat = f_stat, tau2 = tau2, r = r, effect_size = effect_size, n = x$input$n, df1 = x$input$df1, df2 = x$input$df2),
    input       = x$input
  )

  df <- data.frame(
    x         = x_seq,
    prior     = lik.prior,
    posterior = lik.posterior
  )

  return(df)
}

.posterior_plot_data.regression_test <- function(x, prior, x_limit, effect_size){

  x_seq <- seq(x_limit[1], x_limit[2], length.out = 301)

  t_stat <- x$input$t_stat
  r      <- x$r
  omega  <- .posterior_plot_selected_omega(x)

  if(length(t_stat) != 1 || length(x$input$n) != 1 || length(x$input$k) != 1)
    stop("`posterior_plot` for regression-test BFF objects is currently implemented only for a single t statistic.")

  tau2 <- .posterior_plot_selected_tau2(x)
  if(is.null(tau2)){
    tau2 <- get_regression_tau2(n = x$input$n, k = x$input$k, w = omega, r = r)
  }

  if(tau2 <= 0)
    stop("There is no non-local prior distribution that provides more evidence for the null hypothesis than the null prior distribution.")

  alternative <- .posterior_plot_alternative(x)
  one_sided  <- x$input$alternative != "two.sided"
  lik.prior <- .effect_size_density(
    test_type   = "regression_test",
    effect_size = effect_size,
    x           = x_seq,
    density     = function(effect_size) .regression_test.prior(tau2 = tau2, r = r, effect_size = effect_size, n = x$input$n, k = x$input$k, one_sided = one_sided),
    input       = x$input,
    alternative = alternative
  )
  lik.posterior <- .effect_size_density(
    test_type   = "regression_test",
    effect_size = effect_size,
    x           = x_seq,
    density     = function(effect_size) .regression_test.posterior(t_stat = t_stat, tau2 = tau2, r = r, effect_size = effect_size, n = x$input$n, k = x$input$k, one_sided = one_sided),
    input       = x$input,
    alternative = alternative
  )

  df <- data.frame(
    x         = x_seq,
    prior     = lik.prior,
    posterior = lik.posterior
  )

  return(df)
}
