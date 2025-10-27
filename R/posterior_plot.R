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
#'  \item{"x_limit"}{vector defining the plotting range,
#'  defaults to \code{c(-3, 3)}.}
#'  \item{"color"}{vector with color for the posterior and
#'  prior line. Defaults to \code{c("black", "grey")}}
#'  \item{"linetype"}{vector with linetype for the posterior and
#'  prior line. Defaults to \code{c(2, 1)}}
#'  \item{"linewidth"}{vector with linewidth for the posterior and
#'  prior line. Defaults to \code{c(1, 1)}}
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
  x_limit   <- if(is.null(dots[["x_limit"]])) switch(
    x$alternative,
    "two.sided" = c(-3, 3),
    "greater"   = c(0, 3),
    "less"      = c(-3, 0)
  ) else dots[["x_limit"]]
  # TODO: deal with positive/negative only plotting: maybe just implement everything for greater
  # (as in the t_test_BFF function) and then flip the support around y-axis


  # manual dispatching for different tests (could have been done with classes)
  if(x[["test_type"]] == "t_test"){
    plot_data <- .posterior_plot_data.t_test(x, prior, x_limit)
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
    x = if(x$generic_test) expression(tau^2) else .test_effect_size_name(x$test_type),
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

# test specific plotting functions
.posterior_plot_data.t_test <- function(x, prior, x_limit){

  x_seq <- seq(x_limit[1], x_limit[2], length.out = 301)

  # extract fitting information from the model object
  t_stat <- x$input$t_stat
  r      <- x$r
  if(x$input$one_sample){
    tau2 <- get_one_sample_tau2(n = x$input$n, w = x$omega, r = r)
  }else{
    tau2 <- get_two_sample_tau2(n1 = x$input$n1, n2 = x$input$n2, w = x$omega, r = r)
  }

  # terminate if tau2 is equal to 0 -> BFF leads to omega 0 (e.g., opposite direction in one-sided test)
  if(tau2 == 0)
    stop("There is no non-local prior distribution that provides more evidence for the null hypothesis than the null prior distribution.")

  # compute prior and posterior
  lik.prior     <- .t_test.prior(tau2 = tau2, r = r, effect_size = x_seq, n = x$input$n, n1 = x$input$n1, n2 = x$input$n2, one_sample = x$input$one_sample, one_sided = x$alternative != "two.sided")
  lik.posterior <- .t_test.posterior(t_stat = t_stat, tau2 = tau2, r = r, effect_size = x_seq, n = x$input$n, n1 = x$input$n1, n2 = x$input$n2, one_sample = x$input$one_sample, one_sided = x$alternative != "two.sided")

  # create data.frame with values
  posterior = NULL
  df <- data.frame(
    x         = x_seq,
    prior     = lik.prior,
    posterior = lik.posterior
  )

  return(df)
}
