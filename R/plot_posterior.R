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
#' @param ... additional arguments
#'
#' @return
#' @export
#'
#' @examples
posterior_plot <- function(x, prior = FALSE, plot = TRUE, ...){

  # this is the generic function for prior and posterior plot for all tests
  # the individual specific functions are implemented under
  # ".posterior_plot." + test_type

  if(!is.BFF(x))
    stop("`posterior_plot` is implemented only for `BFF` objects.")

  # catch additional arguments (i.e., xlim, col, lwd, etc...)
  dots <- list(...)

  # manual dispatching for different tests (could have been done with classes)
  if(x[["test_type"]] == "t_test"){
    out <- .posterior_plot.t_test(x, prior, plot, dots)
  }


  return(out)
}


# test specific plotting functions
.posterior_plot.t_test <- function(x, prior, plot, dots){

  # TODO: use dots to set up range, colors, etc...
  x_seq <- seq(-3, 3, 0.01)

  # extract fitting information from the model object
  t_stat <- x$input$t_stat
  df     <- x$input$df
  r      <- x$r
  if(x$one_sample){
    tau2 <- get_one_sample_tau2(n = x$input$n, w = x$omega, r = r)
  }else{
    tau2 <- get_two_sample_tau2(n1 = x$input$n1, n2 = x$input$n2, w = x$omega, r = r)
  }

  # compute prior and posterior
  lik.prior     <- .t_test.prior(x_seq, tau2, r, one_sample = x$one_sample, one_sided = x$alternative != "two.sided")
  lik.posterior <- .t_test.posterior(x_seq, tau2, r, t_stat, df, one_sample = x$one_sample, one_sided = x$alternative != "two.sided")

  # create data.frame with values
  df <- data.frame(
    x         = x_seq,
    prior     = lik.prior,
    posterior = lik.posterior
  )

  # return data if no plot requested
  if(!plot){
    return(df)
  }

  # create plot otherwise
  out <- ggplot2::ggplot(df, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = prior), color = "black") +
    ggplot2::labs(x = expression(lambda))

  if(prior){
    out <- out + ggplot2::geom_line(ggplot2::aes(y = posterior), color = "blue")
  }

  return(out)
}
