#' @title Plot Bayes Factor Function
#'
#' @description Creates a Bayes factor function plot
#' of an BFF object. The BFF object needs to be specified
#' with \code{omega = NULL} or .
#'
#' @param x a BFF object
#' @param plot whether plot should be generated.
#' Default to \code{TRUE}. If \code{FALSE} a data
#' frame with the support, prior ordinates, and
#' posterior ordinates is returned instead.
#' @param ... additional arguments to the plotting
#' function. These include: \describe{
#'  \item{"title"}{tile of the figure}
#'  \item{"xlab"}{x-axis label of the figure}
#'  \item{"ylab"}{y-axis label of the figure}
#'  \item{"effect_size"}{effect-size scale for plotting
#'  the Bayes factor function. Defaults to the scale used
#'  to fit the BFF object, or to the package's internal
#'  omega scale.}
#'  \item{"table_dim"}{integer vector \code{c(rows, columns)}
#'  for chi-square transformations that require table dimensions.}
#'  \item{"add_segments"}{whether effect size
#'  segments should be added to the figure. Available only
#'  for standardized effect sizes. Defaults to \code{TRUE}}
#' }
#'
#' @return either a ggplot2 object if \code{plot = TRUE} or a data.frame
#' with a Bayes factor function if \code{plot = FALSE}
#'
#' @seealso [z_test_BFF()], [t_test_BFF()], [chi2_test_BFF()], [f_test_BFF()]
#'
#' @export
plot.BFF = function(x, plot = TRUE,  ...) {

  if (is.null(x$BFF))
    stop("Bayes factor function can be plotted only if a specific omega/tau2 is not user set")

  dots <- list(...)
  effect_size <- if(!x$generic_test) .effect_size_for_object(x, dots[["effect_size"]]) else NULL
  table_dim <- if(is.null(dots[["table_dim"]])) x$input$table_dim else dots[["table_dim"]]

  # extract the BFF object (deal with generic/specific effect sizes)
  log_BF = NULL
  if(x$generic_test){
    df <- data.frame(
      x      = x$BFF$tau2,
      log_BF = x$BFF$log_bf
    )
  }else{
    df <- data.frame(
      x      = .effect_size_from_internal(
        value       = x$BFF$omega,
        test_type   = x$test_type,
        effect_size = effect_size,
        input       = x$input,
        table_dim   = table_dim
      ),
      log_BF = x$BFF$log_bf
    )
  }

  if(!plot)
    return(df)

  # catch and set additional arguments
  title <- if(is.null(dots[["title"]])) "Bayes Factor Function"                else dots[["title"]]
  ylab  <- if(is.null(dots[["ylab"]]))  "Bayes Factor Against Null Hypothesis" else dots[["ylab"]]
  xlab  <- if(is.null(dots[["xlab"]])){
    if(x$generic_test) expression(tau^2) else .effect_size_plot_label(x, dots[["effect_size"]])
  }else dots[["xlab"]]
  add_segments <- if(is.null(dots[["add_segments"]])) TRUE else dots[["add_segments"]]

  out <- ggplot2::ggplot(df) +
    ggplot2::geom_line(ggplot2::aes(x = x, y = log_BF)) +
    ggplot2::ggtitle(title) + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)

  # adding rectangles (these are test specific)
  if(add_segments && !x$generic_test){
    effect_size_cutpoints <- .get_effect_size_cutpoints(x$test_type, effect_size = effect_size)
    effect_size_range     <- .get_effect_size_range(x$test_type)
    effect_size_colors    <- c(
      grDevices::adjustcolor("red"   , 0.1),
      grDevices::adjustcolor("orange", 0.1),
      grDevices::adjustcolor("blue",   0.1),
      grDevices::adjustcolor("green",  0.1)
    )

    if(length(effect_size_cutpoints) > 0){
      for(i in 1:(length(effect_size_cutpoints)+1)){
        out <- out + ggplot2::annotate(
          "rect",
          xmin = if(i == 1) -Inf else effect_size_cutpoints[i-1],
          xmax = if(i == 1) effect_size_cutpoints[1] else if(i == length(effect_size_cutpoints)+1) Inf else effect_size_cutpoints[i],
          ymin = -Inf,
          ymax =  Inf,
          fill = effect_size_colors[i]
        )
      }

      out <- out + ggplot2::geom_vline(
        xintercept = effect_size_cutpoints,
        lwd        = 0.2)
    }
  }

  # add log labels to Bayes factors
  y_ticks  <- sort(do.call(c, lapply(c(3,10),function(x) x*10^(0:nchar(ceiling(exp(max(abs(df$log_BF)))/(10*x)))))))
  # y_ticks1  <- sort(do.call(c, lapply(c(2,10),function(x) x*10^(0:nchar(ceiling(max(abs(df$log_BF)))/(10*x))))))
  #
  # print(max(abs(df$log_BF)))
  # print(df$log_BF)
  # print(y_ticks)
  # print(y_ticks1)
  y_labels <- c(
    paste0("1:",rev(y_ticks)),
    1,
    paste0(y_ticks,":1")
  )
  y_ticks  <- c(1/rev(y_ticks), 1, y_ticks)

  out <- out +
    ggplot2::scale_y_continuous(breaks = log(y_ticks), labels = y_labels) +
    ggplot2::scale_x_continuous(expand = c(0, 0.05)) +
    ggplot2::geom_hline(yintercept = 0, lwd = 0.2) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.grid = ggplot2::element_blank()
    )

  return(out)
}

.get_effect_size_cutpoints <- function(test_type, effect_size = NULL){
  if(is.null(effect_size)){
    return(switch(test_type, # for sources, see each test in https://www.utstat.toronto.edu/~brunner/oldclass/378f16/readings/CohenPower.pdf if not listed below
           "t_test"    = c(0.2, 0.5, 0.8),
           "z_test"    = c(0.2, 0.5, 0.8), # https://pressbooks.bccampus.ca/statspsych/chapter/chapter-11/
           "chi2_test" = c(0.1, 0.3, 0.5),
           "regression_test" = c(0.02, 0.15, 0.35),
           "f_test"    = c(0.1, 0.25, 0.4)))
  }

  effect_size <- .effect_size_normalize(test_type, effect_size)

  switch(
    test_type,
    "t_test" = ,
    "z_test" = c(0.2, 0.5, 0.8),
    "chi2_test" = switch(
      effect_size,
      "omega" = ,
      "cohens_w" = ,
      "phi" = ,
      "cramers_v" = ,
      "tschuprow_t" = c(0.1, 0.3, 0.5),
      numeric(0)
    ),
    "f_test" = switch(
      effect_size,
      "cohens_f" = c(0.1, 0.25, 0.4),
      "cohens_f2" = c(0.01, 0.0625, 0.16),
      "partial_eta2" = ,
      "partial_r2" = c(0.01, 0.06, 0.14),
      numeric(0)
    ),
    "regression_test" = switch(
      effect_size,
      "cohens_f" = sqrt(c(0.02, 0.15, 0.35)),
      "cohens_f2" = ,
      "partial_r2" = c(0.02, 0.15, 0.35),
      "partial_r" = sqrt(c(0.02, 0.15, 0.35)),
      numeric(0)
    ),
    numeric(0)
  )
}
.get_effect_size_range    <-  function(test_type){
  switch(test_type,
         "t_test"    = c(0.0, 1.0),
         "z_test"    = c(0.0, 1.0),
         "chi2_test" = c(0.0, 1.0),
         "regression_test" = c(0.0, 1.0),
         "f_test"    = c(0.0, 1.0))
}
