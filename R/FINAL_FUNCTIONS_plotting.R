require(ggplot2)

substrRight <- function(x, n) {
  substr(x, nchar(x) - n + 1, nchar(x))
}

plot_BFF = function(effect_size,
                    BFF,
                    save = FALSE,
                    savename = NULL,
                    xlab = NULL,
                    ylab = NULL,
                    main = NULL,
                    r = NULL) {
  num_lines = length(BFF)
  plot_legend_names = vector()
  count = 1

  for (i in c(1:length(r))) {
    str_temp = paste("r=", as.character(r[i]), sep = "")
    plot_legend_names[count] = str_temp
    count = count + 1
  }

  # custom plotting
  if (is.null(xlab))
    xlab = expression(paste("RMSE ", tilde(omega)))
  if (is.null(ylab))
    ylab = "Bayes Factor Against Null Hypothesis"
  if (is.null(main))
    main = "BFF"

  if (is.null(savename) && save) {
    print("No savename argument given, plot saving as BFF_plot.pdf")
    savename = "BFF_plot.pdf"
  }
  if (!is.null(savename) && !save) {
    ## TODO string check for .pdf extension
    check_pdf = substrRight(savename, 4)
    if (check_pdf != ".pdf")
      savename = paste(savename, ".pdf", sep = "")
    save = TRUE
    print(savename)
  }
  #Making the plot
  data <- data.frame(effect_size = effect_size, BFF = BFF[[1]])
  p <- ggplot(data, aes(x = effect_size, y = BFF)) +
    geom_line() +
    theme_bw()
  p <- p + labs(x = xlab, y = ylab) +
    ggtitle(main) +
    theme(plot.title = element_text(hjust = 0.5))

  maxval = 1e50

  if (num_lines > 1)
  {
    count = 2
    for (k in c(2:num_lines))
    {
      values = BFF[[k]]
      data2 <- data.frame(effect_size = effect_size, BFF = BFF[[k]])
      p <- ggplot(data2, aes(x = effect_size, y = BFF)) +
        geom_line(linetype = count, colour = k) +
        theme_bw()
      count = count+1
    }
  }


  min_lim_x = min(effect_size)
  max_lim_x = max(effect_size)

  #Adding rectangles
  p <-
    p + annotate(
      "rect",
      xmin = -1,
      xmax = 0.1,
      ymin = -Inf,
      ymax = Inf,
      fill = adjustcolor("red", 0.1)
    )
  p <-
    p + annotate(
      "rect",
      xmin = 0.1,
      xmax = 0.35,
      ymin = -Inf,
      ymax = Inf,
      fill = adjustcolor("orange", 0.1)
    )
  p <-
    p + annotate(
      "rect",
      xmin = 0.35,
      xmax = 0.65,
      ymin = -Inf,
      ymax = Inf,
      fill = adjustcolor("blue", 0.1)
    )
  p <-
    p + annotate(
      "rect",
      xmin = 0.65,
      xmax = 2,
      ymin = -Inf,
      ymax = Inf,
      fill = adjustcolor("green", 0.1)
    )
  p <-
    p + geom_vline(xintercept = c(-1, 0.1, 0.35, 0.65, 2), lwd = 0.2)
  p <- p + coord_cartesian(xlim = c(min_lim_x, max_lim_x))

  #Ading y axis labels
  positive = c(50000,
               10000,
               5000,
               2000,
               1000,
               500,
               200,
               150,
               100,
               75,
               50,
               40,
               30,
               20,
               10,
               5,
               2)
  negative = -1*positive
  positive_labels = vector()
  negative_labels = vector()

  for (k in 1:length(positive)) {
    positive_labels[k] = paste(as.character(format(positive[k],
                                                   scientific = FALSE)),
                               ":1", sep = "")
    negative_labels[k] = paste("1:", as.character(format(positive[k],
                                                         scientific = FALSE)),sep = "")
  }
  axis_position = c(positive, log(1), negative)
  axis_labels = c(positive_labels, "1:1", negative_labels)

  p <-
    p + scale_y_continuous(breaks = axis_position, labels = axis_labels, guide = guide_axis(check.overlap = TRUE))+
    geom_hline(yintercept = log(1), lwd = 0.2)

  p <- p + theme(panel.grid = element_blank())

  if (save) {
    ggsave(savename, plot = p, device = "pdf")
  }

  print(p)
}
