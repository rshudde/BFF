
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

plot_BFF = function(effect_size, BFF,  save = FALSE, savename=NULL, xlab = NULL, ylab=NULL, main=NULL, r = NULL)
{
  num_lines = length(BFF)

  plot_legend_names = vector()
  count = 1
  for (i in c(1:length(r)))
  {
    str_temp = paste("r=", as.character(r[i]), sep="")
    plot_legend_names[count] = str_temp
    count = count + 1
  }

  # custom plotting
  if (is.null(xlab)) xlab = expression(paste("RMSE ", tilde(omega)))
  if (is.null(ylab)) ylab = "Bayes Factor Against Null Hypothesis"
  if (is.null(main)) main = "BFF"
  if (is.null(savename) && save)
  {
    ## TODO string check for .pdf extension
    print("No savename argument given, plot saving as BFF_plot.pdf")
    savename = "BFF_plot.pdf"
  }
  if (!is.null(savename) && !save)
  {
    ## TODO string check for .pdf extension
    check_pdf = substrRight(save, 4)
    if (check_pdf != ".pdf") savename = paste(savename, ".pdf", sep="")
    save=TRUE
    print(savename)
  }

  # if (save) pdf(savename)

  # x11()
  plot(effect_size, BFF[[1]],
       type='l',
       xlab= xlab,
       ylab = ylab,
       main = main,
       lwd=1,
       cex=2,
       lty = 1,
       yaxt="n",
       col=1)

  if (num_lines > 1)
  {
    count = 2
    for (k in c(2:num_lines))
    {
      values = BFF[[k]]
      graphics::lines(effect_size, values, type = "l", lwd = 1, cex = 2, lty = count, col=k)
      count = count+1
    }
  }

  maxval = 1e50

  if (num_lines > 1)
  {
    graphics::legend("topleft", legend=plot_legend_names,
                     col=c(1:num_lines), lty=c(1:length(BFF)), cex=0.6, bty="n")
  }


  max_of_all_BFF = max(unlist(lapply(BFF, max)))


  # add the lines for the colors
  if (min(unlist(lapply(BFF, max)))< 0) abline(h = log(1))
  rect(-1, -log(maxval), 0.1, log(maxval), col=adjustcolor("red", 0.1))
  rect(0.1, -log(maxval), .35, log(maxval), col=adjustcolor("orange", 0.1))
  rect(.35, -log(maxval), .65, log(maxval), col=adjustcolor("blue", 0.1))
  rect(.65, -log(maxval), 2, log(maxval), col=adjustcolor("green", 0.1))

  positive = c(50000, 10000, 5000, 2000, 1000, 500, 200, 150, 100, 75,  50, 40, 30, 20, 10, 5, 2)
  negative = -1 * positive
  positive_labels = vector()
  negative_labels = vector()
  for (k in 1:length(positive)) {
    positive_labels[k] = paste(as.character(format(positive[k],scientific=FALSE)), ":1", sep = "")
    negative_labels[k] = paste(":1", as.character(format(negative[k], scientific=FALSE)), sep = "")
  }
  axis_position = c(positive, log(1), negative)
  axis_labels = c(positive_labels, "1:1", negative_labels)

  axis(2, at=axis_position,labels = axis_labels)
  if(save) {
    print("Will fix this with ggplot2")
    # dev.copy2pdf(file = savename)
    # dev.off()
  }
}
