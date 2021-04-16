setOldClass("ggplot")

#' QcPlot
#'
#' @export
#' @import ggplot2
setClass("QcPlot",
         representation = list(
           plot = "ggplot",
           file = "character",
           width = "numeric",
           height = "numeric",
           data = "data.frame"
         ))
