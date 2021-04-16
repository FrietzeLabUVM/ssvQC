#' QcPlot
#'
#' @export
setClass("QcPlot",
         representation = list(
           plot = "ggplot",
           file = "character",
           width = "numeric",
           height = "numeric",
           data = "data.frame"
         ))
