setOldClass(c( "ggplot"))

#' QcPlot
#'
#' @export
#' @import ggplot2
setClass("QcPlot", contains = "ggplot",
         representation = list(
           plot = "ggplot",
           file = "character",
           width = "numeric",
           height = "numeric",
           data = "data.frame"
         ))


setMethod("initialize","QcPlot", function(.Object,...){
  .Object <- callNextMethod()
  validObject(.Object)
  .Object@plot = ggplot2::ggplot()
  .Object
})

# 
# library(ggplot2)
# isS3method("ggplot") # FALSE ,it means it is an S3 generic for more details: ?isS3method
# 
# # new("QcPlot")
# setGeneric("type", function(x) standardGeneric("type"))
# setMethod("type", signature("matrix"), function(x) "matrix")
# setMethod("type", signature("character"), function(x) "character")
# 
# foo <- structure(list(x = 1), class = "foo")
# type(foo)
# 
# setOldClass("foo")
# setMethod("type", signature("foo"), function(x) "foo")
# 
# type(foo)
# 
# setMethod("+", signature(e1 = "foo", e2 = "numeric"), function(e1, e2) {
#   structure(list(x = e1$x + e2), class = "foo")
# })
# foo + 3
