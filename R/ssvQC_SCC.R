##SCC
#' @export
#' @rdname ssvQC
setGeneric("ssvQC.prepSCC", function(object){standardGeneric("ssvQC.prepSCC")})
setMethod("ssvQC.prepSCC", "ssvQC.complete", function(object){
  object = ssvQC.prepFetch(object)
  feature_names = names(object@features_config$assessment_features)
  names(feature_names) = feature_names
  SCC_data = lapply(feature_names, run_by_match, object = object, callFUN = make_scc_dt.run_by_match)
  object@other_data$SCC = SCC_data
  object
})
setMethod("ssvQC.prepSCC", "ssvQC.featureOnly", function(object){
  stop("Cannot run prepSCC on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.prepSCC", "ssvQC.signalOnly", function(object){
  stop("Cannot run prepSCC on ssvQC with no QcConfigFeature component")
})
#' @export
#' @rdname ssvQC
setGeneric("ssvQC.plotSCC", function(object){standardGeneric("ssvQC.plotSCC")})
setMethod("ssvQC.plotSCC", "ssvQC.complete", function(object){
  if(is.null(object@other_data$SCC)){
    object = ssvQC.prepSCC(object)
  }
  SCC_data = object@other_data$SCC
  
  wrap_plot_scc_dt = function(scc_dt, main_title){
    plot_scc_dt(scc_dt, main_title, name_lev = levels(object@signal_config@meta_data$name_split))
  }
  
  SCC_plots = dbl_lapply(SCC_data, FUN = wrap_plot_scc_dt, FUN_names = dbl_names)
  SCC_dots = dbl_extract(SCC_plots, "scc_dots")
  SCC_curves = dbl_extract(SCC_plots, "scc_curves")
  
  if(is.null(object@plots$SCC)){
    object@plots$SCC = list()
  }
  
  object@plots$SCC$dots = SCC_dots
  object@plots$SCC$curves = SCC_curves
  object
})
setMethod("ssvQC.plotSCC", "ssvQC.featureOnly", function(object){
  stop("Cannot run plotSCC on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.plotSCC", "ssvQC.signalOnly", function(object){
  stop("Cannot run plotSCC on ssvQC with no QcConfigFeature component")
})