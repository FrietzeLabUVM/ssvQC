##SCC
#' @export
#' @rdname ssvQC
setGeneric("ssvQC.prepSCC", function(object){standardGeneric("ssvQC.prepSCC")})
setMethod("ssvQC.prepSCC", "ssvQC.complete", function(object){
  SCC_data = lapply(object@features_config$assessment_features, function(query_gr){
    sig_configs = .make_query_signal_config(object@signal_config)
    lapply(sig_configs, function(sel_sig_config){
      make_scc_dt(as.data.table(sel_sig_config@meta_data), query_gr = query_gr, bfc = object@bfc)
    })
  })
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