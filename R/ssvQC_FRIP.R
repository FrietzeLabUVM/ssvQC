##FRIP
#' @export
#' @rdname ssvQC
setGeneric("ssvQC.prepFRIP", function(object){standardGeneric("ssvQC.prepFRIP")})
setMethod("ssvQC.prepFRIP", "ssvQC.complete", function(object){
  object = ssvQC.prepFetch(object)
  feature_names = names(object@features_config$assessment_features)
  names(feature_names) = feature_names
  FRIP_data = lapply(feature_names, function(name){
    query_gr = object@features_config$assessment_features[[name]]
    sig_configs = .make_query_signal_config(object@signal_config)
    
    must_match = object@matched_only
    if(must_match){
      sig_name = feature_name2signal_name(name)
      if(!sig_name %in% names(sig_configs)){
        must_match = FALSE
      }else{
        out = lapply(sig_configs[feature_name2signal_name(name)], function(sel_sig_config){
          make_frip_dt(as.data.table(sel_sig_config@meta_data), query_gr = query_gr, color_var = sel_sig_config@color_by, bfc = object@bfc, force_overwrite = getOption("SQC_FORCE_CACHE_OVERWRITE", FALSE))
        })    
      }
    }
    if(!must_match){
      out = lapply(sig_configs, function(sel_sig_config){
        make_frip_dt(as.data.table(sel_sig_config@meta_data), query_gr = query_gr, color_var = sel_sig_config@color_by, bfc = object@bfc, getOption("SQC_FORCE_CACHE_OVERWRITE", FALSE))
      })  
    }
    out
  })
  object@other_data$FRIP = FRIP_data
  object
})
setMethod("ssvQC.prepFRIP", "ssvQC.featureOnly", function(object){
  stop("Cannot run prepFRIP on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.prepFRIP", "ssvQC.signalOnly", function(object){
  stop("Cannot run prepFRIP on ssvQC with no QcConfigFeature component")
})

#' @export
#' @rdname ssvQC
setGeneric("ssvQC.plotFRIP", function(object){standardGeneric("ssvQC.plotFRIP")})
setMethod("ssvQC.plotFRIP", "ssvQC.complete", function(object){
  if(is.null(object@other_data$FRIP)){
    object = ssvQC.prepFRIP(object)
  }
  FRIP_data = object@other_data$FRIP
  
  wrap_plot_frip_dt = function(frip_dt, main_title){
    plot_frip_dt(frip_dt, 
                 color_var = object@signal_config@color_by, 
                 color_mapping = object@signal_config@color_mapping,
                 main_title = main_title)
  }
  
  plots = dbl_lapply(FRIP_data, wrap_plot_frip_dt, dbl_names)
  
  if(is.null(object@plots$FRIP)){
    object@plots$FRIP = list()
  }
  
  object@plots$FRIP$reads_per_peak = dbl_extract(plots, "reads_per_peaks")
  object@plots$FRIP$per_peak = dbl_extract(plots, "frip_per_peaks")
  object@plots$FRIP$total = dbl_extract(plots, "frip_total")
  object
})
setMethod("ssvQC.plotFRIP", "ssvQC.featureOnly", function(object){
  stop("Cannot run plotFRIP on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.plotFRIP", "ssvQC.signalOnly", function(object){
  stop("Cannot run plotFRIP on ssvQC with no QcConfigFeature component")
})