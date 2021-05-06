setOldClass(c("theme", "gg"))

setClass("ssvQC",
         representation = list(
           feature_config = "QcConfigFeatures",
           signal_config = "QcConfigSignal",
           signal_data = "list",
           other_data = "list",
           out_dir = "character",
           bfc = "BiocFileCache",
           saving_enabled = "logical",
           plot_post_FUN = "function",
           plots = "list"
         ))

setClass("ssvQC.featureOnly", contains = "ssvQC")
setClass("ssvQC.signalOnly", contains = "ssvQC")
setClass("ssvQC.complete", contains = "ssvQC")

setMethod("initialize","ssvQC", function(.Object,...){
  .Object <- callNextMethod()
  validObject(.Object)
  .Object@plot_post_FUN = function(p)p
  .Object@plots = list()
  .Object@signal_data = list()
  .Object
})

#' ssvQC
#'
#' @param feature_config 
#' @param signal_config 
#' @param out_dir 
#' @param bfc 
#'
#' @return
#' @export
#' @import BiocFileCache
#'
#' @examples
#' options(mc.cores = 10)
#' set.seed(0)
#' feature_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
#' feature_config = QcConfigFeatures.parse(feature_config_file)
#' 
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' bam_config = QcConfigSignal.parse(bam_config_file)
#' 
#' bigwig_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bigwig_config.csv")
#' bigwig_config = QcConfigSignal.parse(bigwig_config_file)
#' 
#' sqc.complete.file = ssvQC(feature_config_file, bam_config_file)
#' 
#' sqc.complete = ssvQC(feature_config, bam_config)
#' 
#' sqc.complete.bw = ssvQC(feature_config, bigwig_config_file)
#' 
#' sqc.signal = ssvQC(signal_config = bam_config)
#' 
#' sqc.feature = ssvQC(feature_config = feature_config)
#' 
#' ssvQC.runAll(sqc.complete)
#' ssvQC.runAll(sqc.signal)
#' ssvQC.runAll(sqc.feature)
#' 
#' sqc = sqc.complete
#' object = sqc
#' 
#' sqc = ssvQC.prepMappedReads(sqc)
#' 
#' p_reads = ssvQC.plotMappedReads(sqc)
#' 
#' sqc = ssvQC.prepSignal(sqc)
#' object = sqc
#' 
#' sqc = ssvQC.prepSCC(sqc)
#' sqc = ssvQC.prepFRIP(sqc)
#' 
#' sqc = ssvQC.plotSCC(sqc)
#' object = sqc
#' sqc@plots
#' 
#' sqc@other_data
ssvQC = function(feature_config = NULL,
                 signal_config = NULL,
                 out_dir = getwd(),
                 bfc = NULL){
  if(is.null(feature_config) & is.null(signal_config)){
    stop("At least one of feature_config or signal_config must be specified.")
  }
  
  if(!is.null(feature_config)){
    if(is.character(feature_config)){
      if(file.exists(feature_config)){
        feature_config = QcConfigFeatures.parse(feature_config)
      }
    }
    if(!"QcConfigFeatures" %in% class(feature_config)){
      stop("feature_config must be either a QcConfigFeatures object or the path to valid configuration file to create one.")
    }  
    stopifnot(file.exists(feature_config@meta_data$file))
  }
  if(!is.null(signal_config)){
    if(is.character(signal_config)){
      if(file.exists(signal_config)){
        signal_config = QcConfigSignal.parse(signal_config)
      }
    }
    if(!"QcConfigSignal" %in% class(signal_config)){
      stop("signal_config must be either a QcConfigSignal object or the path to valid configuration file to create one.")
    }  
    stopifnot(file.exists(signal_config@meta_data$file))
  }
  
  if(is.null(bfc)){
    bfc = BiocFileCache::BiocFileCache()
  }
  
  dir.create(out_dir, showWarnings = FALSE)
  
  if(!is.null(feature_config) & !is.null(signal_config)){
    new("ssvQC.complete",
        feature_config = feature_config,
        signal_config = signal_config,
        signal_data = list(),
        other_data = list(),
        out_dir = out_dir,
        bfc = bfc,
        saving_enabled = TRUE
    )
  }else if(!is.null(feature_config)){
    new("ssvQC.featureOnly",
        feature_config = feature_config,
        signal_config = QcConfigSignal.null(),
        signal_data = list(),
        other_data = list(),
        out_dir = out_dir,
        bfc = bfc,
        saving_enabled = TRUE
    )
  }else if(!is.null(signal_config)){
    new("ssvQC.signalOnly",
        feature_config = QcConfigFeatures.null(),
        signal_config = signal_config,
        signal_data = list(),
        other_data = list(),
        out_dir = out_dir,
        bfc = bfc,
        saving_enabled = TRUE
    )
  }else{
    stop("At least one of feature_config or signal_config must be specified. This should have been caught earlier.")
  }
  
  
}

.make_query_signal_config = function(sc){
  out = split(sc)
  names(out) = paste0(names(out), "_signal")
  out
}

#' Title
#'
#' @param object 
#'
#' @return
#' @export
#'
#' @examples
setGeneric("ssvQC.runAll", function(object){standardGeneric("ssvQC.runAll")})
setMethod("ssvQC.runAll", "ssvQC.complete", function(object){
  message("complete")
})
setMethod("ssvQC.runAll", "ssvQC.featureOnly", function(object){
  message("featureOnly")
})
setMethod("ssvQC.runAll", "ssvQC.signalOnly", function(object){
  message("signalOnly")
})

##MappedReads

#' ssvQC.prepMappedReads
#'
#' @param object 
#'
#' @return
#' @export
#'
#' @examples
setGeneric("ssvQC.prepMappedReads", function(object){standardGeneric("ssvQC.prepMappedReads")})
setMethod("ssvQC.prepMappedReads", "ssvQC.complete", function(object){
  message("complete")
  object@signal_config = ssvQC.prepMappedReads(object@signal_config)
  object
})
setMethod("ssvQC.prepMappedReads", "ssvQC.featureOnly", function(object){
  stop("Cannot run mapped reads on ssvQC with no QcConfigSignal component")
  message("featureOnly")
})
setMethod("ssvQC.prepMappedReads", "ssvQC.signalOnly", function(object){
  message("signalOnly")
  object@signal_config = ssvQC.prepMappedReads(object@signal_config)
  object
})
setMethod("ssvQC.prepMappedReads", c("QcConfigSignal"), function(object){
  message("QcConfigSignal")
  #bam specific independent of peaks
  if(grepl("bam", object@read_mode)){
    object@meta_data$mapped_reads = sapply(object@meta_data$file, get_mapped_reads)
  }else{
    object@meta_data$mapped_reads = NA
  }
  object
})

#' ssvQC.plotMappedReads
#'
#' @param object 
#'
#' @return
#' @export
#'
#' @examples
setGeneric("ssvQC.plotMappedReads", function(object, out_dir){standardGeneric("ssvQC.plotMappedReads")})
setMethod("ssvQC.plotMappedReads", "ssvQC.complete", function(object){
  message("complete")
  ssvQC.plotMappedReads(object@signal_config, ifelse(object@saving_enabled, object@out_dir, "SKIP_SAVING"))
})
setMethod("ssvQC.plotMappedReads", "ssvQC.featureOnly", function(object){
  stop("Cannot run mapped reads on ssvQC with no QcConfigSignal component")
  message("featureOnly")
})
setMethod("ssvQC.plotMappedReads", "ssvQC.signalOnly", function(object){
  message("signalOnly")
  ssvQC.plotMappedReads(object@signal_config, ifelse(object@saving_enabled, object@out_dir, "SKIP_SAVING"))
  
})
setMethod("ssvQC.plotMappedReads", c("QcConfigSignal"), function(object){
  message("QcConfigSignal")
  ssvQC.plotMappedReads(object, "SKIP_SAVING")
})
setMethod("ssvQC.plotMappedReads", c("QcConfigSignal", "character"), function(object, out_dir){
  message("QcConfigSignal, character")
  res_file = function(f)file.path(out_dir, f)
  
  #bam specific independent of peaks
  if(grepl("bam", object@read_mode)){
    
    bam_config_dt = object@meta_data
    if(is.null(bam_config_dt$mapped_reads)){
      stop("ssvQC.prepMappedReads has not been called.")
    }
    
    color_var = object@color_by
    group_var = object@run_by
    color_mapping = object@color_mapping
    
    theme_set(theme(panel.background = element_blank(), axis.text.x = element_text(size = 8)))
    
    p_mapped_reads = ggplot(bam_config_dt, aes_string(x = "mark", y = "mapped_reads", fill = color_var)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = color_mapping) +
      scale_y_continuous(labels = function(x)x/1e6) +
      labs(y = "M mapped reads", fill = color_var, x= "") +
      labs(title = "Mapped reads")
    
    
    if(!out_dir == "SKIP_SAVING") ggsave(res_file("mapped_reads.pdf"), p_mapped_reads, width = 2+.5*nrow(bam_config_dt), height = 3)
    p_mapped_reads
  }else{
    warning("Could not run mapped reads on ")
    NULL
  }
})

##Signal
#' ssvQC.prepSignal
#'
#' @param object 
#'
#' @return
#' @export
#'
#' @examples
setGeneric("ssvQC.prepSignal", function(object){standardGeneric("ssvQC.prepSignal")})
setMethod("ssvQC.prepSignal", "ssvQC.complete", function(object){
  message("complete")
  object@signal_data = lapply(object@feature_config$assessment_features, function(query_gr){
    sig_configs = .make_query_signal_config(object@signal_config)
    lapply(sig_configs, function(sel_sig_config){
      ClusteredSignal.fromConfig(sel_sig_config, 
                                 query_gr, 
                                 signal_var = "y", 
                                 facet_var = "name_split", 
                                 extra_var = union(object@signal_config@color_by, object@signal_config@run_by))
    })
  })
  object
})
setMethod("ssvQC.prepSignal", "ssvQC.featureOnly", function(object){
  message("featureOnly")
  stop("Cannot run prepSignal on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.prepSignal", "ssvQC.signalOnly", function(object){
  message("signalOnly")
  stop("Cannot run prepSignal on ssvQC with no QcConfigFeature component")
})

##FRIP
#' ssvQC.prepFRIP
#'
#' @param object 
#'
#' @return
#' @export
#'
#' @examples
setGeneric("ssvQC.prepFRIP", function(object){standardGeneric("ssvQC.prepFRIP")})
setMethod("ssvQC.prepFRIP", "ssvQC.complete", function(object){
  message("complete")
  FRIP_data = lapply(object@feature_config$assessment_features, function(query_gr){
    sig_configs = .make_query_signal_config(object@signal_config)
    lapply(sig_configs, function(sel_sig_config){
      make_frip_dt(as.data.table(sel_sig_config@meta_data), query_gr = query_gr)
    })
  })
  object@other_data$FRIP = FRIP_data
  object
})
setMethod("ssvQC.prepFRIP", "ssvQC.featureOnly", function(object){
  message("featureOnly")
  stop("Cannot run prepSignal on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.prepFRIP", "ssvQC.signalOnly", function(object){
  message("signalOnly")
  stop("Cannot run prepSignal on ssvQC with no QcConfigFeature component")
})

#' ssvQC.plotFRIP
#'
#' @param object 
#'
#' @return
#' @export
#'
#' @examples
setGeneric("ssvQC.plotFRIP", function(object){standardGeneric("ssvQC.plotFRIP")})
setMethod("ssvQC.plotFRIP", "ssvQC.complete", function(object){
  message("complete")
  if(is.null(object@other_data$FRIP)){
    stop("ssvQC.prepFRIP has not been called.")
  }
  FRIP_data = object@other_data$FRIP
  FRIP_data = lapply(object@feature_config$assessment_features, function(query_gr){
    sig_configs = .make_query_signal_config(object@signal_config)
    lapply(sig_configs, function(sel_sig_config){
      make_frip_dt(as.data.table(sel_sig_config@meta_data), query_gr = query_gr)
    })
  })
  object@other_data$FRIP = FRIP_data
  object
})
setMethod("ssvQC.plotFRIP", "ssvQC.featureOnly", function(object){
  message("featureOnly")
  stop("Cannot run prepSignal on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.plotFRIP", "ssvQC.signalOnly", function(object){
  message("signalOnly")
  stop("Cannot run prepSignal on ssvQC with no QcConfigFeature component")
})

##SCC
#' ssvQC.prepSCC
#'
#' @param object 
#'
#' @return
#' @export
#'
#' @examples
setGeneric("ssvQC.prepSCC", function(object){standardGeneric("ssvQC.prepSCC")})
setMethod("ssvQC.prepSCC", "ssvQC.complete", function(object){
  message("complete")
  SCC_data = lapply(object@feature_config$assessment_features, function(query_gr){
    sig_configs = .make_query_signal_config(object@signal_config)
    lapply(sig_configs, function(sel_sig_config){
      make_scc_dt(as.data.table(sel_sig_config@meta_data), query_gr = query_gr, bfc_corr = object@bfc)
    })
  })
  object@other_data$SCC = SCC_data
  object
})
setMethod("ssvQC.prepSCC", "ssvQC.featureOnly", function(object){
  message("featureOnly")
  stop("Cannot run prepSignal on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.prepSCC", "ssvQC.signalOnly", function(object){
  message("signalOnly")
  stop("Cannot run prepSignal on ssvQC with no QcConfigFeature component")
})

dbl_names = function(name_1, name_2){
  paste(sub("_signal", " signal", name_2), "at", sub("_features", " features", name_1))
}

dbl_lapply = function(in_list, FUN, FUN_names = dbl_names){
  out_list = lapply(names(in_list), function(name_1){
    item_1 = in_list[[name_1]]
    out = lapply(names(item_1), function(name_2){
      item_2 = item_1[[name_2]]
      if(!is.null(FUN_names){
        mt = dbl_names(name_1, name_2)
        FUN(item_2, mt)
      })else{
        FUN(item_2)
      }
    })
    names(out) = names(item_1)
    out
  })
  names(out_list) = names(in_list)
  out_list
}

#' ssvQC.plotSCC
#'
#' @param object 
#'
#' @return
#' @export
#'
#' @examples
setGeneric("ssvQC.plotSCC", function(object){standardGeneric("ssvQC.plotSCC")})
setMethod("ssvQC.plotSCC", "ssvQC.complete", function(object){
  message("complete")
  if(is.null(object@other_data$SCC)){
    stop("ssvQC.prepSCC has not been called.")
  }
  SCC_data = object@other_data$SCC
  SCC_plots = lapply(names(SCC_data), function(scc_features_name){
    scc_features = SCC_data[[scc_features_name]]
    out = lapply(names(scc_features), function(scc_signal_name){
      scc_signal = scc_features[[scc_signal_name]]
      mt = paste(sub("_signal", " signal", scc_signal_name), "at", sub("_features", " features", scc_features_name))
      plot_scc_dt(scc_signal, main_title = mt)
    })
    names(out) = names(scc_features)
    out
  })
  names(SCC_plots) = names(SCC_data)
  SCC_dots = lapply(SCC_plots, function(scc_features){
    lapply(scc_features, function(scc_signal){
      scc_signal$scc_dots
    })
  })
  SCC_curves = lapply(SCC_plots, function(scc_features){
    lapply(scc_features, function(scc_signal){
      scc_signal$scc_curves
    })
  })
  
  object@plots$SCC_dots = SCC_dots
  object@plots$SCC_curves = SCC_curves
  object
})
setMethod("ssvQC.plotSCC", "ssvQC.featureOnly", function(object){
  message("featureOnly")
  stop("Cannot run prepSignal on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.plotSCC", "ssvQC.signalOnly", function(object){
  message("signalOnly")
  stop("Cannot run prepSignal on ssvQC with no QcConfigFeature component")
})



