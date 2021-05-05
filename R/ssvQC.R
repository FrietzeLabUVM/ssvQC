setOldClass(c("theme", "gg"))

setClass("ssvQC",
         representation = list(
           feature_config = "QcConfigFeatures",
           signal_config = "QcConfigSignal",
           signal_data = "ClusteredSignal",
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
  .Object@signal_data = ClusteredSignal.null()
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
        signal_data = ClusteredSignal.null(),
        out_dir = out_dir,
        bfc = bfc,
        saving_enabled = TRUE
    )
  }else if(!is.null(feature_config)){
    new("ssvQC.featureOnly",
        feature_config = feature_config,
        signal_config = QcConfigSignal.null(),
        signal_data = ClusteredSignal.null(),
        out_dir = out_dir,
        bfc = bfc,
        saving_enabled = TRUE
    )
  }else if(!is.null(signal_config)){
    new("ssvQC.signalOnly",
        feature_config = QcConfigFeatures.null(),
        signal_config = signal_config,
        signal_data = ClusteredSignal.null(),
        out_dir = out_dir,
        bfc = bfc,
        saving_enabled = TRUE
    )
  }else{
    stop("At least one of feature_config or signal_config must be specified. This should have been caught earlier.")
  }
  
  
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
  object@signal_data = ClusteredSignal.fromConfig(object@signal_config, object@feature_config$assessment_features, signal_var = "y", facet_var = "name_split", extra_var = union(object@signal_config@color_by, object@signal_config@run_by))
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

