setOldClass(c("theme", "gg"))

setClass("ssvQC",
         representation = list(
           feature_config = "QcConfigFeatures",
           signal_config = "QcConfigSignal",
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
#' 
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
        out_dir = out_dir,
        bfc = bfc,
        saving_enabled = TRUE
    )
  }else if(!is.null(feature_config)){
    new("ssvQC.featureOnly",
        feature_config = feature_config,
        signal_config = QcConfigSignal.null(),
        out_dir = out_dir,
        bfc = bfc,
        saving_enabled = TRUE
    )
  }else if(!is.null(signal_config)){
    new("ssvQC.signalOnly",
        feature_config = QcConfigFeatures.null(),
        signal_config = signal_config,
        out_dir = out_dir,
        bfc = bfc,
        saving_enabled = TRUE
    )
  }else{
    stop("At least one of feature_config or signal_config must be specified. This should have been caught earlier.")
  }
  
  
}

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


.run_mapped_reads = function(sqc){
  out_dir = sqc@out_dir
  res_file = function(f)file.path(out_dir, f)
  
  #bam specific independent of peaks
  if(!grepl("bam", sqc@signal_config@read_mode)){
    stop("Read mode must be bam to assess mapped reads.") 
  }
  
  bam_config_dt = sqc@signal_config@meta_data
  if(is.null(bam_config_dt$mapped_reads)){
    bam_config_dt$mapped_reads = sapply(bam_config_dt$file, get_mapped_reads)
  }
  
  color_var = sqc@signal_config@color_by
  group_var = sqc@signal_config@run_by
  color_mapping = sqc@signal_config@color_mapping
  
  p_mapped_reads = ggplot(bam_config_dt, aes_string(x = "mark", y = "mapped_reads", fill = color_var)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = color_mapping) +
    scale_y_continuous(labels = function(x)x/1e6) +
    labs(y = "M mapped reads", fill = color_var, x= "") +
    labs(title = "Mapped reads") +
    theme(panel.background = element_blank(), axis.text.x = element_text(size = 8))
  
  sqc@plots[["mapped_reads"]] = p_mapped_reads
  
  if(sqc@saving_enabled){
    ggsave(res_file("mapped_reads.pdf"), p_mapped_reads, width = 2+.5*nrow(bam_config_dt), height = 3)  
  }
  
  p_mapped_reads
}



#' Title
#'
#' @param sqc 
#'
#' @return
#' @export
#'
#' @examples
#' feature_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' 
#' sqc = ssvQC(feature_config_file, bam_config_file)
#' 
#' ssvQC(feature_config_file)
ssvQC.runAll = function(sqc){
  out_dir = sqc@out_dir
  res_file = function(f)file.path(out_dir, f)
  
  #bam specific independent of peaks
  if(grepl("bam", sqc@signal_config@read_mode)){
    
    bam_config_dt = sqc@signal_config@meta_data
    if(is.null(bam_config_dt$mapped_reads)){
      bam_config_dt$mapped_reads = sapply(bam_config_dt$file, get_mapped_reads)
    }
    
    color_var = sqc@signal_config@color_by
    group_var = sqc@signal_config@run_by
    color_mapping = sqc@signal_config@color_mapping
    
    theme_set(theme(panel.background = element_blank(), axis.text.x = element_text(size = 8)))
    
    p_mapped_reads = ggplot(bam_config_dt, aes_string(x = "mark", y = "mapped_reads", fill = color_var)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = color_mapping) +
      scale_y_continuous(labels = function(x)x/1e6) +
      labs(y = "M mapped reads", fill = color_var, x= "") +
      labs(title = "Mapped reads")
    
    
    ggsave(res_file("mapped_reads.pdf"), p_mapped_reads, width = 2+.5*nrow(bam_config_dt), height = 3)
  }
  
  
  
  #feature overlaps independendent of signal
  # todo = sqc@feature_config@
}

# qcAssessPeaks = function(){
#   signal_config
# }