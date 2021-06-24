setClass("ssvTSNE", 
         representation = list(
           perplexity = "numeric"
         ),
         contains = "ssvQC.complete")

#' ssvTSNE
#'
#' @param features_config 
#' @param signal_config 
#' @param out_dir 
#' @param bfc 
#' @param matched_only 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
#' options(mc.cores = 1)
#' set.seed(0)
#' features_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
#' features_config = QcConfigFeatures.parse(features_config_file)
#'
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' bam_config = QcConfigSignal.parse(bam_config_file)
#' bam_config$view_size = 600
#' 
#' sts = ssvTSNE(features_config, bam_config)
#' sts = ssvQC.plotFeatures(sts)
#' sts@perplexity = 10
#' sts = ssvQC.prepSignal(sts)
#' sts = ssvQC.plotSignal(sts)
#' sts$plots$TSNE
ssvTSNE = function(features_config = NULL,
                   signal_config = NULL,
                   out_dir = getwd(),
                   bfc = NULL, 
                   matched_only = TRUE,
                   ...){
  if(is.null(features_config) & is.null(signal_config)){
    stop("At least one of features_config or signal_config must be specified.")
  }
  
  features_config = .prep_features_config(features_config)
  signal_config = .prep_signal_config(signal_config)
  
  if(is.null(bfc)){
    bfc = new_cache()
  }
  
  dir.create(out_dir, showWarnings = FALSE)
  
  new("ssvTSNE",
      features_config = features_config,
      signal_config = signal_config,
      signal_data = list(),
      other_data = list(),
      out_dir = out_dir,
      bfc = bfc,
      saving_enabled = TRUE,
      matched_only = matched_only)
}

setMethod("ssvQC.runAll", "ssvTSNE", function(object){
  object = callNextMethod()
  message("NYI")
  object
})

setMethod("ssvQC.prepSignal", "ssvTSNE", function(object){
  object = callNextMethod()
  object@signal_data = lapply(object@signal_data, function(signal_clust_objs){
    lapply(signal_clust_objs, function(clust_obj){
      ClusteredSignal_TSNE.from_ClusteredSignal(clust_obj, object)
    })
  })
  object
})

setMethod("ssvQC.plotSignal", "ssvTSNE", function(object){
  object = callNextMethod()
  browser()
  sts@plots$TSNE = lapply(object@signal_data, function(signal_data_groups){
    lapply(signal_data_groups, function(signal_data){
      tmp = object@signal_data$CTCF_features$CTCF_signal
      if(!"ClusteredSignal_TSNE" %in% class(tmp)){
        stop("Stored signal data is not of class ClusteredSignal_TSNE.")
      }
      tmp@xy_data
      p_summary_profiles = stsPlotSummaryProfiles(tmp@signal_data, tmp@xy_data, 10)
      p_summary_profiles
    })
  })
  
  
  object
})
