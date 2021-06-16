setClass("ssvTSNE", contains = "ssvQC.complete")

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
#' options(mc.cores = 10)
#' set.seed(0)
#' features_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
#' features_config = QcConfigFeatures.parse(features_config_file)
#'
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' bam_config = QcConfigSignal.parse(bam_config_file)
#' 
#' sts = ssvTSNE(features_config, bam_config)
#' sts = ssvQC.plotFeatures(sts)
#' sts = ssvQC.prepSignal(sts)
#' sts$plots$signal
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
  object = ssvQC.plotFeatures(object)
  message("NYI")
  object
})

setMethod("ssvQC.prepSignal", "ssvTSNE", function(object){
  if(length(object@features_config$assessment_features) == 0){
    object = ssvQC.prepFeatures(object)
  }
  if(is.null(object$signal_config$meta_data$mapped_reads) & object$signal_config$read_mode != "bigwig"){
    object = ssvQC.prepMappedReads(object)
  }
  if(is.null(object$signal_config$meta_data$fragLens) & object$signal_config$read_mode == "bam_SE"){
    object = ssvQC.prepFragLens(object)
  }
  if(is.null(object$signal_config$meta_data$cap_value)){
    object = ssvQC.prepCapValue(object)
  }
  
  object@signal_data = lapply(object@features_config$assessment_features, function(query_gr){
    sig_configs = .make_query_signal_config(object@signal_config)
    lapply(sig_configs, function(sel_sig_config){
      ClusteredSignal.fromConfig(sel_sig_config, 
                                 resize(query_gr, object@signal_config@view_size, fix = "center"), 
                                 facet_var = "name_split", 
                                 extra_var = union(object@signal_config@color_by, object@signal_config@run_by), 
                                 bfc = object@bfc)
    })
  })
  object
})

setMethod("ssvQC.plotSignal", "ssvTSNE", function(object){
  if(length(object@signal_data) == 0){
    object = ssvQC.prepSignal(object)
  }
  signal_data = object@signal_data
  clust_sig = signal_data[[1]][[1]]
  sig_config = object@signal_config
})