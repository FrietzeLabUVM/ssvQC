
setClass("ssvQC",
         representation = list(
           feature_config = "QcConfigFeatures",
           signal_config = "QcConfigSignal",
           out_dir = "character",
           bfc = "BiocFileCache",
           saving_enabled = "logical"
         ))

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
#' ssvQC(feature_config_file, bam_config_file)
#' 
#' ssvQC(feature_config, bam_config)
#' 
#' ssvQC(feature_config, bigwig_config_file)
#' 
ssvQC = function(feature_config,
                 signal_config,
                 out_dir = getwd(),
                 bfc = NULL){

  if(is.character(feature_config)){
    if(file.exists(feature_config)){
      feature_config = QcConfigFeatures.parse(feature_config)
    }
  }
  if(!"QcConfigFeatures" %in% class(feature_config)){
    stop("feature_config must be either a QcConfigFeatures object or the path to valid configuration file to create one.")
  }
  if(is.character(signal_config)){
    if(file.exists(signal_config)){
      signal_config = QcConfigSignal.parse(signal_config)
    }
  }
  if(!"QcConfigSignal" %in% class(signal_config)){
    stop("signal_config must be either a QcConfigSignal object or the path to valid configuration file to create one.")
  }
  if(is.null(bfc)){
    bfc = BiocFileCache::BiocFileCache()
  }
  
  stopifnot(file.exists(signal_config@meta_data$file))
  stopifnot(file.exists(signal_config@meta_data$file))
    
  dir.create(out_dir, showWarnings = FALSE)
  new("ssvQC",
      feature_config = feature_config,
      signal_config = signal_config,
      out_dir = out_dir,
      bfc = bfc 
  )
}

#' guess_feature_file_format
#'
#' @param feature_files 
#'
#' @return
#' @export
#'
#' @examples
guess_feature_file_format = function(feature_files){
  .guess_feature_file_format = function(feature_file){
    file_format = "unknown"
    if(grepl("narrowPeak$", feature_file)){
      file_format = "narrowPeak"
    }else if(grepl("broadPeak$", feature_file)){
      file_format = "broadPeak"
    }else if(grepl("bed$", feature_file)){
      file_format = "bed"
    }else{
      warning("Could not guess file format for feature file: ", feature_file)
    }
    file_format
  }
  
  sapply(feature_files, .guess_feature_file_format)
}

#' Title
#'
#' @param signal_file 
#'
#' @return
#' @export
#'
#' @examples
guess_read_mode = function(signal_file){
  if(grepl(".bam$", signal_file[1])){
    mode = "bam_SE"
  }else{
    mode = "bigwig"
  }
  message("read_mode has been guessed as ", mode)
  if(mode == "bam_SE"){
    message("Currently ssvQC cannot guess whether a bam file is SE or PE.  Please manually specify bam_PE if appropriate.")
  }
  mode
}

#' Title
#'
#' @param feature_files 
#'
#' @return
#' @export
#'
#' @examples
get_feature_file_load_function = function(feature_files){
  file_types = guess_feature_file_format(feature_files)
  .get_feature_file_load_function = function(file_type){
    switch (file_type,
            narrowPeak = seqsetvis::easyLoad_narrowPeak,
            broadPeak = seqsetvis::easyLoad_broadPeak,
            seacr = seqsetvis::easyLoad_seacr,
            bed = seqsetvis::easyLoad_bed,
            unknown = {
              warning("Treating unknown file type as bed but if you see errors, check file type.")
              seqsetvis::easyLoad_bed
            },
            stop("'", file_type, "' is not a supported file_type")
    )
  }
  
  sapply(file_types, .get_feature_file_load_function)
  

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
qcBasicMetrics = function(sqc){
  out_dir = sqc@out_dir
  res_file = function(f)file.path(out_dir, f)
  
  if(grepl("bam", sqc@signal_config@read_mode)){

    bam_config_dt = sqc@signal_config@meta_data
    if(is.null(bam_config_dt$mapped_reads)){
      bam_config_dt$mapped_reads = sapply(bam_config_dt$file, get_mapped_reads)
    }
        
    color_var = sqc@signal_config@color_var
    group_var = sqc@signal_config@group_var
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
}

qcAssessPeaks = function(){
  signal_config
}