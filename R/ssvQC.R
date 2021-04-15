
setOldClass("BiocFileCache")
setClass("ssvQC",
         representation = list(
           feature_config = "QcConfigFeatures",
           signal_config = "QcConfigSignal",
           out_dir = "character",
           bfc = "BiocFileCache"
         ))

ssvQC = function(feature_config,
                 signal_config,
                 out_dir = getwd(),
                 bfc = NULL){

  if(is.character(feature_config)){
    
  }
  
  new("ssvQC",
      file_paths =  as.character(file_paths),
      groups = groups,
      group_names = group_names,
      group_colors = group_colors
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



.default_feature_config = function(feature_files){
  QcConfigFeatures(feature_files)
}
