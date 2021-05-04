
#' QcConfigFeatures
#'
#' @slot feature_load_FUN function.
#' @slot n_peaks numeric.
#' @slot consensus_fraction numeric.
#' @slot consensus_n numeric.
#'
#' @export
setClass("QcConfigFeatures", contains = "QcConfig",
         representation = list(
           feature_load_FUN = "function",
           n_peaks = "numeric",
           consensus_fraction = "numeric",
           consensus_n = "numeric",
           loaded_features = "list",
           overlap_gr = "GRanges",
           overlap_extension = "numeric",
           assessment_gr = "GRanges"
           
         ))

setMethod("initialize","QcConfigFeatures", function(.Object,...){
  .Object <- callNextMethod()
  validObject(.Object)
  .Object@loaded_features = list()
  .Object@overlap_gr = GRanges()
  .Object@assessment_gr = GRanges()
  
  .Object
})

setMethod("names", "QcConfigFeatures",
          function(x)
          {
            c("loaded_features", "overlapped_features", "assessment_features", "meta_data", "run_by", "to_run", "to_run_reference", "color_by", "color_mapping")
            
          })


setMethod("$", "QcConfigFeatures",
          function(x, name)
          {
            switch (name,
                    loaded_features = x@loaded_features,
                    overlapped_features = x@overlap_gr,
                    assessment_features = x@assessment_gr,
                    meta_data = x@meta_data,
                    run_by = x@run_by,
                    to_run = x@to_run,
                    to_run_reference = x@to_run_reference,
                    color_by = x@color_by,
                    color_mapping = x@color_mapping
            )
          })

setReplaceMethod("$", "QcConfigFeatures",
                 function(x, name, value)
                 {
                   warn_msg = "This assignment is not supported.  No effect."
                   switch (name,
                           loaded_features = warning(warn_msg),
                           overlapped_features = warning(warn_msg),
                           assessment_features = warning(warn_msg),
                           meta_data = warning(warn_msg),
                           #TODO add checks, call validity?
                           run_by = {
                             x@run_by = value
                           },
                           to_run = {
                             x@to_run = value
                           },
                           to_run_reference = {
                             x@to_run_reference = value
                           },
                           color_by = {
                             x@color_by = value
                             message("Applying option SQC_COLORS for updated color_mapping.")
                             x$color_mapping = getOption("SQC_COLORS")
                           },
                           color_mapping = {
                             col_lev = unique(x@meta_data[[x@color_by]])
                             if(is.null(names(value))){
                               if(length(value) >= length(col_lev)){
                                 value = value[seq_along(col_lev)]
                               }else{
                                 stop("Insufficient colors supplied. Mapping requires at least ", length(col_lev), ".")
                               }
                               names(value) = col_lev
                             }else{
                               if(!all(names(value) %in% col_lev)){
                                 stop(paste(collapse = "\n",
                                            c("Missing name values from color mapping. Required:", 
                                              setdiff(col_lev, names(value)))))
                               }
                             }
                             x@color_mapping = value
                           },
                           {warning(warn_msg)}
                   )
                   x
                 })


setGeneric("featuresList", function(object){standardGeneric("featuresList")})
#' Title
#'
#' @param QcConfigFeatures 
#'
#' @return
#' @export
#'
#' @examples
setMethod("featuresList", "QcConfigFeatures", 
          definition = function(object){
            if(length(object@loaded_features) == 0){
              object@loaded_features = object@feature_load_FUN(object@meta_data$file)
              names(object@loaded_features) = object@meta_data$name_split
            }
            if(length(object@loaded_features) == 0){
              stop("Somehow, no files were loaded. Report this issue at https://github.com/FrietzeLabUVM/ssvQC/issues")   
            }
            object
          })

setGeneric("featuresOverlaps", function(object){standardGeneric("featuresOverlaps")})
#' Title
#'
#' @param QcConfigFeatures 
#'
#' @return
#' @export
#'
#' @examples
setMethod("featuresOverlaps", "QcConfigFeatures", 
          definition = function(object){
            object = featuresList(object)
            if(length(object@overlap_gr) == 0){
              object@overlap_gr = seqsetvis::ssvOverlapIntervalSets(object@loaded_features)
              object@overlap_gr = seqsetvis::prepare_fetch_GRanges_names(object@overlap_gr)
            }
            if(length(object@overlap_gr) == 0){
              stop("No regions in overlap. Check your input or report this issue at https://github.com/FrietzeLabUVM/ssvQC/issues")   
            }
            object
          })

setGeneric("featuresQuery", function(object){standardGeneric("featuresQuery")})
#' Title
#'
#' @param QcConfigFeatures 
#'
#' @return
#' @export
#'
#' @examples
setMethod("featuresQuery", "QcConfigFeatures", 
          definition = function(object){
            # n_peaks = "numeric",
            # consensus_fraction = "numeric",
            # consensus_n = "numeric",
            if(length(object@assessment_gr) == 0){
              object = featuresList(object)
              object = featuresOverlaps(object)
              
              feat_list = object@loaded_features
              olap_gr = object@overlap_gr
              
              f_consensus = floor(object@consensus_fraction * length(feat_list))
              n_consensus = min(length(feat_list), max(object@consensus_n, f_consensus))
              
              if(n_consensus == 1){
                asses_gr.full = olap_gr
              }else{
                asses_gr.full = ssvConsensusIntervalSets(feat_list, min_number = object@consensus_n, min_fraction = object@consensus_fraction)
              }
              object@assessment_gr = sampleCap(asses_gr.full, object@n_peaks)
            }
            if(length(object@assessment_gr) == 0){
              stop("No regions in assessment. Maybe loosen consensus requirements or report this issue at https://github.com/FrietzeLabUVM/ssvQC/issues")   
            }
            object
          })

#' QcConfigFeatures
#'
#' @param config_df 
#' @param run_by 
#' @param color_by 
#' @param color_mapping 
#' @param feature_load_FUN 
#' @param n_peaks 
#' @param consensus_fraction 
#' @param consensus_n 
#'
#' @return
#' @export
#'
#' @examples
#' feature_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
#' config_df = .parse_config_body(feature_config_file)
#' QcConfigFeatures(config_df)
QcConfigFeatures = function(config_df,
                            run_by = "All",
                            to_run = NULL,
                            to_run_reference = NULL,
                            color_by = "file",
                            color_mapping = NULL,
                            feature_load_FUN = NULL,
                            n_peaks = 1e3,
                            overlap_extension = 0,
                            consensus_fraction = getOption("SQC_CONSENSUS_FRACTION", 0),
                            consensus_n = getOption("SQC_CONSENSUS_N", 1),
                            process_features = getOption("SQC_PROCESS_FEATURES", TRUE)){
  .enforce_file_var(config_df)
  if(!run_by %in% colnames(config_df)){
    if(run_by == "All"){
      config_df[[run_by]] = run_by
    }else{
      stop("run_by ", run_by, " was not in column names.")    
    }
  }
  if(!color_by %in% colnames(config_df)){
    stop("color_by ", color_by, " was not in column names.")
  }
  
  if(!is.null(color_mapping)){
    if(!is.null(names(color_mapping))){
      color_names = names(color_mapping)
    }else if(is.factor(config_df[[color_by]])){
      color_names = levels(config_df[[color_by]])
    }else{
      color_names = unique(config_df[[color_by]])
    }
    stopifnot(length(color_names) == length(color_mapping))
    names(color_mapping) = color_names
    
  }else{
    if(is.factor(config_df[[color_by]])){
      color_names = levels(config_df[[color_by]])
    }else{
      color_names = unique(config_df[[color_by]])
    }
    color_mapping = seqsetvis::safeBrew(length(color_names))
    names(color_mapping) = color_names
  }
  
  stopifnot(config_df[[color_by]] %in% names(color_mapping))
  
  if(is.null(feature_load_FUN)){
    file_paths = config_df[["file"]]
    file_formats = guess_feature_file_format(file_paths)
    stopifnot(length(unique(file_formats)) == 1)
    feature_load_FUN = get_feature_file_load_function(file_paths[1])[[1]]
  }
  if(is.null(to_run)){
    to_run = unique(config_df[[run_by]])   
  }
  stopifnot(all(to_run %in% config_df[[run_by]]))
  if(is.null(to_run_reference)){
    to_run_reference = character()
  }
  stopifnot(all(to_run_reference %in% config_df[[run_by]]))
  
  obj = new("QcConfigFeatures",
            meta_data =  config_df,
            run_by = run_by,
            to_run = to_run,
            to_run_reference = to_run_reference,
            color_by = color_by,
            color_mapping = color_mapping,
            feature_load_FUN = feature_load_FUN,
            n_peaks = n_peaks,
            overlap_extension = overlap_extension,
            consensus_fraction = consensus_fraction,
            consensus_n = consensus_n)
  if(process_features){
    obj = featuresQuery(obj)
  }
  obj
  
}

QcConfigFeatures.null = function(){
  qc = suppressWarnings({QcConfigFeatures(data.frame(file = "null"), process_features = FALSE)})
  qc@is_null = TRUE
  qc
}


#' QcConfigFeatures
#'
#' @param file_paths character paths to files
#' @param groups numeric vector of group assignments. 1 is first item in group_names, 2 is second, etc. Default is seq_along(file_path)
#' @param group_names vector of group names to assign from according to groups
#' @param group_colors vector of colors to use per group
#' @param n_peaks number of peaks to subset for
#' @param consensus_n An integer number specifying the absloute minimum of input grs that must overlap for a site to be considered consensus.
#' @param consensus_fraction A numeric between 0 and 1 specifying the fraction of grs that must overlap to be considered consensus.
#'
#' @return a QcConfigFeatures object
#' @export
#'
#' @examples
#' np_files = dir(system.file(package = "ssvQC", "extdata"), pattern = "Peak$", full.names = TRUE)
#' QcConfigFeatures.files(np_files)
QcConfigFeatures.files = function(file_paths,
                                  group_names = NULL,
                                  groups = NULL,
                                  group_colors = NULL,
                                  feature_load_FUN = NULL,
                                  n_peaks = 1e3,
                                  consensus_fraction = getOption("SQC_CONSENSUS_FRACTION", 0),
                                  consensus_n = getOption("SQC_CONSENSUS_N", 1),
                                  process_features = getOption("SQC_PROCESS_FEATURES", TRUE)){
  if(is.null(groups)){
    groups = seq_along(file_paths)
  }
  if(is.null(group_names)){
    group_names = LETTERS[seq_along(unique(groups))]
  }
  if(is.null(group_colors)){
    group_colors = seqsetvis::safeBrew(length(group_names))
  }
  if(is.null(names(group_colors))){
    names(group_colors) = group_names
  }
  if(is.null(feature_load_FUN)){
    feature_load_FUN = get_feature_file_load_function(file_paths[1])[[1]]
  }
  config_df = data.frame(file = as.character(file_paths), group = group_names[groups])
  
  obj = new("QcConfigFeatures",
            meta_data =  config_df,
            run_by = "group",
            color_by = "group",
            color_mapping = group_colors,
            feature_load_FUN = feature_load_FUN,
            n_peaks = n_peaks,
            consensus_fraction = consensus_fraction,
            consensus_n = consensus_n
  )
  if(process_features){
    obj = featuresQuery(obj)
  }
  obj
}


#' QcConfigFeatures.parse
#'
#' @param feature_config_file 
#'
#' @return
#' @export
#'
#' @examples
#' feature_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
#' QcConfigFeatures.parse(feature_config_file)
QcConfigFeatures.parse = function(feature_config_file,
                                  process_features = getOption("SQC_PROCESS_FEATURES", TRUE)){
  peak_config_dt = .parse_config_body(feature_config_file)
  valid_feature_var = c("main_dir", "overlap_extension", "n_peaks", "consensus_n", 
                        "consensus_fraction", "color_by", "color_mapping", "run_by", "to_run")
  cfg_vals = .parse_config_header(feature_config_file, valid_feature_var)
  
  if(!is.null(cfg_vals[["main_dir"]])){
    choose_file_path = function(main_dir, files){
      abs_path = files
      rel_path = file.path(main_dir, files)
      ifelse(file.exists(rel_path), rel_path, abs_path)
    }
    peak_config_dt[, file := choose_file_path(cfg_vals[["main_dir"]], file)]
    cfg_vals[["main_dir"]] = NULL
  }
  if(!all(file.exists(peak_config_dt$file))){
    stop(paste(c("Files specified in config do not exist:", 
                 peak_config_dt$file[!file.exists(peak_config_dt$file)]), collapse = "\n  "))
  }
  
  tfun = function(config_dt, 
                  main_dir = NULL, 
                  n_peaks = 1e3, 
                  overlap_extension = 0,
                  consensus_n = 1, consensus_fraction = 0, 
                  color_by = NULL, color_mapping = NULL, 
                  run_by = NULL, to_run = NULL, to_run_reference = NULL){
    QcConfigFeatures(config_df = config_dt, 
                     run_by = run_by, 
                     to_run = to_run,
                     to_run_reference = to_run_reference,
                     color_by = color_by, 
                     color_mapping = color_mapping, 
                     n_peaks = n_peaks, 
                     overlap_extension = overlap_extension,
                     consensus_fraction = consensus_fraction,
                     consensus_n = consensus_n,
                     process_features = process_features)
  }
  do.call(tfun, c(list(config_dt = peak_config_dt), cfg_vals))
  
}