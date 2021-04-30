
' QcConfigSignal
#'
#' @slot view_size numeric.
#'
#' @export
#'
setClass("QcConfigSignal", contains = "QcConfig",
         representation = list(
           view_size = "numeric",
           read_mode = "character",
           fetch_options = "list"
         ))


#' QcConfigSignal
#'
#' @param config_df 
#' @param run_by 
#' @param color_by 
#' @param color_mapping 
#' @param read_mode 
#' @param view_size 
#'
#' @return
#' @export
#'
#' @examples
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' bam_config_df = .parse_config_body(bam_config_file)
#' QcConfigSignal(bam_config_df)
#' 
#' bigwig_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bigwig_config.csv")
#' bigwig_config_df = .parse_config_body(bigwig_config_file)
#' QcConfigSignal(bigwig_config_df)
QcConfigSignal = function(config_df,
                          run_by = "All",
                          to_run = NULL,
                          to_run_reference = NULL,
                          color_by = "file",
                          color_mapping = NULL,
                          read_mode = NULL,
                          view_size = getOption("SQC_VIEW_SIZE", 3e3)){
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
  
  if(is.null(read_mode)){
    read_mode = guess_read_mode(config_df$file[1])
  }
  
  stopifnot(read_mode %in% c("bam_SE", "bam_PE", "bigwig", "null"))
  
  if(is.null(to_run)){
    to_run = unique(config_df[[run_by]])   
  }
  stopifnot(all(to_run %in% config_df[[run_by]]))
  if(is.null(to_run_reference)){
    to_run_reference = character()
  }
  stopifnot(all(to_run_reference %in% config_df[[run_by]]))
  
  new("QcConfigSignal",
      meta_data =  config_df,
      run_by = run_by,
      to_run = to_run,
      to_run_reference = to_run_reference,
      color_by = color_by,
      color_mapping = color_mapping,
      read_mode = read_mode,
      view_size = view_size)
}

QcConfigSignal.null = function(){
  qc = suppressWarnings({QcConfigSignal(data.frame(file = "null", stringsAsFactors = FALSE))})
  qc@is_null = TRUE
  qc
}

#' Title
#'
#' @param signal_config_file 
#'
#' @return
#' @export
#'
#' @examples
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' QcConfigSignal.parse(bam_config_file)
#' 
#' bigwig_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bigwig_config.csv")
#' QcConfigSignal.parse(bigwig_config_file)
QcConfigSignal.parse = function(signal_config_file){
  signal_config_dt = .parse_config_body(signal_config_file)
  valid_feature_var = c("main_dir", "view_size", "read_mode", 
                        "color_by", "color_mapping", "run_by", "to_run", "to_run_reference", "fetch_options")
  cfg_vals = .parse_config_header(signal_config_file, valid_feature_var)
  
  if(!is.null(cfg_vals[["main_dir"]])){
    choose_file_path = function(main_dir, files){
      abs_path = files
      rel_path = file.path(main_dir, files)
      ifelse(file.exists(rel_path), rel_path, abs_path)
    }
    signal_config_dt[, file := choose_file_path(cfg_vals[["main_dir"]], file)]
    cfg_vals[["main_dir"]] = NULL
  }
  if(!all(file.exists(signal_config_dt$file))){
    stop(paste(c("Files specified in config do not exist:", 
                 signal_config_dt$file[!file.exists(signal_config_dt$file)]), collapse = "\n  "))
  }
  
  tfun = function(config_dt, 
                  main_dir = NULL, 
                  read_mode = NULL,
                  view_size = getOption("SQC_VIEW_SIZE", 3e3),
                  color_by = NULL, color_mapping = NULL, 
                  run_by = NULL, to_run = NULL, to_run_reference = NULL){
    QcConfigSignal(config_df = config_dt, 
                   run_by = run_by, 
                   to_run = to_run,
                   to_run_reference = to_run_reference,
                   color_by = color_by, 
                   color_mapping = color_mapping, 
                   read_mode = read_mode, 
                   view_size = view_size)
  }
  do.call(tfun, c(list(config_dt = signal_config_dt), cfg_vals))
}

#
#' QcConfigSignal.files
#'
#' @param file_paths character paths to files
#' @param groups numeric vector of group assignments. 1 is first item in group_names, 2 is second, etc. Default is seq_along(file_path)
#' @param group_names vector of group names to assign from according to groups
#' @param group_colors vector of colors to use per group
#' @param view_size view size in bp to apply. Defaults to 3000.
#'
#' @return a QcConfigSignal object
#' @export
#'
#' @examples
#' bam_files = dir(system.file(package = "ssvQC", "extdata"), pattern = "CTCF.+bam$", full.names = TRUE)
#' QcConfigSignal.files(bam_files)
QcConfigSignal.files = function(file_paths,
                                groups = NULL,
                                group_names = NULL,
                                group_colors = NULL,
                                view_size = getOption("SQC_VIEW_SIZE", 3e3), 
                                read_mode = NULL){
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
  
  config_df = data.frame(file = as.character(file_paths), group = group_names[groups], stringsAsFactors = FALSE)
  
  QcConfigSignal(config_df, run_by = "group", color_mapping = group_colors)
}

get_fetch_fun = function(read_mode){
  stopifnot(read_mode %in% c("bam_SE", "bam_PE", "bigwig"))
  switch(read_mode,
         bam_SE = {
           seqsetvis::ssvFetchBam
         }, 
         bam_PE = {
           seqsetvis::ssvFetchBamPE
         }, 
         bigwig = {
           seqsetvis::ssvFetchBigwig
         })
}

#' Title
#'
#' @param qc_signal 
#' @param query_gr 
#'
#' @return
#' @export
#'
#' @examples
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' qc_signal = QcConfigSignal.parse(bam_config_file)
#' 
#' feature_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
#' qc_features = QcConfigFeatures.parse(feature_config_file)
#' query_gr = qc_features$assessment_features
#' fetch_signal_at_features(qc_signal, query_gr)
fetch_signal_at_features = function(qc_signal, query_gr){
  extra_args = qc_signal@fetch_options
  call_args = c(list(file_paths = qc_signal@meta_data, qgr = query_gr, return_data.table = TRUE, fragLen = 180), extra_args)
  fetch_FUN = get_fetch_fun(qc_signal@read_mode)
  do.call(fetch_FUN, call_args)
  
}
