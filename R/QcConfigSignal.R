
check_QcConfigSignal = function(object){
  errors <- character()
  
  if(grepl("bam", object@read_mode)){#bam checks
    
  }else{#bigwig checks
    if(object@cluster_value == "RPM"){
      errors = c(errors, "cluster_value cannot be RPM for bigwig files.")
    }
    if(object@sort_value == "RPM"){
      errors = c(errors, "sort_value cannot be RPM for bigwig files.")
    }
  }
  
  if (length(errors) == 0) TRUE else errors
}

' QcConfigSignal
#'
#' @slot view_size numeric.
#' @slot read_mode character.
#' @slot cluster_value 
#' @slot linearQuantile_cutoff 
#' @slot sort_value 
#' @slot sort_method 
#' @slot plot_value 
#' @slot  
#' @slot fetch_options list.
#'
#' @rdname QcConfigSignal
#' @export
#'
setClass("QcConfigSignal", contains = "QcConfig",
         representation = list(
           view_size = "numeric",
           read_mode = "character",
           fetch_options = "list",
           cluster_value = "character",
           linearQuantile_cutoff = "numeric",
           sort_value = "character",
           sort_method = "character",
           plot_value = "character"
         ), validity = check_QcConfigSignal)

valid_cluster_vals = c("raw", "RPM", "linearQuantile")

setMethod("initialize","QcConfigSignal", function(.Object,...){
  .Object <- callNextMethod()
  validObject(.Object)
  .Object
})

#' QcConfigSignal
#'
#' @param config_df 
#' @param run_by 
#' @param color_by 
#' @param color_mapping 
#' @param read_mode 
#' @param view_size 
#' @param to_run 
#' @param to_run_reference 
#' @param fetch_options 
#' @param cluster_value 
#' @param linearQuantile_cutoff 
#' @param sort_value 
#' @param sort_method 
#' @param is_null 
#'
#' @return
#' @export
#' @rdname QcConfigSignal
#' @examples
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' bam_config_df = .parse_config_body(bam_config_file)
#' sig_conf = QcConfigSignal(bam_config_df)
#' 
#' bigwig_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bigwig_config.csv")
#' bigwig_config_df = .parse_config_body(bigwig_config_file)
#' sig_conf.bw = QcConfigSignal(bigwig_config_df)
QcConfigSignal = function(config_df,
                          run_by = "All",
                          to_run = NULL,
                          to_run_reference = NULL,
                          color_by = "file",
                          color_mapping = NULL,
                          read_mode = NULL,
                          view_size = getOption("SQC_VIEW_SIZE", 3e3), 
                          fetch_options = list(),
                          cluster_value = valid_cluster_vals[1],
                          linearQuantile_cutoff = .98,
                          sort_value = valid_cluster_vals[1],
                          sort_method = c("hclust", "sort")[2],
                          plot_value = valid_cluster_vals[1],
                          is_null = FALSE){
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
  if(!cluster_value %in% valid_cluster_vals){
    stop("cluster_value of ", cluster_value, " was not one of : ", paste(valid_cluster_vals, collapse = ", "))
  }
  if(!sort_value %in% valid_cluster_vals){
    stop("sort_value of ", sort_value, " was not one of : ", paste(valid_cluster_vals, collapse = ", "))
  }
  if(linearQuantile_cutoff <= 0 | linearQuantile_cutoff > 1){
    stop("linearQuantile_cutoff must be between 0 and 1. Was ", linearQuantile_cutoff)
  }
  if(!sort_method %in% c("hclust", "sort")){
    stop("sort_method of ", sort_method, " was not one of : ", paste(c("hclust", "sort"), collapse = ", "))
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
      view_size = view_size, 
      fetch_options = fetch_options,
      cluster_value = cluster_value,
      linearQuantile_cutoff = linearQuantile_cutoff,
      sort_value = sort_value,
      sort_method = sort_method,
      plot_value = plot_value,
      is_null = is_null)
}

QcConfigSignal.null = function(){
  qc = suppressWarnings({QcConfigSignal(data.frame(file = "null", stringsAsFactors = FALSE))})
  qc@is_null = TRUE
  qc
}

#' @param signal_config_file Configuration file for signal data.
#'
#' @return
#' @export
#' @rdname QcConfigSignal
#' @examples
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' QcConfigSignal.parse(bam_config_file)
#' 
#' bigwig_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bigwig_config.csv")
#' QcConfigSignal.parse(bigwig_config_file)
QcConfigSignal.parse = function(signal_config_file){
  signal_config_dt = .parse_config_body(signal_config_file)
  valid_signal_var = c("main_dir", 
                       "view_size", 
                       "read_mode", 
                       "color_by", 
                       "color_mapping", 
                       "run_by", 
                       "to_run", 
                       "to_run_reference", 
                       "fetch_options", 
                       "is_null",
                       "cluster_value",
                       "linearQuantile_cutoff",
                       "sort_value",
                       "plot_value",
                       "sort_method")
  cfg_vals = .parse_config_header(signal_config_file, valid_signal_var)
  
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
                  run_by = NULL, 
                  to_run = NULL, 
                  to_run_reference = NULL,
                  cluster_value = valid_cluster_vals[1],
                  linearQuantile_cutoff = .98,
                  sort_value = valid_cluster_vals[1],
                  sort_method = c("hclust", "sort")[2],
                  plot_value = valid_cluster_vals[1],
                  fetch_options = list(), 
                  is_null = FALSE){
    QcConfigSignal(config_df = config_dt, 
                   run_by = run_by, 
                   to_run = to_run,
                   to_run_reference = to_run_reference,
                   color_by = color_by, 
                   color_mapping = color_mapping, 
                   read_mode = read_mode, 
                   view_size = view_size, 
                   fetch_options = fetch_options,
                   cluster_value = cluster_value,
                   linearQuantile_cutoff = linearQuantile_cutoff,
                   sort_value = sort_value,
                   sort_method = sort_method,
                   plot_value = plot_value,
                   is_null = TRUE)
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
#' @rdname QcConfigSignal
#' @examples
#' bam_files = dir(system.file(package = "ssvQC", "extdata"), pattern = "CTCF.+bam$", full.names = TRUE)
#' QcConfigSignal.files(bam_files)
QcConfigSignal.files = function(file_paths,
                                groups = NULL,
                                group_names = NULL,
                                group_colors = NULL,
                                view_size = getOption("SQC_VIEW_SIZE", 3e3), 
                                read_mode = NULL,
                                cluster_value = valid_cluster_vals[1],
                                linearQuantile_cutoff = .98,
                                sort_value = valid_cluster_vals[1],
                                sort_method = c("hclust", "sort")[2]
){
  if(is.null(groups)){
    groups = seq_along(file_paths)
  }
  if(is.null(group_names)){
    if(!any(duplicated(basename(file_paths)))){
      group_names = basename(file_paths)
    }else{
      group_names = LETTERS[seq_along(unique(groups))]  
      message("non-unique file basenames, resorting to A,B,C... names.")
    }
  }
  if(is.null(group_colors)){
    group_colors = seqsetvis::safeBrew(length(group_names))
  }
  if(is.null(names(group_colors))){
    names(group_colors) = group_names
  }
  
  config_df = data.frame(file = as.character(file_paths), group = group_names[groups], stringsAsFactors = FALSE)
  
  config_df$name = basename(config_df$file)
  config_df$name_split = gsub("[_\\.]", "\n", config_df$name)
  
  config_df$name = factor(config_df$name, levels = unique(config_df$name))
  config_df$name_split = factor(config_df$name_split, levels = unique(config_df$name_split))
  
  config_df$All = "All"
  
  QcConfigSignal(config_df, run_by = "All", color_by = "group", color_mapping = group_colors,
                 cluster_value = cluster_value,
                 linearQuantile_cutoff = linearQuantile_cutoff,
                 sort_value = sort_value,
                 sort_method = sort_method)
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

#' @param qc_signal 
#' @param query_gr 
#'
#' @return
#' @export
#' @rdname QcConfigSignal
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
  if(!is.null(qc_signal@meta_data$fragLens)){
    if(!is.null(extra_args$fragLens)){
      if(extra_args$fragLens != "auto"){
        warning("Overwriting configured fragLens with detected fragLens for fetch call.")  
      }
    }
    extra_args$fragLens = qc_signal@meta_data$fragLens
  }
  call_args = c(list(file_paths = qc_signal@meta_data, qgr = query_gr, return_data.table = TRUE), extra_args)
  fetch_FUN = get_fetch_fun(qc_signal@read_mode)
  prof_dt = do.call(fetch_FUN, call_args)
  prof_dt
}

setMethod("split", signature = c("QcConfigSignal"), definition = function(x){
  f = x@meta_data[[x@run_by]]
  split(x, f, FALSE)
})
setMethod("split", signature = c("QcConfigSignal", "factor"), definition = function(x, f){
  split(x, f, FALSE)
})
setMethod("split", signature = c("QcConfigSignal", "factor", "logical"), definition = function(x, f, drop = FALSE){
  config_df = x@meta_data
  meta_split = split(config_df, f)
  refs = config_df[config_df[[x@run_by]] %in% x@to_run_reference,]
  
  meta_split = meta_split[x@to_run]
  meta_split = lapply(meta_split, function(sel_meta_df)rbind(sel_meta_df, refs))
  
  message("split QcConfigSignal")
  
  lapply(meta_split, function(sel_meta_df){
    new("QcConfigSignal",
        meta_data =  sel_meta_df,
        run_by = x@run_by,
        to_run = x@to_run,
        to_run_reference = x@to_run_reference,
        color_by = x@color_by,
        color_mapping = x@color_mapping,
        read_mode = x@read_mode,
        view_size = x@view_size, 
        cluster_value = x@cluster_value,
        linearQuantile_cutoff = x@linearQuantile_cutoff,
        sort_value = x@sort_value,
        sort_method = x@sort_method,
        plot_value = x@plot_value,
        fetch_options = x@fetch_options)
  })
})
setMethod("split", signature = c("QcConfigSignal", "character"), definition = function(x, f){
  split(x, f, FALSE)
})
setMethod("split", signature = c("QcConfigSignal", "character", "logical"), definition = function(x, f, drop = FALSE){
  split(x, factor(f, levels = unique(f)), FALSE)
})

#' @return
#' @export
#' @rdname QcConfigSignal
#' @examples
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' bam_config = QcConfigSignal.parse(bam_config_file)
#' #QcConfigSignal.save_config(bam_config, "bam_config.csv")
#'
#' bigwig_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bigwig_config.csv")
#' bigwig_config = QcConfigSignal.parse(bigwig_config_file)
#' #QcConfigSignal.save_config(bigwig_config, "bigwig_config.csv")
QcConfigSignal.save_config = function(object, file){
  slots_to_save = c(
    "view_size",
    "read_mode",
    "run_by",
    "to_run",
    "to_run_reference",
    "color_by",
    "is_null",
    "cluster_value",
    "linearQuantile_cutoff",
    "sort_value",
    "sort_method"
  )
  kvp_slots = c("color_mapping", "fetch_options")
  # QcConfigSignal.parse(file)
  .save_config(object, file, slots_to_save, kvp_slots)
}
