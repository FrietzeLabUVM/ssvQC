
#' QcConfigFeatures
#'
#' @slot feature_load_FUN function.
#' @slot n_peaks numeric.
#' @slot consensus_fraction numeric.
#' @slot consensus_n numeric.
#' @rdname QcConfigFeatures
#' @export
setClass("QcConfigFeatures", contains = "QcConfig",
         representation = list(
           feature_load_FUN = "function",
           n_peaks = "numeric",
           balance_groups = "logical",
           consensus_fraction = "numeric",
           consensus_n = "numeric",
           loaded_features = "list",
           overlap_gr = "list",
           overlap_extension = "numeric",
           assessment_gr = "list"
           
         ))

setMethod("initialize","QcConfigFeatures", function(.Object,...){
  .Object <- callNextMethod()
  validObject(.Object)
  .Object@loaded_features = list()
  .Object@overlap_gr = list()
  .Object@assessment_gr = list()
  
  .Object
})

#' @export
setMethod("plot", "QcConfigFeatures", definition = function(x).plot_QcConfig(x))

setMethod("names", "QcConfigFeatures",
          function(x)
          {
            c("loaded_features", "overlapped_features", "assessment_features", 
              "n_peaks",
              "meta_data", "run_by", "to_run", "to_run_reference", "color_by", "color_mapping",
              "consensus_fraction", "consensus_n")
            
          })


setMethod("$", "QcConfigFeatures",
          function(x, name)
          {
            switch (name,
                    loaded_features = x@loaded_features,
                    overlapped_features = x@overlap_gr,
                    assessment_features = x@assessment_gr,
                    n_peaks = x@n_peaks,
                    consensus_fraction = x@consensus_fraction,
                    consensus_n = x@consensus_n,
                    meta_data = x@meta_data,
                    run_by = x@run_by,
                    to_run = x@to_run,
                    to_run_reference = x@to_run_reference,
                    color_by = x@color_by,
                    color_mapping = as.list(x@color_mapping)
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
                           meta_data = {
                             x@meta_data = value
                           },
                           n_peaks = {
                             x@assessment_gr = list()
                             x@n_peaks = value
                           },
                           consensus_n = {
                             x@overlap_gr = list()
                             x@assessment_gr = list()
                             x@consensus_n = value
                           },
                           consensus_fraction = {
                             x@overlap_gr = list()
                             x@assessment_gr = list()
                             x@consensus_fraction = value
                           },
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
                             value = unlist(value)
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

.process_features = function(meta_dt, feature_load_FUN, 
                             bfc = new_cache(), 
                             force_overwrite = getOption("SQC_FORCE_CACHE_OVERWRITE", FALSE)){
  bfcif(bfc, digest_args(), force_overwrite = force_overwrite,  function(){
    loaded_features = feature_load_FUN(meta_dt$file)
    names(loaded_features) = as.character(meta_dt$name_split)
    
    if(length(loaded_features) == 0){
      stop("Somehow, no files were loaded. Report this issue at https://github.com/FrietzeLabUVM/ssvQC/issues")   
    }
    loaded_features
  })
}

.process_overlaps = function(loaded_features, 
                             overlap_extension, 
                             bfc = new_cache(), 
                             force_overwrite = getOption("SQC_FORCE_CACHE_OVERWRITE", FALSE)){
  bfcif(bfc, digest_args(), force_overwrite = force_overwrite, function(){
    overlap_gr = seqsetvis::ssvOverlapIntervalSets(loaded_features, ext = overlap_extension)
    overlap_gr = sort(overlap_gr)
    overlap_gr = seqsetvis::prepare_fetch_GRanges_names(overlap_gr)
    if(length(overlap_gr) == 0){
      stop("No regions in overlap. Check your input or report this issue at https://github.com/FrietzeLabUVM/ssvQC/issues")   
    }
    overlap_gr
  })
}

.process_assessment = function(feat_list, 
                               olap_gr, 
                               overlap_extension, 
                               n_peaks, 
                               balance_groups, 
                               consensus_fraction, 
                               consensus_n, bfc = new_cache(), 
                               force_overwrite = getOption("SQC_FORCE_CACHE_OVERWRITE", FALSE)){
  bfcif(bfc, digest_args(), force_overwrite = force_overwrite, function(){
    f_consensus = floor(consensus_fraction * length(feat_list))
    n_consensus = min(length(feat_list), max(consensus_n, f_consensus))
    
    if(n_consensus == 1){
      asses_gr.full = olap_gr
    }else{
      asses_gr.full = ssvConsensusIntervalSets(feat_list, min_number = consensus_n, min_fraction = consensus_fraction, ext = overlap_extension)
    }
    if(n_peaks >= length(asses_gr.full)){
      assessment_gr = asses_gr.full
    }else{
      if(!balance_groups){
        assessment_gr = sort(sampleCap(asses_gr.full, n_peaks))
      }else{
        df = as.data.frame(mcols(asses_gr.full))
        n_peaks.per = floor(n_peaks / ncol(df))
        sel_i = sort(unique(unlist(lapply(seq_len(ncol(df)), function(i){
          x = df[,i]
          sampleCap(which(x), n_peaks.per)
        }))))
        assessment_gr = asses_gr.full[sel_i]
      }  
    }
    if(length(assessment_gr) == 0){
      stop("No regions in assessment. Maybe loosen consensus requirements or report this issue at https://github.com/FrietzeLabUVM/ssvQC/issues")   
    }
    seqsetvis::prepare_fetch_GRanges_names(assessment_gr)
  })
}

#' @param object 
#'
#' @return
#' @rdname QcConfigFeatures
#' @examples
prepFeatures = function(object, bfc = new_cache()){
  object@meta_data[[object@run_by]]
  to_run = object@to_run
  for(tr in to_run){
    tr_name = paste0(tr, "_features")
    rb = object@meta_data[[object@run_by]]
    sel_dt = object@meta_data[rb %in% union(tr, object@to_run_reference),]
    if(is.null(object@loaded_features[[tr_name]])){
      object@loaded_features[[tr_name]] = .process_features(sel_dt, object@feature_load_FUN, bfc = bfc)
    }
    
    if(is.null(object@overlap_gr[[tr_name]])){
      object@overlap_gr[[tr_name]] = .process_overlaps(
        loaded_features = object@loaded_features[[tr_name]], 
        overlap_extension = object@overlap_extension,
        bfc = bfc)
    }
    
    if(is.null(object@assessment_gr[[tr_name]])){
      object@assessment_gr[[tr_name]] = .process_assessment(
        feat_list = object@loaded_features[[tr_name]], 
        olap_gr = object@overlap_gr[[tr_name]], 
        overlap_extension = object@overlap_extension, 
        n_peaks = object@n_peaks, 
        balance_groups = object@balance_groups,
        consensus_fraction = object@consensus_fraction, 
        consensus_n = object@consensus_n, 
        bfc = bfc)
    }
  }
  object
}

#' QcConfigFeatures
#'
#' @param config_df data.frame defining configuration parameters. At a minimum,
#'   paths to valid files in either the first column or a column named "file".
#'   Additional columns defined by color_by and run_by parameters have a big
#'   impact on the configuration.
#' @param run_by character that defines the column of config_df that groups the features.  The default of "All" will simply group all features into a single comparison.
#' @param color_by character that defines the column of config_df that controls color mapping. The default of "file" will assign a unique color to every feature set.
#' @param color_mapping named character vector that maps values of color_by to valid R colors, i.e. "red" or "#FF0000". 
#' @param feature_load_FUN function that takes a vector of file paths and returns list of GRanges.
#' @param n_peaks Number of features to sample from full overlap of feature sets for use in fetching signal.
#' @param balance_groups If TRUE, will attempt to represent imbalanced groups more equally
#' @param consensus_fraction number [0,1] to adjust number of overlap required for consensus dynamically.
#' @param consensus_n number from 1 to number of feature sets to statically set threshold for consensus.
#' @param to_run group name to run
#' @param to_run_reference group name to included in all runs
#' @param overlap_extension bp to extends regions by before calculating overlaps
#' @param process_features if TRUE, features are loaded and overlapped to generate assement regions as part of creating QcConfigFeatures object
#' @param is_null if TRUE this object will be treated as NULL
#'
#' @return QcConfigFeatures object
#' @export
#' @rdname QcConfigFeatures
#' @examples
#' feature_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
#' config_df = .parse_config_body(feature_config_file)
#' config_df$file = file.path(system.file(package = "ssvQC", "extdata"), config_df$file)
#' feature_conf = QcConfigFeatures(config_df, process_features = TRUE)
QcConfigFeatures = function(config_df,
                            run_by = "All",
                            to_run = NULL,
                            to_run_reference = NULL,
                            color_by = "file",
                            color_mapping = NULL,
                            feature_load_FUN = NULL,
                            n_peaks = 1e3,
                            balance_groups = FALSE,
                            overlap_extension = 0,
                            consensus_fraction = getOption("SQC_CONSENSUS_FRACTION", 0),
                            consensus_n = getOption("SQC_CONSENSUS_N", 1),
                            process_features = getOption("SQC_PROCESS_FEATURES", TRUE),
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
            balance_groups = balance_groups,
            overlap_extension = overlap_extension,
            consensus_fraction = consensus_fraction,
            consensus_n = consensus_n,
            is_null = is_null)
  if(process_features){
    obj = prepFeatures(obj)
  }
  obj
  
}

#' QcConfigFeatures null placeholder
#'
#' @return QcConfigFeatures object
#' @export
#' @rdname QcConfigFeatures
#' @examples
#' QcConfigFeatures.null()
QcConfigFeatures.null = function(){
  qc = suppressWarnings({QcConfigFeatures(data.frame(file = "null", stringsAsFactors = FALSE), process_features = FALSE)})
  qc@is_null = TRUE
  qc
}

#' Title
#'
#' @param query_gr either a list of GRanges of a GRanges object.  GRanges object will be handled differently if its mcols are a membership table as generated by seqsetvis::sssvOverlapIntervalSets or seqsetvis::ssvConsensusIntervalSets.
#' @param n_peaks number of peaks to subset for 
#' @param run_separately If TRUE, each item is considered a separate group. Default is TRUE.
#' @param sample_names Unique names to use per GRanges in list or columns in membership table.  Do not set if single GRanges does not have memb table.
#' @param sample_names.split Unique names to use per GRanges in list or columns in membership table.  Do not set if single GRanges does not have memb table.
#' @param group_names vector of group names to assign from according to groups
#' @param group_colors vector of colors to use per group
#' @param balance_groups If TRUE, will attempt to represent imbalanced groups more equally
#' @param overlap_extension bp to extends regions by before calculating overlaps
#' @param consensus_n An integer number specifying the absloute minimum of input grs that must overlap for a site to be considered consensus.
#' @param consensus_fraction A numeric between 0 and 1 specifying the fraction of grs that must overlap to be considered consensus.
#' @param process_features if TRUE, features are loaded and overlapped to generate assement regions as part of creating QcConfigFeatures object
#'
#' @return
#' @export
#'
#' @examples
QcConfigFeatures.GRanges = function(query_gr, 
                                    n_peaks = 2000,
                                    run_separately = TRUE,
                                    sample_names = NULL,
                                    sample_names.split = NULL,
                                    group_names = NULL,
                                    group_colors = NULL,
                                    balance_groups = FALSE,
                                    overlap_extension = 0,
                                    consensus_fraction = getOption("SQC_CONSENSUS_FRACTION", 0),
                                    consensus_n = getOption("SQC_CONSENSUS_N", 1),
                                    process_features = getOption("SQC_PROCESS_FEATURES", TRUE)){
  stopifnot(is(query_gr, "list") | is(query_gr, "GRanges"))
  
  #3 scenarios, list of GRanges to overlap, already overlapped GRanges with memb table in mcols, or single GRanges
  #list of GRanges
  if(is.list(query_gr)){
    if(is.null(sample_names)){
      sample_names = names(query_gr)
    }
    if(is.null(group_names)){
      group_names = rep("query", length(query_gr))
    }
    if(is.null(group_colors)){
      group_colors = get_group_colors(group_names)
    }
    if(is.null(names(group_colors))){
      names(group_colors) = group_names
    }
    config_df = data.frame(file = as.character(seq_along(query_gr)), 
                           group = group_names, 
                           All = "All", 
                           stringsAsFactors = FALSE)
    #memb table  
  }else if(all(sapply(mcols(query_gr), is, class2 = "logical")) & ncol(mcols(query_gr)) > 0){
    if(is.null(sample_names)){
      sample_names = colnames(mcols(query_gr))
    }
    if(is.null(group_names)){
      group_names = rep("query", ncol(mcols(query_gr)))
    }
    if(is.null(group_colors)){
      group_colors = get_group_colors(group_names)
    }
    if(is.null(names(group_colors))){
      names(group_colors) = group_names
    }
    config_df = data.frame(file = as.character(seq_len(ncol(mcols(query_gr)))), 
                           group = group_names, 
                           All = "All", 
                           stringsAsFactors = FALSE)
    #just a GRanges
  }else{
    if(is.null(sample_names)){
      sample_names = "query"
    }
    if(is.null(group_names)){
      group_names = "query"
    }
    if(is.null(group_colors)){
      group_colors = get_group_colors(group_names)
    }
    if(is.null(names(group_colors))){
      names(group_colors) = group_names
    }
    config_df = data.frame(file = "query", 
                           group = group_names, 
                           All = "All", 
                           stringsAsFactors = FALSE)
  }
  
  
  
  if(is.null(sample_names)){
    config_df$name = basename(config_df$file)  
  }else{
    config_df$name = sample_names
  }
  if(is.null(sample_names.split)){
    config_df$name_split = gsub("[_\\. ]", "\n", config_df$name)
  }else{
    config_df$name_split = sample_names.split
  }
  
  config_df$name = factor(config_df$name, levels = unique(config_df$name))
  config_df$name_split = factor(config_df$name_split, levels = unique(config_df$name_split))
  
  run_by = ifelse(run_separately, "group", "All")
  obj = new("QcConfigFeatures",
            meta_data =  config_df,
            run_by = run_by,
            to_run = unique(config_df[[run_by]]),
            to_run_reference = character(),
            color_by = "group",
            color_mapping = group_colors,
            feature_load_FUN = function(...){stop("Peaks should be loaded already")},
            n_peaks = n_peaks,
            balance_groups = balance_groups,
            overlap_extension = overlap_extension,
            consensus_fraction = consensus_fraction,
            consensus_n = consensus_n
  )
  
  #3 scenarios, list of GRanges to overlap, already overlapped GRanges with memb table in mcols, or single GRanges
  #list of GRanges
  if(is.list(query_gr)){
    loaded = query_gr
    loaded = split(loaded, group_names)
    names(loaded) = paste0(names(loaded), "_features")
    obj@loaded_features = loaded
    #memb table  
  }else if(all(sapply(mcols(query_gr), is, class2 = "logical")) & ncol(mcols(query_gr)) > 0){
    loaded = lapply(seq_len(ncol(mcols(query_gr))), function(i){query_gr[mcols(query_gr)[,i]]})
    names(loaded) = colnames(mcols(query_gr))
    loaded = split(loaded, group_names)
    names(loaded) = paste0(names(loaded), "_features")
    obj@loaded_features = loaded
    obj@overlap_gr = list(query_features = query_gr)
    #just a GRanges
  }else{
    loaded = list(list(query_gr))
    names(loaded) = paste0(group_names, "_features")
    names(loaded[[1]]) = sample_names
    
    overlapped = loaded 
    mcols(overlapped[[1]][[1]]) = NULL
    overlapped[[1]][[1]]$query = TRUE
    colnames(mcols(overlapped[[1]][[1]])) = sample_names
    
    overlapped = lapply(overlapped, function(x)x[[1]])
    
    obj@loaded_features = loaded
    obj@overlap_gr = overlapped
    }
  
  if(process_features){
    obj = prepFeatures(obj)
  }
  obj
}


#' QcConfigFeatures for files
#'
#' @param file_paths character paths to files
#' @param run_separately If TRUE, each item is considered a separate group. Default is TRUE.
#' @param groups numeric vector of group assignments. 1 is first item in group_names, 2 is second, etc. Default is seq_along(file_path)
#' @param group_names vector of group names to assign from according to groups
#' @param group_colors vector of colors to use per group
#' @param n_peaks number of peaks to subset for
#' @param balance_groups If TRUE, will attempt to represent imbalanced groups more equally
#' @param consensus_n An integer number specifying the absloute minimum of input grs that must overlap for a site to be considered consensus.
#' @param consensus_fraction A numeric between 0 and 1 specifying the fraction of grs that must overlap to be considered consensus.
#'
#' @return a QcConfigFeatures object
#' @export
#' @rdname QcConfigFeatures
#' @examples
#' np_files = dir(system.file(package = "ssvQC", "extdata"), pattern = "Peak$", full.names = TRUE)
#' object = QcConfigFeatures.files(np_files, balance_groups = TRUE)
#' object = ssvQC.prepFeatures(object)
#' plot(object)
#' 
#' object = QcConfigFeatures.files(np_files, 
#'   group_names = c("10A", "AT1", "CA1"), 
#'   sample_names = c("MCF10A_CTCF", "MCF10AT1_CTCF", "MCF10CA1a_CTCF"))
#' object = ssvQC.prepFeatures(object)
#' plot(object)
QcConfigFeatures.files = function(file_paths,
                                  file_paths.input = character(),
                                  run_separately = TRUE,
                                  sample_names = NULL,
                                  sample_names.split = NULL,
                                  group_names = NULL,
                                  group_name.input = "input",
                                  group_colors = NULL,
                                  feature_load_FUN = NULL,
                                  n_peaks = 1e3,
                                  balance_groups = FALSE,
                                  overlap_extension = 0,
                                  consensus_fraction = getOption("SQC_CONSENSUS_FRACTION", 0),
                                  consensus_n = getOption("SQC_CONSENSUS_N", 1),
                                  process_features = getOption("SQC_PROCESS_FEATURES", TRUE)){
  if(is.null(group_names)){
    group_names = paste(seq_along(file_paths), basename(file_paths))
  }  
  if(is.null(group_colors)){
    group_colors = get_group_colors(group_names)
  }
  if(is.null(names(group_colors))){
    names(group_colors) = group_names
  }
  if(is.null(feature_load_FUN)){
    feature_load_FUN = get_feature_file_load_function(file_paths[1])[[1]]
  }
  config_df = data.frame(file = c(as.character(file_paths), file_paths.input), 
                         group = c(group_names, rep(group_name.input, length(file_paths.input))), 
                         All = "All", 
                         stringsAsFactors = FALSE)
  
  if(is.null(sample_names)){
    config_df$name = basename(config_df$file)  
  }else{
    config_df$name = sample_names
  }
  if(is.null(sample_names.split)){
    config_df$name_split = gsub("[_\\. ]", "\n", config_df$name)
  }else{
    config_df$name_split = sample_names.split
  }
  
  config_df$name = factor(config_df$name, levels = unique(config_df$name))
  config_df$name_split = factor(config_df$name_split, levels = unique(config_df$name_split))
  
  run_by = ifelse(run_separately, "group", "All")
  obj = new("QcConfigFeatures",
            meta_data =  config_df,
            run_by = run_by,
            to_run = unique(config_df[[run_by]]),
            to_run_reference = intersect(group_name.input, config_df$group),
            color_by = "group",
            color_mapping = group_colors,
            feature_load_FUN = feature_load_FUN,
            n_peaks = n_peaks,
            balance_groups = balance_groups,
            overlap_extension = overlap_extension,
            consensus_fraction = consensus_fraction,
            consensus_n = consensus_n
  )
  if(process_features){
    obj = prepFeatures(obj)
  }
  obj
}


#' QcConfigFeatures.parse
#'
#' @param feature_config_file 
#'
#' @return
#' @export
#' @rdname QcConfigFeatures
#' @examples
#' feature_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
#' object = QcConfigFeatures.parse(feature_config_file)
#' plot(object)
QcConfigFeatures.parse = function(feature_config_file,
                                  process_features = getOption("SQC_PROCESS_FEATURES", TRUE)){
  peak_config_dt = .parse_config_body(feature_config_file)
  valid_feature_var = c("main_dir", "overlap_extension", "n_peaks", "balance_groups", "consensus_n", 
                        "consensus_fraction", "color_by", "color_mapping", "run_by", "to_run", "to_run_reference", "is_null")
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
                  balance_groups = FALSE,
                  overlap_extension = 0,
                  consensus_n = 1, consensus_fraction = 0, 
                  color_by = NULL, color_mapping = NULL, 
                  run_by = NULL, to_run = NULL, to_run_reference = NULL, is_null = FALSE){
    QcConfigFeatures(config_df = config_dt, 
                     run_by = run_by, 
                     to_run = to_run,
                     to_run_reference = to_run_reference,
                     color_by = color_by, 
                     color_mapping = color_mapping, 
                     n_peaks = n_peaks, 
                     balance_groups = balance_groups,
                     overlap_extension = overlap_extension,
                     consensus_fraction = consensus_fraction,
                     consensus_n = consensus_n,
                     process_features = process_features,
                     is_null = is_null)
  }
  do.call(tfun, c(list(config_dt = peak_config_dt), cfg_vals))
}

#' Title
#'
#' @param object 
#' @param group_colors 
#' @param n_peaks 
#' @param balance_groups If TRUE, will attempt to represent imbalanced groups more equally
#' @param overlap_extension 
#' @param consensus_n An integer number specifying the absloute minimum of input grs that must overlap for a site to be considered consensus.
#' @param consensus_fraction A numeric between 0 and 1 specifying the fraction of grs that must overlap to be considered consensus.
#' @param process_features if TRUE, features are loaded and overlapped to generate assement regions as part of creating QcConfigFeatures object
#'
#' @return
#' @export
#'
#' @examples
#' np_files = dir(system.file(package = "ssvQC", "extdata"), pattern = "Peak$", full.names = TRUE)
#' object = QcConfigFeatures.files(np_files, run_separately = TRUE)
#' object = ssvQC.prepFeatures(object)
#' plot(object)
QcConfigFeatures.overlap_run_by = function(object,
                                           group_colors = NULL,
                                           n_peaks = 1e3,
                                           balance_groups = FALSE,
                                           overlap_extension = 0,
                                           consensus_fraction = getOption("SQC_CONSENSUS_FRACTION", 0),
                                           consensus_n = getOption("SQC_CONSENSUS_N", 1),
                                           process_features = getOption("SQC_PROCESS_FEATURES", TRUE) ){
  #TODO
  stop("NYI")
  object$loaded_features
  object$overlapped_features
  new_meta = object@meta_data
  new_meta$file = NA
  new_meta[[object@run_by]]
  new_meta[[object@color_by]]
  object@color_mapping
  
  new_obj = new("QcConfigFeatures",
                meta_data =  config_df,
                run_by = run_by,
                to_run = to_run,
                to_run_reference = to_run_reference,
                color_by = color_by,
                color_mapping = color_mapping,
                feature_load_FUN = feature_load_FUN,
                n_peaks = n_peaks,
                balance_groups = balance_groups,
                overlap_extension = overlap_extension,
                consensus_fraction = consensus_fraction,
                consensus_n = consensus_n,
                is_null = is_null)
}

#' @return
#' @export
#' @rdname QcConfigFeatures
#' @examples
#' feature_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
#' feature_config = QcConfigFeatures.parse(feature_config_file)
QcConfigFeatures.save_config = function(object, file){
  slots_to_save = c(
    "n_peaks",
    "balance_groups",
    "consensus_fraction",
    "consensus_n",
    "overlap_extension",
    "run_by",
    "to_run",
    "to_run_reference",
    "color_by",
    "is_null"
  )
  
  kvp_slots = "color_mapping"
  # QcConfigFeatures.parse(file)
  .save_config(object, file, slots_to_save, kvp_slots)
}
