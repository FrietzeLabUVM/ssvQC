#' sync_width
#' 
#' synchronize the element widths (axis, plot, etc.) of a list of ggplots.
#' 
#' To work properly the number and order of elements must match between plots.  In practical terms, this means either all or non must have vertical facets.
#'
#' @param my_plots A
#'
#' @return A list of grobs matching input my_plots list of ggplots.
#' @export
#'
#' @examples
#' library(ggplot2)
#' 
#' theme_update(
#'   plot.title.position = "plot", 
#'   axis.title.y = element_text(angle = 0, vjust = .5, hjust = 0)
#' )
#' 
#' p1 = ggplot(mtcars, aes(x = mpg, y = cyl)) +
#'   geom_point() +
#'   labs(title = "cylinders vs mpg", y = "cylinder") 
#' 
#' p2 = ggplot(mtcars, aes(x = mpg, y = disp)) +
#'   geom_point() +
#'   labs(title = "disp vs mpg") 
#' 
#' p3 = ggplot(mtcars, aes(x = mpg, y = hp)) +
#'   geom_point() +
#'   labs(title = "horsepower vs mpg") 
#' 
#' plots = list(p1, p2, p3)
#' 
#' #without synchronization, x axis positions do not line up between plots.
#' cowplot::plot_grid(plotlist = plots, ncol = 1)
#' 
#' #after synchronization, x axis positions line up perfectly between plots.
#' plots.sync = sync_width(plots)
#' cowplot::plot_grid(plotlist = plots.sync, ncol = 1)
sync_width = function(my_plots){
  stopifnot(class(my_plots) == "list")
  is_ok = sapply(my_plots, function(x){
    "ggplot" %in% class(x) | "grob" %in% class(x)
  })
  stopifnot(all(is_ok))
  my_grobs = lapply(my_plots, function(x){
    if(grid::is.grob(x)){
      x
    }else{
      ggplotGrob(x)
    }
  })
  
  my_widths = lapply(my_grobs, function(gt){
    gt$widths
  })
  maxWidth = my_widths[[1]]
  if(length(my_widths) > 1){
    for(i in 2:length(my_widths)){
      maxWidth = grid::unit.pmax(maxWidth, my_widths[[i]])
    }
  }
  for(j in 1:length(my_grobs)){
    my_grobs[[j]]$widths = maxWidth
  }
  my_grobs
}

#' sync_height
#' 
#' synchronize the element heights (axis, plot, etc.) of a list of ggplots.
#' 
#' To work properly the number and order of elements must match between plots.  In practical terms, this means either all or non must have horizontal facets.
#'
#' @param my_plots A
#'
#' @return A list of grobs matching input my_plots list of ggplots.
#' @export
#'
#' @examples
#' library(ggplot2)
#' 
#' theme_update(
#'   plot.title.position = "plot", 
#'   axis.title.y = element_text(angle = 0, vjust = .5, hjust = 0), 
#'   axis.title.x = element_text(angle = 90, vjust = .5, hjust = 0)
#' )
#' 
#' p1 = ggplot(mtcars, aes(y = mpg, x = cyl)) +
#'   geom_point() +
#'   labs(title = "cylinders vs mpg", x = "cylinder") 
#' 
#' p2 = ggplot(mtcars, aes(y = mpg, x = disp)) +
#'   geom_point() +
#'   labs(title = "disp vs mpg") 
#' 
#' p3 = ggplot(mtcars, aes(y = mpg, x = hp)) +
#'   geom_point() +
#'   labs(title = "horsepower vs mpg") 
#' 
#' plots = list(p1, p2, p3)
#' 
#' #without synchronization, y axis positions do not line up between plots.
#' cowplot::plot_grid(plotlist = plots, nrow = 1)
#' 
#' #after synchronization, y axis positions line up perfectly between plots.
#' plots.sync = sync_height(plots)
#' cowplot::plot_grid(plotlist = plots.sync, nrow = 1)
sync_height = function(my_plots){
  stopifnot(class(my_plots) == "list")
  is_ok = sapply(my_plots, function(x){
    "ggplot" %in% class(x) | "grob" %in% class(x)
  })
  stopifnot(all(is_ok))
  my_grobs = lapply(my_plots, function(x){
    if(grid::is.grob(x)){
      x
    }else{
      ggplotGrob(x)
    }
  })
  
  my_widths = lapply(my_grobs, function(gt){
    gt$heights
  })
  maxWidth = my_widths[[1]]
  if(length(my_widths) > 1){
    for(i in 2:length(my_widths)){
      maxWidth = grid::unit.pmax(maxWidth, my_widths[[i]])
    }
  }
  for(j in 1:length(my_grobs)){
    my_grobs[[j]]$heights = maxWidth
  }
  my_grobs
}

#' sampleCap
#' 
#' A safer version of sample with additional code so that if a sample size of n is greater than length of x, n is reduced and a reordered x is returned.
#'
#' @inheritParams base::sample
#'
#' @return A vector of length size or original length of x if size is too larger, with elements drawn randomly from x.
#' @export
#'
#' @examples
#' x = LETTERS
#' #this would cause an error is normal sample
#' sampleCap(x, 50)
#' #this is equivalent to sample
#' sampleCap(x, 5)
sampleCap = function(x, size = 500){
  size = min(size, length(unique(x)))
  out = sample(unique(x), size)
  # if(is.factor(out)) out = as.character(out) #not sure why this would be necessary
  out
}

#' get_mapped_reads
#'
#' @param bam_file A bam file.  Matching .bai file must exist.
#'
#' @return The number of mapped reads in bam file.
#' @export
#'
#' @examples
#' bam_file = system.file("extdata/MCF10A_CTCF_R1.100peaks.bam", package = "ssvQC")
#' get_mapped_reads(bam_file)
get_mapped_reads = function(bam_file){
  stopifnot(file.exists(bam_file))
  stopifnot(file.exists(paste0(bam_file, ".bai")))
  stats = Rsamtools::idxstatsBam(bam_file)
  sum(stats[,3])
}


#' guess_feature_file_format
#' 
#' Check if feature file has a common format. ie. narrowPeak, broadPeak, or bed.
#'  
#' @param feature_files A feature or interval file path. A vector of multiple files is acceptable.
#'
#' @return A character vector describing each entry in feature_files.
#' @export
#'
#' @examples
#' feature_files = dir(system.file("extdata", package = "ssvQC"), pattern = "bed$|Peak$", full.names = TRUE)
#' guess_feature_file_format(feature_files)
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
      warning("\n", paste0("Could not guess file format for feature file:\n  ", feature_file), "\n")
    }
    file_format
  }
  
  sapply(feature_files, .guess_feature_file_format)
}

#' guess_read_mode
#'
#' @param signal_file A bam of bigwig file path.
#'
#' @return A character describing format of signal_file.
#' @export
#'
#' @examples
#' bam_file = dir(system.file("extdata", package = "ssvQC"), pattern = "bam$", full.names = TRUE)[1]
#' bw_file = dir(system.file("extdata", package = "ssvQC"), pattern = "bw$", full.names = TRUE)[1]
#' guess_read_mode(bam_file)
#' guess_read_mode(bw_file)
guess_read_mode = function(signal_file){
  if(signal_file[1] == "null"){
    return("null")
  }
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

#' get_feature_file_load_function
#'
#' @param feature_files 
#'
#' @return list of appropriate functions for loading feature_files.
#' @export
#'
#' @examples
#' feature_files = dir(system.file("extdata", package = "ssvQC"), pattern = "bed$|Peak$", full.names = TRUE)
#' load_FUNs = get_feature_file_load_function(feature_files)
#' all_loaded = lapply(feature_files, function(f){
#'   get_feature_file_load_function(f)[[1]](f)[[1]]
#' })
get_feature_file_load_function = function(feature_files){
  file_types = guess_feature_file_format(feature_files)
  .get_feature_file_load_function = function(file_type){
    switch (file_type,
            narrowPeak = seqsetvis::easyLoad_narrowPeak,
            broadPeak = seqsetvis::easyLoad_broadPeak,
            seacr = seqsetvis::easyLoad_seacr,
            bed = seqsetvis::easyLoad_bed,
            unknown = {
              warning("Treating unknown file type as bed but if you see errors, check file type and specify appropriate feature_load_FUN when creating QcConfigFeatures.")
              seqsetvis::easyLoad_bed
            },
            stop("'", file_type, "' is not a supported file_type")
    )
  }
  
  sapply(file_types, .get_feature_file_load_function)
}

.enforce_file_var = function(my_df){
  if(!"file" %in% colnames(my_df)){
    file_var = colnames(my_df)[1]
    message("Guessing file paths are in first column, ", file_var)
    colnames(my_df)[colnames(my_df) == file_var] = "file"
  }
  my_df
}

.enforce_name_var = function(my_df){
  if(!"name" %in% colnames(my_df)){
    cn = setdiff(colnames(my_df), c("file", "name", "name_split"))
    message("Creating 'name' from: ", paste(cn, collapse = ", "))
    nams = apply(my_df[, cn], 1, function(x)paste(x, collapse = "_"))
    my_df$name = nams
  }
  if(!"name_split" %in% colnames(my_df)){
    cn = setdiff(colnames(my_df), c("file", "name", "name_split"))
    message("Creating 'name_split' from: ", paste(cn, collapse = ", "))
    nams = apply(my_df[, cn], 1, function(x)paste(x, collapse = "\n"))
    my_df$name_split = nams
  }
  if(any(duplicated(my_df$name))){
    stop("Duplicate entries in 'name', all values must be unique.")
  }
  if(any(duplicated(my_df$name_split))){
    stop("Duplicate entries in 'name_split', all values must be unique.")
  }
  if(!is.factor(my_df$name))
    my_df$name = factor(my_df$name, levels = my_df$name)
  if(!is.factor(my_df$name_split))
    my_df$name_split = factor(my_df$name_split, levels = my_df$name_split)
  my_df
}

.enforce_found_order = function(my_df){
  for(var in colnames(my_df)){
    if(is.character(my_df[[var]])){
      if(var != "file")
        my_df[[var]] = factor(my_df[[var]], levels = unique(my_df[[var]]))
    }
  }
  my_df
}

#' @importFrom utils read.table
.parse_config_body = function(f){
  config_dt = utils::read.table(f, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  config_dt = .enforce_file_var(config_dt)
  config_dt = .enforce_name_var(config_dt)
  config_dt = .enforce_found_order(config_dt)
  #move file to first column
  config_dt = config_dt[, colnames(config_dt)[order(colnames(config_dt) != "file")]]
  config_dt
}


#' parse_fetch_options
#' 
#' @param fop character of options
#'
#' @return named list of options
#' 
#' @examples
#' fop = "win_size:10,win_method:\"summary\",summary_FUN:mean"
#' parse_fetch_options(fop)
parse_fetch_options = function(fop){
  if(is.null(fop)) return(list())
  if(any(is.na(fop))) return(list())
  # fop = strsplit(fop, ",")[[1]]
  fop = strsplit(fop, ":")
  fop_names = sapply(fop, function(x)x[1])
  fop_opts = sapply(fop, function(x)x[2])
  if(any(is.na(fop_names))){
    stop("problem parsing fetch_options, verify fetch_options=key1:val1,key2:val2,... syntax")
  }
  if(any(is.na(fop_opts))){
    stop("problem parsing fetch_options, verify fetch_options=key1:val1,key2:val2,... syntax")
  }
  if("summary_FUN" %in% fop_names){
    message("summary_FUN is not yet fully supported.  Definitions of summary_FUN should work but its value will not be saved via QcConfigSignal.save_config.")
  }
  names(fop_opts) = fop_names
  lapply(fop_opts, function(x){
    #check if number
    if(suppressWarnings({!is.na(as.numeric(x))})){
      x = as.numeric(x)
    }else if(grepl('"', x)){
      x = gsub('"', "", x)
    }else{
      x = get(x)
    }
    x
  })
}

#' .parse_config_header
#'
#' @param f a config file
#' @param valid_feature_var supported variables to attempt to extract
#'
#' @return A named list containing configuration options mapped to values.
#'
#' @examples
#' valid_feature_var = c("main_dir", "overlap_extension", "n_peaks", 
#'   "balance_groups", "consensus_n", 
#'   "consensus_fraction", "color_by", "color_mapping", 
#'   "run_by", "to_run", "to_run_reference", "is_null")
#' cfg_file = system.file("extdata/ssvQC_peak_config.csv", package = "ssvQC")
#' .parse_config_header(cfg_file, valid_feature_var)
.parse_config_header = function(f, valid_feature_var){
  cfg_txt = fread(f, sep = "\n", header = FALSE)[grepl("^#CFG", V1)]
  cfg_txt = paste(sub("#CFG ?", "", cfg_txt$V1), collapse = " ")
  sp = strsplit(cfg_txt, " +")[[1]]
  sp = strsplit(sp, "=")
  sp
  cfg_names = sapply(sp, function(x){
    x[1]
  })
  cfg_vals = sapply(sp, function(x){
    x[2]
  })
  names(cfg_vals) = cfg_names
  cfg_vals = strsplit(cfg_vals, ",")
  
  bad_var = setdiff(cfg_names, valid_feature_var)
  if(length(bad_var) > 0){
    stop("Unrecogized variables in config file: ", paste(bad_var, collapse = ", "))
  }
  
  ### special cases
  if(!is.null(cfg_vals[["main_dir"]])){
    if(cfg_vals[["main_dir"]] == "$SSVQC_DATA"){
      cfg_vals[["main_dir"]] = system.file("extdata", package = "ssvQC")
    }
  }
  #numeric: consensus_n, consensus_fraction, n_peaks, view_size
  for(var in c("consensus_n", "consensus_fraction", "n_peaks", "view_size", "overlap_extension", "heatmap_limit_values", "linearQuantile_cutoff")){
    if(!is.null(cfg_vals[[var]])){
      cfg_vals[[var]] = as.numeric(cfg_vals[[var]])
      if(!is.numeric(cfg_vals[[var]])){
        stop("The variable, '", var, "' is only supported as a numeric currently.")
      }
    }    
  }
  #logical: is_null
  for(var in c("is_null", "lineplot_free_limits", "balance_groups")){
    if(!is.null(cfg_vals[[var]])){
      cfg_vals[[var]] = as.logical(cfg_vals[[var]])
    }    
  }
  
  #parsing the color mapping
  if(!is.null(cfg_vals[["color_mapping"]])){
    cmap = cfg_vals[["color_mapping"]]
    if(any(grepl(":", cmap))){
      cmap = strsplit(cmap, ":")
      color_mapping = sapply(cmap, function(x)x[2])
      names(color_mapping) = sapply(cmap, function(x)x[1])
      cfg_vals[["color_mapping"]] = color_mapping    
    }else{
      color_mapping = cmap
      cfg_vals[["color_mapping"]] = color_mapping
    }
  }
  #read mode synonym
  if(!is.null(cfg_vals[["read_mode"]])){
    read_mode = cfg_vals[["read_mode"]]
    if(read_mode == "SE") read_mode = "bam_SE"
    if(read_mode == "PE") read_mode = "bam_PE"
    if(!read_mode %in% c("bam_SE", "bam_PE", "bigwig")){
      stop("read_mode '", read_mode, "' is not recognized as a valid choice from 'bam_SE', 'bam_PE', or 'bigwig'")
    }
    cfg_vals[["read_mode"]] = read_mode
  }
  #fetch_options
  if(!is.null(cfg_vals[["fetch_options"]])){
    cfg_vals[["fetch_options"]] = parse_fetch_options(cfg_vals[["fetch_options"]])
  }
  
  cfg_vals
}

.test_suff = function(files, suff){
  sapply(files, function(f){
    any(sapply(suff, function(s){
      grepl(paste0(s, "$"), f)
    }))
  })
}

is_feature_file = function(files, suff = getOption("SQC_FEATURE_FILE_SUFF", c("narrowPeak", "broadPeak", "bed", "txt", "tab"))){
  .test_suff(files, suff)
}

is_signal_file = function(files, suff = getOption("SQC_SIGNAL_FILE_SUFF", c("bam", "bigwig", "bw", "bigWig", "BigWig"))){
  .test_suff(files, suff)
}

#' internal function used by QcConfig.save_config QcConfigSignal.save_config and QcConfigFeatures.save_config
.save_config = function(object, file, slots_to_save, kvp_slots, toss_names = "summary_FUN"){
  hdr1 = sapply(slots_to_save, function(x){
    val = slot(object, x)
    ifelse(length(val) > 0,
           paste0("#CFG ", x, "=", paste(val, collapse = ",")),
           character())
  })
  hdr1 = hdr1[!is.na(hdr1)]
  
  hdr2 = sapply(kvp_slots, function(x){
    val = slot(object, x)
    val = val[!names(val) %in% names(toss_names)]
    val = paste(names(val), val, sep = ":", collapse = ",")
    ifelse(length(val) > 0,
           paste0("#CFG ", x, "=", val),
           character())
  })
  hdr2 = hdr2[!is.na(hdr2)]
  
  hdr3 = paste(colnames(object@meta_data), collapse = ",")
  
  
  hdr = c(hdr1, hdr2, hdr3)
  names(hdr) = NULL
  
  writeLines(hdr, file)
  fwrite(object@meta_data, file, append = TRUE)
  invisible(file)
}

get_group_colors = function(group_names){
  cols = getOption("SQC_COLORS", seqsetvis::safeBrew(group_names, "Dark2"))
  cols
}

#' ssv_mclapply
#'
#' @return result of either pblapply or pbmclapply
#'
#' @importFrom pbapply pblapply
#' @importFrom pbmcapply pbmclapply
ssv_mclapply = function(X, FUN, mc.cores = getOption("mc.cores", 1), ...){
  if(.Platform$OS.type == "windows" || mc.cores == 1) {
    pbapply::pblapply(X = X, FUN = FUN, ...)
    
  } else {
    pbmcapply::pbmclapply(X = X, FUN = FUN, mc.cores = mc.cores, ...)
  }
}
