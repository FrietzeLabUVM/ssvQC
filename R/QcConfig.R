DEFAULT_CONSENSUS_N = 1
DEFAULT_CONSENSUS_FRACTION = 0
DEFAULT_VIEW_SIZE = 3e3

#' QcConfig
#'
#' @slot file_paths character.
#' @slot groups numeric.
#' @slot group_names character.
#' @slot group_colors character.
#'
#' @export
setClass("QcConfig",
         representation = list(
             meta_data = "data.frame",
             group_var = "character",
             color_var = "character",
             color_mapping = "character"
         ))

#' QcConfig
#'
#' @param file_paths character paths to files
#' @param groups numeric vector of group assignments. 1 is first item in group_names, 2 is second, etc. Default is seq_along(file_path)
#' @param group_names vector of group names to assign from according to groups
#' @param group_colors vector of colors to use per group
#'
#' @return A QcConfig object
#' @export
#' @importFrom methods new
#' @examples
#' QcConfig(c("A", "B"))
QcConfig.files = function(file_paths,
                          groups = NULL,
                          group_names = NULL,
                          group_colors = NULL){
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
    cfg_df = data.frame(file = as.character(file_paths), group = group_names[groups])
    new("QcConfig",
        meta_data =  cfg_df,
        group_var = "group",
        color_mapping = group_colors
    )
}

#' QcConfigFeatures
#'
#' @slot feature_load_FUN function.
#' @slot n_peaks numeric.
#' @slot min_fraction numeric.
#' @slot min_number numeric.
#'
#' @export
setClass("QcConfigFeatures", contains = "QcConfig",
         representation = list(
             feature_load_FUN = "function",
             n_peaks = "numeric",
             min_fraction = "numeric",
             min_number = "numeric"
         ))


#' QcConfigFeatures
#'
#' @param cfg_df 
#' @param group_var 
#' @param color_var 
#' @param color_mapping 
#' @param feature_load_FUN 
#' @param n_peaks 
#' @param min_fraction 
#' @param min_number 
#'
#' @return
#' @export
#'
#' @examples
#' feature_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
#' cfg_df = .parse_config_body(feature_config_file)
#' QcConfigFeatures(cfg_df)
QcConfigFeatures = function(cfg_df,
                            group_var = "file",
                            color_var = group_var,
                            color_mapping = NULL,
                            feature_load_FUN = NULL,
                            n_peaks = 1e3,
                            min_fraction = DEFAULT_CONSENSUS_FRACTION,
                            min_number = DEFAULT_CONSENSUS_N){
    .enforce_file_var(cfg_df)
    if(!group_var %in% colnames(cfg_df)){
        stop("group_var ", group_var, " was not in column names.")
    }
    if(!color_var %in% colnames(cfg_df)){
        stop("color_var ", color_var, " was not in column names.")
    }
    
    if(!is.null(color_mapping)){
        if(!is.null(names(color_mapping))){
            color_names = names(color_mapping)
        }else if(is.factor(cfg_df[[color_var]])){
            color_names = levels(cfg_df[[color_var]])
        }else{
            color_names = unique(cfg_df[[color_var]])
        }
        stopifnot(length(color_names) == length(color_mapping))
        names(color_mapping) = color_names
        
    }else{
        if(is.factor(cfg_df[[color_var]])){
            color_names = levels(cfg_df[[color_var]])
        }else{
            color_names = unique(cfg_df[[color_var]])
        }
        color_mapping = seqsetvis::safeBrew(length(color_names))
        names(color_mapping) = color_names
    }
    
    stopifnot(cfg_df[[color_var]] %in% names(color_mapping))
    
    if(is.null(feature_load_FUN)){
        file_paths = cfg_df[["file"]]
        file_formats = guess_feature_file_format(file_paths)
        stopifnot(length(unique(file_formats)) == 1)
        feature_load_FUN = get_feature_file_load_function(file_paths[1])[[1]]
    }
    
    new("QcConfigFeatures",
        meta_data =  cfg_df,
        group_var = group_var,
        color_var = color_var,
        color_mapping = color_mapping,
        feature_load_FUN = feature_load_FUN,
        n_peaks = n_peaks,
        min_fraction = min_fraction,
        min_number = min_number)
}



#' QcConfigFeatures
#'
#' @param file_paths character paths to files
#' @param groups numeric vector of group assignments. 1 is first item in group_names, 2 is second, etc. Default is seq_along(file_path)
#' @param group_names vector of group names to assign from according to groups
#' @param group_colors vector of colors to use per group
#' @param n_peaks number of peaks to subset for
#' @param min_number An integer number specifying the absloute minimum of input grs that must overlap for a site to be considered consensus.
#' @param min_fraction A numeric between 0 and 1 specifying the fraction of grs that must overlap to be considered consensus.
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
                                  min_fraction = DEFAULT_CONSENSUS_FRACTION,
                                  min_number = DEFAULT_CONSENSUS_N){
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
    cfg_df = data.frame(file = as.character(file_paths), group = group_names[groups])
    
    new("QcConfigFeatures",
        meta_data =  cfg_df,
        group_var = "group",
        color_var = "group",
        color_mapping = group_colors,
        feature_load_FUN = feature_load_FUN,
        n_peaks = n_peaks,
        min_fraction = min_fraction,
        min_number = min_number
    )
}

.enforce_file_var = function(df){
    if(!"file" %in% colnames(df)){
        file_var = colnames(df)[1]
        message("Guessing file paths are in first column, ", file_var)
        colnames(df)[colnames(df) == file_var] = "file"
    }
    df
}

.enforce_name_var = function(df){
    if(!"name" %in% colnames(df)){
        cn = setdiff(colnames(df), c("file", "name", "name_split"))
        nams = apply(df[, cn], 1, function(x)paste(x, collapse = "_"))
        df$name = nams
    }
    if(!"name_split" %in% colnames(df)){
        cn = setdiff(colnames(df), c("file", "name", "name_split"))
        nams = apply(df[, cn], 1, function(x)paste(x, collapse = "\n"))
        df$name_split = nams
    }
    df
}


#' Title
#'
#' @param f 
#'
#' @return
#' @export
#'
#' @examples
.parse_config_body = function(f){
    config_dt = as.data.table(read.table(f, sep = ",", header = TRUE, stringsAsFactors = FALSE))
    config_dt = .enforce_file_var(config_dt)
    config_dt = .enforce_name_var(config_dt)
    #move file to first column
    config_dt = config_dt[, colnames(config_dt)[order(colnames(config_dt) != "file")], with = FALSE]
    config_dt
}

#' Title
#'
#' @param f 
#' @param valid_feature_var 
#'
#' @return
#' @export
#'
#' @examples
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
    cfg_vals = as.list(cfg_vals)
    
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
    for(var in c("consensus_n", "consensus_fraction", "n_peaks", "view_size")){
        if(!is.null(cfg_vals[[var]])){
            cfg_vals[[var]] = as.numeric(cfg_vals[[var]])
        }    
    }
    #parsing the color mapping
    if(!is.null(cfg_vals[["color_map"]])){
        cmap = cfg_vals[["color_map"]]
        if(grepl(":", cmap)){
            cmap = strsplit(strsplit(cmap, ",")[[1]], ":")
            color_map = sapply(cmap, function(x)x[2])
            names(color_map) = sapply(cmap, function(x)x[1])
            cfg_vals[["color_map"]] = color_map    
        }else{
            color_map = strsplit(cmap, ",")[[1]]
            cfg_vals[["color_map"]] = color_map
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
    
    cfg_vals
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
QcConfigFeatures.parse = function(feature_config_file){
    peak_config_dt = .parse_config_body(feature_config_file)
    valid_feature_var = c("main_dir", "n_peaks", "consensus_n", "consensus_fraction", "color_by", "color_map", "run_by", "to_run")
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
        stop(paste(c("Files specified in config do not exist:", peak_config_dt$file[!file.exists(peak_config_dt$file)]), collapse = "\n  "))
    }
    
    tfun = function(config_dt, 
                    main_dir = NULL, 
                    n_peaks = 1e3, 
                    consensus_n = 1, consensus_fraction = 0, 
                    color_by = NULL, color_map = NULL, 
                    run_by = NULL, to_run = NULL){
        QcConfigFeatures(cfg_df = config_dt, 
                         group_var = run_by, 
                         color_var = color_by, 
                         color_mapping = color_map, 
                         n_peaks = n_peaks, 
                         min_fraction = consensus_fraction,
                         min_number = consensus_n)
    }
    do.call(tfun, c(list(config_dt = peak_config_dt), cfg_vals))
    
}

' QcConfigSignal
#'
#' @slot view_size numeric.
#'
#' @export
#'
setClass("QcConfigSignal", contains = "QcConfig",
         representation = list(
             view_size = "numeric",
             read_mode = "character"
         ))


#' Title
#'
#' @param cfg_df 
#' @param group_var 
#' @param color_var 
#' @param color_mapping 
#' @param read_mode 
#' @param view_size 
#'
#' @return
#' @export
#'
#' @examples
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' bam_cfg_df = .parse_config_body(bam_config_file)
#' QcConfigSignal(bam_cfg_df)
#' 
#' bigwig_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bigwig_config.csv")
#' bigwig_cfg_df = .parse_config_body(bigwig_config_file)
#' QcConfigSignal(bigwig_cfg_df)
QcConfigSignal = function(cfg_df,
                          group_var = "file",
                          color_var = group_var,
                          color_mapping = NULL,
                          read_mode = NULL,
                          view_size = DEFAULT_VIEW_SIZE){
    .enforce_file_var(cfg_df)
    if(!group_var %in% colnames(cfg_df)){
        stop("group_var ", group_var, " was not in column names.")
    }
    if(!color_var %in% colnames(cfg_df)){
        stop("color_var ", color_var, " was not in column names.")
    }
    
    if(!is.null(color_mapping)){
        if(!is.null(names(color_mapping))){
            color_names = names(color_mapping)
        }else if(is.factor(cfg_df[[color_var]])){
            color_names = levels(cfg_df[[color_var]])
        }else{
            color_names = unique(cfg_df[[color_var]])
        }
        stopifnot(length(color_names) == length(color_mapping))
        names(color_mapping) = color_names
        
    }else{
        if(is.factor(cfg_df[[color_var]])){
            color_names = levels(cfg_df[[color_var]])
        }else{
            color_names = unique(cfg_df[[color_var]])
        }
        color_mapping = seqsetvis::safeBrew(length(color_names))
        names(color_mapping) = color_names
    }
    
    stopifnot(cfg_df[[color_var]] %in% names(color_mapping))
    
    if(is.null(read_mode)){
        read_mode = guess_read_mode(cfg_df$file[1])
    }
    
    stopifnot(read_mode %in% c("bam_SE", "bam_PE", "bigwig"))
    
    new("QcConfigSignal",
        meta_data =  cfg_df,
        group_var = group_var,
        color_var = color_var,
        color_mapping = color_mapping,
        read_mode = read_mode,
        view_size = view_size)
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
    valid_feature_var = c("main_dir", "view_size", "read_mode", "color_by", "color_map", "run_by", "to_run")
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
        stop(paste(c("Files specified in config do not exist:", signal_config_dt$file[!file.exists(signal_config_dt$file)]), collapse = "\n  "))
    }
    
    tfun = function(config_dt, 
                    main_dir = NULL, 
                    read_mode = NULL,
                    view_size = DEFAULT_VIEW_SIZE,
                    color_by = NULL, color_map = NULL, 
                    run_by = NULL, to_run = NULL){
        QcConfigSignal(cfg_df = config_dt, 
                       group_var = run_by, 
                       color_var = color_by, 
                       color_mapping = color_map, 
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
#' 
#' np_files = dir(system.file(package = "ssvQC", "extdata"), pattern = "Peak$", full.names = TRUE)
#' QcConfigFeatures.files(np_files)
QcConfigSignal.files = function(file_paths,
                                groups = NULL,
                                group_names = NULL,
                                group_colors = NULL,
                                view_size = DEFAULT_VIEW_SIZE, 
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
    
    cfg_df = data.frame(file = as.character(file_paths), group = group_names[groups])
    
    QcConfigSignal(cfg_df, group_var = "group", color_mapping = group_colors)
}




#' QcColorMapping
#'
#' @param object a QcConfig object
#'
#' @return a named vector of colors
#' @export
#' @rdname QcColorMapping-methods
setGeneric("QcColorMapping", function(object){standardGeneric("QcColorMapping")})

#' QcColorMapping for QcConfig
#'
#' @export
#' @rdname QcColorMapping-methods
#' @aliases QcColorMapping,QcConfig-method
#' @examples
#' QcColorMapping(QcConfig(c("A", "B")))
setMethod("QcColorMapping", c("QcConfig"), function(object){
    cols = object@group_colors
    names(cols) = as.character(object@group_names)
    cols
})


#' QcScaleColor
#'
#' @param object a QcConfig object
#'
#' @return a ggplot2 scale_fill_manual
#' @export
#' @rdname QcScaleColor-methods
setGeneric("QcScaleColor", function(object){standardGeneric("QcScaleColor")})

#' QcScaleColor for QcConfig
#'
#' @export
#' @rdname QcScaleColor-methods
#' @aliases QcScaleColor,QcConfig-method
#' @examples
#' df = data.frame(group = c("A", "B"))
#' df$x = 1:2
#' df$y = 1
#' library(ggplot2)
#' ggplot(df, aes(x = x, y = y, color = group)) +
#'   geom_point(size = 20) +
#'   expand_limits(x = c(0, 3)) +
#'   QcScaleColor(QcConfig(df$group))
setMethod("QcScaleColor", c("QcConfig"), function(object){
    cols = QcColorMapping(object)
    ggplot2::scale_color_manual(values = cols)
})

#' QcScaleFill
#'
#' @param object a QcConfig object
#'
#' @return a ggplot2 scale_fill_manual
#' @export
#' @rdname QcScaleFill-methods
setGeneric("QcScaleFill", function(object){standardGeneric("QcScaleFill")})

#' QcScaleFill for QcConfig
#'
#' @export
#' @rdname QcScaleFill-methods
#' @aliases QcScaleFill,QcConfig-method
#' @examples
#' df = data.frame(group = c("A", "B"))
#' df$x = 1:2
#' df$y = 1
#' library(ggplot2)
#' ggplot(df, aes(x = x, y = y, fill = group)) +
#'   geom_point(size = 20, pch = 22) +
#'   expand_limits(x = c(0, 3)) +
#'   QcScaleFill(QcConfig(df$group))
setMethod("QcScaleFill", c("QcConfig"), function(object){
    cols = object@group_colors
    names(cols) = as.character(object@group_names)
    ggplot2::scale_fill_manual(values = cols)
})

# ob1 = QcConfig(LETTERS[1:5], groups = c(rep(1, 3), rep(2, 2)))
# ob2 = QcConfigFeatures(LETTERS[1:5], groups = c(rep(1, 3), rep(2, 2)))
# ob3 = QcConfigSignal(LETTERS[1:5], groups = c(rep(1, 3), rep(2, 2)))
# QcColorMapping(ob3)
# QcScaleFill(ob3)
