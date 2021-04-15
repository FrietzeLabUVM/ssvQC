

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
        color_var = "group",
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


#' Title
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
#' cfg_df = parse_config_body(feature_config_file)
#' QcConfigFeatures(cfg_df)
QcConfigFeatures = function(cfg_df,
                             group_var = "file",
                             color_var = group_var,
                             color_mapping = NULL,
                             feature_load_FUN = NULL,
                            n_peaks = 1e3,
                            min_fraction = 0,
                            min_number = 1){
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
        names(color_mapping) = color_var
        
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
#' QcConfigFeatures.files(c("A.narrowPeak", "B.narrowPeak"))
QcConfigFeatures.files = function(file_paths,
                                  group_names = NULL,
                                  groups = NULL,
                                  group_colors = NULL,
                                  feature_load_FUN = NULL,
                                  n_peaks = 1e3,
                                  min_fraction = 0,
                                  min_number = 1){
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


#' Title
#'
#' @param f 
#'
#' @return
#' @export
#'
#' @examples
parse_config_body = function(f){
    as.data.table(read.table(f, sep = ",", header = TRUE, stringsAsFactors = FALSE))
}

#' Title
#'
#' @param feature_config_file 
#'
#' @return
#' @export
#' feature_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
#'
#' @examples
parseQcConfigFeatures = function(feature_config_file){
    peak_config_dt = .parse_config_body(feature_config_file)
    peak_config_dt = .enforce_file_var(peak_config_dt)
    #move file to first column
    peak_config_dt = peak_config_dt[, colnames(peak_config_dt)[order(colnames(peak_config_dt) != "file")], with = FALSE]
    
    
    cfg_txt = fread(feature_config_file, sep = "\n", header = FALSE)[grepl("^#CFG", V1)]
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
    
    valid_feature_var = c("main_dir", "consensus_n", "consensus_fraction", "color_by", "color_map", "run_by", "to_run")
    bad_var = setdiff(cfg_names, valid_feature_var)
    if(length(bad_var) > 0){
        stop("Unrecogized variables in config file: ", paste(bad_var, collapse = ", "))
    }
    
    
    tfun = function(config_dt, main_dir, consensus_n, consensus_fraction, color_by, color_map, run_by, to_run){
        message(main_dir)
        message(consensus_n)
    }
    do.call(tfun, as.list(config_dt = peak_config_dt, cfg_vals))
    
}

parseQcConfigSignal = function(signal_config_file){
    
}

#' QcConfigSignal
#'
#' @slot view_size numeric.
#'
#' @export
#'
setClass("QcConfigSignal", contains = "QcConfig",
         representation = list(
             view_size = "numeric"
         ))

#' QcConfigSignal
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
#' QcConfigSignal(c("A", "B"))
QcConfigSignal = function(file_paths,
                          groups = NULL,
                          group_names = NULL,
                          group_colors = NULL,
                          view_size = 3e3){
    if(is.null(groups)){
        groups = seq_along(file_paths)
    }
    if(is.null(group_names)){
        group_names = LETTERS[seq_along(unique(groups))]
    }
    if(is.null(group_colors)){
        group_colors = seqsetvis::safeBrew(length(group_names))
    }
    
    new("QcConfigSignal",
        file_paths = as.character(file_paths),
        groups = groups,
        group_names = group_names,
        group_colors = group_colors,
        view_size = view_size
    )
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
