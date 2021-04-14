

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
             file_paths = "character",
             groups = "numeric",
             group_names = "character",
             group_colors = "character"
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
QcConfig = function(file_paths,
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

    new("QcConfig",
        file_paths =  as.character(file_paths),
        groups = groups,
        group_names = group_names,
        group_colors = group_colors
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
#' QcConfigFeatures(c("A", "B"))
QcConfigFeatures = function(file_paths,
                            groups = NULL,
                            group_names = NULL,
                            group_colors = NULL,
                            n_peaks = 5e3,
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

    new("QcConfigFeatures",
        file_paths =  as.character(file_paths),
        groups = groups,
        group_names = group_names,
        group_colors = group_colors,
        n_peaks = n_peaks,
        min_fraction = min_fraction,
        min_number = min_number
    )
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
