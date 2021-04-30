.onAttach <- function(libname, pkgname) {
    packageStartupMessage("Attaching ssvQC version ",
                          packageDescription("ssvQC")$Version, ".")
    options("SQC_COLORS" = seqsetvis::safeBrew(8, "Dark2"))
    options("SQC_CONSENSUS_N" = 1)
    options("SQC_CONSENSUS_FRACTION" = 0)
    options("SQC_VIEW_SIZE" = 3e3)
    options("SQC_PROCESS_FEATURES" = TRUE)
}
# getOption("SQC_COLORS")

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
             run_by = "character",
             to_run = "character",
             to_run_reference = "character",
             color_by = "character",
             color_mapping = "character",
             is_null = "logical"
         ))

setMethod("initialize","QcConfig", function(.Object,...){
    .Object <- callNextMethod()
    validObject(.Object)
    .Object@is_null = FALSE
    .Object
})

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
#' QcConfig.files(c("A", "B"))
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
    config_df = data.frame(file = as.character(file_paths),  color = group_names[groups], group = "All", stringsAsFactors = FALSE)
    new("QcConfig",
        meta_data =  config_df,
        run_by = "group",
        color_by = "color",
        to_run = unique(config_df$group),
        color_mapping = group_colors
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
#' my_df = data.frame(group = c("A", "B"))
#' my_df$x = 1:2
#' my_df$y = 1
#' library(ggplot2)
#' ggplot(my_df, aes(x = x, y = y, color = group)) +
#'   geom_point(size = 20) +
#'   expand_limits(x = c(0, 3)) +
#'   QcScaleColor(QcConfig(my_df$group))
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
#' my_df = data.frame(group = c("A", "B"))
#' my_df$x = 1:2
#' my_df$y = 1
#' library(ggplot2)
#' ggplot(my_df, aes(x = x, y = y, fill = group)) +
#'   geom_point(size = 20, pch = 22) +
#'   expand_limits(x = c(0, 3)) +
#'   QcScaleFill(QcConfig(my_df$group))
setMethod("QcScaleFill", c("QcConfig"), function(object){
    cols = object@group_colors
    names(cols) = as.character(object@group_names)
    ggplot2::scale_fill_manual(values = cols)
})



#' .show_QcConfig
#' 
#' used by show
#'
#' @param qc A QcConfig object
#'
#' @return
#'
#' @examples
#' feature_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
#' qc = QcConfigFeatures.parse(feature_config_file)
#' qc
.show_QcConfig = function(qc){
    if(qc@is_null){
        msg = "This QcConfig is a NULL placeholder."
    }else{
        msg = paste(sep = "\n",
                    paste("Configuration for", nrow(qc@meta_data), "items."),
                    paste0(length(unique(qc@meta_data[[qc@run_by]])), " potential running group(s) detected for '", qc@run_by, "'."),
                    paste("Will run:", paste(qc@to_run, collapse = ", ")),
                    ifelse(length(qc@to_run_reference) > 0, 
                           paste("Referenced vs:", paste(qc@to_run_reference, collapse = ", ")),
                           "No reference set."),
                    paste0("Use plot() to view color mapping for '", qc@color_by , "'.")
        )
    }
    message(msg)
}

.plot_QcConfig = function(qc){
    plot_dt = as.data.table(qc@meta_data)
    plot_dt[, y := rev(seq_len(.N)), c(qc@run_by)]
    ggplot(plot_dt, aes_string(x = qc@run_by, y = "y", fill = qc@color_by, label = "name_split")) +
        geom_label() +
        scale_fill_manual(values = qc@color_mapping) +
        theme(panel.background = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
        labs(y = "") +
        guides(
            fill = guide_legend(
                override.aes = aes(label = "")
            ))
}

#' Title
#'
#' @param QcConfig 
#'
#' @return
#' @export
#'
#' @examples
setMethod("plot", "QcConfig", definition = function(x).plot_QcConfig(x))


#' Title
#'
#' @param QcConfig 
#'
#' @return
#' @export
#'
#' @examples
setMethod("show", "QcConfig", definition = function(object).show_QcConfig(object))

if(FALSE){
    library(ssvQC)
    feature_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
    object = QcConfigFeatures.parse(feature_config_file)
    plot(object)
    object$loaded_features = "asdf"
    object$run_by = "mark"
    object$run_by
    plot(object)
    object$run_by = "cell"
    object$run_by
    plot(object)
    object$color_by = "mark"
    plot(object)
    object$color_by = "cell"
    plot(object)
    options("SQC_COLORS" = seqsetvis::safeBrew(8, "blues"))
    object$color_by = "cell"
    plot(object)
}
