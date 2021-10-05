
# getOption("SQC_COLORS")

#' QcConfig
#'
#' @slot file_paths character.
#' @slot groups numeric.
#' @slot group_names character.
#' @slot group_colors character.
#' @rdname QcConfig
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
#' @rdname QcConfig
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
#' @rdname QcConfig
setGeneric("QcColorMapping", function(object){standardGeneric("QcColorMapping")})

#' QcColorMapping for QcConfig
#'
#' @export
#' @rdname QcConfig
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
#' @rdname QcConfig
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
#' @rdname QcConfig
setGeneric("QcScaleFill", function(object){standardGeneric("QcScaleFill")})

#' QcScaleFill for QcConfig
#'
#' @export
#' @aliases QcScaleFill,QcConfig-method
#' @rdname QcConfig
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
    if(qc@is_null){
        msg = "This QcConfig is a NULL placeholder."
        ggplot() + labs(title = "This QcConfig is a NULL placeholder.")
    }else{
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
    
}

#' @export
setMethod("plot", "QcConfig", definition = function(x).plot_QcConfig(x))

#' Title
#'
#' @param QcConfig 
#'
#' @return
#' @export
#' @rdname QcConfig
#' @examples
setMethod("show", "QcConfig", definition = function(object).show_QcConfig(object))

#' Title
#'
#' @param object 
#'
#' @return
#' @export
#' @rdname QcConfig
#' @examples
QcConfig.save_config = function(object){
    slots_to_save = c(
        "run_by",
        "to_run",
        "to_run_reference",
        "color_by",
        "is_null"
    )
    kvp_slots = c("color_mapping")
    # QcConfigSignal.parse(file)
    .save_config(object, file, slots_to_save, kvp_slots)
}
