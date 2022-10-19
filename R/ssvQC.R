setOldClass(c("theme", "gg"))

setClass("ssvQC",
         representation = list(
           features_config = "QcConfigFeatures",
           signal_config = "QcConfigSignal",
           signal_data = "list",
           other_data = "list",
           out_dir = "character",
           bfc = "BiocFileCache",
           saving_enabled = "logical",
           plot_post_FUN = "function",
           plots = "list",
           matched_only = "logical"
         ))
#' @export
setClass("ssvQC.featureOnly", contains = "ssvQC")
#' @export
setClass("ssvQC.signalOnly", contains = "ssvQC")
#' @export
setClass("ssvQC.complete", contains = "ssvQC")

#' @export
setClass("ssvQC.signalOnly.bw", contains = "ssvQC.signalOnly")
#' @export
setClass("ssvQC.complete.bw", contains = "ssvQC.complete")

setMethod("initialize","ssvQC", function(.Object,...){
  .Object <- callNextMethod()
  validObject(.Object)
  .Object@plot_post_FUN = function(p)p
  .Object@plots = list()
  .Object@signal_data = list()
  .Object
})

.show_ssvQC = function(qc){
  message("Features configuration:")
  print(qc$features_config)
  message("\n")
  message("Signal configuration:")
  print(qc$signal_config)
  message("\n")
  if(length(qc@signal_data) > 0){
    message("Signal data has been LOADED.")
  }else{
    message("Signal data has NOT been loaded.")
  }
  if(length(qc@other_data) > 0){
    message(paste(names(qc@other_data), collapse = ", "), " have been LOADED.")
  }else{
    message("NO other data have been loaded.")
  }
  if(length(qc@plots) > 0){
    message(paste(names(qc@plots), collapse = ", "), " have been PLOTTED.")
  }else{
    message("NO plots have been made.")
  }
}

.plot_ssvQC = function(qc){
  p1 = plot(qc$signal_config) + labs(title = "Signal configuration")
  p2 = plot(qc$features_config) + labs(title = "Features configuration")
  cowplot::plot_grid(p1, p2)
}

#' @export
setMethod("plot", "ssvQC", definition = function(x).plot_ssvQC(x))

#' ssvQC
#'
#' @param ssvQC 
#'
#' @return
#' @export
#' @rdname ssvQC
#' @examples
setMethod("show", "ssvQC", definition = function(object).show_ssvQC(object))

#' ssvQC.save_config
#'
#' @param object A ssvQC object
#' @param file A .features.csv and .signal.csv file will be written with this prefix.
#'
#' @return invisibly returns config file paths.
#' @export
#'
#' @examples
ssvQC.save_config = function(object, file){
  feature_file = paste0(sub(".csv", "",  file), ".features.csv")
  signal_file = paste0(sub(".csv", "",  file), ".signal.csv")
  QcConfigFeatures.save_config(object@features_config, feature_file)
  QcConfigSignal.save_config(object@signal_config, signal_file)
  invisible(c(features = feature_file, signal = signal_file))
}

#' ssvQC
#'
#' @param features_config Controls features configuration.  May be a:
#'   QcConfigFeatures object, path to a file defining configuration via
#'   QcConfigFeatures.parse, features files to define via
#'   QcConfigFeatures.files, or a data.frame to pass to QcConfigFeatures.
#' @param signal_config Controls signal configuration.  May be a: QcConfigSignal
#'   object, path to a file defining configuration via QcConfigSignal.parse,
#'   features files to define via QcConfigSignal.files, or a data.frame to pass
#'   to QcConfigSignal.
#' @param out_dir NYI
#' @param bfc BiocFileCache object to use for caching. If NULL, default
#'   new_cache() will be used.
#'
#' @return A ssvQC object.  Data needs to be loaded after via ssvQC.runAll or sub-methods ssvQC.plot*.
#' @export
#' @import BiocFileCache
#' 
#' @rdname ssvQC
#' @examples
#' options(mc.cores = 1)
#' set.seed(0)
#' # To make an ssvQC object, confiugration for features (peaks, other genomic 
#' regions) and signal (numeric values on the genome from bam pileups or bigwigs)
#' features_config_file = system.file(
#'   package = "ssvQC", 
#'   "extdata/ssvQC_peak_config.csv"
#' )
#' features_config = QcConfigFeatures.parse(features_config_file)
#'
#' bam_config_file = system.file(
#'   package = "ssvQC", 
#'   "extdata/ssvQC_bam_config.csv"
#' )
#' bam_config = QcConfigSignal.parse(bam_config_file)
#'
#' bigwig_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bigwig_config.csv")
#' bigwig_config = QcConfigSignal.parse(bigwig_config_file)
#' 
#' # Different ways to make ssvQC objects
#' sqc.complete.file = ssvQC(features_config_file, bam_config_file)
#'
#' sqc.complete = ssvQC(features_config, bam_config)
#'
#' sqc.complete.bw = ssvQC(features_config, bigwig_config_file)
#'
#' sqc.signal = ssvQC(signal_config = bam_config)
#'
#' sqc.feature = ssvQC(features_config = features_config)
#'
#' # ssvQC.runAll will run all appropriate QC methods
#' sqc.signal = ssvQC.runAll(sqc.signal)
#' 
#' sqc.feature = ssvQC.runAll(sqc.feature)
#'
#' sqc.complete = ssvQC.runAll(sqc.complete)
#' 
#' sqc.complete$plots$signal$heatmaps
#' sqc.complete$signal_config@plot_value = "RPM"
#' sqc.complete = ssvQC.plotSignal(sqc.complete)
#' sqc.complete$plots$signal$heatmaps
#' 
#' sqc.complete$signal_config@plot_value = "linearQuantile"
#' sqc.complete = ssvQC.plotSignal(sqc.complete)
#' sqc.complete$plots$signal$heatmaps
#'
#' sqc.complete$signal_config@plot_value = SQC_SIGNAL_VALUES$RPM_linearQuantile
#' sqc.complete = ssvQC.plotSignal(sqc.complete)
#' sqc.complete$plots$signal$heatmaps
#' 
#' write_ssvQC.summary(sqc.complete)
#' write_ssvQC.per_peak(sqc.complete)
#' write_ssvQC.correlation(sqc.complete)
#' 
ssvQC = function(features_config = NULL,
                 signal_config = NULL,
                 out_dir = getwd(),
                 bfc = NULL, 
                 matched_only = TRUE){
  if(is.null(features_config) & is.null(signal_config)){
    stop("At least one of features_config or signal_config must be specified.")
  }
  
  features_config = .prep_features_config(features_config)
  signal_config = .prep_signal_config(signal_config)
  
  if(is.null(bfc)){
    bfc = new_cache()
  }
  
  dir.create(out_dir, showWarnings = FALSE)
  
  if(!is.null(features_config) & !is.null(signal_config)){
    if(signal_config@read_mode == SQC_READ_MODES$bigwig){
      new("ssvQC.complete.bw",
          features_config = features_config,
          signal_config = signal_config,
          signal_data = list(),
          other_data = list(),
          out_dir = out_dir,
          bfc = bfc,
          saving_enabled = TRUE,
          matched_only = matched_only
      )
    }else{
      new("ssvQC.complete",
          features_config = features_config,
          signal_config = signal_config,
          signal_data = list(),
          other_data = list(),
          out_dir = out_dir,
          bfc = bfc,
          saving_enabled = TRUE,
          matched_only = matched_only
      )  
    }
  }else if(!is.null(features_config)){
    new("ssvQC.featureOnly",
        features_config = features_config,
        signal_config = QcConfigSignal.null(),
        signal_data = list(),
        other_data = list(),
        out_dir = out_dir,
        bfc = bfc,
        saving_enabled = TRUE,
        matched_only = matched_only
    )
  }else if(!is.null(signal_config)){
    if(signal_config@read_mode == SQC_READ_MODES$bigwig){
      new("ssvQC.signalOnly.bw",
          features_config = QcConfigFeatures.null(),
          signal_config = signal_config,
          signal_data = list(),
          other_data = list(),
          out_dir = out_dir,
          bfc = bfc,
          saving_enabled = TRUE,
          matched_only = matched_only
      )
    }else{
      new("ssvQC.signalOnly",
          features_config = QcConfigFeatures.null(),
          signal_config = signal_config,
          signal_data = list(),
          other_data = list(),
          out_dir = out_dir,
          bfc = bfc,
          saving_enabled = TRUE,
          matched_only = matched_only
      )
    }
  }else{
    stop("At least one of features_config or signal_config must be specified. This should have been caught earlier.")
  }
}

.prep_features_config = function(features_config, ...){
  if(!is.null(features_config)){
    if(is.character(features_config)){
      if(!any(is_feature_file(features_config))){
        features_config = QcConfigFeatures.parse(features_config)
      }else{
        features_config = QcConfigFeatures.files(features_config)
      }
    }else if(is.data.frame(features_config)){
      features_config = QcConfigFeatures(features_config, ...)
    }
    if(!"QcConfigFeatures" %in% class(features_config)){
      stop("features_config must be either a QcConfigFeatures object or the path to valid configuration file to create one.")
    }  
  }
  features_config
}

.prep_signal_config = function(signal_config, ...){
  if(!is.null(signal_config)){
    if(is.character(signal_config)){
      if(!any(is_signal_file(signal_config))){
        signal_config = QcConfigSignal.parse(signal_config)
      }else{
        signal_config = QcConfigSignal.files(signal_config)
      }
    }else if(is.data.frame(signal_config)){
      signal_config = QcConfigSignal(signal_config, ...)
    }
    if(!"QcConfigSignal" %in% class(signal_config)){
      stop("signal_config must be either a QcConfigSignal object or the path to valid configuration file to create one.")
    }  
    f_exists = file.exists(as.character(signal_config@meta_data$file))
    if(!all(f_exists)){
      stop("Not all signal files found! Missing:\n", paste(signal_config@meta_data$file[!f_exists], collapse = "\n"))
    }
  }
  signal_config
}

.make_query_signal_config = function(sc){
  out = split(sc)
  names(out) = paste0(names(out), "_signal")
  out
}

#' @export
#' @rdname ssvQC
setGeneric("ssvQC.runAll", function(object){standardGeneric("ssvQC.runAll")})
setMethod("ssvQC.runAll", "ssvQC.complete", function(object){
  message("run+plot features overlaps")
  object = ssvQC.plotFeatures(object)
  message("run+plot mapped reads")
  # object = ssvQC.prepMappedReads(object)
  if(object$signal_config$read_mode != "bigwig"){
    object = ssvQC.plotMappedReads(object)
  }
  if(object$signal_config$read_mode == "bam_SE"){
    message("run signal fragLens")
    object = ssvQC.prepFragLens(object)  
  }
  message("run signal normalization")
  object = ssvQC.prepCapValue(object)
  message("plot signal")
  object = ssvQC.plotSignal(object)
  message("plot SCC")
  object = ssvQC.plotSCC(object)
  message("plot FRIP")
  object = ssvQC.plotFRIP(object)
  message("plot correlation")
  object = ssvQC.plotCorrelation(object)
  object
})
setMethod("ssvQC.runAll", "ssvQC.complete.bw", function(object){
  message("run+plot features overlaps")
  object = ssvQC.plotFeatures(object)
  message("run signal normalization")
  object = ssvQC.prepCapValue(object)
  message("plot signal")
  object = ssvQC.plotSignal(object)
  message("plot correlation")
  object = ssvQC.plotCorrelation(object)
  object
})

setMethod("ssvQC.runAll", "ssvQC.featureOnly", function(object){
  object = ssvQC.plotFeatures(object)
  object
})
setMethod("ssvQC.runAll", "ssvQC.signalOnly", function(object){
  object = ssvQC.plotMappedReads(object)
  object
})
setMethod("ssvQC.runAll", "ssvQC.signalOnly.bw", function(object){
  message("Nothing to do with signalOnly for bigwigs.")
  object
})

feature_name2signal_name = function(feature_name){
  sub("features$", "signal", feature_name)
}

signal_name2feature_name = function(signal_name){
  sub("signal$", "features", signal_name)
}

dbl_names = function(name_1, name_2){
  paste(sub("_signal", " signal", name_2), "at", sub("_features", " features", name_1))
}

dbl_lapply = function(in_list, FUN, FUN_names = NULL){
  out_list = lapply(names(in_list), function(name_1){
    item_1 = in_list[[name_1]]
    out = lapply(names(item_1), function(name_2){
      item_2 = item_1[[name_2]]
      if(!is.null(FUN_names)){
        mt = dbl_names(name_1, name_2)
        FUN(item_2, mt)
      }else{
        FUN(item_2)
      }
    })
    names(out) = names(item_1)
    out
  })
  names(out_list) = names(in_list)
  out_list
}

dbl_extract = function(in_list, key){
  lapply(in_list, function(item_1){
    lapply(item_1, function(item_2){
      item_2[[key]]
    })
  })
}

single_extract = function(in_list, key){
  lapply(in_list, function(item_1){
    item_1[[key]]
  })
}

#' @export
#' @rdname ssvQC
setGeneric("ssvQC.prepFetch", function(object){standardGeneric("ssvQC.prepFetch")})
setMethod("ssvQC.prepFetch", "ssvQC.complete", function(object){
  if(length(object@features_config$assessment_features) == 0){
    object = ssvQC.prepFeatures(object)
  }
  if(is.null(object$signal_config$meta_data$mapped_reads) & object$signal_config$read_mode != "bigwig"){
    object = ssvQC.prepMappedReads(object)
  }
  if(is.null(object$signal_config$meta_data$fragLens) & object$signal_config$read_mode == "bam_SE"){
    object = ssvQC.prepFragLens(object)
  }
  if(is.null(object$signal_config$meta_data$cap_value)){
    object = ssvQC.prepCapValue(object)
  }
  object
})
setMethod("ssvQC.prepFetch", "ssvQC.featureOnly", function(object){
  stop("Cannot run prepFetch on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.prepFetch", "ssvQC.signalOnly", function(object){
  stop("Cannot run prepFetch on ssvQC with no QcConfigFeature component")
})

#' @export
#' @rdname ssvQC
setGeneric("ssvQC.referenceUsesSameScale", function(object){standardGeneric("ssvQC.referenceUsesSameScale")})
setMethod("ssvQC.referenceUsesSameScale", "ssvQC.complete", function(object){
  object@signal_config = ssvQC.referenceUsesSameScale(object@signal_config)
  object
})
setMethod("ssvQC.referenceUsesSameScale", "ssvQC.featureOnly", function(object){
  stop("Cannot run referenceUsesSameScale on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.referenceUsesSameScale", "ssvQC.signalOnly", function(object){
  object@signal_config = ssvQC.referenceUsesSameScale(object@signal_config)
  object
})
setMethod("ssvQC.referenceUsesSameScale", "QcConfigSignal", function(object){
  reference_uses_same_scale(object)
})

### $ Accessor
setMethod("names", "ssvQC",
          function(x)
          {
            c("plots", "signal_data", "signal_config", "features_config", "SCC", "FRIP", "correlation")
            
          })


setMethod("$", "ssvQC",
          function(x, name)
          {
            switch (name,
                    plots = x@plots,
                    signal_data = x@signal_data,
                    SCC = x@other_data$SCC,
                    FRIP = x@other_data$FRIP,
                    correlation = list(read_count = x@other_data$read_count_correlation, signal_profile = x@other_data$signal_profile_correlation),
                    bfc = x@bfc,
                    features_config = x@features_config,
                    signal_config = x@signal_config
                    
            )
          })

setReplaceMethod("$", "ssvQC",
                 function(x, name, value)
                 {
                   warn_msg = "This assignment is not supported.  No effect."
                   switch (name,
                           features_config = {
                             x@features_config = value
                           },
                           signal_config = {
                             x@signal_config = value
                           },
                           {warning(warn_msg)}
                           
                   )
                   
                   #TODO, some assignments may be appropriate
                   x
                 })


