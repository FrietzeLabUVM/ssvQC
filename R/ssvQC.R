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

setClass("ssvQC.featureOnly", contains = "ssvQC")
setClass("ssvQC.signalOnly", contains = "ssvQC")
setClass("ssvQC.complete", contains = "ssvQC")

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
    message(paste(names(qc@other_data), collapse = ", "), "have been LOADED.")
  }else{
    message("NO other data have been loaded.")
  }
  if(length(qc@plots) > 0){
    message(paste(names(qc@plots), collapse = ", "), "have been PLOTTED.")
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
#' features_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
#' features_config = QcConfigFeatures.parse(features_config_file)
#'
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' bam_config = QcConfigSignal.parse(bam_config_file)
#'
#' bigwig_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bigwig_config.csv")
#' bigwig_config = QcConfigSignal.parse(bigwig_config_file)
#'
#' sqc.complete.file = ssvQC(features_config_file, bam_config_file)
#'
#' sqc.complete = ssvQC(features_config, bam_config)
#' object = sqc.complete
#'
#' sqc.complete.bw = ssvQC(features_config, bigwig_config_file)
#'
#' sqc.signal = ssvQC(signal_config = bam_config)
#'
#' sqc.feature = ssvQC(features_config = features_config)
#'
#' sqc.signal = ssvQC.runAll(sqc.signal)
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
                 matched_only = TRUE,
                 ...){
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
    stopifnot(file.exists(signal_config@meta_data$file))
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
setMethod("ssvQC.runAll", "ssvQC.featureOnly", function(object){
  object = ssvQC.plotFeatures(object)
  object
})
setMethod("ssvQC.runAll", "ssvQC.signalOnly", function(object){
  object = ssvQC.plotMappedReads(object)
  object
})

##FragLens for SE bams
.prepFragLens = function(bam_f, peak_gr, name, n_regions, bfc){
  fl = bfcif(bfc, digest::digest(list(bam_f, peak_gr, n_regions)), function(){
    seqsetvis::fragLen_calcStranded(bam_f, peak_gr, n_regions = n_regions)
  })
  dt = data.table(fragLens = fl)
  dt$name = name
  dt
}


#' @export
#' @rdname ssvQC
setGeneric("ssvQC.prepFragLens", function(object, query, bfc, use_matched){standardGeneric("ssvQC.prepFragLens")})
setMethod("ssvQC.prepFragLens", "ssvQC.complete", function(object){
  object@signal_config = ssvQC.prepFragLens(object@signal_config, object@features_config, object@bfc, object@matched_only)
  object
})
setMethod("ssvQC.prepFragLens", "ssvQC.featureOnly", function(object){
  stop("Cannot run prepCapValue on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.prepFragLens", "ssvQC.signalOnly", function(object){
  stop("Cannot run prepCapValue on ssvQC with no QcConfigFeatures component")
})
setMethod("ssvQC.prepFragLens", c("QcConfigSignal", "QcConfigFeatures", "BiocFileCache", "logical"), function(object, query, bfc, use_matched){
  if(object@read_mode != "bam_SE"){
    stop("ssvQC.prepFragLens only appropriate for read_mode bam_SE")
  }
  
  #bam specific independent of peaks
  
  sig_dt = as.data.table(object@meta_data)[, .(file, name)][order(file)]
  peak_dt = as.data.table(query@meta_data)[, .(file, name)][order(file)]
  query = ssvQC.prepFeatures(query)
  
  if(use_matched){
    matched_dt = merge(sig_dt[, .(bam_file = file, name)], peak_dt[, .(peak_file = file, name)], by = "name")
    
    if(nrow(matched_dt) == 0){
      use_matched = FALSE
    }else if(any(!file.exists(matched_dt$peak_file))){
      use_matched = FALSE
    }else{
      matched_peaks_gr = query@feature_load_FUN(matched_dt$peak_file)
      matched_dt = matched_dt[lengths(matched_peaks_gr) > 0,]
      matched_peaks_gr = matched_peaks_gr[lengths(matched_peaks_gr) > 0]
      
      unmatched_dt = sig_dt[!name %in% matched_dt$name]    
    }
  }
  if(!use_matched){
    matched_dt = data.table()
    unmatched_dt = sig_dt
  }
  
  
  fl_dt = bfcif(bfc, digest::digest(list(sig_dt, peak_dt, "ssvQC.prepFragLens")), function(){
    if(nrow(matched_dt) > 0){
      fl_dt.matched = rbindlist(lapply(seq_len(nrow(matched_dt)), function(i){
        peak_f = matched_dt[i,]$peak_file
        bam_f = matched_dt[i,]$bam_file
        name = matched_dt[i,]$name
        peak_gr = matched_peaks_gr[[1]]
        # rname = digest::digest(list(bam_f, name, peak_gr, "ssvQC.prepFragLens"))
        .prepFragLens(bam_f, peak_gr, name, 500, bfc)
      }))
    }else{
      fl_dt.matched = NULL
    }
    
    if(nrow(unmatched_dt) > 0){
      peak_gr = unlist(GRangesList(query$assessment_features))
      
      fl_dt.unmatched = rbindlist(lapply(seq_len(nrow(unmatched_dt)), function(i){
        bam_f = unmatched_dt[i,]$file
        name = unmatched_dt[i,]$name
        # rname = digest::digest(list(bam_f, name, peak_gr, "ssvQC.prepFragLens"))
        .prepFragLens(bam_f, peak_gr, name, 500, bfc)
      }))
    }else{
      fl_dt.unmatched = NULL
    }
    fl_dt = rbind(fl_dt.matched, fl_dt.unmatched)
    
    if(!setequal(fl_dt$name, object@meta_data$name)){
      stop("something has gone wrong assigning fragLens, pease report this issue at https://github.com/FrietzeLabUVM/ssvQC/issues") 
    }
    setkey(fl_dt, "name")
    fl_dt
  })
  
  object@meta_data$fragLens = fl_dt[.(object@meta_data$name)]$fragLens
  object
})

setMethod("ssvQC.prepFragLens", c("QcConfigSignal", "QcConfigFeatures"), function(object, query){
  ssvQC.prepFragLens(object, query, new_cache())
})

##MappedReads

#' @export
#' @rdname ssvQC
setGeneric("ssvQC.prepMappedReads", function(object){standardGeneric("ssvQC.prepMappedReads")})
setMethod("ssvQC.prepMappedReads", "ssvQC.complete", function(object){
  object@signal_config = ssvQC.prepMappedReads(object@signal_config)
  object
})
setMethod("ssvQC.prepMappedReads", "ssvQC.featureOnly", function(object){
  stop("Cannot run mapped reads on ssvQC with no QcConfigSignal component")
  message("featureOnly")
})
setMethod("ssvQC.prepMappedReads", "ssvQC.signalOnly", function(object){
  object@signal_config = ssvQC.prepMappedReads(object@signal_config)
  object
})
setMethod("ssvQC.prepMappedReads", c("QcConfigSignal"), function(object){
  #bam specific independent of peaks
  if(grepl("bam", object@read_mode)){
    object@meta_data$mapped_reads = sapply(object@meta_data$file, get_mapped_reads)
  }else{
    warning("ssvQC.prepMappedReads called on bigwig read_mode QcConfigSignal.")
    object@meta_data$mapped_reads = NA
  }
  object
})

##CapValue used by linearQuantile

#' @export
#' @rdname ssvQC
setGeneric("ssvQC.prepCapValue", function(object, query, bfc, use_matched){standardGeneric("ssvQC.prepCapValue")})
setMethod("ssvQC.prepCapValue", "ssvQC.complete", function(object){
  object@signal_config = ssvQC.prepCapValue(object@signal_config, object@features_config, object@bfc, object@matched_only)
  object
})
setMethod("ssvQC.prepCapValue", "ssvQC.featureOnly", function(object){
  stop("Cannot run prepCapValue on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.prepCapValue", "ssvQC.signalOnly", function(object){
  stop("Cannot run prepCapValue on ssvQC with no QcConfigFeatures component")
})
setMethod("ssvQC.prepCapValue", c("QcConfigSignal", "QcConfigFeatures", "BiocFileCache"), function(object, query, bfc, use_matched){
  #bam specific independent of peaks
  if(is.null(object@meta_data$fragLens)){
    sig_dt = as.data.table(object@meta_data)[, .(file, name)][order(file)]
  }else{
    sig_dt = as.data.table(object@meta_data)[, .(file, name, fragLens)][order(file)]  
  }
  object = ssvQC.prepMappedReads(object)
  
  setkey(sig_dt, "name")
  peak_dt = as.data.table(query@meta_data)[, .(file, name)][order(file)]
  query = ssvQC.prepFeatures(query)
  
  if(use_matched){
    matched_dt = merge(sig_dt[, .(bam_file = file, name)], peak_dt[, .(peak_file = file, name)], by = "name")
    
    if(nrow(matched_dt) == 0){
      use_matched = FALSE
    }else if(any(!file.exists(matched_dt$peak_file))){
      use_matched = FALSE
    }else{
      matched_peaks_gr = query@feature_load_FUN(matched_dt$peak_file)
      matched_dt = matched_dt[lengths(matched_peaks_gr) > 0,]
      matched_peaks_gr = matched_peaks_gr[lengths(matched_peaks_gr) > 0]
      
      unmatched_dt = sig_dt[!name %in% matched_dt$name]  
    }
  }
  if(!use_matched){
    matched_dt = data.table()
    unmatched_dt = sig_dt
  }
  
  cap_dt = bfcif(bfc, digest::digest(list(sig_dt, peak_dt, "ssvQC.prepCapValue")), function(){
    if(nrow(matched_dt) > 0){
      cap_dt.matched = rbindlist(lapply(seq_len(nrow(matched_dt)), function(i){
        peak_f = matched_dt[i,]$peak_file
        bam_f = matched_dt[i,]$bam_file
        # if(grepl("bam", object@read_mode)){
        #   bam_f = data.table(file = bam_f, mapped_reads = get_mapped_reads(bam_f))
        # }
        name_i = matched_dt[i,]$name
        peak_gr = query@feature_load_FUN(peak_f)[[1]]
        peak_gr = peak_gr[sample(min(5e3, length(peak_gr)))]
        
        if(object@read_mode == "bam_SE"){
          fragLens = sig_dt[.(name_i)]$fragLens  
          fragLens = ifelse(is.null(fragLens), "auto", fragLens)
          max_dt = get_fetch_fun(object@read_mode)(bam_f, 
                                                   peak_gr, 
                                                   fragLens = fragLens,
                                                   win_size = 1, 
                                                   win_method = "summary", 
                                                   summary_FUN = function(x,w)max(x), 
                                                   return_data.table = TRUE, 
                                                   n_region_splits = getOption("mc.cores", 1))
        }else{
          max_dt = get_fetch_fun(object@read_mode)(bam_f, 
                                                   peak_gr, 
                                                   win_size = 1, 
                                                   win_method = "summary", 
                                                   summary_FUN = function(x,w)max(x), 
                                                   return_data.table = TRUE, 
                                                   n_region_splits = getOption("mc.cores", 1))
        }
        cap_dt = seqsetvis::calc_norm_factors(max_dt)
        cap_dt$name = name_i
        cap_dt
      }))
    }else{
      cap_dt.matched = NULL
    }
    
    if(nrow(unmatched_dt) > 0){
      query = ssvQC.prepFeatures(query)
      peak_gr = unlist(GRangesList(query$assessment_features))
      names(peak_gr) = NULL
      
      cap_dt.unmatched = rbindlist(lapply(seq_len(nrow(unmatched_dt)), function(i){
        bam_f = unmatched_dt[i,]$file
        name_i = unmatched_dt[i,]$name
        # rname = digest::digest(list(bam_f, name, peak_gr, "ssvQC.prepFragLens"))
        if(object@read_mode == "bam_SE"){
          fragLens = sig_dt[.(name_i)]$fragLens  
          fragLens = ifelse(is.null(fragLens), "auto", fragLens)
          max_dt = get_fetch_fun(object@read_mode)(bam_f, 
                                                   peak_gr, 
                                                   fragLens = fragLens,
                                                   win_size = 1, 
                                                   win_method = "summary", 
                                                   summary_FUN = function(x,w)max(x), 
                                                   return_data.table = TRUE, 
                                                   n_region_splits = getOption("mc.cores", 1))
        }else{
          max_dt = get_fetch_fun(object@read_mode)(bam_f, 
                                                   peak_gr, 
                                                   win_size = 1, 
                                                   win_method = "summary", 
                                                   summary_FUN = function(x,w)max(x), 
                                                   return_data.table = TRUE, 
                                                   n_region_splits = getOption("mc.cores", 1))
        }
        cap_dt = seqsetvis::calc_norm_factors(max_dt)
        cap_dt$name = name_i
        cap_dt
      }))
    }else{
      cap_dt.unmatched = NULL
    }
    cap_dt = rbind(cap_dt.matched, cap_dt.unmatched)
    
    if(!setequal(cap_dt$name, object@meta_data$name)){
      stop("something has gone wrong assigning fragLens, pease report this issue at https://github.com/FrietzeLabUVM/ssvQC/issues") 
    }
    setkey(cap_dt, "name")
    cap_dt
  })
  
  cap_dt[y_cap_value <= 1, y_cap_value := 1]
  
  if(grepl("bam", object@read_mode)){
    cap_dt[, mapped_reads := get_mapped_reads(as.character(sample)), .(sample) ]
    cap_dt[, y_RPM_cap_value := y_cap_value / max(mapped_reads, 1) * 1e6]
  }
  
  object@meta_data$cap_value = cap_dt[.(object@meta_data$name)]$y_cap_value
  if(!is.null(cap_dt$y_RPM_cap_value)){
    object@meta_data$RPM_cap_value = cap_dt[.(object@meta_data$name)]$y_RPM_cap_value  
  }
  object
})

.ssvQC.plotMappedReads = function(object){
  full_bam_config = object@signal_config
  full_meta = full_bam_config@meta_data
  if(is.null(full_meta$mapped_reads)){
    full_bam_config = ssvQC.prepMappedReads(full_bam_config)
    object@signal_config = full_bam_config
  }
  
  todo = .make_query_signal_config(full_bam_config)
  plots = lapply(todo, function(sig_config){
    if(grepl("bam", sig_config@read_mode)){
      bam_config_dt = sig_config@meta_data
      
      color_var = sig_config@color_by
      group_var = sig_config@run_by
      color_mapping = sig_config@color_mapping
      
      p_mapped_reads = ggplot(bam_config_dt, aes_string(x = "name_split", y = "mapped_reads", fill = color_var)) +
        geom_bar(stat = "identity", position = "dodge", color = "black") +
        scale_fill_manual(values = color_mapping) +
        scale_y_continuous(labels = function(x)x/1e6) +
        labs(y = "M mapped reads", fill = color_var, x= "") +
        labs(title = "Mapped reads")
      
      if(is.null(object@plots$reads)) object@plots$reads = list()
      
    }else{
      p_mapped_reads = ggplot() + theme_void() + labs(title = "Could not run mapped reads on non-bam")
    }
    p_mapped_reads
  })
  object@plots$reads = plots
  
  object
}

#' @export
#' @rdname ssvQC
setGeneric("ssvQC.plotMappedReads", function(object){standardGeneric("ssvQC.plotMappedReads")})
setMethod("ssvQC.plotMappedReads", "ssvQC.complete", .ssvQC.plotMappedReads)
setMethod("ssvQC.plotMappedReads", "ssvQC.featureOnly", function(object){
  stop("Cannot run mapped reads on ssvQC with no QcConfigSignal component")
  message("featureOnly")
})
setMethod("ssvQC.plotMappedReads", "ssvQC.signalOnly", .ssvQC.plotMappedReads)

feature_name2signal_name = function(feature_name){
  sub("features$", "signal", feature_name)
}

signal_name2feature_name = function(signal_name){
  sub("signal$", "features", signal_name)
}



##FRIP
#' @export
#' @rdname ssvQC
setGeneric("ssvQC.prepFRIP", function(object){standardGeneric("ssvQC.prepFRIP")})
setMethod("ssvQC.prepFRIP", "ssvQC.complete", function(object){
  if(length(object@features_config$assessment_features) == 0){
    object = ssvQC.prepFeatures(object)
  }
  
  feature_names = names(object@features_config$assessment_features)
  names(feature_names) = feature_names
  FRIP_data = lapply(feature_names, function(name){
    query_gr = object@features_config$assessment_features[[name]]
    sig_configs = .make_query_signal_config(object@signal_config)
    
    must_match = object@matched_only
    if(must_match){
      sig_name = feature_name2signal_name(name)
      if(!sig_name %in% names(sig_configs)){
        must_match = FALSE
      }else{
        out = lapply(sig_configs[feature_name2signal_name(name)], function(sel_sig_config){
          make_frip_dt(as.data.table(sel_sig_config@meta_data), query_gr = query_gr, color_var = sel_sig_config@color_by)
        })    
      }
    }
    if(!must_match){
      out = lapply(sig_configs, function(sel_sig_config){
        make_frip_dt(as.data.table(sel_sig_config@meta_data), query_gr = query_gr, color_var = sel_sig_config@color_by)
      })  
    }
    out
  })
  object@other_data$FRIP = FRIP_data
  object
})
setMethod("ssvQC.prepFRIP", "ssvQC.featureOnly", function(object){
  stop("Cannot run prepFRIP on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.prepFRIP", "ssvQC.signalOnly", function(object){
  stop("Cannot run prepFRIP on ssvQC with no QcConfigFeature component")
})

#' @export
#' @rdname ssvQC
setGeneric("ssvQC.plotFRIP", function(object){standardGeneric("ssvQC.plotFRIP")})
setMethod("ssvQC.plotFRIP", "ssvQC.complete", function(object){
  if(is.null(object@other_data$FRIP)){
    object = ssvQC.prepFRIP(object)
  }
  FRIP_data = object@other_data$FRIP
  
  wrap_plot_frip_dt = function(frip_dt, main_title){
    plot_frip_dt(frip_dt, 
                 color_var = object@signal_config@color_by, 
                 color_mapping = object@signal_config@color_mapping,
                 main_title = main_title)
  }
  
  plots = dbl_lapply(FRIP_data, wrap_plot_frip_dt, dbl_names)
  
  if(is.null(object@plots$FRIP)){
    object@plots$FRIP = list()
  }
  
  object@plots$FRIP$reads_per_peak = dbl_extract(plots, "reads_per_peaks")
  object@plots$FRIP$per_peak = dbl_extract(plots, "frip_per_peaks")
  object@plots$FRIP$total = dbl_extract(plots, "frip_total")
  object
})
setMethod("ssvQC.plotFRIP", "ssvQC.featureOnly", function(object){
  stop("Cannot run plotFRIP on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.plotFRIP", "ssvQC.signalOnly", function(object){
  stop("Cannot run plotFRIP on ssvQC with no QcConfigFeature component")
})

##Correlation
#' @export
#' @rdname ssvQC
setGeneric("ssvQC.prepCorrelation", function(object){standardGeneric("ssvQC.prepCorrelation")})
setMethod("ssvQC.prepCorrelation", "ssvQC.complete", function(object){
  if(is.null(object@other_data$FRIP)){
    object = ssvQC.prepFRIP(object)
  }
  if(length(object@signal_data) == 0){
    object = ssvQC.prepSignal(object)
  }
  Correlation_data = lapply(object@other_data$FRIP, function(feature_FRIP_data){
    lapply(feature_FRIP_data, function(frip_dt){
      reads_dt = dcast(frip_dt, id~name_split, value.var = "reads_in_peak")
      mat = as.matrix(reads_dt[,-1])
      rownames(mat) = reads_dt$id
      cor_mat = cor(mat)
      row_clust = hclust(dist(cor_mat))
      col_clust = hclust(dist(t(cor_mat)))
      mat_dt = as.data.table(reshape2::melt(cor_mat))
      setnames(mat_dt, c("Var1", "Var2"), c("row_name_split", "col_name_split"))
      # mat_dt$row_name_split = factor(mat_dt$row_name_split, levels = rownames(cor_mat)[row_clust$order])
      # mat_dt$col_name_split = factor(mat_dt$col_name_split, levels = colnames(cor_mat)[col_clust$order])
      return(list(mat = cor_mat, dt = mat_dt, row_hclust = row_clust, col_hclust = col_clust))
    })
  })
  
  Correlation_data.profile = lapply(object@signal_data, function(feature_signal_data){
    lapply(feature_signal_data, function(sig_obj){
      sig_obj = object@signal_data[[1]][[1]]
      signal_dt = sig_obj@signal_data
      
      reads_dt = dcast(signal_dt, id+x~name_split, value.var = "y")
      mat = as.matrix(reads_dt[,-1:-2])
      rownames(mat) = paste(reads_dt$id, reads_dt$x)
      cor_mat = cor(mat)
      row_clust = hclust(dist(cor_mat))
      col_clust = hclust(dist(t(cor_mat)))
      mat_dt = as.data.table(reshape2::melt(cor_mat))
      setnames(mat_dt, c("Var1", "Var2"), c("row_name_split", "col_name_split"))
      # mat_dt$row_name_split = factor(mat_dt$row_name_split, levels = rownames(cor_mat)[row_clust$order])
      # mat_dt$col_name_split = factor(mat_dt$col_name_split, levels = colnames(cor_mat)[col_clust$order])
      return(list(mat = cor_mat, dt = mat_dt, row_hclust = row_clust, col_hclust = col_clust))
    })
  })
  object@other_data$read_count_correlation = Correlation_data
  object@other_data$signal_profile_correlation = Correlation_data.profile
  object
})
setMethod("ssvQC.prepCorrelation", "ssvQC.featureOnly", function(object){
  stop("Cannot run prepCorrelation on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.prepCorrelation", "ssvQC.signalOnly", function(object){
  stop("Cannot run prepCorrelation on ssvQC with no QcConfigFeature component")
})

#' @export
#' @rdname ssvQC
#' @import pheatmap
setGeneric("ssvQC.plotCorrelation", function(object){standardGeneric("ssvQC.plotCorrelation")})
setMethod("ssvQC.plotCorrelation", "ssvQC.complete", function(object){
  if(is.null(object@other_data$read_count_correlation) || is.null(object@other_data$signal_profile_correlation)){
    object = ssvQC.prepCorrelation(object)
  }
  
  wrap_corr_plots = function(corr_res, main_title){
    brks = seq(-1, 1, by = .05)
    
    cr = colorRamp(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))
    
    color = rgb(cr(scales::rescale(brks, c(0,1)))/255)
    p_pheat = pheatmap::pheatmap(corr_res$mat, main = main_title,
                       cluster_rows = corr_res$row_hclust, 
                       cluster_cols = corr_res$col_hclust, 
                       scale = "none", 
                       breaks = brks, 
                       color = color, 
                       silent = TRUE)[[4]]
    p_pheat = ggplotify::as.ggplot(p_pheat)
    
    p_dt = corr_res$dt
    p_dt
    p_dt[, label := format(round(value, digits = 2))]
    if(is.factor(p_dt$row_name_split)){
      p_dt$row_name_split = factor(p_dt$row_name_split, rev(levels(p_dt$row_name_split)))  
    }
    
    p_gg = ggplot(corr_res$dt, aes(x = col_name_split, y = row_name_split, fill = value, label = label)) +
      geom_tile() +
      geom_text() +
      scale_fill_gradientn(colors = color, limits = range(brks)) +
      labs(title = main_title, fill = "pearson correlation", subtitle = "correlation of read count at assessed peaks") +
      theme(panel.background = element_blank(), panel.grid = element_blank())
    
    # cowplot::plot_grid(p_pheat, p_gg)
    list(pheatmap = p_pheat, ggplot_heatmap = p_gg)
  }
  
  if(is.null(object@plots$correlation)){
    object@plots$correlation = list()
  }
  
  Correlation_data = object@other_data$read_count_correlation
  plots = dbl_lapply(Correlation_data, wrap_corr_plots, dbl_names)
  object@plots$correlation$read_count_pheatmap = dbl_extract(plots, "pheatmap")
  object@plots$correlation$read_count_ggplot_heatmap = dbl_extract(plots, "ggplot_heatmap")
 
  Correlation_data = object@other_data$signal_profile_correlation
  plots = dbl_lapply(Correlation_data, wrap_corr_plots, dbl_names)
  object@plots$correlation$signal_profile_pheatmap = dbl_extract(plots, "pheatmap")
  object@plots$correlation$signal_profile_ggplot_heatmap = dbl_extract(plots, "ggplot_heatmap")
  
  object
})
setMethod("ssvQC.plotCorrelation", "ssvQC.featureOnly", function(object){
  stop("Cannot run plotCorrelation on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.plotCorrelation", "ssvQC.signalOnly", function(object){
  stop("Cannot run plotCorrelation on ssvQC with no QcConfigFeature component")
})

##SCC
#' @export
#' @rdname ssvQC
setGeneric("ssvQC.prepSCC", function(object){standardGeneric("ssvQC.prepSCC")})
setMethod("ssvQC.prepSCC", "ssvQC.complete", function(object){
  SCC_data = lapply(object@features_config$assessment_features, function(query_gr){
    sig_configs = .make_query_signal_config(object@signal_config)
    lapply(sig_configs, function(sel_sig_config){
      make_scc_dt(as.data.table(sel_sig_config@meta_data), query_gr = query_gr, bfc = object@bfc)
    })
  })
  object@other_data$SCC = SCC_data
  object
})
setMethod("ssvQC.prepSCC", "ssvQC.featureOnly", function(object){
  stop("Cannot run prepSCC on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.prepSCC", "ssvQC.signalOnly", function(object){
  stop("Cannot run prepSCC on ssvQC with no QcConfigFeature component")
})

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
setGeneric("ssvQC.plotSCC", function(object){standardGeneric("ssvQC.plotSCC")})
setMethod("ssvQC.plotSCC", "ssvQC.complete", function(object){
  if(is.null(object@other_data$SCC)){
    object = ssvQC.prepSCC(object)
  }
  SCC_data = object@other_data$SCC
  
  wrap_plot_scc_dt = function(scc_dt, main_title){
    plot_scc_dt(scc_dt, main_title, name_lev = levels(object@signal_config@meta_data$name_split))
  }
  
  SCC_plots = dbl_lapply(SCC_data, FUN = wrap_plot_scc_dt, FUN_names = dbl_names)
  SCC_dots = dbl_extract(SCC_plots, "scc_dots")
  SCC_curves = dbl_extract(SCC_plots, "scc_curves")
  
  if(is.null(object@plots$SCC)){
    object@plots$SCC = list()
  }
  
  object@plots$SCC$dots = SCC_dots
  object@plots$SCC$curves = SCC_curves
  object
})
setMethod("ssvQC.plotSCC", "ssvQC.featureOnly", function(object){
  stop("Cannot run plotSCC on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.plotSCC", "ssvQC.signalOnly", function(object){
  stop("Cannot run plotSCC on ssvQC with no QcConfigFeature component")
})

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

##Signal
#' @export
#' @rdname ssvQC
setGeneric("ssvQC.prepSignal", function(object){standardGeneric("ssvQC.prepSignal")})
setMethod("ssvQC.prepSignal", "ssvQC.complete", function(object){
  object = ssvQC.prepFetch(object)
  object@signal_data = lapply(object@features_config$assessment_features, function(query_gr){
    sig_configs = .make_query_signal_config(object@signal_config)
    lapply(sig_configs, function(sel_sig_config){
      ClusteredSignal.fromConfig(sel_sig_config, 
                                 resize(query_gr, object@signal_config@view_size, fix = "center"), 
                                 facet_var = "name_split", 
                                 extra_var = union(object@signal_config@color_by, object@signal_config@run_by), 
                                 bfc = object@bfc)
    })
  })
  object
})
setMethod("ssvQC.prepSignal", "ssvQC.featureOnly", function(object){
  stop("Cannot run prepSignal on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.prepSignal", "ssvQC.signalOnly", function(object){
  stop("Cannot run prepSignal on ssvQC with no QcConfigFeature component")
})

#' @export
#' @rdname ssvQC
setGeneric("ssvQC.plotSignal", function(object){standardGeneric("ssvQC.plotSignal")})
setMethod("ssvQC.plotSignal", "ssvQC.complete", function(object){
  if(length(object@signal_data) == 0){
    object = ssvQC.prepSignal(object)
  }
  signal_data = object@signal_data
  sig_config = object@signal_config
  
  wrap_plot_signal_dt = function(clust_sig, main_title = NULL){
    is_bam = grepl("bam", sig_config@read_mode)
    value_label = ifelse(is_bam, 
                         val2lab[sig_config@plot_value], 
                         val2bwlab[sig_config@plot_value])
    x_label =  paste(sig_config@view_size, "bp view")
    extra_vars = unique(c(sig_config@color_by, sig_config@run_by, "name_split"))
    extra_vars = intersect(extra_vars, colnames(clust_sig@signal_data))
    p_heatmap = seqsetvis::ssvSignalHeatmap.ClusterBars(clust_sig@signal_data, 
                                                        fill_ = val2var[sig_config@plot_value],
                                                        facet_ = "name_split", 
                                                        fill_limits = sig_config@heatmap_limit_values,
                                                        rel_widths = c(1, 20),
                                                        FUN_format_heatmap = function(p){
                                                          p + labs(x = x_label, fill = value_label, title = main_title)
                                                        })
    clust_sig.agg = clust_sig@signal_data[, .(y = mean(y), y_RPM = mean(y_RPM), y_linQ = mean(y_linQ), y_RPM_linQ = mean(y_RPM_linQ)), 
                                          c("x", extra_vars)]
    clust_sig.agg_per_cluster = clust_sig@signal_data[, .(y = mean(y), y_RPM = mean(y_RPM), y_linQ = mean(y_linQ), y_RPM_linQ = mean(y_RPM_linQ)), 
                                                      c("x", "cluster_id", extra_vars)]
    
    p_line = ggplot(clust_sig.agg, 
                    aes_string(x = "x", 
                               y = val2var[sig_config@plot_value], 
                               color = sig_config@color_by, 
                               group = "name_split")) +
      geom_path() +
      facet_grid(paste0(".~", sig_config@run_by)) +
      labs(x = x_label, 
           y = value_label, 
           subtitle = "mean at assessed features", 
           title = main_title) +
      scale_color_manual(values = sig_config@color_mapping)
    
    p_heatmap.line = ggplot(clust_sig.agg_per_cluster, 
                            aes_string(x = "x", 
                                       y = val2var[sig_config@plot_value], 
                                       color = sig_config@color_by, 
                                       group = "name_split")) +
      geom_path() +
      facet_grid(paste0("cluster_id~", sig_config@run_by), 
                 scales = ifelse(sig_config@lineplot_free_limits, "free_y", "fixed")) +
      labs(x = x_label, 
           y = value_label,
           subtitle = "mean per cluster", 
           title = main_title)+
      scale_color_manual(values = sig_config@color_mapping)
    
    return(list(
      heatmap = p_heatmap,
      heatmap.lines = p_heatmap.line,
      lines = p_line
    ))
  }
  signal_plots = dbl_lapply(signal_data, FUN = wrap_plot_signal_dt, FUN_names = dbl_names)
  signal_heatmaps = dbl_extract(signal_plots, "heatmap")
  signal_heatmaps.lines = dbl_extract(signal_plots, "heatmap.lines")
  signal_lines = dbl_extract(signal_plots, "lines")
  
  if(is.null(object@plots$signal)){
    object@plots$signal = list()
  }
  
  object@plots$signal$heatmaps = signal_heatmaps
  object@plots$signal$heatmaps.lines = signal_heatmaps.lines
  object@plots$signal$lines = signal_lines
  object
})
setMethod("ssvQC.plotSignal", "ssvQC.featureOnly", function(object){
  stop("Cannot run plotSignal on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.plotSignal", "ssvQC.signalOnly", function(object){
  stop("Cannot run plotSignal on ssvQC with no QcConfigFeature component")
})

### Features
#' @export
#' @rdname ssvQC
setGeneric("ssvQC.prepFeatures", function(object, bfc){standardGeneric("ssvQC.prepFeatures")})
setMethod("ssvQC.prepFeatures", "ssvQC.complete", function(object){
  object@features_config = ssvQC.prepFeatures(object@features_config, object@bfc)
  object
})
setMethod("ssvQC.prepFeatures", "ssvQC.featureOnly", function(object){
  object@features_config = prepFeatures(object@features_config)
  object
})
setMethod("ssvQC.prepFeatures", "ssvQC.signalOnly", function(object){
  stop("Cannot run prepFeatures on ssvQC with no QcConfigFeature component")
})
setMethod("ssvQC.prepFeatures", c("QcConfigFeatures", "BiocFileCache"), function(object, bfc){
  prepFeatures(object, bfc)
})
setMethod("ssvQC.prepFeatures", "QcConfigFeatures", function(object){
  prepFeatures(object)
})

.plotFeatures = function(object, force_euler = FALSE){
  feat_config = object@features_config
  need_prep_features = 
    length(feat_config$loaded_features) == 0 |
    length(feat_config$overlapped_features) == 0 |
    length(feat_config$assessment_features) == 0
  if(need_prep_features){
    object = ssvQC.prepFeatures(object)
    feat_config = object@features_config
  }
  
  feature_plots =lapply(names(feat_config$loaded_features), function(feat_nam){
    feat_label = sub("_features", " features", feat_nam)
    peak_grs = feat_config$loaded_features[[feat_nam]]
    peak_dt = data.table(N = lengths(peak_grs), name_split = names(peak_grs))
    peak_dt = merge(feat_config@meta_data, peak_dt, by = "name_split")
    
    out = list()
    
    p_peak_count = ggplot(peak_dt, aes_string(x = "name_split", y = "N", fill = feat_config@color_by)) +
      geom_bar(stat = "identity", color = "black") +
      scale_fill_manual(values = feat_config@color_mapping) +
      labs(x = "", y = "feature count", title = feat_label)
    
    olap_gr = feat_config$overlapped_features[[feat_nam]]
    n_feature_sets = ncol(mcols(olap_gr))
    
    p_binary_heatmap = seqsetvis::ssvFeatureBinaryHeatmap(olap_gr, raster_approximation = TRUE)
    p_upset = seqsetvis::ssvFeatureUpset(olap_gr)
    
    if(n_feature_sets < 4){
      p_venn = seqsetvis::ssvFeatureVenn(olap_gr, circle_colors = feat_config@color_mapping[peak_dt[[feat_config@color_by]]])  
      
    }else{
      venn_msg = "Venn diagrams not supported for more than 3 feature sets."
      p_venn = ggplot() + theme_void() + labs(title = venn_msg)
      message(venn_msg)
    }
    
    if(n_feature_sets < 9 | force_euler){
      p_euler = seqsetvis::ssvFeatureEuler(olap_gr, circle_colors = feat_config@color_mapping[peak_dt[[feat_config@color_by]]])  
      
    }else{
      euler_msg = "Euler diagrams not generated for more than 8 feature sets by default due to slow computation speed.  You may override this behavior by calling ssvQC.plotFeatures with force_euler = TRUE."
      p_euler = ggplot() + theme_void() + labs(title = euler_msg)
      message(euler_msg)
    }
    
    out$peak_count = p_peak_count
    out$binary_heatmap = p_binary_heatmap
    out$UpSet = p_upset
    out$venn = p_venn
    out$euler = p_euler
    
    out
  })
  names(feature_plots) = names(feat_config$loaded_features)
  
  if(is.null(object@plots$features)){
    object@plots$features = list()
  }
  object@plots$features$count = single_extract(feature_plots, "peak_count")
  object@plots$features$binary_heatmap = single_extract(feature_plots, "binary_heatmap")
  object@plots$features$UpSet = single_extract(feature_plots, "UpSet")
  object@plots$features$venn = single_extract(feature_plots, "venn")
  object@plots$features$euler = single_extract(feature_plots, "euler")
  
  object
}


#' @param force_euler If TRUE forces Euler plots to be generated for a list of feature sets longer than 8.  Euler plots can take quite a long time to generate as more feature sets are generated.
#' @export
#' @rdname ssvQC
setGeneric("ssvQC.plotFeatures", function(object, force_euler){standardGeneric("ssvQC.plotFeatures")})

setMethod("ssvQC.plotFeatures", c("ssvQC.complete"), function(object){
  ssvQC.plotFeatures(object, force_euler = FALSE)
})
setMethod("ssvQC.plotFeatures", c("ssvQC.complete", "logical"), .plotFeatures)
setMethod("ssvQC.plotFeatures", "ssvQC.signalOnly", function(object){
  stop("Cannot run plotFeatures on ssvQC with no QcConfigFeature component")
})
setMethod("ssvQC.plotFeatures", c("ssvQC.featureOnly"), function(object){
  ssvQC.plotFeatures(object, force_euler = FALSE)
})
setMethod("ssvQC.plotFeatures", c("ssvQC.featureOnly", "logical"), .plotFeatures)



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