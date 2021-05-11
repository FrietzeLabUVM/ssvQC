setOldClass(c("theme", "gg"))

setClass("ssvQC",
         representation = list(
           feature_config = "QcConfigFeatures",
           signal_config = "QcConfigSignal",
           signal_data = "list",
           other_data = "list",
           out_dir = "character",
           bfc = "BiocFileCache",
           saving_enabled = "logical",
           plot_post_FUN = "function",
           plots = "list"
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


ssvQC.save_config = function(object, file){
  feature_file = paste0(sub(".csv", "",  file), ".features.csv")
  signal_file = paste0(sub(".csv", "",  file), ".signal.csv")
  QcConfigFeatures.save_config(object@feature_config, feature_file)
  QcConfigSignal.save_config(object@signal_config, signal_file)
  invisible(c(features = feature_file, signal = signal_file))
}

#' ssvQC
#'
#' @param feature_config Controls features configuration.  May be a:
#'   QcConfigFeatures object, path to a file defining configuration via
#'   QcConfigFeatures.parse, features files to define via
#'   QcConfigFeatures.files, or a data.frame to pass to QcConfigFeatures.
#' @param signal_config Controls signal configuration.  May be a: QcConfigSignal
#'   object, path to a file defining configuration via QcConfigSignal.parse,
#'   features files to define via QcConfigSignal.files, or a data.frame to pass
#'   to QcConfigSignal.
#' @param out_dir NYI
#' @param bfc BiocFileCache object to use for caching. If NULL, default
#'   BiocFileCache::BiocFileCache() will be used.
#'
#' @return A ssvQC object.  Data needs to be loaded after via ssvQC.runAll or sub-methods ssvQC.plot*.
#' @export
#' @import BiocFileCache
#' 
#' @rdname ssvQC
#' @examples
#' options(mc.cores = 10)
#' set.seed(0)
#' feature_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
#' feature_config = QcConfigFeatures.parse(feature_config_file)
#'
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' bam_config = QcConfigSignal.parse(bam_config_file)
#'
#' bigwig_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bigwig_config.csv")
#' bigwig_config = QcConfigSignal.parse(bigwig_config_file)
#'
#' sqc.complete.file = ssvQC(feature_config_file, bam_config_file)
#'
#' sqc.complete = ssvQC(feature_config, bam_config)
#' object = sqc.complete
#'
#' sqc.complete.bw = ssvQC(feature_config, bigwig_config_file)
#'
#' sqc.signal = ssvQC(signal_config = bam_config)
#'
#' sqc.feature = ssvQC(feature_config = feature_config)
#'
#' sqc.signal = ssvQC.runAll(sqc.signal)
#' sqc.feature = ssvQC.runAll(sqc.feature)
#'
#' sqc.complete = ssvQC.runAll(sqc.complete)
ssvQC = function(feature_config = NULL,
                 signal_config = NULL,
                 out_dir = getwd(),
                 bfc = NULL, 
                 ...){
  if(is.null(feature_config) & is.null(signal_config)){
    stop("At least one of feature_config or signal_config must be specified.")
  }
  
  if(!is.null(feature_config)){
    if(is.character(feature_config)){
      if(!any(is_feature_file(feature_config))){
        feature_config = QcConfigFeatures.parse(feature_config)
      }else{
        feature_config = QcConfigFeatures.files(feature_config)
      }
    }else if(is.data.frame(feature_config)){
      feature_config = QcConfigFeatures(feature_config, ...)
    }
    if(!"QcConfigFeatures" %in% class(feature_config)){
      stop("feature_config must be either a QcConfigFeatures object or the path to valid configuration file to create one.")
    }  
    stopifnot(file.exists(feature_config@meta_data$file))
  }
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
  
  if(is.null(bfc)){
    bfc = BiocFileCache::BiocFileCache()
  }
  
  dir.create(out_dir, showWarnings = FALSE)
  
  if(!is.null(feature_config) & !is.null(signal_config)){
    new("ssvQC.complete",
        feature_config = feature_config,
        signal_config = signal_config,
        signal_data = list(),
        other_data = list(),
        out_dir = out_dir,
        bfc = bfc,
        saving_enabled = TRUE
    )
  }else if(!is.null(feature_config)){
    new("ssvQC.featureOnly",
        feature_config = feature_config,
        signal_config = QcConfigSignal.null(),
        signal_data = list(),
        other_data = list(),
        out_dir = out_dir,
        bfc = bfc,
        saving_enabled = TRUE
    )
  }else if(!is.null(signal_config)){
    new("ssvQC.signalOnly",
        feature_config = QcConfigFeatures.null(),
        signal_config = signal_config,
        signal_data = list(),
        other_data = list(),
        out_dir = out_dir,
        bfc = bfc,
        saving_enabled = TRUE
    )
  }else{
    stop("At least one of feature_config or signal_config must be specified. This should have been caught earlier.")
  }
  
  
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
  object = ssvQC.plotFeatures(object)
  
  object = ssvQC.prepMappedReads(object)
  object = ssvQC.prepFragLens(object)
  object = ssvQC.prepCapValue(object)
  
  object = ssvQC.plotMappedReads(object)
  object = ssvQC.plotSignal(object)
  object = ssvQC.plotSCC(object)
  object = ssvQC.plotFRIP(object)
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

#' @export
#' @rdname ssvQC
setGeneric("ssvQC.prepFragLens", function(object, query){standardGeneric("ssvQC.prepFragLens")})
setMethod("ssvQC.prepFragLens", "ssvQC.complete", function(object){
  message("complete")
  object@signal_config = ssvQC.prepFragLens(object@signal_config, object@feature_config)
  object
})
setMethod("ssvQC.prepFragLens", "ssvQC.featureOnly", function(object){
  stop("Cannot run prepCapValue on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.prepFragLens", "ssvQC.signalOnly", function(object){
  stop("Cannot run prepCapValue on ssvQC with no QcConfigFeatures component")
})
setMethod("ssvQC.prepFragLens", c("QcConfigSignal", "QcConfigFeatures"), function(object, query){
  if(object@read_mode != "bam_SE"){
    stop("ssvQC.prepFragLens only appropriate for read_mode bam_SE")
  }
  
  #bam specific independent of peaks
  
  sig_dt = as.data.table(object@meta_data)
  peak_dt = as.data.table(query@meta_data)
  
  matched_dt = merge(sig_dt[, .(bam_file = file, name)], peak_dt[, .(peak_file = file, name)], by = "name")
  unmatched_dt = sig_dt[!name %in% matched_dt$name]
 
  if(nrow(matched_dt) > 0){
    fl_dt.matched = rbindlist(lapply(seq_len(nrow(matched_dt)), function(i){
      peak_f = matched_dt[i,]$peak_file
      bam_f = matched_dt[i,]$bam_file
      name = matched_dt[i,]$name
      peak_gr = query@feature_load_FUN(peak_f)[[1]]
      # rname = digest::digest(list(bam_f, name, peak_gr, "ssvQC.prepFragLens"))
      fl = seqsetvis::fragLen_calcStranded(bam_f, peak_gr, n_regions = 500)
      
      dt = data.table(fragLens = fl)
      dt$name = name
      dt
    }))
  }else{
    fl_dt.matched = NULL
  }
  
  if(nrow(unmatched_dt) > 0){
    query = ssvQC.prepFeatures(query)
    peak_gr = unlist(GRangesList(query$assessment_features))
    
    fl_dt.unmatched = rbindlist(lapply(seq_len(nrow(unmatched_dt)), function(i){
      bam_f = unmatched_dt[i,]$file
      name = unmatched_dt[i,]$name
      # rname = digest::digest(list(bam_f, name, peak_gr, "ssvQC.prepFragLens"))
      fl = seqsetvis::fragLen_calcStranded(bam_f, peak_gr, n_regions = 500)
      
      dt = data.table(fragLens = fl)
      dt$name = name
      dt
    }))
  }else{
    fl_dt.unmatched = NULL
  }
  fl_dt = rbind(fl_dt.matched, fl_dt.unmatched)
  
  if(!setequal(fl_dt$name, object@meta_data$name)){
    stop("something has gone wrong assigning fragLens, pease report this issue at https://github.com/FrietzeLabUVM/ssvQC/issues") 
  }
  
  setkey(fl_dt, "name")
  object@meta_data$fragLens = fl_dt[.(object@meta_data$name)]$fragLens
  object
})

##MappedReads

#' @export
#' @rdname ssvQC
setGeneric("ssvQC.prepMappedReads", function(object){standardGeneric("ssvQC.prepMappedReads")})
setMethod("ssvQC.prepMappedReads", "ssvQC.complete", function(object){
  message("complete")
  object@signal_config = ssvQC.prepMappedReads(object@signal_config)
  object
})
setMethod("ssvQC.prepMappedReads", "ssvQC.featureOnly", function(object){
  stop("Cannot run mapped reads on ssvQC with no QcConfigSignal component")
  message("featureOnly")
})
setMethod("ssvQC.prepMappedReads", "ssvQC.signalOnly", function(object){
  message("signalOnly")
  object@signal_config = ssvQC.prepMappedReads(object@signal_config)
  object
})
setMethod("ssvQC.prepMappedReads", c("QcConfigSignal"), function(object){
  message("QcConfigSignal")
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
setGeneric("ssvQC.prepCapValue", function(object, query){standardGeneric("ssvQC.prepCapValue")})
setMethod("ssvQC.prepCapValue", "ssvQC.complete", function(object){
  message("complete")
  object@signal_config = ssvQC.prepCapValue(object@signal_config, object@feature_config)
  object
})
setMethod("ssvQC.prepCapValue", "ssvQC.featureOnly", function(object){
  stop("Cannot run prepCapValue on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.prepCapValue", "ssvQC.signalOnly", function(object){
  stop("Cannot run prepCapValue on ssvQC with no QcConfigFeatures component")
})
setMethod("ssvQC.prepCapValue", c("QcConfigSignal", "QcConfigFeatures"), function(object, query){
  message("QcConfigSignal")
  #bam specific independent of peaks
  
  sig_dt = as.data.table(object@meta_data)
  setkey(sig_dt, "name")
  peak_dt = as.data.table(query@meta_data)
  
  matched_dt = merge(sig_dt[, .(bam_file = file, name)], peak_dt[, .(peak_file = file, name)], by = "name")
  unmatched_dt = sig_dt[!name %in% matched_dt$name]
  
  if(nrow(matched_dt) > 0){
    cap_dt.matched = rbindlist(lapply(seq_len(nrow(matched_dt)), function(i){
      peak_f = matched_dt[i,]$peak_file
      bam_f = matched_dt[i,]$bam_file
      name_i = matched_dt[i,]$name
      peak_gr = query@feature_load_FUN(peak_f)[[1]]
      peak_gr = peak_gr[sample(min(5e3, length(peak_gr)))]
      
      if(!is.null(matched_dt$fragLens))
      fragLens = match_
      
      if(object@read_mode == "bam_SE"){
        fragLens = sig_dt[.(name_i)]$fragLens  
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
      if(!is.null(matched_dt$fragLens))
        fragLens = match_
      
      if(object@read_mode == "bam_SE"){
        fragLens = sig_dt[.(name_i)]$fragLens  
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
  object@meta_data$cap_value = cap_dt[.(object@meta_data$name)]$y_cap_value
  object
})

.ssvQC.plotMappedReads = function(object){
  todo = .make_query_signal_config(object@signal_config)
  plots = lapply(todo, function(sig_config){
    if(grepl("bam", sig_config@read_mode)){
      bam_config_dt = sig_config@meta_data
      
      if(is.null(bam_config_dt$mapped_reads)){
        sig_config = ssvQC.prepMappedReads(sig_config)
        object@signal_config = sig_config
        bam_config_dt = sig_config@meta_data
      }
      
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

##FRIP
#' @export
#' @rdname ssvQC
setGeneric("ssvQC.prepFRIP", function(object){standardGeneric("ssvQC.prepFRIP")})
setMethod("ssvQC.prepFRIP", "ssvQC.complete", function(object){
  FRIP_data = lapply(object@feature_config$assessment_features, function(query_gr){
    sig_configs = .make_query_signal_config(object@signal_config)
    lapply(sig_configs, function(sel_sig_config){
      make_frip_dt(as.data.table(sel_sig_config@meta_data), query_gr = query_gr, color_var = sel_sig_config@color_by)
    })
  })
  object@other_data$FRIP = FRIP_data
  object
})
setMethod("ssvQC.prepFRIP", "ssvQC.featureOnly", function(object){
  stop("Cannot run prepSignal on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.prepFRIP", "ssvQC.signalOnly", function(object){
  stop("Cannot run prepSignal on ssvQC with no QcConfigFeature component")
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
  stop("Cannot run prepSignal on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.plotFRIP", "ssvQC.signalOnly", function(object){
  stop("Cannot run prepSignal on ssvQC with no QcConfigFeature component")
})

##SCC
#' @export
#' @rdname ssvQC
setGeneric("ssvQC.prepSCC", function(object){standardGeneric("ssvQC.prepSCC")})
setMethod("ssvQC.prepSCC", "ssvQC.complete", function(object){
  SCC_data = lapply(object@feature_config$assessment_features, function(query_gr){
    sig_configs = .make_query_signal_config(object@signal_config)
    lapply(sig_configs, function(sel_sig_config){
      make_scc_dt(as.data.table(sel_sig_config@meta_data), query_gr = query_gr, bfc_corr = object@bfc)
    })
  })
  object@other_data$SCC = SCC_data
  object
})
setMethod("ssvQC.prepSCC", "ssvQC.featureOnly", function(object){
  stop("Cannot run prepSignal on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.prepSCC", "ssvQC.signalOnly", function(object){
  stop("Cannot run prepSignal on ssvQC with no QcConfigFeature component")
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
  stop("Cannot run prepSignal on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.plotSCC", "ssvQC.signalOnly", function(object){
  stop("Cannot run prepSignal on ssvQC with no QcConfigFeature component")
})

##Signal
#' @export
#' @rdname ssvQC
setGeneric("ssvQC.prepSignal", function(object){standardGeneric("ssvQC.prepSignal")})
setMethod("ssvQC.prepSignal", "ssvQC.complete", function(object){
  object@signal_data = lapply(object@feature_config$assessment_features, function(query_gr){
    sig_configs = .make_query_signal_config(object@signal_config)
    lapply(sig_configs, function(sel_sig_config){
      ClusteredSignal.fromConfig(sel_sig_config, 
                                 query_gr, 
                                 facet_var = "name_split", 
                                 extra_var = union(object@signal_config@color_by, object@signal_config@run_by))
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
  sig_dt = signal_data[[1]][[1]]
  sig_config = object@signal_config
  
  wrap_plot_signal_dt = function(sig_dt, main_title = NULL){
    is_bam = grepl("bam", sig_config@read_mode)
    value_label = ifelse(is_bam, "read\npileup", "bigWig\nsignal")
    x_label =  paste(sig_config@view_size, "bp view")
    extra_vars = unique(c(sig_config@color_by, sig_config@run_by, "name_split"))
    
    p_heatmap = seqsetvis::ssvSignalHeatmap.ClusterBars(sig_dt@signal_data, facet_ = "name_split", rel_widths = c(1, 20),
                                                        FUN_format_heatmap = function(p){
                                                          p + labs(x = x_label, fill = value_label, title = main_title)
                                                        })
    sig_dt.agg = sig_dt@signal_data[, .(y = mean(y)), c("x", extra_vars)]
    sig_dt.agg_per_cluster = sig_dt@signal_data[, .(y = mean(y)), c("x", "cluster_id", extra_vars)]
    
    p_line = ggplot(sig_dt.agg, aes_string(x = "x", y = "y", color = sig_config@color_by, group = "name_split")) +
      geom_path() +
      facet_grid(paste0(".~", sig_config@run_by)) +
      labs(x = x_label, y = value_label, subtitle = "mean at assessed features", title = main_title) +
      scale_color_manual(values = sig_config@color_mapping)
    
    p_heatmap.line = ggplot(sig_dt.agg_per_cluster, aes_string(x = "x", y = "y", color = sig_config@color_by, group = "name_split")) +
      geom_path() +
      facet_grid(paste0("cluster_id~", sig_config@run_by)) +
      labs(x = x_label, y = value_label, subtitle = "mean per cluster", title = main_title)+
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
  stop("Cannot run prepSignal on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.plotSignal", "ssvQC.signalOnly", function(object){
  stop("Cannot run prepSignal on ssvQC with no QcConfigFeature component")
})

### Features
#' @export
#' @rdname ssvQC
setGeneric("ssvQC.prepFeatures", function(object){standardGeneric("ssvQC.prepFeatures")})
setMethod("ssvQC.prepFeatures", "ssvQC.complete", function(object){
  object@feature_config = ssvQC.prepFeatures(object@feature_config)
  object
})
setMethod("ssvQC.prepFeatures", "ssvQC.featureOnly", function(object){
  object@feature_config = prepFeatures(object@feature_config)
  object
})
setMethod("ssvQC.prepFeatures", "ssvQC.signalOnly", function(object){
  stop("Cannot run prepSignal on ssvQC with no QcConfigFeature component")
})
setMethod("ssvQC.prepFeatures", "QcConfigFeatures", function(object){
  prepFeatures(object)
})

.plotFeatures = function(object, force_euler = FALSE){
  feat_config = object@feature_config
  if(length(feat_config$loaded_features) == 0 |
     length(feat_config$overlapped_features) == 0 |
     length(feat_config$assessment_features) == 0){
    object = ssvQC.prepFeatures(object)
    feat_config = object@feature_config
  }
  
  
  feature_plots =lapply(names(feat_config$loaded_features), function(feat_nam){
    feat_nam = names(feat_config$loaded_features)
    feat_label = sub("_features", " features", feat_nam)
    peak_grs = feat_config$loaded_features[[feat_nam]]
    peak_dt = data.table(N = lengths(peak_grs), name_split = names(peak_grs))
    peak_dt = merge(peak_dt, feat_config@meta_data, by = "name_split")
    
    out = list()
    
    p_peak_count = ggplot(peak_dt, aes_string(x = "name_split", y = "N", fill = feat_config@color_by)) +
      geom_bar(stat = "identity", color = "black") +
      scale_fill_manual(values = feat_config@color_mapping) +
      labs(x = "", y = "feature count", title = feat_label)
    
    olap_gr = feat_config$overlapped_features[[feat_nam]]
    n_feature_sets = ncol(mcols(olap_gr))
    
    p_binary_heatmap = seqsetvis::ssvFeatureBinaryHeatmap(olap_gr)
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
  stop("Cannot run prepSignal on ssvQC with no QcConfigFeature component")
})
setMethod("ssvQC.plotFeatures", c("ssvQC.featureOnly"), function(object){
  ssvQC.plotFeatures(object, force_euler = FALSE)
})
setMethod("ssvQC.plotFeatures", c("ssvQC.featureOnly", "logical"), .plotFeatures)


### $ Accessor
setMethod("names", "ssvQC",
          function(x)
          {
            c("plots", "signal_data", "signal_config", "features_config", "SCC", "FRIP")
            
          })


setMethod("$", "ssvQC",
          function(x, name)
          {
            switch (name,
                    plots = x@plots,
                    signal_data = x@signal_data,
                    SCC = x@other_data$SCC,
                    FRIP = x@other_data$FRIP,
                    bfc = x@bfc,
                    feature_config = x@feature_config,
                    signal_config = x@signal_config
                    
            )
          })

setReplaceMethod("$", "ssvQC",
                 function(x, name, value)
                 {
                   warn_msg = "This assignment is not supported.  No effect."
                   switch (name,
                           feature_config = {
                             x@feature_config = value
                           },
                           signal_config = {
                             x@signal_config = value
                           },
                           {warning(warn_msg)}
                           
                   )
                   
                   #TODO, some assignments may be appropriate
                   x
                 })