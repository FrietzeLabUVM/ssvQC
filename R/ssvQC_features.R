### Features
#' @export
#' @rdname ssvQC
setGeneric("ssvQC.prepFeatures", function(object, bfc){standardGeneric("ssvQC.prepFeatures")})
setMethod("ssvQC.prepFeatures", "ssvQC.complete", function(object){
  object@features_config = ssvQC.prepFeatures(object@features_config, object@bfc)
  object
})
setMethod("ssvQC.prepFeatures", "ssvQC.featureOnly", function(object){
  object@features_config = .prepFeatures(object@features_config)
  object
})
setMethod("ssvQC.prepFeatures", "ssvQC.signalOnly", function(object){
  stop("Cannot run prepFeatures on ssvQC with no QcConfigFeature component")
})
setMethod("ssvQC.prepFeatures", c("QcConfigFeatures", "BiocFileCache"), function(object, bfc){
  .prepFeatures(object, bfc)
})
setMethod("ssvQC.prepFeatures", "QcConfigFeatures", function(object){
  .prepFeatures(object)
})


#' @param object 
#'
#' @return
#' @rdname QcConfigFeatures
#' @examples
.prepFeatures = function(object, bfc = new_cache()){
  if(object@is_null) return(object)
  
  if(length(object$assessment_features) > 0){
    # message("Not prepping features.  Call resetFeatures() on object to reset before prepping.")
    return(object)
  }
  
  to_run = object@to_run
  for(tr in to_run){
    tr_name = paste0(tr, "_features")
    rb = object@meta_data[[object@run_by]]
    sel_dt = object@meta_data[rb %in% union(tr, object@to_run_reference),]
    if(is.null(object@loaded_features[[tr_name]])){
      object@loaded_features[[tr_name]] = .process_features(sel_dt, object@feature_load_FUN, bfc = bfc)
    }
    is_error = sapply(object@loaded_features[[tr_name]], is, class2 = "try-error")
    if(any(is_error)){
      stop("Failed to load some feature files:\n  ",
           paste(sel_dt$file[is_error], collapse = "\n  "),
           "\nPlease supply appropriate feature_load_FUN when creating QcConfigFeatures.")
    }
    
    if(is.null(object@overlap_gr[[tr_name]])){
      object@overlap_gr[[tr_name]] = .process_overlaps(
        loaded_features = object@loaded_features[[tr_name]], 
        overlap_extension = object@overlap_extension,
        bfc = bfc)
    }
    
    if(length(object@overlap_gr[[tr_name]]) < 5){
      stop("Currently, fewer than 5 regions are not allowed. This will be addressed in future versions.")
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
    
    if(length(object@assessment_gr[[tr_name]]) < 5){
      stop("Currently, fewer than 5 regions are not allowed. This will be addressed in future versions.")
    }
  }
  object
}

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
    peak_dt = merge(as.data.table(feat_config@meta_data), peak_dt, by = "name_split")
    peak_dt$name_split = factor(peak_dt$name_split, levels = levels(feat_config@meta_data$name_split))
    
    setkey(peak_dt, "name_split")
    out = list()
    
    p_peak_count = ggplot(peak_dt, aes_string(x = "name_split", y = "N", fill = feat_config@color_by)) +
      geom_bar(stat = "identity", color = "black") +
      scale_fill_manual(values = feat_config@color_mapping) +
      labs(x = "", y = "feature count", title = feat_label) +
      scale_y_continuous(labels = function(x)paste(x/1e3, " k"))
    
    olap_gr = feat_config$overlapped_features[[feat_nam]]
    #colors assigned by color_by in order of GRanges names
    name_o_cols = feat_config@color_mapping[peak_dt[.(colnames(mcols(olap_gr)))][[feat_config@color_by]]]
    
    n_feature_sets = ncol(mcols(olap_gr))
    
    p_binary_heatmap = seqsetvis::ssvFeatureBinaryHeatmap(olap_gr, raster_approximation = TRUE)
    p_upset = seqsetvis::ssvFeatureUpset(olap_gr)
    
    if(n_feature_sets < 4){
      p_venn = seqsetvis::ssvFeatureVenn(olap_gr, circle_colors = name_o_cols)  
      
    }else{
      venn_msg = "Venn diagrams not supported for more than 3 feature sets."
      p_venn = ggplot() + theme_void() + labs(title = venn_msg)
      message(venn_msg)
    }
    
    if(n_feature_sets < 9 | force_euler){
      colnames(mcols(olap_gr))
      
      p_euler = seqsetvis::ssvFeatureEuler(olap_gr, circle_colors = name_o_cols)  
      
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

