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
