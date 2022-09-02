setClass("ssvTSNE", 
         representation = list(
           perplexity = "numeric",
           n_glyphs_x = "numeric",
           n_glyphs_y = "numeric",
           n_heatmap_pixels_x = "numeric",
           n_heatmap_pixels_y = "numeric"
         ),
         contains = "ssvQC.complete")

#' ssvTSNE
#'
#' @param features_config 
#' @param signal_config 
#' @param out_dir 
#' @param bfc 
#' @param matched_only 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
#' options(mc.cores = 1)
#' set.seed(0)
#' features_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
#' features_config = QcConfigFeatures.parse(features_config_file)
#'
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' bam_config = QcConfigSignal.parse(bam_config_file)
#' bam_config$view_size = 600
#' 
#' # TSNE plot paramters can be controlled at ssvTSNE object creation or set later
#' # if the y settings aren't specified they are the same as x
#' sts = ssvTSNE(
#'   features_config, 
#'   bam_config, 
#'   n_glyphs_x = 3, 
#'   n_heatmap_pixels_x = 5)
#' sts = ssvQC.plotFeatures(sts)
#' sts$perplexity = 10
#' sts = ssvQC.prepFetch(sts)
#' sts = ssvQC.referenceUsesSameScale(sts)
#' sts = ssvQC.prepSignal(sts)
#' sts = ssvQC.plotSignal(sts)
#' sts$plots$TSNE
#' 
#' sts$n_glyphs_x = 7 
#' sts$n_glyphs_y = 5
#' 
#' sts$n_heatmap_pixels_x = 30
#' sts$n_heatmap_pixels_y = 60
#' sts.replot = ssvQC.plotSignal(sts)
#' 
ssvTSNE = function(features_config = NULL,
                   signal_config = NULL,
                   out_dir = getwd(),
                   bfc = NULL, 
                   matched_only = TRUE,
                   n_glyphs_x = 10,
                   n_glyphs_y = n_glyphs_x,
                   n_heatmap_pixels_x = 25,
                   n_heatmap_pixels_y = n_heatmap_pixels_x,
                   ...){
  if(is.null(features_config) & is.null(signal_config)){
    stop("At least one of features_config or signal_config must be specified.")
  }
  
  features_config = .prep_features_config(features_config)
  signal_config = .prep_signal_config(signal_config)
  
  signal_config@cluster_value = SQC_SIGNAL_VALUES$linearQuantile
  signal_config@sort_value = SQC_SIGNAL_VALUES$linearQuantile
  signal_config@plot_value = SQC_SIGNAL_VALUES$linearQuantile
  
  if(is.null(bfc)){
    bfc = new_cache()
  }
  
  dir.create(out_dir, showWarnings = FALSE)
  
  new("ssvTSNE",
      features_config = features_config,
      signal_config = signal_config,
      signal_data = list(),
      other_data = list(),
      out_dir = out_dir,
      bfc = bfc,
      saving_enabled = TRUE,
      n_glyphs_x = n_glyphs_x,
      n_glyphs_y = n_glyphs_y,
      n_heatmap_pixels_x = n_heatmap_pixels_x,
      n_heatmap_pixels_y = n_heatmap_pixels_y,
      matched_only = matched_only)
}

setMethod("ssvQC.runAll", "ssvTSNE", function(object){
  object = callNextMethod()
  message("NYI")
  object
})

setMethod("ssvQC.prepSignal", "ssvTSNE", function(object){
  object = callNextMethod()
  object@signal_data = lapply(object@signal_data, function(signal_clust_objs){
    lapply(signal_clust_objs, function(clust_obj){
      ClusteredSignal_TSNE.from_ClusteredSignal(clust_obj, object)
    })
  })
  object
})

setMethod("ssvQC.plotSignal", "ssvTSNE", function(object){
  object = callNextMethod()
  
  n_glyphs_x = object@n_glyphs_x
  n_glyphs_y = object@n_glyphs_y
  
  n_heatmap_pixels_x = object@n_heatmap_pixels_x
  n_heatmap_pixels_y = object@n_heatmap_pixels_y
  
  object@plots$TSNE$regional_glyphs = lapply(object@signal_data, function(signal_data_groups){
    lapply(signal_data_groups, function(signal_data){
      if(!"ClusteredSignal_TSNE" %in% class(signal_data)){
        stop("Stored signal data is not of class ClusteredSignal_TSNE.")
      }
      p_summary_profiles = stsPlotSummaryProfiles(signal_data@signal_data,
                                                  signal_data@xy_data,
                                                  x_points = n_glyphs_x,
                                                  y_points = n_glyphs_y,
                                                  y_var = val2var[object@signal_config@plot_value])
      p_summary_profiles
    })
  })
  
  object@plots$TSNE$cluster_glyphs = lapply(object@signal_data, function(signal_data_groups){
    lapply(signal_data_groups, function(signal_data){
      if(!"ClusteredSignal_TSNE" %in% class(signal_data)){
        stop("Stored signal data is not of class ClusteredSignal_TSNE.")
      }
      
      use_tsne_clust = TRUE
      if(use_tsne_clust){
        cluster_result = ssvQC:::nn_clust(signal_data@xy_data)
        p_cluster_profiles = stsPlotClusterProfiles(signal_data@signal_data, cluster_result$data)
      }else{
        p_cluster_profiles = stsPlotClusterProfiles(signal_data@signal_data,
                                                    merge(
                                                      signal_data@signal_data,
                                                      signal_data@xy_data, by = "id"))
      }
      p_cluster_profiles
    })
  })
  
  object@plots$TSNE$regional_heatmap = lapply(object@signal_data, function(signal_data_groups){
    lapply(signal_data_groups, function(signal_data){
      if(!"ClusteredSignal_TSNE" %in% class(signal_data)){
        stop("Stored signal data is not of class ClusteredSignal_TSNE.")
      }
      prof_dt = signal_data@signal_data
      xy_dt = signal_data@xy_data
      
      x_var = "x" 
      y_var = "y" 
      id_var = "id"
      wide_var = "name" 
      tall_var = "tall_none"
      agg_FUN = max
      
      y_var = val2var[object@signal_config@plot_value]
      y_lab = names(y_var)
      
      xy_dt[, xbin := bin_values(tx, n_heatmap_pixels_x, c(-.5, .5))]
      xy_dt[, ybin := bin_values(ty, n_heatmap_pixels_y, c(-.5, .5))]
      
      heat_dt = merge(prof_dt, xy_dt[, .(id, xbin, ybin)], by = "id")[, .(y = agg_FUN(get(y_var))), c(id_var, wide_var, "xbin", "ybin")]
      
      
      ggplot(heat_dt, aes(x = xbin, y = ybin, fill = y)) + 
        geom_raster() +
        facet_wrap(~name) +
        scale_fill_viridis_c()
    })
  })
  
  
  saveRDS(object, "dev_object.ssvQC.plotSignal.Rds")
  object
})


### $ Accessor
setMethod("names", "ssvTSNE",
          function(x)
          {
            c("plots", "signal_data", "signal_config", "features_config", "SCC", "FRIP", "correlation",
              "perplexity",
              "n_glyphs_x",
              "n_glyphs_y",
              "n_heatmap_pixels_x",
              "n_heatmap_pixels_y")
            
          })


setMethod("$", "ssvTSNE",
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
                    signal_config = x@signal_config,
                    perplexity = x@perplexity,
                    n_glyphs_x = x@n_glyphs_x,
                    n_glyphs_y = x@n_glyphs_y,
                    n_heatmap_pixels_x = x@n_heatmap_pixels_x,
                    n_heatmap_pixels_y = x@n_heatmap_pixels_y
                    
            )
          })

setReplaceMethod("$", "ssvTSNE",
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
                           perplexity = {
                             x@perplexity = value
                           },
                           n_glyphs_x = {
                             x@n_glyphs_x = value
                           },
                           n_glyphs_y = {
                             x@n_glyphs_y = value
                           },
                           n_heatmap_pixels_x = {
                             x@n_heatmap_pixels_x = value
                           },
                           n_heatmap_pixels_y = {
                             x@n_heatmap_pixels_y = value
                           },
                           {warning(warn_msg)}
                           
                   )
                   
                   #TODO, some assignments may be appropriate
                   x
                 })
