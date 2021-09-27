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
    }else{
      p_dt$row_name_split = factor(p_dt$row_name_split, rev(sort(unique(p_dt$row_name_split))))  
    }
    p_gg = ggplot(p_dt, aes(x = col_name_split, y = row_name_split, fill = value, label = label)) +
      geom_tile() +
      geom_text() +
      scale_fill_gradientn(colors = color, limits = range(brks)) +
      labs(title = main_title, 
           fill = "pearson correlation", 
           subtitle = "correlation of read count at assessed peaks", 
           x = "", y = "") +
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