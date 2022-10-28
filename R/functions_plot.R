#' plot_fq_dt
#'
#' @param fq_dt output of \code{\link{make_fq_dt}}
#'
#' @return a ggplot shoing read counts from fastq files
#' @export
#'
#' @import ggplot2
#' @rawNamespace import(data.table, except = c(shift, first, second, last))
#'
#' @examples
#' fq_files = dir("inst/extdata",
#'   pattern = "(fq$)|(fq.gz$)|(fastq$)|(fastq.gz$)",
#'   full.names = TRUE)
#' #no idea why this make_fq_dt example won't run
#' #fq_dt = make_fq_dt(fq_files,
#' #   fastq_names = c("4_reads_fq", "4_reads_gz", "5_reads_fq", "5_reads_gz"))
#' #plot_fq_dt(fq_dt)
plot_fq_dt = function(fq_dt){
  name = count = treatment = NULL #global binding for data.table
  p_fq1 = ggplot(fq_dt,
                 aes(x = name, y = count, fill = treatment)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.margin = margin(.01, .01, .01, .1, unit = "npc")) +
    scale_y_continuous(labels = function(x)x/1e6) +
    labs(y = "FASTQ reads (M)", x = "")
  p_fq1
}

#' plot_feature_comparison
#'
#' @param peak_grs list of GRanges peak sets
#' @param min_fraction numeric fraction from 0 to 1. min_fraction of peaks for
#'   consensus. Default is 0.
#' @param min_number numeric from 2 to length of peak_gres. min_number of peaks
#'   for consesensus. default is 2.
#' @param force_euler logical. If TRUE, Euler plots will generate for more than
#'   8 peak sets, may be very slow! Default is FALSE
#'
#' @return Nested list of plots for peak overlaps. Outer list contains "all" and
#'   "consensus" sets of plots.  Inner list has a barplot, binary heatmap and
#'   upset
#' @export
#'
#' @examples
#' peak_files = dir(system.file("extdata", package = "seqqc"), pattern = "Peak$", full.names = TRUE)
#' peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files)
#' names(peak_grs) = sub("_rand.+", "", names(peak_grs))
#' peak_plots = plot_feature_comparison(peak_grs)
#'
#' #plotting looks like this:
#' peak_plots$all$venn
#'
#' #possible assembly of plots
#' cowplot::plot_grid(ncol = 1,
#'   cowplot::plot_grid(plotlist = peak_plots$all[1:3], nrow = 1),
#'   cowplot::plot_grid(plotlist = peak_plots$all[4:5], nrow = 1)
#' )
#'
#' cowplot::plot_grid(ncol = 1,
#'   cowplot::plot_grid(plotlist = peak_plots$consensus[1:3], nrow = 1),
#'   cowplot::plot_grid(plotlist = peak_plots$consensus[4:5], nrow = 1)
#' )
#'
plot_feature_comparison = function(peak_grs, min_fraction = 0, min_number = 2, force_euler = FALSE, peak_colors = NULL){
  olaps_all = seqsetvis::ssvOverlapIntervalSets(c(peak_grs))

  p_peak_counts_all = seqsetvis::ssvFeatureBars(peak_grs, show_counts = FALSE, counts_text_colors = "gray60", bar_colors = peak_colors) +
    guides(fill = "none") +
    scale_y_continuous(labels = function(x)x/1e3) +
    labs(y = "peak count (k)", title = "all peaks", subtitle = paste("counts:", formatC(max(lengths(peak_grs)), big.mark = ",", format = "d"), "max"), x = "")

  p_peak_overlaps_all = seqsetvis::ssvFeatureBinaryHeatmap(olaps_all, raster_approximation = TRUE) +
    labs(title = "all peaks", subtitle = paste("overlaps:", formatC(length(olaps_all), big.mark = ",", format = "d"), "regions"))

  p_peak_upset_all = seqsetvis::ssvFeatureUpset(olaps_all) +
    labs(title = "all peaks", subtitle = paste("overlaps:", formatC(length(olaps_all), big.mark = ",", format = "d"), "regions"))

  if(length(peak_grs) < 4){
    p_peak_venn_all = seqsetvis::ssvFeatureVenn(olaps_all, circle_colors = peak_colors) +
      labs(title = "all peaks", subtitle = paste("overlaps:", formatC(length(olaps_all), big.mark = ",", format = "d"), "regions"))
  }else{
    p_peak_venn_all = ggplot() + theme_void() + labs(title = "Can't run venn for more than 3 groups")
  }

  if(length(peak_grs) < 9 || force_euler){
    p_peak_euler_all = seqsetvis::ssvFeatureEuler(olaps_all, circle_colors = peak_colors) +
      labs(title = "all peaks", subtitle = paste("overlaps:", formatC(length(olaps_all), big.mark = ",", format = "d"), "regions"))
  }else{
    p_peak_euler_all = ggplot() + theme_void() + labs(title = "Can't run Euler for more than 8 groups.\nCan force with force_euler = TRUE.")
  }

  message("consensus peaks...")
  olaps_consensus = seqsetvis::ssvConsensusIntervalSets(peak_grs, min_fraction = min_fraction, min_number = max(2, min(min_number, length(peak_grs))))
  p_peak_counts_consenus = seqsetvis::ssvFeatureBars(olaps_consensus, bar_colors = "black", show_counts = FALSE, counts_text_colors = "gray60" ) +
    guides(fill = "none") +
    scale_y_continuous(labels = function(x)x/1e3) +
    labs(y = "peak count (k)", title = "consensus peaks only", subtitle = paste("counts:", formatC(max(colSums(as.data.frame(mcols(olaps_consensus)))), big.mark = ",", format = "d"), "max"))

  p_peak_overlaps_consensus = seqsetvis::ssvFeatureBinaryHeatmap(olaps_consensus, raster_approximation = TRUE) +
    labs(title = "consensus peaks only", subtitle = paste("overlaps:", formatC(length(olaps_consensus), big.mark = ",", format = "d"), "regions"))

  p_peak_upset_consensus = seqsetvis::ssvFeatureUpset(olaps_consensus) +
    labs(title = "consensus peaks only", subtitle = paste("overlaps:", formatC(length(olaps_consensus), big.mark = ",", format = "d"), "regions"))

  if(length(peak_grs) < 4){
    p_peak_venn_consensus = seqsetvis::ssvFeatureVenn(olaps_consensus) +
      labs(title = "consensus peaks only", subtitle = paste("overlaps:", formatC(length(olaps_consensus), big.mark = ",", format = "d"), "regions"))
  }else{
    p_peak_venn_consensus = ggplot() + theme_void() + labs(title = "Can't run venn for more than 3 groups")
  }

  if(length(peak_grs) < 9 || force_euler){
    p_peak_euler_consensus = seqsetvis::ssvFeatureEuler(olaps_consensus) +
      labs(title = "consensus peaks only", subtitle = paste("overlaps:", formatC(length(olaps_consensus), big.mark = ",", format = "d"), "regions"))
  }else{
    p_peak_euler_consensus = ggplot() + theme_void() + labs(title = "Can't run Euler for more than 8 groups.\nCan force with force_euler = TRUE.")
  }

  return(list(
    all = list(
      count = p_peak_counts_all,
      overlap = p_peak_overlaps_all,
      upset = p_peak_upset_all,
      venn = p_peak_venn_all,
      euler = p_peak_euler_all),
    consensus = list(
      count = p_peak_counts_consenus,
      overlap = p_peak_overlaps_consensus,
      upset = p_peak_upset_consensus,
      venn = p_peak_venn_consensus,
      euler = p_peak_euler_consensus)
  ))
}

#' plot_frip_dt
#' 
#' internal function used by ssvQC.plotFRIP
#'
#' @param frip_dt output from \code{\link{make_frip_dt}}
#' @param sort_by character. One of "frip" or "reads_in_peak". Should plots be
#'   sorted by decreasing FRIP ("frip") or total reads ("reads_in_peak")?
#'   Default is "frip".
#' @param name_lev character. Name levels to impose manual ordering.  sort_by is
#'   not ignored if this is not NULL. Default is NULL.
#'
#' @return list of ggplot plots relevant to FRIP.
#' @importFrom stats median quantile
#'
#' @examples
#' peak_files = dir(system.file("extdata", package = "seqqc"), pattern = "Peak$", full.names = TRUE)
#' peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files)
#' query_gr = resize(seqsetvis::ssvOverlapIntervalSets(peak_grs), 6e2, fix = "center")
#'
#' bam_files = dir(system.file("extdata", package = "seqqc"), pattern = "^M.+bam$", full.names = TRUE)
#' query_dt.bam = make_dt(bam_files)
#'
#' frip_dt = make_frip_dt(query_dt.bam, query_gr)
#' frip_plots = plot_frip_dt(frip_dt)
#' frip_plots$frip_per_peaks
#' frip_plots$frip_total
plot_frip_dt = function(frip_dt, sort_by = c("none", "frip", "reads_in_peak")[1], name_lev = NULL, name_var = "name_split", color_var = name_var, color_mapping = NULL, main_title = NULL){
  reads_in_peak = name = V1 = frip = fastq_count = treatment = mapped_reads = peak_count = NULL #global binding for data.table
  if(is.null(name_lev)){
    if(sort_by == "none"){
      if(is.factor(frip_dt[[name_var]])){
        name_lev = levels(frip_dt[[name_var]])  
      }else{
        name_lev = sort(frip_dt[[name_var]])
      }
    }else if(sort_by == "reads_in_peak"){
      name_lev = as.character(frip_dt[, stats::median(reads_in_peak) , c(name_var), with = FALSE][rev(order(V1))][[name_var]])
    }else if(sort_by == "frip"){
      name_lev = as.character(frip_dt[, stats::median(frip) , c(name_var), with = FALSE][rev(order(V1))][[name_var]])
    }else{
      stop("sort_by must be one of frip or reads_in_peak")
    }
  }
  if(!sort_by == "none"){
    stopifnot(all(frip_dt[[name_var]] %in% name_lev))
    # stopifnot(all(name_lev %in% frip_dt[[name_var]]))
    frip_dt[[name_var]] = factor(frip_dt[[name_var]], levels = name_lev)  
  }
  p_reads1 = ggplot(frip_dt,
                    aes_string(x = name_var, y = "reads_in_peak", color = color_var)) +
    geom_boxplot(outlier.shape = NA) +
    coord_cartesian(ylim = stats::quantile(frip_dt$reads_in_peak, c(0.1, 0.96))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), plot.margin = margin(.01, .01, .01, .1, unit = "npc")) +
    labs(y = "Read count per peak", x = "")

  p_frip1 = ggplot(frip_dt,
                   aes_string(x = name_var, y = "frip", color = color_var)) +
    geom_boxplot(outlier.shape = NA) +
    coord_cartesian(ylim = stats::quantile(frip_dt$frip, c(0.1, 0.96))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), plot.margin = margin(.01, .01, .01, .1, unit = "npc")) +
    labs(y = "FRIP per peak", x = "")

  tmp = frip_dt[, list(reads_in_peak = sum(reads_in_peak), mapped_reads = unique(mapped_reads)), c(union(color_var, name_var))]
  tmp[, frip := reads_in_peak / mapped_reads]
  p_fripSum1 = ggplot(tmp,
                      aes_string(x = name_var, y = "frip", fill = color_var)) +
    geom_bar(stat = "identity", color = "black") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), plot.margin = margin(.01, .01, .01, .1, unit = "npc")) +
    labs(x = "", y = "FRIP")

  if(!is.null(color_mapping)){
    p_reads1 = p_reads1 + scale_color_manual(values = color_mapping)
    p_frip1 = p_frip1 + scale_color_manual(values = color_mapping)
    p_fripSum1 = p_fripSum1 + scale_fill_manual(values = color_mapping)
  }
  if(!is.null(main_title)){
    p_reads1 = p_reads1 + labs(title = main_title)
    p_frip1 = p_frip1 + labs(title = main_title)
    p_fripSum1 = p_fripSum1 + labs(title = main_title)
  }
  
  return(list(
    reads_per_peaks = p_reads1,
    frip_per_peaks = p_frip1,
    frip_total = p_fripSum1,
    levels = name_lev
  ))

}

#' plot_scc_dt
#'
#' plot strand cross correlation (SCC) info
#'
#' @param scc_dt output from \code{\link{make_scc_dt}}
#'
#' @return list of ggplots relevant to SCC
#' @export
#'
#' @examples
#' bam_files = dir(system.file("extdata", package = "seqqc"), pattern = "^M.+bam$", full.names = TRUE)
#' query_dt.bam = make_dt(bam_files)
#'
#' peak_files = dir(system.file("extdata", package = "seqqc"), pattern = "Peak$", full.names = TRUE)
#' peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files)
#' query_gr = resize(seqsetvis::ssvOverlapIntervalSets(peak_grs), 6e2, fix = "center")
#'
#' scc_dt = make_scc_dt(query_dt.bam, query_gr)
#' plot_scc_dt(scc_dt)
plot_scc_dt = function(scc_dt, main_title = NULL, name_var = "name_split", name_lev = NULL){
  correlation = read_length = fragment_length = name = id = read_correlation = fragment_correlation = NULL #global bindings for data.table
  scc_dt_agg = scc_dt$average_correlation

  if(!is.null(name_lev)){
    stopifnot(unique(scc_dt_agg[[name_var]]) %in% name_lev)
    scc_dt_agg[[name_var]] = factor(scc_dt_agg[[name_var]], levels = name_lev)
    
    scc_dt$read_length[[name_var]] = factor(scc_dt$read_length[[name_var]], levels = name_lev)
    scc_dt$fragment_length[[name_var]] = factor(scc_dt$fragment_length[[name_var]], levels = name_lev)
  }
  
  p_scc_correlation = ggplot(scc_dt_agg, aes(x = shift, y = correlation)) +
    geom_path() +
    facet_wrap(paste0("~", name_var)) +
    geom_vline(data = scc_dt$read_length, aes(xintercept = read_length), color = "red", linetype = 2) +
    geom_vline(data = scc_dt$fragment_length, aes(xintercept = fragment_length), color = "blue", linetype = 2) +
    labs(subtitle = "Average Strand Cross Correlation (SCC)", caption = "estimated fragment size in blue, read length in red")

  tmp = scc_dt$read_correlation[, c(name_var, "id", "correlation"), with = FALSE]
  setnames(tmp, "correlation", "read_correlation")
  tmp2 = scc_dt$stable_fragment_correlation[, c(name_var, "id", "correlation"), with = FALSE]
  setnames(tmp2, "correlation", "fragment_correlation")
  scc_dt_p = merge(tmp, 
                   tmp2,
                   by = c(name_var, "id")
                   )
  if(!is.null(name_lev)){
    stopifnot(unique(scc_dt_p[[name_var]]) %in% name_lev)
    scc_dt_p[[name_var]] = factor(scc_dt_p[[name_var]], levels = name_lev)
  }
  p_scc_frag_vs_read = ggplot(scc_dt_p, aes(x = read_correlation, y = fragment_correlation)) +
    annotate("rect", xmin = .9, xmax = 1.03,
             ymin = min(scc_dt_p$fragment_correlation),
             ymax = max(scc_dt_p$fragment_correlation), fill = "#FF000015") +
    geom_point(alpha = .1) +
    expand_limits(x = c(0, 1), y = c(0, 1)) +
    facet_wrap(paste0("~", name_var)) +
    theme(panel.background = element_blank(), panel.grid = element_blank()) +
    labs(subtitle = "Strand Cross Correlation (SCC) per Feature", caption = "Features in the red zone have quite high\ncorrelation at read length and are likely artifacts",
         x = "read length SCC",
         y = "fragment length estimate SCC")
  
  if(!is.null(main_title)){
    p_scc_correlation = p_scc_correlation + labs(title = main_title)
    p_scc_frag_vs_read = p_scc_frag_vs_read + labs(title = main_title)
  }
  return(list(
    scc_curves = p_scc_correlation,
    scc_dots = p_scc_frag_vs_read
  ))
}

#'plot_signals
#'
#'Performs clustering on signal profiles in prof_dt and produces various plots.
#'Optionally outputs annotation (with anno_grs) and/or FRIP (with frip_dt) plots
#'that incorporate signal clustering information.
#'
#'@param prof_dt data.table of signal profiles from
#'  \link[seqsetvis]{ssvFetchBam} or \link[seqsetvis]{ssvFetchBigwig}.
#'@param query_gr The GRanges used to produce prof_dt.
#'@param assign_dt (optional) Precalculated clustering.  From running
#'  \link[seqsetvis]{ssvSignalClustering} followed by
#'  \code{\link{make_assign_dt}}
#'@param n_to_plot numeric. Number of items plotted in heatmap.  Default is 500.
#'@param fill_var character. Variable in prof_dt to cluster and plot on. Default
#'  of y_relative is the fraction of max y per id.
#'@param anno_grs (optional) output from \code{\link{make_anno_grs}}
#'@param frip_dt (optional) output from \code{\link{make_frip_dt}}
#'@param scc_dt (optional) output from \code{\link{make_scc_dt}}
#'@param name_lev (optional) manual specification for levels of "name" to
#'  control ordering.
#'@param nclust Number of clusters for heatmap. Defaults to 6.
#'
#'@return named list of ggplot plots and relevant data
#'@export
#'
#' @examples
#' bw_files = dir(system.file("extdata", package = "seqqc"), pattern = "^M.+bw$", full.names = TRUE)
#' query_dt = make_dt(bw_files)
#' query_dt[, sample := sub("_FE_random100.A", "", name)]
#'
#' peak_files = dir(system.file("extdata", package = "seqqc"), pattern = "Peak$", full.names = TRUE)
#' peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files)
#' query_gr = resize(seqsetvis::ssvOverlapIntervalSets(peak_grs), 6e2, fix = "center")
#'
#' prof_dt = seqsetvis::ssvFetchBigwig(query_dt, query_gr, return_data.table = TRUE)
#'
#' sig_res = plot_signals(prof_dt, query_gr)
#' sig_res$heatmap
#' sig_res$heatmap_sidebar
#' sig_res$cluster_assignment
#'
#' gtf_file = system.file(package = "seqqc", "extdata/gencode.v35.annotation.at_peaks.gtf")
#' anno_grs = make_anno_grs(gtf_file)
#'
#' sig_res.anno = plot_signals(prof_dt, query_gr, anno_grs = anno_grs)
#' sig_res.anno$heatmap
#' sig_res.anno$annotation_heatmap
#'
#' bam_files = dir(system.file("extdata", package = "seqqc"), pattern = "^M.+bam$", full.names = TRUE)
#' query_dt.bam = make_dt(bam_files)
#'
#' frip_dt = make_frip_dt(query_dt.bam, query_gr)
#'
#' sig_res.frip = plot_signals(prof_dt, query_gr, frip_dt = frip_dt)
#' sig_res.frip$heatmap
#' sig_res.frip$frip_bars_per_cluster
#' sig_res.frip$frip_text_per_cluster
#'
#' scc_dt = make_scc_dt(query_dt.bam, query_gr)
#' sig_res.scc = plot_signals(prof_dt, query_gr, scc_dt = scc_dt)
#'
#'sig_res.scc$scc_curves_per_cluster
#'sig_res.scc$scc_dots_per_cluster
#'
plot_signals = function(prof_dt, query_gr, assign_dt = NULL, n_to_plot = 500, fill_var = "y_relative",
                        anno_grs = NULL, frip_dt = NULL, scc_dt = NULL, name_lev = NULL, nclust = 6){
  y_relative = y = id = cluster_id = xmin = xmax = N = ymax = x = reads_in_peak = mapped_reads =
    treatment = name = frip = cluster_label = txt_mean = median_frip = txt = correlation =
    read_length = fragment_length = read_correlation = fragment_correlation = ymin = mean_frip = txt_median = NULL #global binding for data.table
  if(is.null(prof_dt$name)){
    prof_dt$name = prof_dt$sample
  }
  if(is.null(name_lev)){
    if(is.factor(prof_dt$name)){
      name_lev = levels(prof_dt$name)
    }else if(is.character(prof_dt$name)){
      name_lev = unique(prof_dt$name)
    }else{
      stop("prof_dt$name not character or factor")
    }
  }
  stopifnot(all(unique(prof_dt$name) %in% name_lev))
  prof_dt$name = factor(prof_dt$name, levels = name_lev)
  prof_dt$facet = prof_dt$name
  levels(prof_dt$facet) = gsub("_", "\n", levels(prof_dt$facet))

  prof_dt[, y_relative := y / max(y), list(id)]
  if(!is.null(prof_dt$cluster_id)){
    clust_dt = prof_dt
    assign_dt = unique(clust_dt[, list(id, cluster_id)])
  }else if(is.null(assign_dt)){
    set.seed(0)
    clust_dt = seqsetvis::ssvSignalClustering(prof_dt, fill_ = fill_var, max_cols = Inf, facet_ = "facet", max_rows = Inf, nclust = nclust)
    assign_dt = unique(clust_dt[, list(id, cluster_id)])
  }else{
    # assign_dt = assign_dt[order(id)]
    clust_dt = merge(prof_dt, assign_dt, by = "id")
    clust_dt$id = factor(clust_dt$id, levels = levels(assign_dt$id))
  }

  toplot_id = sampleCap(unique(clust_dt$id), n_to_plot)

  x_vals = unique(prof_dt$x)
  view_size = diff(range(x_vals)) + abs(x_vals[2] - x_vals[1])

  p_heat = seqsetvis::ssvSignalHeatmap(clust_dt[id %in% toplot_id],
                            fill_ = fill_var,
                            max_cols = Inf,
                            facet_ = "facet", show_cluster_bars = FALSE)

  p_heat_sb = seqsetvis::ssvSignalHeatmap.ClusterBars(clust_dt[id %in% toplot_id],
                                           fill_ = fill_var,
                                           max_cols = Inf,
                                           facet_ = "facet", 
                                           show_cluster_bars = FALSE, 
                                           return_unassembled_plots = TRUE)
  p_heat_sb$heatmap = p_heat_sb$heatmap +
    labs(x = paste(view_size, "bp view size"), fill = "relative pileup") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")

  p_heat_sb_assembled = seqsetvis::assemble_heatmap_cluster_bars(p_heat_sb, rel_widths = c(1, 9))


  p_heat = p_heat +
    labs(x = paste(view_size, "bp view size"), fill = "relative pileup") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")

  clust_sizes = unique(clust_dt[, list(cluster_id, id)])[, .N, list(cluster_id)][rev(order(cluster_id))]
  clust_sizes[, xmin := 0]
  clust_sizes[, xmax := 1]
  clust_sizes[, ymin := c(0, cumsum(N)[-length(N)])]
  clust_sizes[, ymax := cumsum(N)]
  clust_sizes[, col := as.character(as.numeric(cluster_id)%%2)]
  clust_sizes$facet = ""
  p_clust_boxes = ggplot(clust_sizes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
    geom_rect(aes(fill = col), color = "black") +
    geom_text(aes(label = cluster_id, x = (xmin+xmax)/2, y = (ymin+ymax)/2)) +
    scale_fill_manual(values = c("0" = "gray70", "1" = 'gray90')) +
    guides(fill = 'none') +
    coord_cartesian(expand = FALSE) +
    theme_void() +
    facet_grid(~facet)

  if(is.null(anno_grs)){
    p_heat_anno = ggplot() + theme_void() + labs(title = "anno_grs not provided")
  }else{
    anno_signal_dt = make_feature_as_signal_dt(anno_grs, query_gr)
    anno_clust_dt = anno_signal_dt[id %in% toplot_id]
    anno_clust_dt$id = factor(anno_clust_dt$id, levels = levels(clust_dt$id))
    anno_clust_dt$feature_type = factor(anno_clust_dt$feature_type, levels = rev(names(anno_grs)))
    p_heat_anno = ggplot(anno_clust_dt, aes(x = x, y = id, fill = y>0)) +
      geom_raster() +
      facet_wrap(~feature_type, nrow = 1) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      scale_fill_manual(values = c("FALSE" = "gray80", "TRUE" = "gray20")) +
      scale_x_continuous(labels = function(x)x/1e3) +
      labs(fill = "feature overlap", x= "kbp", y = "peak region") +
      theme(panel.background = element_blank(), panel.grid = element_blank())
  }

  if(is.null(frip_dt)){
    p_fripHeat = ggplot() + theme_void() + labs(title = "frip_dt no provided")
    p_clust_text = ggplot() + theme_void() + labs(title = "frip_dt no provided")
  }else{
    stopifnot(any(frip_dt$id %in% assign_dt$id))
    dt.heat = copy(frip_dt)
    dt.heat = merge(dt.heat, assign_dt, by = "id")
    tmp = dt.heat[, list(reads_in_peak = sum(reads_in_peak), mapped_reads = unique(mapped_reads)), list(treatment, name, cluster_id)]
    tmp[, frip := reads_in_peak / mapped_reads]

    .genome_fraction = function(x){
      sum(width(query_gr[x]))  /3.2e9
    }

    w_frac = lapply(split(as.character(assign_dt$id), assign_dt$cluster_id), .genome_fraction)
    w_frac = unlist(w_frac)
    tmp[, cluster_label := paste0("cluster", cluster_id, "\n", w_frac[cluster_id])]


    p_fripHeat = ggplot(tmp,
                        aes(x = name, y = frip, fill = treatment)) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), plot.margin = margin(.01, .01, .01, .1, unit = "npc")) +
      labs(subtitle = paste(sum(width((query_gr)))/3.2e9, "of genome covered by peaks"), x = "", y = "FRIP") +
      facet_grid(cluster_label~., scales = "free_y", switch = "y") +
      theme(strip.text.y = element_text(angle = 0), strip.placement = "outside")

    dt.heat_summary = dt.heat[, list(mean_frip = mean(frip), median_frip = median(frip)), list(name, cluster_id)]
    dt.heat_summary[, txt_mean := round(mean_frip*1e6, 1)]
    dt.heat_summary[, txt_median := round(median_frip*1e6, 1)]
    dt.heat_summary[, txt := paste0(round(mean_frip*1e6, 1), "\n", round(median_frip*1e6, 1))]

    p_clust_text = ggplot(dt.heat_summary, aes(x = name, y = factor(cluster_id), label = txt)) +
      geom_text(size = 5, color = NA) +
      geom_text(data = dt.heat_summary, aes(x = name, as.numeric(cluster_id)+.1, label = txt_mean), color = "red", vjust = 0) +
      geom_text(data = dt.heat_summary, aes(x = name, as.numeric(cluster_id)-.1, label = txt_median), color = "blue", vjust = 1) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), panel.background = element_blank(), panel.grid = element_blank()) +
      labs(y = "cluster", x= "", title = "mean FRIP e6", subtitle = "median FRIP e6") +
      theme(plot.title = element_text(size = 14, color = "red"), plot.subtitle = element_text(size = 14, color = "blue"))
  }

  if(is.null(scc_dt)){
    p_scc_correlation = ggplot() + theme_void() + labs(title = "scc_dt no provided")
    p_scc_frag_vs_read = ggplot() + theme_void() + labs(title = "scc_dt no provided")
  }else{

    scc_dt_agg = scc_dt$full_correlation_results
    stopifnot(any(scc_dt_agg$id %in% assign_dt$id))
    scc_dt_agg = merge(scc_dt_agg, assign_dt, by = "id")
    scc_dt_agg = scc_dt_agg[, list(correlation = mean(correlation)), list(shift, name, cluster_id)]

    p_scc_correlation = ggplot(scc_dt_agg, aes(x = shift, y = correlation)) +
      geom_path() +
      facet_grid(cluster_id~name) +
      geom_vline(data = scc_dt$read_length, aes(xintercept = read_length), color = "red", linetype = 2) +
      geom_vline(data = scc_dt$fragment_length, aes(xintercept = fragment_length), color = "blue", linetype = 2) +
      labs(title = "Strand Cross Correlation (SCC)", subtitle = "estimated fragment size in blue, read length in red")


    scc_dt_p = merge(scc_dt$read_correlation[, list(name, id, read_correlation = correlation)],
                     scc_dt$stable_fragment_correlation[, list(name, id, fragment_correlation = correlation)], by = c("name", "id"))
    scc_dt_p = merge(scc_dt_p, assign_dt, by = "id")
    p_scc_frag_vs_read = ggplot(scc_dt_p, aes(x = read_correlation, y = fragment_correlation)) +
      annotate("rect", xmin = .9, xmax = 1.03,
               ymin = min(scc_dt_p$fragment_correlation),
               ymax = max(scc_dt_p$fragment_correlation), fill = "#FF000015") +
      annotate("point", x = scc_dt_p$read_correlation, y = scc_dt_p$fragment_correlation, color = 'gray80', size = .4) +
      geom_point() +
      expand_limits(x = c(0, 1), y = c(0, 1)) +
      facet_grid(cluster_id~name) +
      theme(panel.background = element_blank(), panel.grid = element_blank()) +
      labs(title = "Peaks in the red zone have quite high\ncorrelation at read length and are likely artifacts",
           x = "read length SCC",
           y = "fragment length estimate SCC")
  }

  return(list(
    # assembled = pg_heat,
    heatmap = p_heat,
    heatmap_sidebar = p_heat_sb_assembled,
    cluster_assignment = assign_dt,
    annotation_heatmap = p_heat_anno,
    frip_bars_per_cluster = p_fripHeat,
    frip_text_per_cluster = p_clust_text,
    scc_curves_per_cluster = p_scc_correlation,
    scc_dots_per_cluster = p_scc_frag_vs_read
  ))
}

#' plot_anno_overlap
#'
#' @param anno_dt Output from \code{\link{make_anno_dt}}
#' @param name_lev Optional name levels top apply.
#'
#' @return bar plot of annotation overlap frequency
#' @export
#'
#' @examples
#' gtf_file = system.file(package = "seqqc", "extdata/gencode.v35.annotation.at_peaks.gtf")
#' anno_grs = make_anno_grs(gtf_file)
#'
#' peak_files = dir(system.file("extdata", package = "seqqc"), pattern = "Peak$", full.names = TRUE)
#' peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files)
#' names(peak_grs) = sub("_rand.+", "", names(peak_grs))
#'
#' anno_dt = make_anno_dt(peak_grs, anno_grs)
#' plot_anno_overlap(anno_dt)
plot_anno_overlap = function(anno_dt, name_lev = NULL){
  sample_cnt = fraction = feature = NULL #global binding for data.table
  if(!is.null(name_lev)){
    stopifnot(all(anno_dt$sample %in% name_lev))
    anno_dt$sample = factor(anno_dt$sample, levels = name_lev)
  }
  anno_dt = anno_dt[order(sample)]
  anno_dt$sample_cnt = factor(anno_dt$sample_cnt, levels = unique(anno_dt$sample_cnt))

  p_features = ggplot(anno_dt, aes(x = sample_cnt, y = fraction, fill = feature)) +
    geom_bar(stat = "identity") +
    labs(x = "sample\npeak count") +
    theme(panel.background = element_blank(),
          panel.grid.major.y = element_line(color = "black"),
          panel.grid.minor.y = element_line(color = "black"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    scale_fill_brewer(palette = "Set1") +
    labs(title = "Feature overlaps for peaks") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  return(list(plot = p_features, data = anno_dt))
}


#' plot_feature_overlap_signal_profiles
#'
#' @param grouped_prof_dt Output of \code{make_feature_overlap_signal_profiles}
#' @param group_var character. Must already be present in grouped_prof_dt.  Use \code{make_feature_overlap_signal_profiles} to add properly. Default is "overlap_group".
#' @param rank_var character. Must already be present in grouped_prof_dt.  Use \code{make_feature_overlap_signal_profiles} to add properly. Default is "rnk".
#' @param fill_limits Limits of heatmap fill scale.
#' @param signal_var value to plot
#'
#' @return list of ggplot plots that use feature overlaps to plot signal
#' @export
#'
#' @examples
#' bw_files = dir(system.file("extdata", package = "seqqc"), pattern = "^M.+bw$", full.names = TRUE)
#' query_dt = make_dt(bw_files)
#' query_dt[, sample := sub("_FE_random100.A", "", name)]
#' query_dt[, name_split := gsub("_", "\n", name)]
#'
#' peak_files = dir(system.file("extdata", package = "seqqc"), pattern = "Peak$", full.names = TRUE)
#' names(peak_files) = sub("_CTCF_rand.+", "", basename(peak_files))
#' peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files, )
#' overlaps_gr = seqsetvis::ssvOverlapIntervalSets(peak_grs)
#' query_gr = resize(overlaps_gr, 6e2, fix = "center")
#'
#' group_prof_dt = make_feature_overlap_signal_profiles(query_dt, overlaps_gr)
#' 
#' plots = plot_feature_overlap_signal_profiles(group_prof_dt)
#' cowplot::plot_grid(plotlist = plots)
plot_feature_overlap_signal_profiles = function(grouped_prof_dt, 
                                                group_var = "overlap_group", 
                                                rank_var = "rnk", 
                                                fill_limits = NULL, 
                                                signal_var = "y", 
                                                color_var = "name", 
                                                color_mapping = NULL,
                                                facet_var = "name_split",
                                                heatmap_free_y = TRUE){
  plot_fill_ = y = NULL #global binding for data.table
  stopifnot(group_var %in% colnames(grouped_prof_dt))
  if(!is.factor(grouped_prof_dt[[group_var]])){
    grouped_prof_dt[[group_var]] = factor(grouped_prof_dt[[group_var]])
  }
  levels(grouped_prof_dt[[group_var]]) = gsub(" ", "\n", levels(grouped_prof_dt[[group_var]]))
  stopifnot(rank_var %in% colnames(grouped_prof_dt))
  grouped_prof_dt[[rank_var]] = factor(grouped_prof_dt[[rank_var]])
  grouped_prof_dt[[rank_var]] = factor(grouped_prof_dt[[rank_var]], levels = rev(levels(grouped_prof_dt[[rank_var]])))
  if(is.null(fill_limits)){
    fill_limits = range(grouped_prof_dt[[signal_var]])
  }
  grouped_prof_dt[, plot_fill_ := get(signal_var)]

  grouped_prof_dt[plot_fill_ > max(fill_limits), plot_fill_ := max(fill_limits)]
  grouped_prof_dt[plot_fill_ < min(fill_limits), plot_fill_ := min(fill_limits)]

  p2_heat_overlaps = ggplot(grouped_prof_dt, aes_string(x = "x", y = rank_var, fill = "plot_fill_")) +
    geom_raster() +
    scale_fill_viridis_c(limits = fill_limits) +
    facet_grid(paste0(group_var, "~name_split"), scales = ifelse(heatmap_free_y, "free_y", "fixed")) +
    labs(fill = "read pileup", y = "", x = "bp", title = "Signal at peak overlap sets") +
    theme(panel.background = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          strip.text.y = element_text(angle = 0), 
          legend.position = "bottom")
  # p2_heat_overlaps

  agg_dt = grouped_prof_dt[, list(plot_fill_ = mean(plot_fill_)), c("name", "name_split", group_var, "x")]

  p2_line_facets = ggplot(agg_dt, aes_string(x = "x", y = "plot_fill_", color = color_var)) +
    geom_path() +
    facet_grid(paste0(group_var, "~name_split"), scales = "free_y") +
    labs(y = "mean read pileup", x = "bp", color = "") +
    theme(strip.text.y = element_text(angle = 0),
          legend.position = "bottom")
  if(!is.null(color_mapping)){
    p2_line_facets = p2_line_facets + 
      scale_color_manual(values = color_mapping)
  }
  list(heatmap = p2_heat_overlaps, lineplot = p2_line_facets)
}

# plot_signal_cluster_artifacts = function(artifact_profiles, signal_profiles){
#   xy_dt = prof_dt[abs(x) < 200, list(y = mean(y)), list(group, name, id)]
#   xy_dt = dcast(xy_dt, id+group~name, value.var = "y")
#
#   xy_dt[, signal := (TR_V_1 + TR_V_2 + TR_T3_1 + TR_T3_2)/4]
#
#   cutoff = 4
#
#   ggplot(xy_dt, aes(x = log10(signal+1), y = log10(control_IgG_1+1), color = group)) +
#     geom_point() +
#     annotate("line", x = range(log10(xy_dt$signal)), y = log10(cutoff+1))
#
#   # ggplot(xy_dt, aes(x = signal, y = control_IgG_1-signal, color = group)) +
#   #   geom_point()
#
#   xy_dt[, artifact := ifelse(control_IgG_1 > signal | control_IgG_1 > cutoff, "artifact", "signal")]
#   table(xy_dt$artifact)
#
#   prof_dt$artifact = NULL
#   prof_dt = merge(prof_dt, xy_dt[, list(id, artifact)], by = "id")
#
#   prof_dt[, .N, list(id, sample, x)][order(N)]
#   prof_dt[id == 1003 & sample == "TR_T3_1" & x == -2975]
#
#
#   clust_dt.artifact = ssvSignalClustering(prof_dt[artifact == "artifact"], facet_ = "sample",
#                                           max_cols = Inf, max_rows = Inf)
#   ssvSignalHeatmap(clust_dt.artifact, facet_ = "sample", fill_limits = c(0, 20)) +
#     labs(title = "artifact sites")
#
#   clust_dt.not_artifact = ssvSignalClustering(prof_dt[artifact != "artifact"], facet_ = "sample",
#                                               max_cols = Inf, max_rows = Inf)
#   ssvSignalHeatmap(clust_dt.not_artifact, facet_ = "sample", fill_limits = c(0, 20)) +
#     labs(title = "artifact sites")
# }
