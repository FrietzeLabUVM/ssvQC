

#' ClusteredSignal
#'
#' @slot signal_data 
#' @slot query_gr 
#' @slot signal_var 
#' @slot facet_var 
#' @slot extra_var 
#' @export
setClass("ClusteredSignal",
         representation = list(
           signal_data = "data.table",
           query_gr = "GRanges",
           signal_var = "character",
           facet_var = "character",
           extra_var = "character"
         ))

#' Title
#'
#' @param signal_profile_dt 
#' @param query_gr 
#' @param manual_assigned 
#' @param nclust 
#' @param signal_var 
#' @param facet_var 
#' @param extra_var 
#'
#' @return
#' @export
#'
#' @examples
ClusteredSignal = function(signal_profile_dt,
                           query_gr,
                           manual_assigned = list(),
                           nclust = 6,
                           signal_var = "y",
                           signal_var.within = "y",
                           facet_var = "name_split",
                           extra_var = character()){
  query_gr = seqsetvis::prepare_fetch_GRanges_names(query_gr)
  
  if(is.null(signal_profile_dt[[signal_var]])) stop("signal_var \"", signal_var, "\" not found in signal_data." )
  clust_dt = seqsetvis::ssvSignalClustering(signal_profile_dt, nclust = nclust, facet_ = "name_split", max_cols = Inf, max_rows = Inf, fill_ = signal_var)
  if(signal_var != signal_var.within){
    clust_dt = within_clust_sort(clust_dt = clust_dt, facet_ = "name_split", fill_ = signal_var.within)
  }
  
  if(length(manual_assigned) > 0){
    stop("manual_assigned NYI")
  }
  
  new("ClusteredSignal",
      signal_data =  clust_dt,
      query_gr = query_gr,
      signal_var = signal_var,
      facet_var = facet_var,
      extra_var = extra_var)
}

#' ClusteredSignal.fromConfig
#'
#' @return
#' @export
#'
#' @examples
#' options(mc.cores = 10)
#' set.seed(0)
#' feature_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
#' feature_config = QcConfigFeatures.parse(feature_config_file, process_features= TRUE)
#'
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' bam_config = QcConfigSignal.parse(bam_config_file)
#' 
#' sqc = ssvQC(feature_config, bam_config)
#' 
#' prepS
#' query_gr = feature_config$assessment_features$CTCF_features
#' sclust = ClusteredSignal(sqc@signal_profile, sqc@feature_config$assessment_features)
#' sclust$
ClusteredSignal.fromConfig = function(signal_config,
                                      query_gr,
                                      manual_assigned = list(),
                                      nclust = 6,
                                      facet_var = "name_split",
                                      extra_var = character(),
                                      bfc = new_cache()){
  if(signal_config@cluster_value == "RPM" | signal_config@sort_value == "RPM"){
    if(is.null(signal_config@meta_data$mapped_reads)){
      stop("Call ssvQC.prepMappedReads() on signal_config first.")
    }
  }
  
  if(signal_config@cluster_value == "linearQuantile" | signal_config@sort_value == "linearQuantile"){
    if(is.null(signal_config@meta_data$cap_value)){
      stop("Call ssvQC.prepCapValue() on signal_config first.")
    }
  }
  
  query_gr = seqsetvis::prepare_fetch_GRanges_names(query_gr)
  
  prof_dt = fetch_signal_at_features(signal_config, query_gr, bfc)
  # if(signal_config@cluster_value == "RPM" | signal_config@sort_value == "RPM"){
  prof_dt[, y_RPM := y / mapped_reads * 1e6]
  # }
  # if(signal_config@cluster_value == "linearQuantile" | signal_config@sort_value == "linearQuantile"){
  prof_dt[, y_linQ := y / cap_value]
  prof_dt[y_linQ > 1, y_linQ := 1]
  # }
  clust_dt = ClusteredSignal(prof_dt, query_gr, 
                             manual_assigned = manual_assigned,
                             nclust = nclust,
                             signal_var = val2var[signal_config@cluster_value],
                             signal_var.within = val2var[signal_config@sort_value],
                             facet_var = facet_var,
                             extra_var = extra_var)
  clust_dt
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
ClusteredSignal.null = function(){
  new("ClusteredSignal",
      signal_data =  data.table(),
      query_gr = GRanges(),
      signal_var = character(),
      facet_var = character(),
      extra_var = character())
}

setMethod("names", "ClusteredSignal",
          function(x)
          {
            c("signal_data", "query_gr", "assignment_data", "query_gr.cluster_list", "signal_data.mean_per_cluster")
            
          })


setMethod("$", "ClusteredSignal",
          function(x, name)
          {
            switch (name,
                    signal_data = x@signal_data,
                    query_gr = x@query_gr,
                    assignment_data = unique(x@signal_data[, .(id, cluster_id)]),
                    query_gr.cluster_list = {
                      assign_dt = x$assignment_data
                      qgr = x@query_gr[assign_dt$id]
                      split(qgr, assign_dt$cluster_id)
                    },
                    signal_data.mean_per_cluster = {
                      sig_var = x@signal_var
                      agg_dt = x@signal_data[, .(y = mean(get(sig_var))), c("x", "cluster_id", union(x@facet_var, x@extra_var))]
                      setnames(agg_dt, "y", sig_var)
                      agg_dt
                    }
            )
          })

setReplaceMethod("$", "ClusteredSignal",
                 function(x, name, value)
                 {
                   warn_msg = "This assignment is not supported.  No effect."
                   warning(warn_msg)
                   NULL
                 })
