
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
                           facet_var = "name_split",
                           extra_var = character()){
  query_gr = seqsetvis::prepare_fetch_GRanges_names(query_gr)
  
  if(is.null(signal_profile_dt[[signal_var]])) stop("signal_var \"", signal_var, "\" not found in signal_data." )
  clust_dt = seqsetvis::ssvSignalClustering(signal_profile_dt, nclust = nclust, facet_ = "name_split", max_cols = Inf, max_rows = Inf, fill_ = signal_var)
  
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
#' signal_data = sqc@signal_profile
#' query_gr = sqc@feature_config$assessment_features
#' sclust = ClusteredSignal(sqc@signal_profile, sqc@feature_config$assessment_features)
#' sclust$
ClusteredSignal.fromConfig = function(signal_config,
                                      query_gr,
                                      manual_assigned = list(),
                                      nclust = 6,
                                      signal_var = "y",
                                      facet_var = "name_split",
                                      extra_var = character()){
  query_gr = seqsetvis::prepare_fetch_GRanges_names(query_gr)
  
  prof_dt = fetch_signal_at_features(signal_config, query_gr)
  
  ClusteredSignal(prof_dt, query_gr, 
                  manual_assigned = manual_assigned,
                  nclust = nclust,
                  signal_var = signal_var,
                  facet_var = facet_var,
                  extra_var = extra_var)
  
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
