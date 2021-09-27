.remove_signal_ids = function(x, keep_ids){
  new_sig_dt = x$signal_data[id %in% keep_ids]
  
  new_nclust = length(unique(new_sig_dt$cluster_id))
  
  new_sig_dt$cluster_id = NULL
  new_sig_dt$id = as.character(new_sig_dt$id)
  
  new_sig_dt
  new_query_gr = x$query_gr[keep_ids]
  
  new_man = lapply(x@manual_assigned, function(x){
    intersect(x, keep_ids)
  })
  
  ClusteredSignal(signal_profile_dt = new_sig_dt, 
                  query_gr = new_query_gr, 
                  manual_assigned = new_man,
                  nclust = new_nclust,
                  signal_var = x@signal_var,
                  signal_var.within = x@signal_var.within,
                  facet_var = x@facet_var,
                  extra_var = x@extra_var)
}

.remove_SCC_ids = function(x, keep_ids, meta_dt){
  new_corr_res = x$full_correlation_results[id %in% keep_ids]  
  new_corr_res_l = split(new_corr_res, new_corr_res$name)
  
  rls = sapply(split(x$read_length, x$read_length$name), function(x)x$read_length)
  todo = names(new_corr_res_l)
  names(todo) = todo
  scc_res_l = lapply(todo, function(nam){
    corr_res = gather_metrics(peak_strand_corr = new_corr_res_l[[nam]], read_length = rls[nam])
    vnames = names(corr_res)
    corr_res = lapply(vnames, function(nam){
      xv = corr_res[[nam]]
      if(!is.data.table(xv)){
        xv = data.table(xv)
        setnames(xv, nam)
      }
      xv
    })
    names(corr_res) = vnames
    corr_res
  })
  
  vnames = names(scc_res_l[[1]])
  
  nam = "read_correlation"
  scc_res = lapply(vnames, function(nam){
    
    part = lapply(scc_res_l, function(x){
      xv = x[[nam]]
      # if(!is.data.table(xv)){
      #   xv = data.table(xv)
      #   setnames(xv, nam)
      # }
      xv
    })
    dt = rbindlist(part, idcol = "name")
    mdt = meta_dt
    mdt = mdt[, union("name", setdiff(colnames(mdt), colnames(dt))), with = FALSE]
    if(ncol(mdt) > 1){
      merge(dt, mdt, by = "name")  
    }else{
      dt
    }
    
  })
  names(scc_res) = vnames
  scc_res
}

.remove_FRIP_ids = function(x, keep_ids){
  x[id %in% keep_ids]
}


#' ssvQC.selectFeatures
#'
#' @param sqc 
#' @param ids 
#' @param grs 
#' @param features_name 
#' @param invert 
#'
#' @return
#' @export
#'
#' @examples
ssvQC.selectFeatures = function(sqc, ids = NULL, grs = NULL, features_name = NULL, invert = FALSE){
  ssvQC.removeFeatures(sqc, ids = ids, grs = grs, features_name = features_name, invert = !invert)
}

#' ssvQC.removeFeatures
#' 
#' @param sqc 
#' @param ids 
#' @param grs 
#' @param invert 
#' 
#' @return 
#' @export
#' 
#' @example 
#' library(ssvQC)
#' options(mc.cores = 1)
#' set.seed(0)
#' features_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
#' features_config = QcConfigFeatures.parse(features_config_file)
#'
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' bam_config = QcConfigSignal.parse(bam_config_file)
#'
#' sqc.complete.file = ssvQC(features_config_file, bam_config_file)
#'
#' sqc.complete = ssvQC(features_config, bam_config)
#' sqc = ssvQC.runAll(sqc.complete)
#' ids = c("region_1", "region_2")
#' grs = NULL
#' invert = FALSE
#' features_name = NULL
#' 
#' sqc$plots$SCC$dots
#' 
#' sqc.filtered = ssvQC.removeFeatures(sqc, ids = c("region_1", "region_2"))
#' 
#' sqc.filtered2 = ssvQC.removeFeatures(sqc.filtered, ids = c("region_1", "region_2"))
#' 
#' sqc.filtered.inv = ssvQC.removeFeatures(sqc, ids = c("region_1", "region_2", "region_3", "region_4", "region_5"), invert = TRUE)
#' 
#' sqc.filtered.gr = ssvQC.removeFeatures(sqc, grs = GRanges(paste0("chr", 1:16), IRanges(1, 1e9)))
#' sqc.filtered.gr$features_config$assessment_features
#' sqc.filtered.gr$plots$signal$heatmaps
#' sqc.filtered.gr$plots$SCC$dots
#' 
#' sqc.filtered.gr_inv = ssvQC.removeFeatures(sqc, grs = GRanges(paste0("chr", 1:16), IRanges(1, 1e9)), invert = TRUE)
#' sqc.filtered.gr_inv$features_config$assessment_features
#' sqc.filtered.gr_inv$plots$signal$heatmaps
#' sqc.filtered.gr_inv$plots$SCC$dots
ssvQC.removeFeatures = function(sqc, ids = NULL, grs = NULL, features_name = NULL, invert = FALSE){
  stopifnot(is(sqc, "ssvQC"))
  stopifnot(is.logical(invert))
  
  if(length(sqc$features_config$assessment_features) == 0) stop("Call ssvQC.prepFeatures first.")
  
  poss_feature_names = names(sqc$features_config$assessment_features)
  if(is.null(features_name)){
    if(length(poss_feature_names) > 1){
      stop("features_name is required when there is more than one feature assessment group. Should be one of: ", paste(poss_feature_names, collapse = ", ") )
    }else{
      features_name = poss_feature_names[1]
    }
  }
  stopifnot(is.character(features_name))
  stopifnot(length(features_name) == 1)
  stopifnot(features_name %in% poss_feature_names)
  
  if(is.null(ids) && is.null(grs)){
    stop("One of ids or grs is required.")
  }
  if(!is.null(ids) && !is.null(grs)){
    stop("Only one of ids or grs is allowed.")
  }
  
  assessed_grs = sqc$features_config$assessment_features[[features_name]]
  
  if(!is.null(grs)){
    stopifnot(is(grs, "GRanges"))
    ids = names(subsetByOverlaps(assessed_grs, grs))
    #select ids based on grs
  }
  if(is.factor(ids)) ids = as.character(ids)
  stopifnot(is(ids, "character"))
  
  
  if(invert){
    keep_ids = intersect(ids, names(assessed_grs))
    remv_ids = setdiff(names(assessed_grs), ids)
  }else{
    keep_ids = setdiff(names(assessed_grs), ids)
    remv_ids = intersect(ids, names(assessed_grs))
    if(!setequal(ids, remv_ids)){
      message("Not all ids present in current assessment set.")
    }
  }
  
  keep_grs = assessed_grs[keep_ids]
  remv_grs = assessed_grs[remv_ids]
  
  if(length(sqc$features_config@assessment_gr) == 0){
    stop("Call ssvQC.prepFeatures first.")
  }
  
  sqc$features_config@assessment_gr[[features_name]] = keep_grs
  sqc$features_config@overlap_gr[[features_name]] = subsetByOverlaps(sqc$features_config$overlapped_features[[features_name]], remv_grs, invert = TRUE)
  
  sqc$features_config@loaded_features[[features_name]] = lapply(sqc$features_config@loaded_features[[features_name]], function(gr){
    subsetByOverlaps(gr, remv_grs, invert = TRUE)
  })
  
  if(length(sqc$signal_data) > 0){
    sqc@signal_data[[features_name]] = lapply(sqc$signal_data[[features_name]], function(feat_sig){
      .remove_signal_ids(feat_sig, keep_ids)
    })  
  }
  
  if(!is.null(sqc@other_data$SCC)){
    sqc@other_data$SCC[[features_name]] = lapply(sqc@other_data$SCC[[features_name]], function(x){
      .remove_SCC_ids(x, keep_ids, as.data.table(sqc$signal_config$meta_data))
    })  
  }
  if(!is.null(sqc@other_data$FRIP)){
    sqc@other_data$FRIP[[features_name]] = lapply(sqc@other_data$FRIP[[features_name]], function(x){
      .remove_FRIP_ids(x, keep_ids)
    })  
  }
  
  if(!is.null(sqc@other_data$read_count_correlation) | !is.null(sqc@other_data$signal_profile_correlation)){
    sqc@other_data$read_count_correlation = NULL
    sqc@other_data$signal_profile_correlation = NULL
    sqc = ssvQC.prepCorrelation(sqc)
  }
  
  if(!is.null(sqc$plots$features)){
    sqc = ssvQC.plotFeatures(sqc)
  }
  if(!is.null(sqc$plots$signal)){
    sqc = ssvQC.plotSignal(sqc)
  }
  if(!is.null(sqc$plots$SCC)){
    sqc = ssvQC.plotSCC(sqc)
  }
  if(!is.null(sqc$plots$FRIP)){
    sqc = ssvQC.plotFRIP(sqc)
  }
  if(!is.null(sqc$plots$correlation)){
    sqc = ssvQC.plotCorrelation(sqc)
  }
  
  sqc
}

#' ssvQC.removeClusters
#'
#' @param sqc 
#' @param cluster_numbers 
#' @param features_name 
#' @param signals_name 
#' @param invert 
#'
#' @return
#' @export
#'
#' @examples
#' library(ssvQC)
#' options(mc.cores = 1)
#' set.seed(0)
#' features_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
#' features_config = QcConfigFeatures.parse(features_config_file)
#'
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' bam_config = QcConfigSignal.parse(bam_config_file)
#'
#' sqc.complete.file = ssvQC(features_config_file, bam_config_file)
#'
#' sqc.complete = ssvQC(features_config, bam_config)
#' sqc = ssvQC.runAll(sqc.complete)
#' sqc$plots$signal$heatmaps$CTCF_features$CTCF_signal
#' 
#' sqc.rmclust = ssvQC.removeClusters(sqc, cluster_numbers = 1:4)
#' sqc.rmclust$plots$signal$heatmaps$CTCF_features$CTCF_signal
#' 
#' sqc.rmclust_inv = ssvQC.removeClusters(sqc, cluster_numbers = 1:4, invert = TRUE)
#' sqc.rmclust_inv$plots$signal$heatmaps$CTCF_features$CTCF_signal
ssvQC.removeClusters = function(sqc, cluster_numbers, features_name = NULL, signals_name = NULL, invert = FALSE){
  stopifnot(is(sqc, "ssvQC"))
  stopifnot(is.logical(invert))
  
  if(length(sqc$features_config$assessment_features) == 0) stop("Call ssvQC.prepFeatures first.")
  
  poss_feature_names = names(sqc$features_config$assessment_features)
  if(is.null(features_name)){
    if(length(poss_feature_names) > 1){
      stop("features_name is required when there is more than one feature assessment group. Should be one of: ", paste(poss_feature_names, collapse = ", ") )
    }else{
      features_name = poss_feature_names[1]
    }
  }
  stopifnot(is.character(features_name))
  stopifnot(length(features_name) == 1)
  stopifnot(features_name %in% poss_feature_names)
  
  if(length(sqc$signal_data) == 0) stop("Call ssvQC.prepSignal first.")
  
  poss_signal_names = names(sqc$signal_data[[features_name]])
  if(is.null(signals_name)){
    if(length(poss_signal_names) > 1){
      stop("signals_name is required when there is more than one feature assessment group. Should be one of: ", paste(poss_signal_names, collapse = ", ") )
    }else{
      signals_name = poss_signal_names[1]
    }
  }
  stopifnot(is.character(signals_name))
  stopifnot(length(signals_name) == 1)
  stopifnot(signals_name %in% poss_signal_names)
  
  
  
  x = sqc$signal_data[[features_name]][[signals_name]]
  sel = x$signal_data[cluster_id %in% cluster_numbers]
  table(sel$cluster_id)
  ids = as.character(unique(sel$id))
  
  ssvQC.removeFeatures(sqc, ids = ids, features_name = features_name, invert = invert)
}

#' Title
#'
#' @param sqc 
#' @param cluster_numbers 
#' @param features_name 
#' @param signals_name 
#' @param invert 
#'
#' @return
#' @export
#'
#' @examples
ssvQC.selectClusters = function(sqc, cluster_numbers, features_name = NULL, signals_name = NULL, invert = FALSE){
  ssvQC.removeClusters(sqc, cluster_numbers, features_name = features_name, signals_name = signals_name, invert = !invert)
}