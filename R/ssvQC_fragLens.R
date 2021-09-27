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