##CapValue used by linearQuantile

.prepCapValue = function(object, query, bfc, use_matched){
  #bam specific independent of peaks
  if(is.null(object@meta_data$fragLens)){
    sig_dt = as.data.table(object@meta_data)[, .(file, name)][order(file)]
  }else{
    sig_dt = as.data.table(object@meta_data)[, .(file, name, fragLens)][order(file)]  
  }
  if(object@read_mode != SQC_READ_MODES$bigwig){
    object = ssvQC.prepMappedReads(object)  
  }
  
  setkey(sig_dt, "name")
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
  
  cap_dt = bfcif(bfc, digest::digest(list(sig_dt, peak_dt, "ssvQC.prepCapValue")), function(){
    if(nrow(matched_dt) > 0){
      cap_dt.matched = rbindlist(lapply(seq_len(nrow(matched_dt)), function(i){
        peak_f = matched_dt[i,]$peak_file
        bam_f = matched_dt[i,]$bam_file
        # if(grepl("bam", object@read_mode)){
        #   bam_f = data.table(file = bam_f, mapped_reads = get_mapped_reads(bam_f))
        # }
        name_i = matched_dt[i,]$name
        peak_gr = query@feature_load_FUN(peak_f)[[1]]
        peak_gr = peak_gr[sample(min(5e3, length(peak_gr)))]
        
        if(object@read_mode == "bam_SE"){
          fragLens = sig_dt[.(name_i)]$fragLens  
          fragLens = ifelse(is.null(fragLens), "auto", fragLens)
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
        if(object@read_mode == "bam_SE"){
          fragLens = sig_dt[.(name_i)]$fragLens  
          fragLens = ifelse(is.null(fragLens), "auto", fragLens)
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
    cap_dt
  })
  
  cap_dt[y_cap_value <= 1, y_cap_value := 1]
  
  if(grepl("bam", object@read_mode)){
    cap_dt[, mapped_reads := get_mapped_reads(as.character(sample)), .(sample) ]
    cap_dt[, y_RPM_cap_value := y_cap_value / max(mapped_reads, 1) * 1e6]
  }
  
  object@meta_data$cap_value = cap_dt[.(object@meta_data$name)]$y_cap_value
  if(!is.null(cap_dt$y_RPM_cap_value)){
    object@meta_data$RPM_cap_value = cap_dt[.(object@meta_data$name)]$y_RPM_cap_value  
  }
  object
}

#' @export
#' @rdname ssvQC
setGeneric("ssvQC.prepCapValue", function(object, query, bfc, use_matched){standardGeneric("ssvQC.prepCapValue")})
setMethod("ssvQC.prepCapValue", "ssvQC.complete", function(object){
  object@signal_config = ssvQC.prepCapValue(object@signal_config, object@features_config, object@bfc, object@matched_only)
  object
})
setMethod("ssvQC.prepCapValue", "ssvQC.featureOnly", function(object){
  stop("Cannot run prepCapValue on ssvQC with no QcConfigSignal component")
})
setMethod("ssvQC.prepCapValue", "ssvQC.signalOnly", function(object){
  stop("Cannot run prepCapValue on ssvQC with no QcConfigFeatures component")
})
setMethod("ssvQC.prepCapValue", c("QcConfigSignal", "QcConfigFeatures", "BiocFileCache"), function(object, query, bfc, use_matched){
  .prepCapValue(object, query, bfc, use_matched)
})
