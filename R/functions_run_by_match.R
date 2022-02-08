run_by_match = function(name, object = NULL, callFUN = NULL){
  must_match = object@matched_only
  query_gr = object@features_config$assessment_features[[name]]
  sig_configs = .make_query_signal_config(object@signal_config)
  if(must_match){
    sig_name = feature_name2signal_name(name)
    if(!sig_name %in% names(sig_configs)){
      must_match = FALSE
    }else{
      out = lapply(sig_configs[feature_name2signal_name(name)], 
                   callFUN, 
                   query_gr = query_gr, 
                   bfc = object@bfc)    
    }
  }
  if(!must_match){
    out = lapply(sig_configs, 
                 callFUN, 
                 query_gr = query_gr, 
                 bfc = object@bfc)      
  }
  out
}

make_frip_dt.run_by_match = function(sel_sig_config, query_gr, bfc){
  make_frip_dt(as.data.table(sel_sig_config@meta_data), 
               query_gr = query_gr, 
               color_var = sel_sig_config@color_by, 
               bfc = bfc, 
               force_overwrite = getOption("SQC_FORCE_CACHE_OVERWRITE", FALSE)
  )
}

make_scc_dt.run_by_match = function(query_gr, sel_sig_config, bfc){
  
}