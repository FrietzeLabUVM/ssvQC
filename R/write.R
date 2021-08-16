#' write_ssvQC.per_peak
#'
#' @param sqc 
#' @param f 
#'
#' @return
#' @export
#'
#' @examples
write_ssvQC.per_peak = function(sqc, out_dir = getwd()){
  make_qc_table = function(x, value.var, value.var.final = value.var){
    x[, name := gsub("\n", "_", name_split)]
    xx = x[, c("id", "name", value.var), with = FALSE]
    setnames(xx, value.var, value.var.final)
    xx = dcast(xx, id~name, value.var = value.var.final)
    xx = as.data.frame(xx)
    rownames(xx) = xx$id
    xx$id = NULL
    xx
  }
  
  
  lapply(names(sqc$FRIP), function(f_name){
    lapply(names(sqc$FRIP[[f_name]]), function(s_name){
      qc_gr = sqc$features_config$assessment_features[[f_name]]
      qc_gr.memb = as.data.frame(mcols(qc_gr))
      mcols(qc_gr) = NULL
      qc_gr.df = as.data.frame(qc_gr)
      qc_gr.df$width = NULL
      qc_gr.df$id = rownames(qc_gr.df)
      
      colnames(qc_gr.memb) = gsub("\\.", "_", colnames(qc_gr.memb))
      colnames(qc_gr.memb) = paste0("peakcall_", colnames(qc_gr.memb))
      
      qc_gr.df = cbind(qc_gr.df, qc_gr.memb[rownames(qc_gr.df),])
      
      qc_gr = as.data.frame(sqc$features_config$assessment_features[[f_name]])
      qc_gr$width = NULL
      
      qc_frip = make_qc_table(sqc$FRIP[[f_name]][[s_name]], c("reads_in_peak", "frip"), c("reads", "FRIP"))
      head(qc_frip)
      
      qc_gr.df = cbind(qc_gr.df, qc_frip[rownames(qc_gr.df),])
      head(qc_gr.df)
      
      qcc_scc.stable = make_qc_table(sqc$SCC[[f_name]][[s_name]]$stable_fragment_correlation, 
                                     c("shift", "correlation"),
                                     c("stableShift", "stableCorrelation"))
      
      qcc_scc.flex = make_qc_table(sqc$SCC[[f_name]][[s_name]]$flex_fragment_correlation, 
                                   c("shift", "correlation"),
                                   c("flexShift", "flexCorrelation"))
      
      qcc_scc.read = make_qc_table(sqc$SCC[[f_name]][[s_name]]$read_correlation, 
                                   c("shift", "correlation"),
                                   c("readShift", "readCorrelation"))
      
      qc_gr.df = cbind(qc_gr.df, 
                       qcc_scc.stable[rownames(qc_gr.df),], 
                       qcc_scc.flex[rownames(qc_gr.df),], 
                       qcc_scc.read[rownames(qc_gr.df),])
      out_f = file.path(out_dir, paste0("qc_perpeak.", f_name, ".", s_name, ".bedlike.txt"))
      fwrite(qc_gr.df, out_f, sep = "\t")
      message("wrote ", out_f)
    })
  })
  invisible()
}

#' write_ssvQC.summary
#'
#' @param sqc 
#' @param out_dir 
#'
#' @return
#' @export
#'
#' @examples
write_ssvQC.summary = function(sqc, out_dir = getwd()){
  if(is.null(sqc@other_data$FRIP)){
    stop("missing FRIP not yet supported")
  }
  if(is.null(sqc$SCC)){
    stop("missing SCC not yet supported") 
  }
  lapply(names(sqc$FRIP), function(f_name){
    lapply(names(sqc$FRIP[[f_name]]), function(s_name){
      qc_dt = sqc$FRIP[[f_name]][[s_name]][, list(FRIP_assessed = sum(reads_in_peak) / mapped_reads), list(name_split, mapped_reads)]
      qc_dt$mapped_reads = NULL
      qc_dt = merge(qc_dt, sqc$SCC[[f_name]][[s_name]]$fragment_length, by = "name_split")
      qc_dt$name_split = NULL
      sel_cn = colnames(sqc$signal_config$meta_data)
      setcolorder(qc_dt, intersect(sel_cn, colnames(qc_dt)))
      out_f = file.path(out_dir, paste0("qc_summary.", f_name, ".", s_name, ".csv"))
      fwrite(qc_dt, out_f, sep = ",")
      message("wrote ", out_f)
    })
  })
  invisible()
}

#' write_ssvQC.correlation
#'
#' @param sqc 
#' @param out_dir 
#'
#' @return
#' @export
#'
#' @examples
write_ssvQC.correlation = function(sqc, out_dir = getwd()){
  if(is.null(sqc@other_data$read_count_correlation)){
    stop("missing read_count_correlation not yet supported")
  }
  if(is.null(sqc@other_data$signal_profile_correlation)){
    stop("missing signal_profile_correlation not yet supported") 
  }
  todo = sqc@other_data$signal_profile_correlation
  hidden = lapply(names(todo), function(f_name){
    lapply(names(todo[[f_name]]), function(s_name){
      cor_mat = todo[[f_name]][[s_name]]$mat
      rownames(cor_mat) = gsub("\n", "_", rownames(cor_mat))
      colnames(cor_mat) = gsub("\n", "_", colnames(cor_mat))
      out_f = file.path(out_dir, paste0("qc_signal_profile_correlation.", f_name, ".", s_name, ".csv"))
      write.csv(cor_mat, out_f)
      message("wrote ", out_f)
    })
  })
  
  todo = sqc@other_data$read_count_correlation
  hidden = lapply(names(todo), function(f_name){
    lapply(names(todo[[f_name]]), function(s_name){
      cor_mat = todo[[f_name]][[s_name]]$mat
      rownames(cor_mat) = gsub("\n", "_", rownames(cor_mat))
      colnames(cor_mat) = gsub("\n", "_", colnames(cor_mat))
      out_f = file.path(out_dir, paste0("qc_read_count_correlation.", f_name, ".", s_name, ".csv"))
      write.csv(cor_mat, out_f)
      message("wrote ", out_f)
    })
  })
  invisible()
}