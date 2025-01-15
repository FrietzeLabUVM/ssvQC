##MappedReads

#' @export
#' @rdname ssvQC
setGeneric("ssvQC.prepMappedReads", function(object){standardGeneric("ssvQC.prepMappedReads")})
setMethod("ssvQC.prepMappedReads", "ssvQC.complete", function(object){
  object@signal_config = ssvQC.prepMappedReads(object@signal_config)
  object
})
setMethod("ssvQC.prepMappedReads", "ssvQC.featureOnly", function(object){
  stop("Cannot run mapped reads on ssvQC with no QcConfigSignal component")
  message("featureOnly")
})
setMethod("ssvQC.prepMappedReads", "ssvQC.signalOnly", function(object){
  object@signal_config = ssvQC.prepMappedReads(object@signal_config)
  object
})
setMethod("ssvQC.prepMappedReads", c("QcConfigSignal"), function(object){
  #bam specific independent of peaks
  if(grepl("bam", object@read_mode)){
    object@meta_data$mapped_reads = sapply(object@meta_data$file, get_mapped_reads)
  }else{
    # warning("ssvQC.prepMappedReads called on bigwig read_mode QcConfigSignal.")
    # object@meta_data$mapped_reads = NA
  }
  object
})

.ssvQC.plotMappedReads = function(object){
  full_bam_config = object@signal_config
  full_meta = full_bam_config@meta_data
  if(is.null(full_meta$mapped_reads)){
    full_bam_config = ssvQC.prepMappedReads(full_bam_config)
    object@signal_config = full_bam_config
  }
  
  todo = .make_query_signal_config(full_bam_config)
  plots = lapply(todo, function(sig_config){
    if(grepl("bam", sig_config@read_mode)){
      bam_config_dt = sig_config@meta_data
      
      color_var = sig_config@color_by
      group_var = sig_config@run_by
      color_mapping = sig_config@color_mapping
      
      p_mapped_reads = ggplot(bam_config_dt, aes_string(x = "name", y = "mapped_reads", fill = color_var)) +
        geom_bar(stat = "identity", position = "dodge", color = "black") +
        scale_fill_manual(values = color_mapping) +
        scale_y_continuous(labels = function(x)paste(x/1e6, "M")) +
        labs(y = "mapped reads", fill = color_var, x= "") +
        labs(title = "Mapped reads")
      
      if(is.null(object@plots$reads)) object@plots$reads = list()
      
    }else{
      p_mapped_reads = ggplot() + theme_void() + labs(title = "Could not run mapped reads on non-bam")
    }
    p_mapped_reads
  })
  object@plots$reads = plots
  
  object
}

#' @export
#' @rdname ssvQC
setGeneric("ssvQC.plotMappedReads", function(object){standardGeneric("ssvQC.plotMappedReads")})
setMethod("ssvQC.plotMappedReads", "ssvQC.complete", .ssvQC.plotMappedReads)
setMethod("ssvQC.plotMappedReads", "ssvQC.featureOnly", function(object){
  stop("Cannot run mapped reads on ssvQC with no QcConfigSignal component")
  message("featureOnly")
})
setMethod("ssvQC.plotMappedReads", "ssvQC.signalOnly", .ssvQC.plotMappedReads)
