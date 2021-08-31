signal_vars= c(raw = "raw", RPM = "RPM", linearQuantile = "linearQuantile", RPM_linearQuantile = "RPM_linearQuantile")
sqc_read_modes = list(bam_SE = "bam_SE", bam_PE = "bam_PE", bigwig = "bigwig", null = "null")

val2var = c(
  raw = "y",
  RPM = "y_RPM",
  linearQuantile = "y_linQ",
  RPM_linearQuantile = "y_RPM_linQ"  
)

val2lab = c(
  raw = "read\npileup",
  RPM = "RPM\npileup",
  linearQuantile = "normalized\npileup",
  RPM_linearQuantile = "RPM normalized\npileup"  
)

val2bwlab = c(
  raw = "bigWig\nsignal",
  RPM = "***ERROR***",
  linearQuantile = "normalized\nbigWig\nsignal",
  RPM_linearQuantile = "***ERROR***"    
)

stopifnot(setequal(signal_vars, names(val2var)))
stopifnot(setequal(signal_vars, names(val2lab)))
stopifnot(setequal(signal_vars, names(val2bwlab)))

sqc_signal_values = as.list(signal_vars)
names(sqc_signal_values) = signal_vars

get_default_signal_var = function(read_mode){
  if(read_mode %in% c(sqc_read_modes$bigwig, sqc_read_modes$null)){
    sqc_signal_values$raw
  }else if (read_mode %in% c(sqc_read_modes$bam_SE, sqc_read_modes$bam_PE)){
    sqc_signal_values$RPM
  }else{
    stop("Unrecognized read_mode: ", read_mode)
  }
}

