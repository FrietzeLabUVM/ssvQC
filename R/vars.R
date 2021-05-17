signal_vars= c("raw", "RPM", "linearQuantile")

val2var = c(
  raw = "y",
  RPM = "y_RPM",
  linearQuantile = "y_linQ"  
)

val2lab = c(
  raw = "read\npileup",
  RPM = "RPM\npileup",
  linearQuantile = "normalized\npileup"  
)

val2bwlab = c(
  raw = "bigWig\nsignal",
  RPM = "***ERROR***",
  linearQuantile = "normalized\nbigWig\nsignal"  
)

stopifnot(setequal(signal_vars, names(val2var)))
stopifnot(setequal(signal_vars, names(val2lab)))
stopifnot(setequal(signal_vars, names(val2bwlab)))

sqc_signal_values = as.list(signal_vars)
names(sqc_signal_values) = signal_vars
