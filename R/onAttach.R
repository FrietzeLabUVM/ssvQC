.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Attaching ssvQC version ",
                        packageDescription("ssvQC")$Version, ".")
  options("SQC_COLORS" = seqsetvis::safeBrew(8, "Dark2"))
  options("SQC_CONSENSUS_N" = 1)
  options("SQC_CONSENSUS_FRACTION" = 0)
  options("SQC_VIEW_SIZE" = 3e3)
  options("SQC_PROCESS_FEATURES" = FALSE)
  options("SQC_FEATURE_FILE_SUFF" = c("narrowPeak", "broadPeak", "bed", "txt", "tab"))
  options("SQC_SIGNAL_FILE_SUFF" = c("bam", "bigwig", "bw", "bigWig", "BigWig"))
}