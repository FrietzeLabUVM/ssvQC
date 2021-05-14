.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Attaching ssvQC version ",
                        packageDescription("ssvQC")$Version, ".")
  #When adding new options here, also add them to the "names" setMethod below
  options("SQC_COLORS" = seqsetvis::safeBrew(8, "Dark2"))
  options("SQC_CONSENSUS_N" = 1)
  options("SQC_CONSENSUS_FRACTION" = 0)
  options("SQC_VIEW_SIZE" = 3e3)
  options("SQC_PROCESS_FEATURES" = FALSE)
  options("SQC_FEATURE_FILE_SUFF" = c("narrowPeak", "broadPeak", "bed", "txt", "tab"))
  options("SQC_SIGNAL_FILE_SUFF" = c("bam", "bigwig", "bw", "bigWig", "BigWig"))
  options("SQC_FORCE_CACHE_OVERWRITE" = FALSE)
  options("SQC_CACHE_VERSION" = "v1")
  options("SQC_CACHE_PATH" = "~/.cache")
  SQC_OPTIONS <<- new("SQC_OPTIONS")
}

setClass("SQC_OPTIONS", representation = list(
  is_valid = "logical"
  ))

setMethod("names", "SQC_OPTIONS",
          function(x)
          {
            c(
              "SQC_COLORS",
              "SQC_CONSENSUS_N",
              "SQC_CONSENSUS_FRACTION",
              "SQC_VIEW_SIZE",
              "SQC_PROCESS_FEATURES",
              "SQC_FEATURE_FILE_SUFF",
              "SQC_SIGNAL_FILE_SUFF",
              "SQC_FORCE_CACHE_OVERWRITE",
              "SQC_CACHE_VERSION",
              "SQC_CACHE_PATH",
              "mc.cores"
            )
          })


setMethod("$", "SQC_OPTIONS",
          function(x, name)
          {
            getOption(name)
          })

setReplaceMethod("$", "SQC_OPTIONS",
                 function(x, name, value)
                 {
                   warn_msg = "This assignment is not supported.  No effect."
                   value = list(value)
                   names(value) = name
                   do.call("options", value)
                   x
                 })



setMethod("show", "SQC_OPTIONS",
          function(object)
          {
            message("Use the $ accessor (i.e. SQC_OPTIONS$SQC_COLORS) to get/set SQC relevant options.")
            message("Use names(SQC_OPTIONS) to view all options.")
          })
