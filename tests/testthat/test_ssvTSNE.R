library(ssvQC)
library(testthat)
options(mc.cores = 1)
SQC_OPTIONS$SQC_FORCE_CACHE_OVERWRITE = TRUE
set.seed(0)
features_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
features_config = QcConfigFeatures.parse(features_config_file)

bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
bam_config = QcConfigSignal.parse(bam_config_file)
bam_config$view_size = 600

# TSNE plot paramters can be controlled at ssvTSNE object creation or set later
# if the y settings aren't specified they are the same as x
sts = ssvTSNE(
  features_config,
  bam_config,
  n_glyphs_x = 3,
  n_heatmap_pixels_x = 5)
# 
# 
# sts$perplexity
# 
# # sts = ssvQC.plotFeatures(sts)
# # sts$perplexity = 10
# sts = ssvQC.prepFetch(sts)
# # debug(ssvQC:::ClusteredSignal_TSNE.from_ClusteredSignal)
# sts = ssvQC.prepSignal(sts)
# sts$signal_data
# sts = ssvQC.referenceUsesSameScale(sts)
# sts = ssvQC.prepSignal(sts)
# sts = ssvQC.plotSignal(sts)

suppressWarnings({
  sts = ssvQC.plotSignal(sts)  
})


test_that("QcConfigFeatures.parse", {
  expect_is(sts$plots$TSNE, "list")
  expect_is(sts$plots$TSNE$regional_glyphs$CTCF_features$CTCF_signal, "ggplot")
  expect_is(sts$plots$TSNE$cluster_glyphs$CTCF_features$CTCF_signal, "ggplot")
  expect_is(sts$plots$TSNE$regional_heatmap$CTCF_features$CTCF_signal, "ggplot")
})


