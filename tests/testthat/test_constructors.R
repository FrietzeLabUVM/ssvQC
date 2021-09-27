testthat::context("run_ssvQC")
library(ssvQC)

#QcConfigFeatures
test_that("QcConfigFeatures.parse", {
  features_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
  features_config = QcConfigFeatures.parse(features_config_file)
  expect_s4_class(features_config, "QcConfigFeatures")
})

test_that("QcConfigFeatures.files", {
  np_files = dir(system.file(package = "ssvQC", "extdata"), pattern = "narrowPeak", full.names = TRUE)
  features_config = QcConfigFeatures.files(np_files)
  expect_s4_class(features_config, "QcConfigFeatures")
})

test_that("QcConfigFeatures.GRanges", {
  np_files = dir(system.file(package = "ssvQC", "extdata"), pattern = "narrowPeak", full.names = TRUE)
  np_grs = seqsetvis::easyLoad_narrowPeak(np_files)
  features_config = QcConfigFeatures.GRanges(np_grs)
  expect_s4_class(features_config, "QcConfigFeatures")
})

test_that("QcConfigFeatures.null", {
  features_config = QcConfigFeatures.null()
  expect_s4_class(features_config, "QcConfigFeatures")
})

#QcConfigSignal
test_that("QcConfigSignal.parse", {
  bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
  bam_config = QcConfigSignal.parse(bam_config_file)
  expect_s4_class(bam_config, "QcConfigSignal")
})

test_that("QcConfigSignal.files", {
  bam_files = dir(system.file(package = "ssvQC", "extdata"), pattern = "100peaks.bam$", full.names = TRUE)
  bam_config = QcConfigSignal.files(bam_files)
  expect_s4_class(bam_config, "QcConfigSignal")
})

test_that("QcConfigSignal.null", {
  bam_config = QcConfigSignal.null()
  expect_s4_class(bam_config, "QcConfigSignal")
})


#ssvQC
test_that("ssvQC with created configs", {
  features_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
  features_config = QcConfigFeatures.parse(features_config_file)
  
  bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
  bam_config = QcConfigSignal.parse(bam_config_file)
  
  sqc.complete = ssvQC(features_config, bam_config)
  expect_s4_class(sqc.complete, "ssvQC.complete")
  
  sqc.signal = ssvQC(signal_config = bam_config)
  expect_s4_class(sqc.signal, "ssvQC.signalOnly")
  
  sqc.feature = ssvQC(features_config = features_config)
  expect_is(sqc.feature, "ssvQC.featureOnly")
})