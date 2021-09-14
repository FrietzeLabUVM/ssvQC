testthat::context("unique_names")
#all names and names_split must be unique
library(ssvQC)
library(testthat)

options(mc.cores = 1)
set.seed(0)


features_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
features_config = QcConfigFeatures.parse(features_config_file)

bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
bam_config = QcConfigSignal.parse(bam_config_file)

sqc.complete = ssvQC(features_config, bam_config)
sqc.signal = ssvQC(signal_config = bam_config)
sqc.feature = ssvQC(features_config = features_config)

#QcConfigFeatures
#required attributes not NULL
test_that("get error when setting all name to NULL", {
  to_test = features_config
  expect_error({to_test$meta_data$name = NULL}, regexp = 'invalid class') 
})
test_that("get error when setting all name_split to NULL", {
  to_test = features_config
  expect_error({to_test$meta_data$name_split = NULL}, regexp = 'invalid class') 
})
#no duplicates
test_that("get error when setting all name to A", {
  to_test = features_config
  expect_error({to_test$meta_data$name = "A"}, regexp = 'invalid class') 
})

test_that("get error when setting all name to A", {
  to_test = features_config
  expect_error({to_test$meta_data$name_split = "A"}, regexp = 'invalid class') 
})
