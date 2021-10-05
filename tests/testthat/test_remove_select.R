testthat::context("remove_select")
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

sqc.complete.prepFeatures = suppressWarnings({ssvQC.prepFeatures(sqc.complete)})
sqc.complete.prepSignal = suppressWarnings({ssvQC.prepSignal(sqc.complete)})

sqc.complete.ran = ssvQC.runAll(sqc.complete)
sqc.signal.ran = ssvQC.runAll(sqc.signal)
sqc.feature.ran = ssvQC.runAll(sqc.feature)

plot_names = c("features",    "reads",       "signal",      "SCC",         "FRIP",        "correlation")
SCC_names = c("read_length", "fragment_length", "read_correlation", "flex_fragment_correlation", "stable_fragment_correlation", "full_correlation_results", "average_correlation")
FRIP_names = c("name_split", "id", "reads_in_peak", "cell", "sample", "mapped_reads", "frip")
SCC_partnames = list(
    read_length = c("name_split", "read_length", "file", "cell", "mark", "rep", "name", "mapped_reads", "fragLens", "cap_value", "RPM_cap_value"),
    fragment_length = c("name_split", "fragment_length", "file", "cell", "mark", "rep", "name", "mapped_reads", "fragLens", "cap_value", "RPM_cap_value"),
    read_correlation = c("name_split", "shift", "id", "correlation", "file", "cell", "mark", "rep", "name", "mapped_reads", "fragLens", "cap_value", "RPM_cap_value"),
    flex_fragment_correlation = c("name_split", "shift", "id", "correlation", "file", "cell", "mark", "rep", "name", "mapped_reads", "fragLens", "cap_value", "RPM_cap_value"),
    stable_fragment_correlation = c("name_split", "shift", "id", "correlation", "file", "cell", "mark", "rep", "name", "mapped_reads", "fragLens", "cap_value", "RPM_cap_value"),
    full_correlation_results = c("name_split", "shift", "id", "correlation", "file", "cell", "mark", "rep", "name", "mapped_reads", "fragLens", "cap_value", "RPM_cap_value"),
    average_correlation = c("name_split", "shift", "correlation", "file", "cell", "mark", "rep", "name", "mapped_reads", "fragLens", "cap_value", "RPM_cap_value")
)


print_names = function(x){
    (paste0('c("', paste(x, collapse = '", "'), '")'))
}
names(sqc.complete.ran$SCC[[1]][[1]])
hidden = lapply(sqc.complete.ran$SCC[[1]][[1]], function(x){
    message(print_names(colnames(x)))
})

test_that("error if not yet ssvQC.prepFeatures", {
  expect_error(ssvQC.removeFeatures(sqc.complete, "region_1", features_name = "CTCF_features"), regexp = "Call ssvQC.prepFeatures first.")
  expect_error(ssvQC.removeClusters(sqc.complete, "1", features_name = "CTCF_features"), regexp = "Call ssvQC.prepFeatures first.")
})

test_that("can remove by id if not yet ssvQC.prepSignal, but not cluster", {
  sqc.rm1 = ssvQC.removeFeatures(sqc.complete.prepFeatures, "region_45", features_name = "CTCF_features")
  expect_equal(length(sqc.rm1$features_config$assessment_features[[1]]), 49)
  expect_true("region_45" %in% names(sqc.complete.prepFeatures$features_config$assessment_features[[1]]))
  expect_false("region_45" %in% names(sqc.rm1$features_config$assessment_features[[1]]))
  
  sqc.complete.ran.rm1 = ssvQC.removeFeatures(sqc.complete.ran, "region_45", features_name = "CTCF_features")
  expect_equal(length(sqc.complete.ran.rm1$features_config$assessment_features[[1]]), 49)
  expect_false("region_45" %in% names(sqc.complete.ran.rm1$features_config$assessment_features[[1]]))
  
  expect_failure(expect_equivalent(sqc.complete.ran$plots$signal$heatmaps, sqc.complete.ran.rm1$plots$signal$heatmaps))
  expect_failure(expect_equivalent(sqc.complete.ran$plots$SCC$dots, sqc.complete.ran.rm1$plots$SCC$dots))
  expect_failure(expect_equivalent(sqc.complete.ran$plots$FRIP$total, sqc.complete.ran.rm1$plots$FRIP$total))
  
  expect_error(ssvQC.removeClusters(sqc.complete.prepFeatures, "1", features_name = "CTCF_features"), regexp = "Call ssvQC.prepSignal first.")
})

test_that("can remove by cluster after running ssvQC.prepSignal", {
  id_removed.ran = names(sqc.complete.ran$signal_data$CTCF_features$CTCF_signal$query_gr.cluster_list$`1`)
  
  sqc.complete.prepSignal.rm1 = ssvQC.removeClusters(sqc.complete.prepSignal, "1", features_name = "CTCF_features")
  sqc.complete.ran.rm1 = ssvQC.removeClusters(sqc.complete.ran, "1", features_name = "CTCF_features")
  
  expect_false(any(id_removed.ran %in% names(sqc.complete.ran.rm1$features_config$assessment_features$CTCF_features)))
  expect_false(any(id_removed.ran %in% sqc.complete.ran.rm1$signal_data$CTCF_features$CTCF_signal$signal_data$id))
  
  expect_true(all(id_removed.ran %in% names(sqc.complete.ran$features_config$assessment_features$CTCF_features)))
  expect_true(all(id_removed.ran %in% sqc.complete.ran$signal_data$CTCF_features$CTCF_signal$signal_data$id))
})

###selection
test_that("error if not yet ssvQC.prepFeatures", {
  expect_error(ssvQC.selectClusters(sqc.complete, "region_1", features_name = "CTCF_features"), regexp = "Call ssvQC.prepFeatures first.")
  expect_error(ssvQC.selectClusters(sqc.complete, "1", features_name = "CTCF_features"), regexp = "Call ssvQC.prepFeatures first.")
})

test_that("can select by id if not yet ssvQC.prepSignal, but not cluster", {
  expect_error(ssvQC.selectFeatures(sqc.complete.prepFeatures, "region_45", features_name = "CTCF_features"), "Currently, fewer than 5 regions are not allowed. This will be addressed in future versions.")
  sqc.rm1 = ssvQC.selectFeatures(sqc.complete.prepFeatures, paste0("region_", 45:50), features_name = "CTCF_features")
  expect_equal(length(sqc.rm1$features_config$assessment_features[[1]]), 6)
  expect_true(all(paste0("region_", 45:50) %in% names(sqc.complete.prepFeatures$features_config$assessment_features[[1]])))
  #can't select fewer than 5
  sqc.complete.ran.rm1 = ssvQC.selectFeatures(sqc.complete.ran, paste0("region_", 45:50), features_name = "CTCF_features")
  expect_equal(length(sqc.complete.ran.rm1$features_config$assessment_features[[1]]), 6)
  expect_true(all(paste0("region_", 45:50) %in% names(sqc.complete.ran.rm1$features_config$assessment_features[[1]])))
  
  expect_failure(expect_equivalent(sqc.complete.ran$plots$signal$heatmaps, sqc.complete.ran.rm1$plots$signal$heatmaps))
  expect_failure(expect_equivalent(sqc.complete.ran$plots$SCC$dots, sqc.complete.ran.rm1$plots$SCC$dots))
  expect_failure(expect_equivalent(sqc.complete.ran$plots$FRIP$total, sqc.complete.ran.rm1$plots$FRIP$total))
  
  expect_error(ssvQC.selectClusters(sqc.complete.prepFeatures, "1", features_name = "CTCF_features"), regexp = "Call ssvQC.prepSignal first.")
})

test_that("can select by cluster after running ssvQC.prepSignal", {
  id_selectd.ran = names(sqc.complete.ran$signal_data$CTCF_features$CTCF_signal$query_gr.cluster_list$`1`)
  
  sqc.complete.prepSignal.rm1 = ssvQC.selectClusters(sqc.complete.prepSignal, "1", features_name = "CTCF_features")
  sqc.complete.ran.rm1 = ssvQC.selectClusters(sqc.complete.ran, "1", features_name = "CTCF_features")
  
  expect_true(all(id_selectd.ran %in% names(sqc.complete.ran.rm1$features_config$assessment_features$CTCF_features)))
  expect_true(all(id_selectd.ran %in% sqc.complete.ran.rm1$signal_data$CTCF_features$CTCF_signal$signal_data$id))
})