testthat::context("run_ssvQC")
library(ssvQC)

options(mc.cores = 1)
set.seed(0)


features_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
features_config = QcConfigFeatures.parse(features_config_file)

bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
bam_config = QcConfigSignal.parse(bam_config_file)

sqc.complete = ssvQC(features_config, bam_config)
sqc.signal = ssvQC(signal_config = bam_config)
sqc.feature = ssvQC(features_config = features_config)

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

test_that("complete plots", {
    expect_equal(length(sqc.complete$plots), 0)
    expect_setequal(names(sqc.complete.ran$plots), plot_names)
    
    
    expect_null(sqc.complete$FRIP)
    expect_null(sqc.complete$correlation[[1]])
    
    expect_true(all(c("SCC", "FRIP", "correlation") %in% names(sqc.complete.ran)))
})

test_that("complete SCC", {
    expect_null(sqc.complete$SCC)
    
    expect_equal(names(sqc.complete.ran$SCC), "CTCF_features")
    expect_equal(names(sqc.complete.ran$SCC[[1]]), "CTCF_signal")
    expect_is(sqc.complete.ran$SCC[[1]][[1]], "list")
    expect_setequal(names(sqc.complete.ran$SCC[[1]][[1]]), SCC_names)
    
    lapply(names(sqc.complete.ran$SCC[[1]][[1]]), function(nam){
      expect_equal(colnames(sqc.complete.ran$SCC[[1]][[1]][[nam]]), SCC_partnames[[nam]])  
    })
})

test_that("complete FRIP", {
    expect_null(sqc.complete$FRIP)
    
    expect_equal(names(sqc.complete.ran$FRIP), "CTCF_features")
    expect_equal(names(sqc.complete.ran$FRIP[[1]]), "CTCF_signal")
    expect_is(sqc.complete.ran$FRIP[[1]][[1]], "data.table")
    expect_setequal(colnames(sqc.complete.ran$FRIP[[1]][[1]]), FRIP_names)
    
})

test_that("complete correlation", {
    expect_null(sqc.complete$correlation[[1]])
    expect_null(sqc.complete$correlation[[2]])
    expect_equal(length(sqc.complete$correlation), 2)
    expect_setequal(names(sqc.complete$correlation), c("read_count", "signal_profile"))
    
    expect_equal(names(sqc.complete.ran$correlation$read_count), "CTCF_features")
    expect_equal(names(sqc.complete.ran$correlation$read_count[[1]]), "CTCF_signal")
    expect_setequal(names(sqc.complete.ran$correlation$read_count[[1]][[1]]), c("mat",        "dt",         "row_hclust", "col_hclust"))
    
    expect_equal(names(sqc.complete.ran$correlation$signal_profile), "CTCF_features")
    expect_equal(names(sqc.complete.ran$correlation$signal_profile[[1]]), "CTCF_signal")
    expect_setequal(names(sqc.complete.ran$correlation$signal_profile[[1]][[1]]), c("mat",        "dt",         "row_hclust", "col_hclust"))
    
})

###TODO test for incomplete

