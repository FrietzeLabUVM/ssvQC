testthat::context("bug_fixes")
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
test_that("QcConfigFeatures file can be factor", {
  np_files = dir(system.file(package = "ssvQC", "extdata"), pattern = "Peak$", full.names = TRUE)
  np_conf = QcConfigFeatures.files(factor(np_files), 
                                  group_names = c("10A", "AT1", "CA1"), 
                                  sample_names = c("MCF10A_CTCF", "MCF10AT1_CTCF", "MCF10CA1_CTCF"))
  np_conf$meta_data$file = factor(np_conf$meta_data$file)
  
  bam_files = dir(system.file(package = "ssvQC", "extdata"), 
                  pattern = "CTCF.+bam$", full.names = TRUE)
  bam_conf = QcConfigSignal.files(bam_files)
  bam_conf$meta_data$file = factor(bam_conf$meta_data$file)
  
  sqc = ssvQC(np_conf, bam_conf)
  testthat::expect_is(sqc, "ssvQC")
})

test_that("QcConfigFeatures to_run and reference updates", {
  np_files = dir(system.file(package = "ssvQC", "extdata"), pattern = "Peak$", full.names = TRUE)
  np_conf = QcConfigFeatures.files(factor(np_files), 
                                   group_names = c("MCF10A", "MCF10AT1", "MCF10CA1"),
                                   sample_names = c("MCF10A_CTCF", "MCF10AT1_CTCF", "MCF10CA1_CTCF"))
  
  np_conf$meta_data$mark = "CTCF"
  
  to_run_orig = np_conf$to_run
  to_run_reference_orig = np_conf$to_run_reference
  
  np_conf$to_run_reference = "MCF10A"
  
  to_run_2 = np_conf$to_run
  to_run_reference_2 = np_conf$to_run_reference
  np_conf$run_by = "mark"
  to_run_new = np_conf$to_run
  to_run_reference_new = np_conf$to_run_reference
  
  expect_equal(to_run_orig, c("MCF10A", "MCF10AT1", "MCF10CA1"))
  expect_equal(to_run_reference_orig, character())
  
  expect_equal(to_run_2, c("MCF10AT1", "MCF10CA1"))
  expect_equal(to_run_reference_2, c("MCF10A"))
  
  expect_equal(to_run_new, c("CTCF"))
  expect_equal(to_run_reference_new, character())
  
  expect_error({np_conf$to_run_reference = "CTCF"}, regexp = "1 or nothing will be run")
})

test_that("QcConfigSignal to_run and reference updates", {
  bigwig_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bigwig_config.csv")
  bw_conf = QcConfigSignal.parse(bigwig_config_file)
  bw_conf$meta_data

  bw_conf$run_by = "cell"  
  
  to_run_orig = bw_conf$to_run
  to_run_reference_orig = bw_conf$to_run_reference
  
  bw_conf$to_run_reference = "MCF10A"
  
  to_run_2 = bw_conf$to_run
  to_run_reference_2 = bw_conf$to_run_reference
  bw_conf$run_by = "mark"
  to_run_new = bw_conf$to_run
  to_run_reference_new = bw_conf$to_run_reference
  
  expect_equal(to_run_orig, c("MCF10A", "MCF10AT1", "MCF10CA1"))
  expect_equal(to_run_reference_orig, character())
  
  expect_equal(to_run_2, c("MCF10AT1", "MCF10CA1"))
  expect_equal(to_run_reference_2, c("MCF10A"))
  
  expect_equal(to_run_new, c("CTCF"))
  expect_equal(to_run_reference_new, character())
  
  expect_error({bw_conf$to_run_reference = "CTCF"}, regexp = "1 or nothing will be run")
})
