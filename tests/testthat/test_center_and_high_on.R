testthat::context("center_and_high_on")
library(ssvQC)
library(testthat)
options(mc.cores = 1)
SQC_OPTIONS$SQC_FORCE_CACHE_OVERWRITE = TRUE
set.seed(0)

test_that("ssvQC center_at_max and high_on_right", {
  #### setup configs ####
  features_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
  features_config = QcConfigFeatures.parse(features_config_file)
  features_config$n_peaks = 10
  
  bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
  bam_config = QcConfigSignal.parse(bam_config_file)
  bam_config$view_size = 600
  bam_config$sort_value = SQC_SIGNAL_VALUES$raw
  bam_config$cluster_value = SQC_SIGNAL_VALUES$raw
  bam_config$plot_value = SQC_SIGNAL_VALUES$raw
  bam_config$heatmap_limit_values = c(0, 100)
  
  sqc.default = ssvQC(
    features_config,
    bam_config
  )
  expect_equal(sqc.default$signal_config$center_signal_at_max, FALSE)
  expect_equal(sqc.default$signal_config$flip_signal_mode, "none")
  
  bam_config$center_signal_at_max = TRUE
  sqc.default.centered = ssvQC(
    features_config,
    bam_config,
  )
  expect_equal(sqc.default.centered$signal_config$center_signal_at_max, TRUE)
  expect_equal(sqc.default.centered$signal_config$flip_signal_mode, "none")
  
  bam_config$center_signal_at_max = FALSE
  bam_config$flip_signal_mode = SQC_FLIP_SIGNAL_MODES$high_on_left
  sqc.high_on_left = ssvQC(
    features_config,
    bam_config,
  )
  expect_equal(sqc.high_on_left$signal_config$center_signal_at_max, FALSE)
  expect_equal(sqc.high_on_left$signal_config$flip_signal_mode, "high_on_left")
  
  bam_config$flip_signal_mode = SQC_FLIP_SIGNAL_MODES$high_on_right
  sqc.high_on_right = ssvQC(
    features_config,
    bam_config,
  )
  expect_equal(sqc.high_on_right$signal_config$center_signal_at_max, FALSE)
  expect_equal(sqc.high_on_right$signal_config$flip_signal_mode, "high_on_right")
  
  bam_config$center_signal_at_max = TRUE
  sqc.high_on_right.centered = ssvQC(
    features_config,
    bam_config,
  )
  expect_equal(sqc.high_on_right.centered$signal_config$center_signal_at_max, TRUE)
  expect_equal(sqc.high_on_right.centered$signal_config$flip_signal_mode, "high_on_right")
  
  set.seed(0)
  sqc.default = 
    ssvQC.plotSignal(sqc.default)
  set.seed(0)
  sqc.default.centered = 
    ssvQC.plotSignal(sqc.default.centered)
  set.seed(0)
  sqc.high_on_left = 
    ssvQC.plotSignal(sqc.high_on_left)
  set.seed(0)
  sqc.high_on_right = 
    ssvQC.plotSignal(sqc.high_on_right)
  set.seed(0)
  sqc.high_on_right.centered = 
    ssvQC.plotSignal(sqc.high_on_right.centered)
  #### test for   query_gr ####
  #test start change properly
  expect_false(all(  
    start(sqc.default$signal_data$CTCF_features$CTCF_signal$query_gr) == 
      start(sqc.default.centered$features_config$assessment_features$CTCF_features)
  ))
  
  expect_true(all(  
    start(sqc.default$signal_data$CTCF_features$CTCF_signal$query_gr) == 
      start(sqc.high_on_left$features_config$assessment_features$CTCF_features)
  ))
  
  expect_true(all(  
    start(sqc.default$signal_data$CTCF_features$CTCF_signal$query_gr) == 
      start(sqc.high_on_right$features_config$assessment_features$CTCF_features)
  ))
  
  expect_false(all(  
    start(sqc.default$signal_data$CTCF_features$CTCF_signal$query_gr) == 
      start(sqc.high_on_right.centered$features_config$assessment_features$CTCF_features)
  ))
  #test strand change properly
  expect_true(all(  
    strand(sqc.default$signal_data$CTCF_features$CTCF_signal$query_gr) == 
      strand(sqc.default.centered$features_config$assessment_features$CTCF_features)
  ))
  
  expect_false(all(  
    strand(sqc.default$signal_data$CTCF_features$CTCF_signal$query_gr) == 
      strand(sqc.high_on_left$features_config$assessment_features$CTCF_features)
  ))
  
  expect_false(all(  
    strand(sqc.default$signal_data$CTCF_features$CTCF_signal$query_gr) == 
      strand(sqc.high_on_right$features_config$assessment_features$CTCF_features)
  ))
  
  expect_false(all(  
    strand(sqc.default$signal_data$CTCF_features$CTCF_signal$query_gr) == 
      strand(sqc.high_on_right.centered$features_config$assessment_features$CTCF_features)
  ))
  #high_on_right and high_on_left should have no match
  expect_false(all(  
    strand(sqc.high_on_right$features_config$assessment_features$CTCF_features) == 
      strand(sqc.high_on_left$features_config$assessment_features$CTCF_features)
  ))
  #### test for   signal_data ####
  prof_dt.default = sqc.default$signal_data$CTCF_features$CTCF_signal$signal_data
  prof_dt.default$group = "default"
  prof_dt.default.centered = sqc.default.centered$signal_data$CTCF_features$CTCF_signal$signal_data
  prof_dt.default.centered$group = "centered"
  prof_dt.high_on_left = sqc.high_on_left$signal_data$CTCF_features$CTCF_signal$signal_data
  prof_dt.high_on_left$group = "high_on_left"
  prof_dt.high_on_right = sqc.high_on_right$signal_data$CTCF_features$CTCF_signal$signal_data
  prof_dt.high_on_right$group = "high_on_right"
  prof_dt.high_on_right.centered = sqc.high_on_right.centered$signal_data$CTCF_features$CTCF_signal$signal_data
  prof_dt.high_on_right.centered$group = "high_on_right.centered"
  
  sel_id = "region_1"
  sel_name = "MCF10A_CTCF_rep1"
  plot_dt = rbind(
    prof_dt.default[id == sel_id & name == sel_name],
    prof_dt.default.centered[id == sel_id & name == sel_name],
    prof_dt.high_on_left[id == sel_id & name == sel_name],
    prof_dt.high_on_right[id == sel_id & name == sel_name],
    prof_dt.high_on_right.centered[id == sel_id & name == sel_name]
  )
  plot_dt$group = factor(plot_dt$group)
  plot_dt[, y_shift := y + 4*as.numeric(group)]
  ggplot(plot_dt, aes(x = x, y = y_shift, color = group, group = group)) +
    geom_path() +
    facet_wrap(~name)
  
  expect_true(setequal(prof_dt.default$y, prof_dt.high_on_left$y))
  expect_true(setequal(prof_dt.default$y, prof_dt.high_on_right$y))
  
  expect_true(setequal(prof_dt.default.centered$y, prof_dt.high_on_right.centered$y))
  
  expect_false(setequal(prof_dt.default$y, prof_dt.default.centered$y))
  expect_false(setequal(prof_dt.default$start, prof_dt.default.centered$start))

  expect_equal(
    sqc.default$signal_data$CTCF_features$CTCF_signal$query_gr,
    resize(sqc.default$features_config$assessment_features$CTCF_features, 600, fix = "center")
  )  
  expect_equal(
    sqc.default.centered$signal_data$CTCF_features$CTCF_signal$query_gr,
    sqc.default.centered$features_config$assessment_features$CTCF_features
  )
  expect_equal(
    sqc.high_on_left$signal_data$CTCF_features$CTCF_signal$query_gr,
    sqc.high_on_left$features_config$assessment_features$CTCF_features
  )
  expect_equal(
    sqc.high_on_right$signal_data$CTCF_features$CTCF_signal$query_gr,
    sqc.high_on_right$features_config$assessment_features$CTCF_features
  )
  expect_equal(
    sqc.high_on_right.centered$signal_data$CTCF_features$CTCF_signal$query_gr,
    sqc.high_on_right.centered$features_config$assessment_features$CTCF_features
  )
})

# cowplot::plot_grid(
#   sqc.default$plots$signal$heatmaps$CTCF_features$CTCF_signal + cowplot::draw_text(x = 1, y = 1, hjust = 1, vjust = 1, text = "default"),
#   sqc.default.centered$plots$signal$heatmaps$CTCF_features$CTCF_signal  + cowplot::draw_text(x = 1, y = 1, hjust = 1, vjust = 1, text = "centered"),
#   sqc.high_on_left$plots$signal$heatmaps$CTCF_features$CTCF_signal  + cowplot::draw_text(x = 1, y = 1, hjust = 1, vjust = 1, text = "high on left"),
#   sqc.high_on_right$plots$signal$heatmaps$CTCF_features$CTCF_signal  + cowplot::draw_text(x = 1, y = 1, hjust = 1, vjust = 1, text = "high on right")
# )


