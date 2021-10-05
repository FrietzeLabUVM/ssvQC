?QcConfigSignal

bigwig_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bigwig_config.csv")
bw_conf = QcConfigSignal.parse(bigwig_config_file)

np_files = dir(system.file(package = "ssvQC", "extdata"), pattern = "Peak$", full.names = TRUE)
np_conf = QcConfigFeatures.files(factor(np_files),
                                 group_names = c("MCF10A", "MCF10AT1", "MCF10CA1a"),
                                 sample_names = c("MCF10A_CTCF", "MCF10AT1_CTCF", "MCF10CA1a_CTCF"))
np_conf$meta_data$mark = "CTCF"
np_conf$run_by = "mark"
#
# to_run_orig = np_conf$to_run
# to_run_reference_orig = np_conf$to_run_reference
# np_conf$run_by = "mark"
# to_run_new = np_conf$to_run
# to_run_reference_new = np_conf$to_run_reference

bw_conf$to_run
np_conf$to_run

sqc = ssvQC(np_conf, bw_conf)

sqc = ssvQC.plotMappedReads(sqc)

sqc = ssvQC.prepFeatures(sqc)
sqc = ssvQC.prepCapValue(sqc)
sqc = ssvQC.prepSignal(sqc)
sqc = ssvQC.plotSignal(sqc)

sqc = ssvQC.runAll(sqc)

sqc$plots$reads$CTCF_signal
