# testthat::context("Test1")
# library(seqqc)
# library(testthat)
# library(data.table)
#
# #sampling
# test_that("test 1", {
# })
#informal tests
if(FALSE){
    anno_grs = make_anno_grs("~/gencode.v35.annotation.gtf.gz")
    stopifnot(is.list(anno_grs))

    bw_files = dir(system.file("extdata", package = "seqsetvis"), pattern = "M.+bw$", full.names = TRUE)
    qdt = data.table(file = bw_files)
    qdt[, cell := tstrsplit(basename(file), "_", keep = 1)]

    peak_grs = CTCF_in_10a_narrowPeak_grs
    qgr = CTCF_in_10a_overlaps_gr



    qgr.c = make_centered_query_gr(qdt, qgr = qgr)
    stopifnot(!all(start(qgr.c) == start(qgr)))

    make_dt(bw_files, basename(bw_files))
    make_dt(bw_files, "extdata")

    fq_files = dir("/slipstream/home/joeboyd/rausch_atac_pipeline/ATACseq/output_MCF10A_progression.hg38/MCF10A_rep1_001",
                   pattern = "fq", full.names = TRUE)
    fq_files = fq_files[!grepl(".cnt$", fq_files)]

    make_fq_dt(fq_files, fastq_names = LETTERS[seq_along(fq_files)], fastq_treatments = 1:3)

    bam_files = c(
        "/slipstream/home/joeboyd/rausch_atac_pipeline/ATACseq/output_MCF10A_progression.hg38/MCF10A_rep1_001/MCF10A_rep1_001.final.bam",
        "/slipstream/home/joeboyd/rausch_atac_pipeline/ATACseq/output_MCF10A_progression.hg38/MCF10A_rep2/MCF10A_rep2.final.bam")

    # qdt_frip = data.table(
    #     file = bam_files
    # )
    qdt_frip = make_dt(bam_files, group_lev = c("MCF10A_rep1_001", "MCF10A_rep2"))
    # qdt_frip[, name := sub(".bam", "", basename(file))]
    # qdt_frip[, name := basename(file)]
    # qdt_frip$treatment = LETTERS[1:2]

    peak_files = c(
        "/slipstream/home/joeboyd/rausch_atac_pipeline/ATACseq/output_MCF10A_progression.hg38/MCF10A_rep1_001/MCF10A_rep1_001.final.peaks",
        "/slipstream/home/joeboyd/rausch_atac_pipeline/ATACseq/output_MCF10A_progression.hg38/MCF10A_rep2/MCF10A_rep2.final.peaks"
    )

    qdt_peaks = make_dt(peak_files, group_lev = c("MCF10A_rep1_001", "MCF10A_rep2"))
    peaks_frip = seqsetvis::easyLoad_narrowPeak(qdt_peaks$file, file_names = qdt_peaks$name)
    qgr_frip = seqsetvis::ssvConsensusIntervalSets(peaks_frip)

    frip_dt = make_frip_dt(qdt_frip, qgr_frip, is_PE = FALSE)
    peak_dt = make_peak_dt(peaks_frip)
    # peak_dt$name = qdt_frip$name

    frip_plots = plot_frip_metrics(frip_dt, peak_dt, name_lev = peak_dt$name)

    frip_plots$aligned_reads
    frip_plots$peaks
    frip_plots$reads_per_peaks
    frip_plots$frip_per_peaks
    frip_plots$frip_total

    peak_overlap_plots = plot_features_overlap(peaks_frip, anno_grs)
    peak_comparison_plots = plot_feature_comparison(peaks_frip)

    sel_qgr = sample(qgr_frip, 100)

    prof_dt = seqsetvis::ssvFetchBam(qdt_frip, sel_qgr, fragLens = 210, return_data.table = TRUE)
    anno_dt = make_anno_dt(anno_grs, qgr_frip)

    singal_plots = plot_signals(prof_dt, sel_qgr, frip_dt = frip_dt, anno_dt = anno_dt)

    write_bed_frip()
    write_bed_overlaps()
}
# fwrite(meta_dt, paste0("sample_info.", file_tag, ".csv"), sep = ",")
# ggsave(paste0("frip_heatmap.", file_tag, ".pdf"), pg_heat, width = 8+length(peak_grs), height = 14)
# ggsave(paste0("feature_overlap.", file_tag, ".pdf"), p_features, width = 2+length(peak_grs)*.6, height = 6)

# fwrite(anno_cnt,
#        file = paste0("feature_counts.", file_tag, ".csv"),
#        sep = ",")
