#' write_bed_frip
#'
#' @param query_gr GRanges of regions used to make frip_dt and assign_dt
#' @param assign_dt data.table with "id" and "cluster_id" that maps region ids to clusters. In output from \code{\link{plot_signals}}.
#' @param frip_dt output from \code{\link{make_frip_dt}}
#' @param file File to write to
#'
#' @return invisibly returns data that is written to file
#' @export
#'
#' @examples
#' bw_files = dir(system.file("extdata", package = "seqqc"), pattern = "^M.+bw$", full.names = TRUE)
#' query_dt = make_dt(bw_files)
#' query_dt[, sample := sub("_FE_random100.A", "", name)]
#'
#' peak_files = dir(system.file("extdata", package = "seqqc"), pattern = "Peak$", full.names = TRUE)
#' peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files)
#' query_gr = resize(seqsetvis::ssvOverlapIntervalSets(peak_grs), 6e2, fix = "center")
#'
#' prof_dt = seqsetvis::ssvFetchBigwig(query_dt, query_gr, return_data.table = TRUE)
#'
#' bam_files = dir(system.file("extdata", package = "seqqc"), pattern = "^M.+bam$", full.names = TRUE)
#' query_dt.bam = make_dt(bam_files)
#' frip_dt = make_frip_dt(query_dt.bam, query_gr)
#'
#' sig_res = plot_signals(prof_dt, query_gr)
#'
#' assign_dt = sig_res$cluster_assignment
#'
#' outdir = system.file("extdata", package = "seqqc")
#' write_bed_frip(query_gr, assign_dt, frip_dt,
#'   file = file.path(outdir, "regions_with_FRIP.txt"))
write_bed_frip = function(query_gr, assign_dt, frip_dt, file = "regions_with_FRIP.txt"){
  id = name = reads_in_peak = treatment = frip = NULL #global binding for data.table
  region_frip_dt = dcast(frip_dt[, list(id, name, reads_in_peak, treatment, frip)], id~name, value.var = "frip")
  setkey(region_frip_dt, id)
  setkey(assign_dt, id)
  bed_frip_towrite = query_gr

  mcols(bed_frip_towrite) = region_frip_dt[list(names(bed_frip_towrite))]
  # bed_frip_towrite$id = NULL

  message("write files...")
  bed_frip_towrite$cluster_id = assign_dt[list(names(bed_frip_towrite))]$cluster_id
  # rtracklayer::export.bed(bed_frip_towrite, paste0("consensus_regions_with_FRIP.", file_tag, ".bed"))

  bed_frip_towrite_dt = as.data.table(bed_frip_towrite)
  bed_frip_towrite_dt$strand = "."
  bed_frip_towrite_dt$score = 0
  extra_cols = setdiff(colnames(mcols(bed_frip_towrite)), "id")
  extra_cols = gsub("-", ".", extra_cols)
  fwrite(bed_frip_towrite_dt[, c("seqnames", "start", "end", "id", "score", "strand", extra_cols), with = FALSE],
         file = file,
         sep = "\t")
  invisible(bed_frip_towrite_dt)
}

#' write_bed_overlaps
#'
#' @param overlaps_gr GRanges of regions used to make assign_dt with memb table in mcols.
#' @param assign_dt data.table with "id" and "cluster_id" that maps region ids to clusters. In output from \code{\link{plot_signals}}.
#' @param file File to write to
#'
#' @return invisibly returns data that is written to file
#' @export
#'
#' @examples
#' bw_files = dir(system.file("extdata", package = "seqqc"), pattern = "^M.+bw$", full.names = TRUE)
#' query_dt = make_dt(bw_files)
#' query_dt[, sample := sub("_FE_random100.A", "", name)]
#'
#' peak_files = dir(system.file("extdata", package = "seqqc"), pattern = "Peak$", full.names = TRUE)
#' peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files)
#' overlaps_gr = seqsetvis::ssvOverlapIntervalSets(peak_grs)
#' query_gr = resize(overlaps_gr, 6e2, fix = "center")
#'
#' prof_dt = seqsetvis::ssvFetchBigwig(query_dt, query_gr, return_data.table = TRUE)
#'
#' sig_res = plot_signals(prof_dt, query_gr)
#' assign_dt = sig_res$cluster_assignment
#'
#' outdir = system.file("extdata", package = "seqqc")
#' write_bed_overlaps(overlaps_gr, assign_dt,
#'   file = file.path(outdir, "regions_with_overlaps.txt"))
write_bed_overlaps = function(overlaps_gr, assign_dt, file = "regions_with_overlaps.txt"){
  id = NULL #global binding for data.table
  setkey(assign_dt, id)
  bed_peak_towrite = overlaps_gr
  bed_peak_towrite$cluster_id = assign_dt[list(names(bed_peak_towrite))]$cluster_id
  bed_peak_towrite_dt = as.data.table(bed_peak_towrite)
  bed_peak_towrite_dt$id = names(bed_peak_towrite)
  bed_peak_towrite_dt$strand = "."
  bed_peak_towrite_dt$score = 0
  extra_cols = setdiff(colnames(mcols(bed_peak_towrite)), "id")
  extra_cols = gsub("-", ".", extra_cols)
  fwrite(bed_peak_towrite_dt[, c("seqnames", "start", "end", "id", "score", "strand", extra_cols), with = FALSE],
         file = file,
         sep = "\t")
}
