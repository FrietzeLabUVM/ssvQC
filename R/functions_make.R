
#' make_dt
#'
#' Creates a data.table from a vector of files.
#'
#' @param files character vector of file paths.
#' @param group_lev If provided, should be directory names at equivalent depth
#'   in file paths. Default of NULL disables this functionality and group will
#'   be none and batch will be A.
#' @param max_name_len Names are truncated to this length. Default is 30.
#'
#' @return A data.table with attribute "file" and "name", "group", and "batch" attributes.
#' @export
#'
#' @examples
#' files = c("exp1/file1", "exp1/file2", "exp2/file3")
#' #default no group levels
#' make_dt(files)
#' dt1 = make_dt(files, group_lev = c("exp1", "exp2"))
#' levels(dt1$group)
#' dt1$batch
#' #reversed group_lev
#' dt2 = make_dt(files, group_lev = rev(c("exp1", "exp2")))
#' levels(dt2$group)
#' dt2$batch
make_dt = function(files, group_lev = NULL, max_name_len = 30){
  file_dir = group = batch = name = NULL#global data.table bindings
  p_dt = data.table(file = files)
  p_dt[, sample := sub("\\..+", "", basename(file))]
  p_dt[, sample := sub("(?<=rep[0-9])_.+", "", sample, perl = TRUE)]
  
  p_dt[, file_dir := file]
  
  if(is.null(group_lev)){
    p_dt$group = "none"
    p_dt$batch = factor("A")
  }else{
    while(!any(p_dt$file_dir == "/")){
      if(all(basename(p_dt$file_dir) %in% group_lev)){
        p_dt[, group := basename(file_dir)]
        break
      }else{
        p_dt[, file_dir := dirname(file_dir)]
      }
    }
    p_dt$group = factor(p_dt$group, levels = group_lev)
    p_dt[, batch := LETTERS[as.numeric(group)]]
  }
  p_dt[, name := paste0(substr(sample, 0, max_name_len), ".", batch)]
  dupe_num = 1
  if(any(duplicated(p_dt$name))){
    p_dt[, dupe_num := seq(.N)-1, list(name)]
    p_dt[dupe_num > 0, name := paste0(name, ".", dupe_num) ]
    p_dt$dupe_num = NULL
  }
  p_dt[]
}

#' make_anno_grs
#'
#' Loads a GTF file and processes it to list of GRanges suitable for serial
#' annotation by \code{\link{make_feature_as_signal_dt}}.
#'
#' Output is intended for use with \code{\link{make_feature_as_signal_dt}}
#'
#' @param gtf_file A gtf or gtf.gz file from GENCODE  Other sources may work but
#'   have not been tested.
#'
#' @return a named list of GRanges where each corresponds to an annotated
#'   feature type
#' @export
#'
#' @import rtracklayer GenomicRanges
#'
#' @examples
#' gtf_file = system.file(package = "seqqc", "extdata/gencode.v35.annotation.at_peaks.gtf")
#' make_anno_grs(gtf_file)
make_anno_grs = function(gtf_file){
  type = tag = NULL # global binding for data.table
  ref_gr = rtracklayer::import.gff(gtf_file, format = "gtf")
  gene_gr = reduce(subset(ref_gr, type == "gene"))
  exon_gr = reduce(subset(ref_gr, type == "exon" & tag == "basic"))
  intron_gr = setdiff(gene_gr, exon_gr)
  tx_gr = subset(ref_gr, type == "transcript" & tag == "basic")
  tss_gr = flank(tx_gr, 1e3, start = TRUE, both = TRUE)
  tts_gr = flank(tx_gr, 1e3, start = FALSE, both = TRUE)
  artifact_gr = rtracklayer::import.bed("~/R/qc_cutnrun/reference/blacklist_hg38.bed")
  anno_grs = rev(list(artifact = artifact_gr, tss = tss_gr, tts = tts_gr, exon = exon_gr, intron = intron_gr, genebody = gene_gr))
  anno_grs
}

#' make_anno_dt
#'
#' @param peak_grs a named list of GRanges where each feature set should be
#'   annotated. Use output from \link[seqsetvis]{easyLoad_narrowPeak} or
#'   similar.
#' @param anno_grs a named list of GRanges where each corresponds to an
#'   annotated feature type. Use output from \code{\link{make_anno_grs}}
#' @param name_lev Optional levels to impose order in feature sets.  Default of
#'   NULL uses input order of peak_grs.
#'
#' @return data.table with annotation overlaps for inerval sets in peak_grs.
#' @export
#'
#' @examples
#' gtf_file = system.file(package = "seqqc", "extdata/gencode.v35.annotation.at_peaks.gtf")
#' anno_grs = make_anno_grs(gtf_file)
#'
#' peak_files = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "Peak$", full.names = TRUE)
#' peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files)
#'
#' make_anno_dt(peak_grs, anno_grs)
make_anno_dt =  function(peak_grs, anno_grs, name_lev = NULL){
  if(!is.list(peak_grs)){
    if("GRanges" %in% class(peak_grs)){
      peak_grs = list(peak_grs)
      names(peak_grs) = "peaks"
    }else{
      stop("peak_grs must be a list of GRanges or a single GRanges.")
    }
  }else{
    if(!all(sapply(peak_grs, function(x)"GRanges" %in% class(x)))){
      stop("non-GRanges detected in peak_grs. Must be a list of GRanges or a single GRanges.")
    }
  }
  sample_cnt = count = fraction = type = tag = NULL #data.table global declaration
  if(is.null(names(peak_grs))){
    names(peak_grs) = paste("peaks", LETTERS[seq_along(peak_grs)])
  }
  if(is.null(name_lev)){
    name_lev = names(peak_grs)
  }
  
  message("features overlap")
  .apply_annotation = function(x){
    x$anno = "intergenic"
    for(i in seq_along(anno_grs)){
      overlaps_gr = findOverlaps(x, anno_grs[[i]])
      x$anno[S4Vectors::queryHits(overlaps_gr)] = names(anno_grs)[i]
      
    }
    x
  }
  
  peak_grs.anno = lapply(peak_grs, .apply_annotation)
  
  .dt_count_anno = function(x){
    tab = table(x$anno)
    data.table(feature = names(tab), count = as.numeric(tab))
  }
  
  peak_grs.anno_cnt = lapply(peak_grs.anno, .dt_count_anno)
  
  anno_cnt = rbindlist(peak_grs.anno_cnt, idcol = "sample")
  anno_cnt[, sample_cnt := paste0(sample, "\n", sum(count)), list(sample)]
  anno_cnt[, fraction := count / sum(count), list(sample)]
  stopifnot(setequal(anno_cnt$sample, name_lev))
  anno_cnt$sample = factor(anno_cnt$sample, levels = name_lev)
  anno_cnt = anno_cnt[order(sample)]
  anno_cnt$sample_cnt = factor(anno_cnt$sample_cnt, levels = unique(anno_cnt$sample_cnt))
  anno_cnt
}

#' make_assign_dt
#'
#' @param clust_dt Output from \link[seqsetvis]{ssvSignalClustering}
#' @param cluster_var variable name with cluster assignment info. Default is "cluster_id".
#' @param id_var variable name with id info. Default is "id".
#'
#' @return data.table with id and cluster_id assignment info.
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
#' clust_dt = seqsetvis::ssvSignalClustering(prof_dt, nclust = 3)
#'
#' assign_dt = make_assign_dt(clust_dt)
#' assign_dt
#'
make_assign_dt = function(clust_dt, cluster_var = "cluster_id", id_var = "id"){
  assign_dt = unique(clust_dt[, list(get(cluster_var), get(id_var))])
  setnames(assign_dt, c(cluster_var, id_var))
  assign_dt
}

#' make_feature_as_signal_dt
#'
#' Create a data.table of overlaps between query_gr and anno_grs (from
#' \code{\link{make_anno_grs}} )
#'
#' @param anno_grs named list of GRanges objects where each is an annotation
#'   class.  Use \code{\link{make_anno_grs}} to generate.
#' @param query_gr GRanges of regions to characterize via anno_grs.
#'
#' @return a data.table of counts for each annotation class.
#' @export
#'
#' @import seqsetvis
#'
#' @examples
#' gtf_file = system.file(package = "seqqc", "extdata/gencode.v35.annotation.at_peaks.gtf")
#' anno_grs = make_anno_grs(gtf_file)
#'
#' peak_files = dir(system.file("extdata", package = "seqqc"), pattern = "Peak$", full.names = TRUE)
#' peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files)
#' query_gr = resize(seqsetvis::ssvOverlapIntervalSets(peak_grs), 2e4, fix = "center")
#' anno_dt = make_feature_as_signal_dt(anno_grs, query_gr)
#' anno_dt$id = factor(anno_dt$id, levels = rev(unique(anno_dt$id)))
#' ggplot(anno_dt, aes(x = x, fill = y, y = id)) + geom_tile() + facet_wrap(~feature_type)
make_feature_as_signal_dt = function(anno_grs, query_gr){
  f_dt = seqsetvis::ssvFetchGRanges(anno_grs, query_gr, return_data.table = TRUE)
  setnames(f_dt, "sample", "feature_type")
  f_dt
}


#' make_fq_dt
#'
#' @param fastq_files paths to fast files (can be gzipped with .gz extension)
#' @param fastq_names optional parallel vector of names for fastq files.
#'   Defaults to basename of fastq_files. Should be unique.
#' @param fastq_treatments optional parallel vector of treatments. Defaults to
#'   fastq_names. May be duplicated.
#' @param n_cores number of cores to use to count lines in fastq files. Defaults
#'   to mc.cores if set or 1.
#' @param cache_counts logical. Should the counts be saved to *.cnt files
#'   alongside the fastq_files?  Default is TRUE
#'
#' @return a data.table countaining fastq, count, name, and treatment attributes
#' @export
#'
#' @examples
#' fq_files = dir("inst/extdata", pattern = "(fq$)|(fq.gz$)|(fastq$)|(fastq.gz$)", full.names = TRUE)
#'  #no idea why this make_fq_dt example won't run
#'  #make_fq_dt(fq_files,
#'  #  fastq_names = c("4_reads_fq", "4_reads_gz", "5_reads_fq", "5_reads_gz"),
#'  #  cache_counts = FALSE)
make_fq_dt = function(fastq_files, fastq_names = basename(fastq_files), fastq_treatments = fastq_names, n_cores = getOption("mc.cores", 1), cache_counts = TRUE){
  message("count fastq reads...")
  dt_mode = FALSE
  if(is.data.table(fastq_files)){
    dt_mode = TRUE
    .config_dt = fastq_files
    fastq_files = .get_files(.config_dt)
    .config_dt$file = fastq_files
  }
  
  stopifnot(length(fastq_files) == length(fastq_names))
  stopifnot(length(fastq_files) == length(fastq_treatments))
  .cnt_fq = function(f){
    cnt_f = paste0(f, ".cnt")
    if(!file.exists(cnt_f)){
      if(grepl(".gz$", f)){
        cnt = system(paste0("gunzip -c ", f, "| wc -l"), intern = TRUE)
      }else{
        cnt = system(paste0("cat ", f, "| wc -l"), intern = TRUE)
      }
      
      cnt_dt = data.table(f, as.numeric(cnt)/4)
      if(cache_counts){
        fwrite(cnt_dt, cnt_f, sep = "\t", col.names = FALSE)
      }
    }else{
      cnt_dt = fread(cnt_f, sep = "\t")
    }
    # message(class(cnt_dt))
    cnt_dt
  }
  
  fq_dt = data.table::rbindlist(pbmcapply::pbmclapply(fastq_files, mc.cores = n_cores, .cnt_fq))
  setnames(fq_dt, c("file", "count"))
  
  if(dt_mode){
    fq_dt = merge(fq_dt, .config_dt, by = "file")
  }else{
    fq_dt$name = fastq_names
    fq_dt$treatment = fastq_treatments
  }
  fq_dt
}

.get_fetch_fun = function(files){
  is_bam = grepl(".bam$", files)
  if(all(is_bam)){
    fetch_fun = seqsetvis::ssvFetchBam
  }else if(all(!is_bam)){
    fetch_fun = seqsetvis::ssvFetchBigwig
  }else{
    stop("Files should be either all bams or all bigwigs. No mixing!")
  }
  fetch_fun
}

#' make_centered_query_gr
#'
#' Returns version of query_gr GRanges centered on the maximum signal in
#' bams/bigwigs supplied by query_dt.
#'
#' @param query_dt data.table with query information. Only really needs file as
#'   first column.
#' @param query_gr GRanges to be centered.
#' @param view_size Size of final regions.
#' @param n_cores Number of cores to use. Defaults to mc.cores if set or 1.
#' @param ... passed on the ssvFechBam or ssvFetchBigwig. Do not use, n_cores,
#'   win_size, win_method, return_data.table or n_region_splits.
#'
#' @return GRanges centered on the maximum signal
#' @export
#'
#' @examples
#' library(seqsetvis)
#' #bigwig example with 3 bigwigs
#' bw_files = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "bw$", full.names = TRUE)
#' query_dt = make_dt(bw_files)
#' query_dt[, sample := sub("_FE_random100.A", "", name)]
#'
#' peak_files = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "Peak$", full.names = TRUE)
#' peak_grs = easyLoad_narrowPeak(peak_files)
#' query_gr = resize(ssvOverlapIntervalSets(peak_grs), 700, fix = "center")
#'
#'
#'
#' #heatmap before centering
#' set.seed(0)
#' ssvSignalHeatmap(ssvFetchBigwig(query_dt, query_gr), nclust = 4) +
#'   labs(title = "Before centering")
#'
#' query_gr.centered = make_centered_query_gr(query_dt, query_gr)
#' set.seed(0)
#' ssvSignalHeatmap(ssvFetchBigwig(query_dt, query_gr.centered), nclust = 4)+
#'   labs(title = "After centering")
#'
#' #bam example with 1 bam, query_dt can just be file paths
#' peak_file = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "test_peaks.bed$", full.names = TRUE)
#' bam_file = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "test_peaks.bam$", full.names = TRUE)
#'
#' query_gr2 = easyLoad_bed(peak_file)[[1]]
#' set.seed(0)
#' ssvSignalHeatmap(ssvFetchBam(bam_file, query_gr2), nclust = 4) +
#'   labs(title = "Before centering")
#'
#' query_gr2.centered = make_centered_query_gr(bam_file, query_gr2, fragLens = 180)
#' set.seed(0)
#' ssvSignalHeatmap(ssvFetchBam(bam_file, query_gr2.centered), nclust = 4) +
#'   labs(title = "After centering")
make_centered_query_gr = function(query_dt, query_gr, view_size = NULL, fetch_fun = NULL,
                                  n_cores = getOption("mc.cores", 1), ...){
  if(is.character(query_dt)) query_dt = data.table(file = query_dt)
  stopifnot(is.data.table(query_dt))
  if(!is.null(query_dt$file)){
    files = query_dt$file
  }else{
    files = query_dt[[1]]
  }
  if(is.null(fetch_fun)){
    fetch_fun = .get_fetch_fun(files)
  }
  
  stopifnot("GRanges" %in% class(query_gr))
  if(is.null(view_size)){
    view_size = median(width(query_gr))
  }
  n_region_splits = max(1, floor(length(query_gr) / 1e3))
  if(n_cores == 1) n_region_splits = 1
  if(max(width(query_gr)) > 100){
    message("seeking coarse maxima...")
    ws = ceiling(max(width(query_gr)) / 100)
    if(ws > min(width(query_gr))){
      stop("Too much variation in width of input query_gr.  Please increase minimum or decrease maximum width and retry.")
    }
    prof_dt.coarse = fetch_fun(query_dt, query_gr, n_cores = n_cores,
                               win_size = ws, win_method = "summary",
                               return_data.table = TRUE,
                               n_region_splits = n_region_splits, ...)
    query_gr = unique(resize(centerGRangesAtMax(prof_dt.coarse, query_gr), 100, fix = 'center'))
  }
  message("seeking fine maxima...")
  prof_dt.fine = fetch_fun(query_dt, query_gr, n_cores = n_cores,
                           win_size = 1, win_method = "sample",
                           return_data.table = TRUE,
                           n_region_splits = n_region_splits, ...)
  
  centered_gr = unique(resize(centerGRangesAtMax(prof_dt.fine, query_gr), view_size, fix = 'center'))
  centered_gr
}

#' make_frip_dt
#'
#' create a data.table with FRIP data.  To be used with
#' \code{\link{plot_frip_dt}}
#'
#' @param query_dt data.table with query information. Only really needs file as
#'   first column.
#' @param query_gr GRanges to be centered.
#' @param n_cores Number of cores to use. Defaults to mc.cores if set or 1.
#' @param name_lev Optional levels to impose order in feature sets.  Default of
#'   NULL uses decreasing FRIP score.
#' @return a data.table with reads_in_peak, mapped_reads, and frip data for all
#'   bam files in query_dt at each region in query_gr.
#' @export
#' @import Rsamtools
#' @importFrom stats median
#'
#' @examples
#' peak_file = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "test_peaks.bed$", full.names = TRUE)
#' bam_file = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "test_peaks.bam$", full.names = TRUE)
#'
#' query_gr = seqsetvis::easyLoad_bed(peak_file)[[1]]
#'
#' make_frip_dt(bam_file, query_gr)
make_frip_dt = function(query_dt, query_gr, n_cores = getOption("mc.cores", 1), name_lev = NULL, name_var = "name_split", color_var = name_var){
  treatment = qname = id = frip = N = mapped_reads = V1 = NULL#global data.table bindings
  aes_vars = union(name_var, color_var)
  if(is.character(query_dt)){
    query_dt = data.table(file = query_dt)
  }
  if(is.null(query_dt$file)){
    stop("query_dt must contain file variable")
  }
  if(is.null(query_dt[[name_var]])){
    query_dt[[name_var]] = basename(query_dt$file)
  }
  if(!is.null(name_lev)) stopifnot(all(query_dt[[name_var]] %in% name_lev))
  if(is.null(query_dt[[color_var]])){
    query_dt[[color_var]] = query_dt[[name_var]]
  }
  
  n_region_splits = max(1, floor(length(query_gr) / 1e3))
  message("fetch read counts...")
  # reads_dt = seqsetvis::ssvFetchBam(query_dt, query_gr, fragLens = NA, return_unprocessed = TRUE, n_region_splits = n_region_splits, n_cores = n_cores)
  # frip_dt = reads_dt[, list(N = length(unique(qname))), list(id, name, treatment, sample)]
  frip_dt = seqsetvis::ssvFetchBam(query_dt, query_gr, fragLens = 1, win_size = 1, win_method = "summary", summary_FUN = function(x,w){sum(x)}, n_region_splits = n_region_splits, n_cores = n_cores, return_data.table = TRUE)
  setnames(frip_dt, "y", "N")
  
  frip_dt_filled = melt(dcast(frip_dt, paste0("id~", name_var), value.var = "N", fill = 0), id.vars = "id", value.name = "N", variable.name = name_var)
  frip_dt = merge(frip_dt_filled, unique(frip_dt[, c(aes_vars, "sample"), with = FALSE]), by = name_var)
  
  message("fetch total mapped reads...")
  mapped_counts = sapply(query_dt$file, function(f){
    stats = Rsamtools::idxstatsBam(f)
    stats = subset(stats, grepl("chr[0-9XY]+$", seqnames ))
    sum(stats[,3])
  })
  frip_dt$mapped_reads = mapped_counts[frip_dt$sample]
  frip_dt[, frip := N/mapped_reads]
  
  if(!is.null(name_lev)){
    # name_lev = frip_dt[, stats::median(N) , list(name)][rev(order(V1))]$name
    stopifnot(all(frip_dt[[name_var]] %in% name_lev))
    frip_dt[[name_var]] = factor(frip_dt[[name_var]], levels = name_lev)
  }
  
  
  setnames(frip_dt, "N", "reads_in_peak")
  frip_dt
  
}

#' make_peak_dt
#'
#' @param peak_grs list of GRanges for peak sets
#' @param treatments opional character vector of treatments for each peak set.
#'
#' @return a data.table with peak_count data for each GRanges in peak_grs.
#' @export
#' @import seqsetvis
#' @examples
#' peak_files = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "Peak$", full.names = TRUE)
#' peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files)
#'
#' make_peak_dt(peak_grs)
make_peak_dt = function(peak_grs, treatments = NULL){
  treatment = NULL#global binding for data.table
  peak_dt = seqsetvis::ssvFeatureBars(peak_grs, return_data = TRUE)
  if(is.null(treatments)){
    peak_dt$treatment = peak_dt$group
  }else{
    stopifnot(nrow(peak_dt) == length(treatments))
    peak_dt$treatment = treatment
  }
  if(is.null(names(peak_grs))){
    peak_dt$name = peak_dt$treatment
  }else{
    peak_dt$name = names(peak_grs)
  }
  setnames(peak_dt, "count", "peak_count")
  peak_dt
}

.get_files = function(query_dt, name_var = "name_split"){
  if(!is.null(query_dt$file)){
    files = query_dt$file
  }else{
    files = query_dt[[1]]
  }
  if(!is.null(query_dt[[name_var]])){
    uniq_names = query_dt[[name_var]]
  }else if(!is.null(query_dt$name)){
    uniq_names = query_dt$name
  }else if(!is.null(query_dt$sample)){
    uniq_names = query_dt$sample
  }else{
    uniq_names = basename(files)
  }
  stopifnot(!any(duplicated(uniq_names)))
  names(files) = uniq_names
  files
}

#' make_scc_dt
#'
#' Calculate Strand Cross Correlation (SCC) for several bam files defined in
#' query_dt at regions defined by query_gr and return tidy data.table.
#'
#' @param query_dt data.table with query information. Only really needs file as
#'   first column.
#' @param query_gr GRanges of regions to calculate SCC for
#' @param frag_sizes optional numeric. Fragment sizes to calculate correlation
#'   at.  The higher the resolution the longer calculation will take.  The
#'   default is to count by 10 from 50 to 350.
#' @param fetch_size optional numeric. Size in bp centered around each interval
#'   in query_gr to retrieve.  Should be greater than max frag_size. The default
#'   is 3*max(frag_sizes).
#' @param bfc_corr BiocFileCache object to use.
#' @param cache_version Modifying the cache version will force recalulation of
#'   all results going forward. Default is v1.
#' @param force_overwrite Logical, if TRUE, cache contents will be overwritten.
#' @param name_var Character. Variable where name information is stored.
#' @param n_cores Number of cores to use. Defaults to mc.cores if set or 1.
#' @param ... passed to Rsamtools::ScanBamParam()
#'
#' @return list fo tidy data.table of SCC data for every bam file in query_dt
#' @export
#'
#' @examples
#' peak_file = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "test_peaks.bed$", full.names = TRUE)
#' bam_file = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "test_peaks.bam$", full.names = TRUE)
#'
#' query_gr = seqsetvis::easyLoad_bed(peak_file)[[1]]
#' query_dt = data.table(file = rep(bam_file, 2))
#' #attributes added here will be carried forward to the final data.tables
#' query_dt$name = c("A", "B")
#' query_dt$value = c(.3, .7)
#'
#' scc_res = make_scc_dt(query_dt, query_gr)
#'
#' #making the data.table can be skipped for conveinence
#' #also run at higher resolution if a more precise estimate of fragment size is required
#' scc_res = make_scc_dt(bam_file, query_gr, frag_sizes = seq(150, 210, 1))
make_scc_dt = function(query_dt,
                       query_gr,
                       frag_sizes = seq(50, 350, 10),
                       fetch_size = 3*max(frag_sizes),
                       bfc_corr = BiocFileCache("~/.cache_peakrefine"),
                       cache_version = "v1",
                       force_overwrite = FALSE,
                       name_var = "name_split",
                       n_cores = getOption("mc.cores", 1L),
                       ...){
  
  if(is.character(query_dt)) query_dt = data.table(file = query_dt, name = basename(query_dt))
  if(!name_var %in% colnames(query_dt)){
    stop(
      paste(collapse = "\n",
            c(paste0("name_var \"", name_var, "\" was not found in colnames of query_dt!"),
              colnames(query_dt))
      )
    )
  }
  stopifnot()
  files = .get_files(query_dt)
  scc_res_l = lapply(files, function(f){
    message("SCC for ", f)
    make_scc_dt.single(f, query_gr,
                       frag_sizes = frag_sizes,
                       fetch_size = fetch_size,
                       bfc_corr = bfc_corr,
                       cache_version = cache_version,
                       force_overwrite = force_overwrite,
                       n_cores = n_cores,
                       ...)
  })
  
  vnames = names(scc_res_l[[1]])
  
  scc_res = lapply(vnames, function(nam){
    part = lapply(scc_res_l, function(x){
      xv = x[[nam]]
      # if(!is.data.table(xv)){
      #   xv = data.table(xv)
      #   setnames(xv, nam)
      # }
      xv
    })
    dt = rbindlist(part, idcol = name_var)
    merge(dt, query_dt, by = name_var)
  })
  names(scc_res) = vnames
  scc_res
}

#' make_scc_dt.single
#'
#' based on peakrefine::calcSCCMetrics
#'
#' @param bam_file a single path to a bam file
#' @param query_gr GRanges of regions to calculate SCC for
#' @param frag_sizes optional numeric. Fragment sizes to calculate correlation
#'   at.  The higher the resolution the longer calculation will take.  The
#'   default is to count by 10 from 50 to 350.
#' @param fetch_size optional numeric. Size in bp centered around each interval
#'   in query_gr to retrieve.  Should be greater than max frag_size. The default
#'   is 3*max(frag_sizes).
#' @param bfc_corr BiocFileCache object to use.
#' @param cache_version Modifying the cache version will force recalulation of
#'   all results going forward. Default is v1.
#' @param force_overwrite Logical, if TRUE, cache contents will be overwritten.
#' @param n_cores Number of cores to use. Defaults to mc.cores if set or 1.
#' @param ... passed to Rsamtools::ScanBamParam()
#'
#' @return list fo tidy data.table of SCC data for bam_file
#'
#' @import tools digest pbmcapply
#'
#' @examples
#' peak_file = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "test_peaks.bed$", full.names = TRUE)
#' bam_file = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "test_peaks.bam$", full.names = TRUE)
#'
#' query_gr = seqsetvis::easyLoad_bed(peak_file)[[1]]
#'
#' scc_res = seqqc:::make_scc_dt.single(bam_file, query_gr, seq(50,300, by = 10))
make_scc_dt.single = function(bam_file,
                              query_gr,
                              frag_sizes,
                              fetch_size = 3*max(frag_sizes),
                              bfc_corr = BiocFileCache("~/.cache_peakrefine"),
                              cache_version = "v1",
                              force_overwrite = FALSE,
                              n_cores = getOption("mc.cores", 1L),
                              ...){
  bam_md5 = NULL
  qgr_md5 = NULL
  if(is.null(bam_md5)){
    bam_md5 = tools::md5sum(bam_file)
  }
  if(is.null(qgr_md5)){
    qgr_md5 = digest::digest(query_gr)
  }
  if(fetch_size <= max(frag_sizes)){
    stop("fetch_size (", fetch_size, ") must exceed max of frag_sizes (", max(frag_sizes), ").")
  }
  stopifnot(file.exists(bam_file))
  stopifnot(class(query_gr) == "GRanges")
  if(!file.exists(paste0(bam_file, ".bai"))){
    stop("bam_file index not found. expected at ", paste0(bam_file, ".bai"),
         "\ntry running:\nsamtools index ", bam_file)
  }
  stopifnot(n_cores >= 1)
  
  query_gr = resize(query_gr, fetch_size, fix = 'center')
  if(is.null(query_gr$name)){
    if(is.null(names(query_gr))){
      query_gr$name = paste0("peak_", seq_along(query_gr))
    }else{
      query_gr$name = names(query_gr)
    }
  }else{
    if(is.null(names(query_gr))){
      names(query_gr) = query_gr$name
    }else{
      #both names() and $name are set, leave it alone
    }
  }
  stopifnot(!any(duplicated(names(query_gr))))
  
  corr_key = paste(qgr_md5, bam_md5, digest::digest(frag_sizes), fetch_size, cache_version, sep = "_")
  corr_res = bfcif(bfc_corr, corr_key, function(){
    message("cached results not found, gathering correlation info.")
    nper = ceiling(length(query_gr) / n_cores)
    grps = ceiling(seq_along(query_gr)/ nper)
    table(grps)
    # browser()
    rl = getReadLength(bam_file, query_gr)
    lres = pbmcapply::pbmclapply(unique(grps), function(g){
      k = grps == g
      crossCorrByRle(bam_file, query_gr[k], fragment_sizes = frag_sizes, read_length = rl, ...)
    })
    peak_strand_corr = rbindlist(lres)
    gather_metrics(peak_strand_corr, rl)
  }, force_overwrite = force_overwrite)
  
  #make everything a data.table
  vnames = names(corr_res)
  corr_res = lapply(vnames, function(nam){
    xv = corr_res[[nam]]
    if(!is.data.table(xv)){
      xv = data.table(xv)
      setnames(xv, nam)
    }
    xv
  })
  names(corr_res) = vnames
  corr_res
}

#' make_feature_overlap_signal_profiles
#'
#' Creates a data.table with grouping for profiles based on overlaps in overlap_gr and ranking of signal within each.
#'
#' @param query_dt data.table that specifies file paths with $file.  Suitable for use with \link[seqsetvis]{ssvFetchBam} or \link[seqsetvis]{ssvFetchBigwig}
#' @param overlaps_gr GRanges with membership table in mcols.  Create with \link[seqsetvis]{ssvOverlapIntervalSets} or \link[seqsetvis]{ssvConsensusIntervalSets}.
#' @param feature_groups Optional data.frame defining "id" and "group". Default of NULL uses \link[seqsetvis]{ssvFactorizeMembTable} on overlaps_gr
#' @param view_size Size of regions to fetch.  Will be median of widths of overlaps_gr if NULL. Default is NULL.
#' @param group_var character. Must not already be present in query_dt  Will be added. Default is "overlap_group".
#' @param rank_var character. Must not already be present in query_dt  Will be added. Default is "rnk".
#' @param min_group_size Groups smaller than this are allowed but larger groups won't be shrunk below this length. Default is 50.
#' @param max_group_size Even if the smallest group is larger than this, no group will exceed this size. Default is 500.
#' @param discard_group_below Groups smaller than this will be discarded entirely. Default is 1.
#' @param ... Passed to \link[seqsetvis]{ssvFetchBam} or \link[seqsetvis]{ssvFetchBigwig} as appropriate to file type.
#'
#' @return A profile data.table with overlap grouping information add and independent ranking added to each group.
#' @export
#'
#' @examples
#' bw_files = dir(system.file("extdata", package = "seqqc"), pattern = "^M.+bw$", full.names = TRUE)
#' query_dt = make_dt(bw_files)
#' query_dt[, sample := sub("_FE_random100.A", "", name)]
#'
#' peak_files = dir(system.file("extdata", package = "seqqc"), pattern = "Peak$", full.names = TRUE)
#' names(peak_files) = sub("_CTCF_rand.+", "", basename(peak_files))
#' peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files, )
#' overlaps_gr = seqsetvis::ssvOverlapIntervalSets(peak_grs)
#' query_gr = resize(overlaps_gr, 6e2, fix = "center")
#'
#' group_prof_dt = make_feature_overlap_signal_profiles(query_dt, overlaps_gr)
#' plot_feature_overlap_signal_profiles(group_prof_dt)
make_feature_overlap_signal_profiles = function(query_dt,
                                                overlaps_gr,
                                                feature_groups = NULL,
                                                view_size = NULL,
                                                group_var = "overlap_group",
                                                rank_var = "rnk",
                                                min_group_size = 50,
                                                max_group_size = 500,
                                                discard_group_below = 1,
                                                ...){
  y = rnk_ = id = NULL #global binding for data.table
  if(is.character(query_dt)) query_dt = data.table(file = query_dt)
  if(is.null(view_size)) view_size = median(width(overlaps_gr))
  
  if(group_var %in% colnames(query_dt)) stop("group_var \"", group_var, "\" cannot be in colnames of query_dt.")
  if(rank_var %in% colnames(query_dt)) stop("rank_var \"", rank_var, "\" cannot be in colnames of query_dt.")
  
  if(is.null(feature_groups)){
    feature_groups = as.data.table(seqsetvis::ssvFactorizeMembTable(overlaps_gr))
  }else{
    stopifnot(c("id", "group") %in% colnames(feature_groups))
  }
  
  # hist(width(overlaps_gr))
  qgr = resize(overlaps_gr, view_size, fix = "center")
  qgr_l = split(qgr[feature_groups$id], feature_groups$group)
  qgr_l$none = NULL
  
  sm_len = max(min(lengths(qgr_l)), min_group_size)
  qgr_l.sample = lapply(qgr_l, function(x){sample(x, min(max_group_size, sm_len, length(x)))})
  qgr_l.sample = qgr_l.sample[!lengths(qgr_l.sample) < discard_group_below]
  
  files = .get_files(query_dt)
  fetch_fun = .get_fetch_fun(files)
  
  fetch_list_profiles = function(qgr_list){
    prof_dt_l = lapply(qgr_list, function(x){
      prof_dt = suppressMessages({
        fetch_fun(query_dt,
                  unique_names = query_dt$name,
                  resize(x, view_size, fix = "center"),
                  return_data.table = TRUE,
                  ...)
      })
      
      prof_dt
    })
    
    prof_dt = rbindlist(prof_dt_l, idcol = group_var)
    rnk_dt = prof_dt[, list(y = max(y)), c("id", group_var)]
    rnk_dt[, rnk_ := frank(-y, ties.method = "first"), c(group_var)]
    prof_dt = merge(prof_dt, rnk_dt[, list(id, rnk_)], by = "id")
    setnames(prof_dt, "rnk_", rank_var)
    
    prof_dt
  }
  prof_dt = fetch_list_profiles(qgr_l.sample)
  prof_dt
}
