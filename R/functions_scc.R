#functions from peakrefine

crossCorrByRle = function(bam_file,
                          query_gr,
                          max_dupes = 1,
                          fragment_sizes = 50:300,
                          read_length = NULL,
                          include_plots = TRUE,
                          ...){
  rn = cor = NULL # reserve for data.table
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
  names(query_gr) = query_gr$name
  # query_gr = resize(query_gr, 500, fix = "center")

  query_gr = harmonize_seqlengths(query_gr, bam_file)

  Param <- Rsamtools::ScanBamParam(which=query_gr,
                        what=c("flag","mapq"),
                        ...)
  temp <- GenomicAlignments::readGAlignments(bam_file,param=Param)
  dt = as.data.table(temp)
  # browser()
  if(is.null(read_length)){
    read_length = getReadLength(bam_file, query_gr)
  }
  if(is.na(read_length)){
    read_length = numeric()
  }
  fragment_sizes = sort(union(read_length, fragment_sizes))

  PosCoverage <- coverage(GenomicRanges::shift(GRanges(temp[strand(temp)=="+"])), -read_length)
  PosCoverage = PosCoverage[query_gr]
  names(PosCoverage) = query_gr$name

  NegCoverage <- coverage(GRanges(temp[strand(temp)=="-"]))
  NegCoverage = NegCoverage[query_gr]
  names(NegCoverage) = query_gr$name
  ShiftMatCor = pbapply::pbsapply(seq_along(query_gr), simplify = FALSE, function(i){
    ShiftsCorTemp <- S4Vectors::shiftApply(fragment_sizes,
                                           PosCoverage[[i]],
                                           NegCoverage[[i]],
                                           cor, simplify = FALSE,
                                           verbose = FALSE)
  })
  #necessary due to singleton query_gr or shift not resulting in matrix
  ShiftMatCor = matrix(unlist(ShiftMatCor),
                       byrow = FALSE,
                       nrow = length(fragment_sizes),
                       ncol = length(query_gr))
  ShiftMatCor[is.nan(ShiftMatCor)] = 0

  colnames(ShiftMatCor) = query_gr$name
  rownames(ShiftMatCor) = fragment_sizes
  shift_dt = as.data.table(ShiftMatCor, keep.rownames = TRUE)
  shift_dt[, shift := as.numeric(rn)]
  shift_dt$rn = NULL
  shift_dt = melt(shift_dt, id.vars = "shift",
                  variable.name = "id", value.name = "correlation")
  return(shift_dt)
}

getReadLength = function(bam_file,
                         query_gr){
  Param <- Rsamtools::ScanBamParam(which=sample(query_gr, min(10, length(query_gr))),
                                   what=c("flag","mapq"))
  temp <- GenomicAlignments::readGAlignments(bam_file,param=Param)
  readlength=as.numeric(names(sort(table(width(temp)), decreasing = TRUE))[1])
  readlength
}

#' Forces seqlengths in GRanges to be that in header of bam_file
#'
#' This is important to avoid errors while using the gr to query the bam_file.
#'
#' @param gr GRanges
#' @param bam_file bam_file
#'
#' @return GRanges with seqlengths matching bam_file
#' @import Rsamtools GenomeInfoDb
#' @examples
#' peak_files = dir(system.file("extdata", package = "seqqc"), pattern = "Peak$", full.names = TRUE)
#' peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files)
#' query_gr = resize(seqsetvis::ssvOverlapIntervalSets(peak_grs), 6e2, fix = "center")
#'
#' bam_files = dir(system.file("extdata", package = "seqqc"), pattern = "^M.+bam$", full.names = TRUE)
#'
#' query_gr.harmonized = seqqc:::harmonize_seqlengths(query_gr, bam_files[1])
#' seqlengths(query_gr.harmonized)
#'
harmonize_seqlengths = function(gr, bam_file){
  chr_lengths = Rsamtools::scanBamHeader(bam_file)[[1]]$targets
  GenomeInfoDb::seqlengths(gr) = chr_lengths[names(GenomeInfoDb::seqlengths(gr))]
  too_long = end(gr) > GenomeInfoDb::seqlengths(gr)[as.character(seqnames(gr))]
  if(any(too_long)){
    message(sum(too_long), " region shifted for extending beyond seqlengths")
    fix_gr = gr[too_long]
    shift_by = -(end(fix_gr) - GenomeInfoDb::seqlengths(fix_gr)[as.character(seqnames(fix_gr))])
    gr[too_long] = GenomicRanges::shift(fix_gr, shift_by)
  }
  too_short = start(gr) < 1
  if(any(too_short)){
    message(sum(too_short), " region shifted for starting before seqlengths")
    fix_gr = gr[too_short]
    shift_by = 1 - start(fix_gr)
    gr[too_short] = GenomicRanges::shift(fix_gr, shift_by)
  }
  gr
}

bfcif = function(bfc, rname, FUN, force_overwrite = FALSE){
  # is rname in cache?
  if(nrow(BiocFileCache::bfcquery(bfc, query = rname, field = "rname")) == 0){
    cache_path = BiocFileCache::bfcnew(bfc, rname = rname)

  }else{
    cache_path = BiocFileCache::bfcrpath(bfc, rname)
  }
  # does cached file exist?
  if(file.exists(cache_path) && !force_overwrite){
    load(BiocFileCache::bfcrpath(bfc, rname))
  }else{
    res = FUN()
    save(res, file = cache_path)
  }
  # return either new results or cached results
  res
}

gather_metrics = function(peak_strand_corr, read_length = NULL){
  correlation = id = NULL
  max_dt = peak_strand_corr[, list(shift = shift[which.max(correlation)], correlation = max(correlation)), by = list(id)]
  # if(is.null(read_length)){
  #     fl = round(.my_mode(max_dt[shift != min(shift, na.rm = TRUE) & shift != max(shift, na.rm = TRUE)]$shift, na.rm = TRUE))
  # }else{
  #     fl = round(.my_mode(max_dt[shift != read_length][shift != min(shift, na.rm = TRUE) & shift != max(shift, na.rm = TRUE)]$shift, na.rm = TRUE))
  # }
  flex_frag_corrs = max_dt[, list(shift, id, correlation)]


  average_corr = peak_strand_corr[, list(correlation = mean(correlation)), list(shift)]

  fl = average_corr[, shift[which.max(correlation)[1]]]


  stable_frag_corrs = peak_strand_corr[shift == fl]


  if(!is.null(read_length)){
    read_corrs = peak_strand_corr[shift == read_length]
    out = list(read_length = read_length,
               fragment_length = fl,
               read_correlation = read_corrs,
               flex_fragment_correlation = flex_frag_corrs,
               stable_fragment_correlation = stable_frag_corrs,
               full_correlation_results = peak_strand_corr,
               average_correlation = average_corr)
  }else{
    out = list(fragment_length = fl,
               flex_fragment_correlation = flex_frag_corrs,
               stable_fragment_correlation = stable_frag_corrs,
               full_correlation_results = peak_strand_corr,
               average_correlation = average_corr)
  }
  out
}
