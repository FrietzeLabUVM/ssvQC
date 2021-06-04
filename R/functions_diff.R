# wds = c(
#   dir("/slipstream/galaxy/uploads/working/qc_framework/output", pattern = "MCF7.+K4ME3.+R[12]", full.names = TRUE),
#   dir("/slipstream/galaxy/uploads/working/qc_framework/output", pattern = "MCF10A.+K4ME3.+R[12]", full.names = TRUE)
# )
# wds.pool = c(
#   dir("/slipstream/galaxy/uploads/working/qc_framework/output", pattern = "MCF7.+K4ME3.+pool", full.names = TRUE),
#   dir("/slipstream/galaxy/uploads/working/qc_framework/output", pattern = "MCF10A.+K4ME3.+pool", full.names = TRUE)
# )
# parse_n = function(txt, n){
#   sapply(strsplit(basename(txt), "_"), function(x)x[n])
# }
# 
# parse_n(wds, 1)
# parse_n(wds, 2)
# parse_n(wds, 3)
# 
# pool_names = parse_n(wds, 1)
# rep_bams = dir(wds, pattern = "bam$", full.names = TRUE)
# peaks = dir(wds.pool, pattern = "peaks.narrowPeak$", full.names = TRUE)
# 
# setup_diff = function(pool_names, rep_bams, peaks){
#   conditional_split = function(txt){
#     if(length(txt) == 1)
#       if(grepl(",", txt)){
#         txt = strsplit(txt, ",")[[1]]
#       }
#     txt
#   }
#   pool_names = conditional_split(pool_names)
#   rep_bams = conditional_split(rep_bams)
#   peaks = conditional_split(peaks)
#   peaks = unique(peaks)
#   
#   stopifnot(length(pool_names) == length(rep_bams))
#   stopifnot(is.character(pool_names))
#   stopifnot(is.character(rep_bams))
#   stopifnot(is.character(peaks))
#   
#   load_funs = ssvQC::get_feature_file_load_function(ssvQC::guess_feature_file_format(peaks))
#   peaks_gr = list()
#   for(i in seq_along(peaks)){
#     peaks_gr[[basename(peaks[i])]] = load_funs[[i]](peaks[i])[[1]]
#   }
#   
#   peaks_gr = seqsetvis::ssvOverlapIntervalSets(seqsetvis::easyLoad_narrowPeak(peaks))
#   
#   samp_names = sub(".bam", "", basename(rep_bams))
#   names(rep_bams) = samp_names
#   names(pool_names) = samp_names
#   return(list(pool_names = pool_names, rep_bams = rep_bams, peaks_gr = peaks_gr))
# }
# 
# run_diffbind.rscript = function(pool_names, rep_bams, peaks, bfc = BiocFileCache::BiocFileCache()){
#   library(DiffBind) #following error without DiffBind attached:
#   #Error in as.environment(where) : 
#   #no item called "package:parallel" on the search list
# 
#   setup = setup_diff(pool_names, rep_bams, peaks)
#   pool_names = setup$pool_names
#   rep_bams = setup$rep_bams
#   peaks_gr = setup$peaks_gr
#   n_cores = getOption("mc.cores", 1)#DiffBind is doing weird stuff with mc.cores creating nested lists
#   dba = bfcif(bfc, digest::digest(list(setup, "run_diffbind.rscript full")), function(){
#     dba = DiffBind::dba.peakset(peaks = peaks_gr, bamReads = rep_bams[1], sampID = samp_names[1], treatment = pool_names[1])
#     for(sn in samp_names){
#       dba = DiffBind::dba.peakset(dba, peaks = peaks_gr, bamReads = rep_bams[sn], sampID = sn, treatment = pool_names[sn])
#     }
#     dba = DiffBind::dba.count(dba)
#     dba = DiffBind::dba.contrast(dba, categories = DiffBind::DBA_TREATMENT, minMembers = 2) 
#     dba = DiffBind::dba.analyze(dba)
#     dba
#   })
#   if(is.list(getOption("mc.cores", 1))){
#     options(mc.cores = n_cores)
#   }
#   rpt = DiffBind::dba.report(dba)
#   plot(dba)
#   subset(rpt, FDR < .01 & abs(Fold) > 2)
# }
# 
# library(csaw)W
# 
# run_csaw.rscript = function(pool_names, rep_bams, peaks, bfc = BiocFileCache::BiocFileCache()){
#   setup = setup_diff(pool_names, rep_bams, peaks)
#   pool_names = setup$pool_names
#   rep_bams = setup$rep_bams
#   peaks_gr = setup$peaks_gr
#   
# }