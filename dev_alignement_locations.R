#curious about read alignment per chromosome
# 1) distribution of total reads per chromosome
# 2) chrM fraction specifically
# 3) distribution per hetero/euchromatin
# 4) tracks per chromosome

library(ssvRecipes)
library(magrittr)
library(GenomicRanges)/
?ssvRecipes::ssvR_plot_ideogogram_data
library(data.table)

gen = "hg38"
bfc = BiocFileCache::BiocFileCache()
grIdeo = ssvRecipes::bfcif(bfc, rname = paste("ideo", gen), function(){
  biovizBase::getIdeogram(gen)
})   

f = system.file("extdata/pcaOut.PC1.chr2x.txt", package = "ssvRecipes")
dat = fread(f)
dat = dat[, .(x = (start+end)/2, y = PC1, seqnames = chr)]
ssvR_plot_ideogogram_data(dat, grIdeo = grIdeo, facet_cols = 2, chr_to_show = c("chr20", "chr21", "chr22"))

bam_files.1 = dir("/slipstream/home/conggao/ChIP_seq/MCF10_H4Kac/", full.names = TRUE) %>%
  dir(full.names = TRUE, pattern = "pool.+bam$")

bam_files.2 = dir("/slipstream/galaxy/uploads/working/qc_framework/output_AF_MCF10_CTCF", full.names = TRUE) %>%
  dir(full.names = TRUE, pattern = "input_R.+.bam$")

bam_files.3 = dir("data", "MCF10DCIS.+input.+bam", full.names = TRUE)

bam_files = c(bam_files.1, bam_files.2, bam_files.3)

bam_files = bam_files[!grepl("DCIS_input", bam_files)]

paste(bam_files, collapse = " ")

dt = data.table(file = bam_files)
dt = dt[file.exists(file),]

dt[, c("cell", "mark", "rep") := tstrsplit(basename(file), "[_\\.]", keep = 1:3)]
dt[, mapped_reads_total := ssvQC::get_mapped_reads(file), .(file)]
dt 

f = dt$file[1]

#based on ssvQC::get_mapped_reads
get_perchrm_reads = function(f){
  stats = Rsamtools::idxstatsBam(f)
  pc_dt = as.data.table(stats)
  pc_dt$file = f
  pc_dt
}

#1)
pc_dt = lapply(dt$file, get_perchrm_reads) %>% rbindlist()
pc_dt = merge(pc_dt, dt, by = "file")                                                    
ggplot(pc_dt, aes(x = seqnames, y = mapped)) +
  geom_bar(stat = "identity") +
  facet_wrap(~cell+mark+rep, scales = "free_y", ncol = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8))

ggplot(pc_dt, aes(x = seqlength, y = mapped, label = seqnames)) +
  geom_text(stat = "identity", size = 3) +
  facet_wrap(~cell+mark+rep, scales = "free_y", ncol = 4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8))


ggplot(pc_dt[mark == "input"], aes(x = seqlength, y = mapped, label = seqnames)) +
  geom_text(stat = "identity", size = 3) +
  facet_wrap(~cell+mark+rep, scales = "free_y", ncol = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8))

#2)
frac_dt = pc_dt[, .(fraction = mapped[seqnames == "chrM"] / sum(mapped)), .(file, cell, mark, rep)]
ggplot(frac_dt, aes(x = paste(cell, mark, rep), y = fraction, fill = mark == "input")) +
  geom_bar(stat = "identity")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8))

ggplot(frac_dt[mark != "input"], aes(x = paste(cell, mark, rep), y = fraction, fill = mark == "input")) +
  geom_bar(stat = "identity")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8))


#3)
ssvRecipes::ssvR_plot_ideogogram
gen = "hg38"
library(BiocFileCache)
bfc = BiocFileCache()
grIdeo = bfcif(bfc, rname = paste("ideo", gen), function(){
  biovizBase::getIdeogram(gen)
})

library(seqsetvis)
sl_dt = unique(pc_dt[, .(seqnames, seqlength)])
qgr = GRanges(seqnames = sl_dt$seqnames, IRanges(1, sl_dt$seqlength))
qgr_l = split(qgr, seqnames(qgr))
qgr_l = qgr_l[names(qgr_l) != "chrM"]

dt[cell == "MCF10CA1", cell := "MCF10CA1a"]
dt[, name := paste(cell, mark, rep, sep = "_")]

bsize = 5e5

res_l = bfcif(bfc, digest::digest(list(dt, qgr_l, bsize)), function(){
  lapply(qgr_l, function(qgr){
    seqname = as.character(seqnames(qgr))
    message(seqname)
    nbins = ceiling(width(qgr)/bsize)
    ssvRecipes::bfcif(bfc, paste(seqname, digest::digest(list(dt, bsize))), function(){
      ssvFetchBam(dt, qgr, fragLens = NA, win_size = nbins, win_method = "summary", return_data.table = TRUE, n_cores = 36, n_region_splits = 100)
    })
  })
})

res_dt = rbindlist(res_l)
dim(res_dt)

res_dt[, x := (start + end)/2]
res_dt[, y_norm := y / mapped_reads_total * 1e6]

ggplot(res_dt[seqnames %in% paste0("chr", 1:12)], aes(x = x, y = y_norm, group = name, color = cell)) + 
  geom_path() +
  facet_grid(seqnames~mark) +
  coord_cartesian(ylim = c(0, .2))

cells = res_dt$cell %>% unique
cols = RColorBrewer::brewer.pal(length(cells), "Dark2")
names(cols) = cells

my_plot = function(sel_chr, sel_cell, sel_mark = c("H4K5ac", "H4K8ac", "input"), ymax = .2){
  ggplot(res_dt[seqnames %in% sel_chr & cell %in% sel_cell & mark %in% sel_mark], aes(x = x, y = y_norm, group = name, color = cell)) + 
    geom_path() +
    facet_grid(seqnames~mark) +
    coord_cartesian(ylim = c(0, ymax)) +
    scale_color_manual(values = cols) +
    scale_x_continuous(labels = function(x)x/1e6) +
    labs(x = "Mbp", y = "binned RPM") +
    theme(legend.position = "bottom", panel.background = element_blank(), panel.grid = element_blank())
}

p_mcf10a = ggplot(res_dt[seqnames %in% paste0("chr", 1:9) & cell == "MCF10A"], aes(x = x, y = y_norm, group = name, color = cell)) + 
  geom_path() +
  facet_grid(seqnames~mark+name) +
  coord_cartesian(ylim = c(0, .2)) +
  scale_color_manual(values = cols) +
  scale_x_continuous(labels = function(x)x/1e6) +
  labs(x = "Mbp", y = "binned RPM") +
  theme(legend.position = "bottom", panel.background = element_blank(), panel.grid = element_blank())

p_mcf10ca1a = ggplot(res_dt[seqnames %in% paste0("chr", 1:9) & cell == "MCF10CA1a"], aes(x = x, y = y_norm, group = name, color = cell)) + 
  geom_path() +
  facet_grid(seqnames~mark+name) +
  coord_cartesian(ylim = c(0, .2)) +
  scale_color_manual(values = cols) +
  scale_x_continuous(labels = function(x)x/1e6) +
  labs(x = "Mbp", y = "binned RPM") +
  theme(legend.position = "bottom", panel.background = element_blank(), panel.grid = element_blank())

my_plot("chrX", c("MCF10A", "MCF10CA1a"))
cols
p = my_plot(paste0("chr", c(1:22, "X", "Y")), c("MCF10A", "MCF10AT1",  "MCF10CA1a", "MCFDCIS"), sel_mark = "H4K8ac")
p

my_plot(sel_chr = "chr9", c("MCF10A", "MCF10AT1",  "MCF10CA1a", "MCFDCIS"), sel_mark = "H4K8ac")

sel_chr = paste0("chr", c(1:22, "X", "Y"))
sel_chr = c("chr9", "chr10", "chr11", "chr8", "chr7")
my_plot(sel_chr, c("MCF10A", "MCF10AT1",  "MCF10CA1a", "MCFDCIS"), sel_mark = "H4K8ac", ymax = .3)
my_plot(sel_chr, c("MCF10A", "MCF10AT1",  "MCF10CA1a", "MCFDCIS"), sel_mark = "H4K5ac", ymax = .3)
my_plot(sel_chr, c("MCF10A", "MCF10AT1",  "MCF10CA1a", "MCFDCIS"), sel_mark = "input", ymax = .3)



ggsave("tmp.pdf", p, width = 9, height = 20)
my_plot(c("chrX", "chr1"), c("MCF10A", "MCF10CA1a"), sel_mark = "H4K5ac")
my_plot(c("chrX", "chr1"), c("MCF10A", "MCF10CA1a"), sel_mark = "input")

cowplot::plot_grid(p_mcf10a, p_mcf10ca1a)

my_plot_ideogram_data = function(data_dt,
                                 gen = "hg38",
                                 chr_to_show = paste0("chr", c(1:22, "X", "Y")),
                                 gr_highlights = NULL,
                                 bfc = BiocFileCache::BiocFileCache(),
                                 grIdeo = NULL,
                                 ideo_ymin = -1,
                                 ideo_ymax = 0,
                                 highlight_fill = "green",
                                 highlight_color = NA,
                                 highlight_alpha = .5,
                                 facet_cols = 1,
                                 facet_by_row = FALSE,
                                 data_ymin = 0, data_ymax = 1,
                                 print_plot = TRUE,
                                 color_var = "cell", 
                                 color_values = NULL){
  ideo_res = ssvR_plot_ideogogram(gen = gen,
                                  chr_to_show = chr_to_show,
                                  gr_highlights = gr_highlights,
                                  bfc = bfc,
                                  grIdeo = grIdeo,
                                  ideo_ymin = ideo_ymin,
                                  ideo_ymax = ideo_ymax,
                                  highlight_fill = highlight_fill,
                                  highlight_color = highlight_color,
                                  highlight_alpha = highlight_alpha,
                                  facet_cols = facet_cols,
                                  facet_by_row = facet_by_row, print_plot = FALSE)
  
  if(is.null(color_values)){
    if(!is.null(data_dt[[color_var]])){
      cells = data_dt[[color_var]] %>% unique
      color_values = seqsetvis::safeBrew(length(cells), "Dark2")
      names(color_values) = cells  
    }
  }
  ideo_res$plot
  dt = data.table::copy(data_dt)
  dt = dt[seqnames %in% chr_to_show]
  dt[, y := y - min(y)]
  dt[, y := y / max(y)]
  dt[, y := y * (data_ymax - data_ymin) + data_ymin]
  
  pdt = dt[order(x)][order(seqnames)]
  pdt$seqnames = factor(pdt$seqnames, chr_to_show)
  
  if(!is.null(data_dt[[color_var]]) & color_var != "none"){
    ideo_res$plot = ideo_res$plot + geom_path(data = pdt, aes_string(x = "x", y = "y", color = color_var)) +
      scale_color_manual(values = color_values)
  }else{
    ideo_res$plot = ideo_res$plot + geom_path(data = pdt, aes(x = x, y = y))  
  }
  
  
  if(print_plot) plot(ideo_res$plot)
  ideo_res$data$line_data = pdt
  invisible(ideo_res)
}

dat
data_dt = res_dt[name == "MCF10A_H4K5ac_pooled"][, .(x = x, y = y_norm, seqnames, cell)]

cowplot::plot_grid(ncol = 1,
                   ggplot(res_dt[name == "MCF10A_H4K5ac_pooled"][, .(x = x, y = y_norm, seqnames)][seqnames == "chr9"], aes(x = x, y = y)) +
                     geom_path(),
                   my_plot_ideogram_data(data_dt = data_dt, 
                                         gen = "hg38", 
                                         grIdeo = grIdeo, 
                                         chr_to_show = c("chr9"))$plot
)

data_dt = res_dt[mark == "input"][, .(x = x, y = y_norm, seqnames, cell)]
data_dt = res_dt[mark == "H4K5ac"][, .(x = x, y = y_norm, seqnames, cell)]
data_dt = res_dt[mark == "H4K8ac"][, .(x = x, y = y_norm, seqnames, cell)]

cols

my_plot_ideogram_data(data_dt = data_dt, gen = "hg38", grIdeo = grIdeo, chr_to_show = c("chr9"),)$plot

seq_todo = paste0("chr", c(1:22, "X", "Y"))
names(seq_todo) = seq_todo
mark_todo = c("input", "H4K5ac", "H4K8ac")
names(mark_todo) = mark_todo

if(FALSE){
  plots = lapply(mark_todo, function(m){
    lapply(seq_todo, function(seq){
      data_dt = res_dt[mark == m][, .(x = x, y = y_norm, seqnames, cell)]
      p = my_plot_ideogram_data(data_dt = data_dt, 
                                gen = "hg38", 
                                grIdeo = grIdeo, 
                                chr_to_show = seq, 
                                color_var = "cell", 
                                color_values = cols)$plot
      p + labs(title = paste(m, "on", seq))
    })
  })
  
  pg = cowplot::plot_grid(
    plotlist = lapply(plots$H4K5ac, function(p){
      p +
        labs(title = "") + 
        guides(color = "none")}), 
    ncol = 4)
  ggsave("tmp.H4K5ac.pdf", pg, width = 10, height = 16)  
  
  pg = cowplot::plot_grid(
    plotlist = lapply(plots$H4K8ac, function(p){
      p +
        labs(title = "") + 
        guides(color = "none")}), 
    ncol = 4)
  ggsave("tmp.H4K8ac.pdf", pg, width = 10, height = 16)  
  
  pg = cowplot::plot_grid(
    plotlist = lapply(plots$input, function(p){
      p +
        labs(title = "") + 
        guides(color = "none")}), 
    ncol = 4)
  ggsave("tmp.input.pdf", pg, width = 10, height = 16)  
}
dev.size()

pdf("tmp.all.pdf", width = 10, height = 16)
ggplot(res_dt[mark == "H4K8ac"], aes(x = x, y = y_norm, color = cell)) + 
  geom_path() + 
  facet_grid(seqnames~cell, scales = "free_y") +
  scale_color_manual(values = cols) +
  labs(title = "H4K8ac")

ggplot(res_dt[mark == "H4K5ac"], aes(x = x, y = y_norm, color = cell)) + 
  geom_path() + 
  facet_grid(seqnames~cell, scales = "free_y") +
  scale_color_manual(values = cols) +
  labs(title = "H4K5ac")

ggplot(res_dt[mark == "input"], aes(x = x, y = y_norm, color = cell)) + 
  geom_path() + 
  facet_grid(seqnames~cell, scales = "free_y") +
  scale_color_manual(values = cols) +
  labs(title = "input")
dev.off()



m_dt = ssvFetchBam(dt, 
                   subset(qgr, seqnames == "chrM"), 
                   fragLens = NA, 
                   win_size = 1e3, 
                   win_method = "summary", 
                   return_data.table = TRUE, 
                   n_cores = 36, 
                   n_region_splits = 100)

m_dt[, y_norm := y / mapped_reads_total * 1e6]

ggplot(m_dt, aes(x = x, y = y_norm, color = cell)) +
  geom_path() +
  facet_grid(mark~cell, scales = "free_y") +
  scale_color_manual(values = cols) +
  labs(title = "chrM")

grIdeo$gieStain %>% table

gr_by_gieStain = split(subset(grIdeo, seqnames %in% paste0("chr", c(1:22, "X", "Y"))), grIdeo$gieStain)
sapply(gr_by_gieStain, function(x)sum(width(x))/1e6)

nam = names(gr_by_gieStain)[1]

message("fetch tiled gieStain regions")
nam = names(gr_by_gieStain)[2]
# debug(ssvFetchSignal)
# undebug(ssvFetchSignal)
tile_file = "bam_tile_dt.1bin.Rds"
if(!file.exists(tile_file)){
  bam_tile_dtl = lapply(names(gr_by_gieStain), function(nam){
    message(nam)
    gr = gr_by_gieStain[[nam]]
    gr_tiles = tile(gr, width = 5e4)
    gr_tiles = unlist(gr_tiles)
    hist(width(gr_tiles))
    gr_tiles = gr_tiles[width(gr_tiles) > 4.5e4]
    gr_tiles = prepare_fetch_GRanges_names(gr_tiles)
    
    bam_tile_dt = ssvFetchBam(dt,
                              gr_tiles,
                              fragLens = NA, 
                              # win_size = 1e3,
                              win_size = 1, 
                              win_method = "summary", 
                              return_data.table = TRUE, 
                              n_cores = 36, 
                              n_region_splits = 100)
    bam_tile_dt$gieStain = nam
    bam_tile_dt 
  })
  
  bam_tile_dt = rbindlist(bam_tile_dtl)
  saveRDS(bam_tile_dt, file = tile_file)
}else{
  bam_tile_dt = readRDS(tile_file)
}


format(object.size(bam_tile_dt), units = "GB")
format(object.size(bam_tile_dt), units = "MB")

ggplot(bam_tile_dt, aes(x = gieStain, y = y+1)) + geom_boxplot() +
  facet_wrap(~name) +
  scale_y_log10()
