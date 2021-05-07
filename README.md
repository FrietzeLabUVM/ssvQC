# ssvQC
R package for QC of enrichment based NGS assays.  ChIP-seq, cut&amp;run, ATAC-seq, etc.

# Installation
```
devtools::install_github("FrietzeLabUVM/ssvQC")
```

# Quick Start

## Simple 
The simplest (though not necessarily recommended) way to start using ssvQC is by generating 2 character vectors of file paths: one of signal files (bam or bigwig), the second for feature files (bed, narrowPeak, etc.)

ssvQC comes with example files that consist of 100 random peaks from a CTCF ChIP-seq in breast cancer cell lines.

```
pkg_path = system.file("extdata", package = "ssvQC")
np_files = dir(pkg_path, pattern = "random100.narrowPeak", full.names = TRUE)
bam_files = dir(pkg_path, pattern = "100peaks.bam$", full.names = TRUE)
```

Next we'll create an ssvQC object using these files.

```
options(mc.cores = 8)
set.seed(0)
sqc = ssvQC(np_files, bam_files)
sqc = ssvQC.runAll(sqc)
```

ssvQC.runAll does all the data processing and generates all plots.  

ssvQC objects are intended to be accessed via the $ accessor.  

For instance, looking at the peak overlap is done like so.

```
sqc$plots$features$count
sqc$plots$features$venn
```

Another example, the strand-cross correlation (SCC) plot is accessed like so.

```
sqc$plots$SCC$curves
```

Here we can look at the read pileup in each bam at the peaks.

```
cowplot::plot_grid(
  sqc$plots$signal$heatmaps$All_features$All_signal,
  sqc$plots$signal$heatmaps.lines$All_features$All_signal
)
```

We can also fiddle with the colors

```
sqc$feature_config$color_mapping$MCF10A_CTCF_random100.narrowPeak = "red"
sqc$feature_config$color_mapping$MCF10AT1_CTCF_random100.narrowPeak = "blue"
sqc$feature_config$color_mapping$MCF10CA1_CTCF_random100.narrowPeak = "green"
#plots need to be regenerated after a change like this, unfortunately
sqc = ssvQC.runAll(sqc)
sqc$plots$features$count
sqc$plots$features$venn
```

## Recommended

The most powerful way to control ssvQC is by the use of configuration files.

```
bam_config = file.path(pkg_path, "ssvQC_bam_config.csv")
peak_config = file.path(pkg_path, "ssvQC_peak_config.csv")

sqc.config = ssvQC(peak_config, bam_config)
sqc.config = ssvQC.runAll(sqc.config)

sqc.config$plots$features$count
sqc.config$plots$features$venn
```