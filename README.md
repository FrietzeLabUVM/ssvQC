# ssvQC
R package for QC of enrichment based NGS assays.  ChIP-seq, cut&amp;run, ATAC-seq, etc.

# Installation
```
#unless you're using the most recent dev version of BioConductor, you'll need the install seqsetvis from github.
devtools::install_github("jrboyd/seqsetvis")
devtools::install_github("FrietzeLabUVM/ssvQC")
```

For older versions of R (prior to R 4.0) it is tricky to get the requried version of seqsetvis installed alongside older BiocConductor packages.

The following works in R 3.6.3 with BioConductor 3.10, though you will see various deprecation messages once per session.

```
#this creates a designated personal library so that you don't muck up your existing R setup
dev_lib = paste0(.libPaths()[1], ".dev")
dir.create(dev_lib)
.libPaths(dev_lib)
.libPaths()

install.packages("devtools")

#dplyr_1.0.0  dbplyr_1.3.0 works
#Do not update existing packages if asked!
devtools::install_version("dplyr", "1.0.0")
devtools::install_version("dbplyr", "1.3.0")

devtools::install_github("jrboyd/seqsetvis")
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
