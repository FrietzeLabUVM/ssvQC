% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClusteredSignal.R
\name{ClusteredSignal}
\alias{ClusteredSignal}
\title{ClusteredSignal}
\usage{
ClusteredSignal(
  signal_profile_dt,
  query_gr,
  manual_assigned = list(),
  nclust = 6,
  signal_var = "y",
  signal_var.within = signal_var,
  facet_var = "name_split",
  extra_var = character(),
  bfc = new_cache()
)
}
\arguments{
\item{signal_profile_dt}{Tidy data.table containing profile information.  See output of seqsetvis::ssvFetchBam.}

\item{query_gr}{A GRanges containing regions to retrieve signal data at.}

\item{manual_assigned}{NYI but should allow manual cluster assignment.}

\item{nclust}{Number of k-means clusters to calculate. Default is 6.}

\item{signal_var}{Variable name for signal information to cluster upon in signal_profile_dt. Default is "y".}

\item{signal_var.within}{Variable name for ranking items within clusters.  The Default is the same as signal_var.}

\item{facet_var}{Variable that will eventually be used in heatmap facets.  Ensures it is preserved and not aggregated away.  Default is "name_split".}

\item{extra_var}{Any extra variables to preserve and avaoid aggregating away.}

\item{bfc}{BiocFileCache to use, uses default location otherwise.}
}
\value{
A ClusteredSignal object containing clustering information.
}
\description{
ClusteredSignal
}
\examples{
signal_profile_dt = seqsetvis::CTCF_in_10a_profiles_dt
setnames(signal_profile_dt, "sample", "name_split")
query_gr = seqsetvis::CTCF_in_10a_overlaps_gr
clust_sig = ClusteredSignal(signal_profile_dt, query_gr)
}
