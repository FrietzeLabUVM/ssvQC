% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_bfcif.R
\name{bfcif}
\alias{bfcif}
\title{bfcif}
\usage{
bfcif(
  bfc,
  rname,
  FUN,
  version = getOption("SQC_CACHE_VERSION", "v4"),
  force_overwrite = getOption("SQC_FORCE_CACHE_OVERWRITE", FALSE),
  return_path_only = FALSE,
  verbose = getOption("SQC_CACHE_VERBOSE", FALSE)
)
}
\arguments{
\item{bfc}{A BiocFileCache object, typically from BiocFileCache::BiocFileCache().}

\item{rname}{The unique identifier for the results in the cache. The recommendation is to use either a unique and meaningful description or digest::digest() on a list containing FUN and it's parameters.}

\item{FUN}{A function that takes zero arguments.  This function can be a wrapper around other functions; ie. FUN = function(){mean(x)}.}

\item{version}{A version indicator string to further distinguish cache entries.  Typically, you want to iterate this value when code outside of FUN or input parameters have changed and you want to force FUN to reevaluate. Default is SQC_CACHE_VERSION option or v3 if option is not set.}

\item{force_overwrite}{If TRUE, FUN will be rerun regardless of cache state. If it exists, current cache contents will be overwritten. Default is SQC_FORCE_CACHE_OVERWRITE option or FALSE if option is not set.}

\item{return_path_only}{If TRUE, FUN will not be run and instead the cache path is returned. Default is FALSE.}

\item{verbose}{If TRUE, status is reported via messages. Default is value of SQC_CACHE_VERBOSE option or FALSE if option not set.}
}
\value{
Result of FUN, from cache if available.
}
\description{
Conditionally runs a function if it's results don't exist in the cache.  Cache is controlled by BiocFileCache.
}
\examples{
bfc = BiocFileCache::BiocFileCache()
bfcif(bfc, "test1", function(x)mean(seq(10)), verbose = TRUE)
bfcif(bfc, "test1", function(x)mean(seq(10)), verbose = TRUE, 
  return_path_only = TRUE)
bfcif(bfc, "test1", function(x)mean(seq(10)), verbose = TRUE, 
  force_overwrite = TRUE)
}
