
#' bfcif
#' 
#' Conditionally runs a function if it's results don't exist in the cache.  Cache is controlled by BiocFileCache.
#'
#' @param bfc A BiocFileCache object, typically from BiocFileCache::BiocFileCache().
#' @param rname The unique identifier for the results in the cache. The recommendation is to use either a unique and meaningful description or digest::digest() on a list containing FUN and it's parameters.
#' @param FUN A function that takes zero arguments.  This function can be a wrapper around other functions; ie. FUN = function(){mean(x)}.
#' @param version A version indicator string to further distinguish cache entries.  Typically, you want to iterate this value when code outside of FUN or input parameters have changed and you want to force FUN to reevaluate. Default is SQC_CACHE_VERSION option or v3 if option is not set.
#' @param force_overwrite If TRUE, FUN will be rerun regardless of cache state. If it exists, current cache contents will be overwritten. Default is SQC_FORCE_CACHE_OVERWRITE option or FALSE if option is not set.
#' @param return_path_only If TRUE, FUN will not be run and instead the cache path is returned. Default is FALSE.
#' @param verbose If TRUE, status is reported via messages. Default is value of SQC_CACHE_VERBOSE option or FALSE if option not set.
#'
#' @return Result of FUN, from cache if available.
#' @export
#'
#' @examples
#' bfc = BiocFileCache::BiocFileCache()
#' bfcif(bfc, "test1", function(x)mean(seq(10)), verbose = TRUE)
#' bfcif(bfc, "test1", function(x)mean(seq(10)), verbose = TRUE, 
#'   return_path_only = TRUE)
#' bfcif(bfc, "test1", function(x)mean(seq(10)), verbose = TRUE, 
#'   force_overwrite = TRUE)
bfcif = function(bfc, rname, FUN, 
                 version = getOption("SQC_CACHE_VERSION", "v4"),
                 force_overwrite = getOption("SQC_FORCE_CACHE_OVERWRITE", FALSE),
                 return_path_only = FALSE, 
                 verbose = getOption("SQC_CACHE_VERBOSE", FALSE)){
  # is rname in cache?
  vrname = paste0(rname, "_", version)
  if(nrow(BiocFileCache::bfcquery(bfc, query = vrname, field = "rname")) == 0){
    if(verbose) message("results not in cache. ")
    cache_path = BiocFileCache::bfcnew(bfc, rname = vrname)
  }else{
    if(verbose) message("previous cache results found. ")
    cache_path = BiocFileCache::bfcrpath(bfc, vrname)
  }
  if(return_path_only){
    if(verbose) message("returning cache path.")
    if(verbose) message(cache_path)
    return(cache_path)
  }
  # does cached file exist?
  if(file.exists(cache_path) && !force_overwrite){
    if(verbose) message("loading previous cache results...")
    if(verbose) message(cache_path)
    load(cache_path)
  }else{
    if(verbose) message("running function...")
    res = FUN()
    if(verbose) message("caching results...")
    if(verbose) message(cache_path)
    save(res, file = cache_path)
  }
  # return either new results or cached results
  if(is(res, "data.table")){
    data.table::setalloccol(res)
  }
  if(is.list(res)){
    if(any(sapply(res, is, class2 = "data.table"))){
      res = lapply(res, function(x){
        if(is(x, "data.table")){
          data.table::setalloccol(x)
        }
      })
    }
  }
  res
}

new_cache = function(cache_path = getOption("SQC_CACHE_PATH", "~/.cache")){
  BiocFileCache::BiocFileCache(cache_path)
}

#' get_args
#'
#' returns parameters of calling function as a named list.
#'
#' @param env Environment to collect variable from. By default this is the
#'   calling environment, i.e. the function containing a call to get_args. This
#'   may also be a function that returns an environment.
#' @param to_ignore character names of variables to ignore.
#' @param ... Additional variables that should be considered as part of
#'   environment.
#'
#' @return parameters of calling function as a named list.
#' @export
#'
#' @examples
#' #The most common usage is to simply collect all local variables in a function
#' test_fun = function(x = 1, y = 2){
#'   get_args()
#' }
#' test_fun()
#'
#' #Specified variables may be ignored
#' test_fun2 = function(x = 1, y = 2){
#'   get_args(to_ignore = "x")
#' }
#' test_fun2()
#'
#' #Additional variables can also be added from higher environments
#' global_z = 3
#' test_fun3 = function(x = 1, y = 2){
#'   get_args(env = parent.frame, to_ignore = character(), z = global_z)
#' }
#' test_fun3()
get_args = function(env = parent.frame(), to_ignore = character(), ...){
  if(is.function(env)) env = env()
  stopifnot(is.environment(env))
  args = c(as.list(env), list(...))
  args = args[!names(args) %in% to_ignore]
  args[order(names(args))]
}
#' digest_args
#'
#' returns digest::digest results of named list of arguments to calling function
#'
#' @inheritParams get_args
#'
#' @return A md5 digest of all variabes in calling environment. For a function
#'   is the same as calling arguments.
#' @export
#'
#' @examples
#' #The most common usage is to simply collect all local variables in a function
#' test_fun = function(x = 1, y = 2){
#'   digest_args()
#' }
#' test_fun()
#'
#' #Specified variables may be ignored
#' test_fun2 = function(x = 1, y = 2){
#'   digest_args(to_ignore = "x")
#' }
#' test_fun2()
#'
#' #Additional variables can also be added from higher environments
#' global_z = 3
#' test_fun3 = function(x = 1, y = 2){
#'   digest_args(env = parent.frame, to_ignore = character(), z = global_z)
#' }
#' test_fun3()
digest_args = function(env = parent.frame(), to_ignore = character(), ...){
  digest::digest(get_args(env, to_ignore, ...))
}