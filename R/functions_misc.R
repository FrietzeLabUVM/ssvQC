#' Title
#'
#' @param my_plots 
#'
#' @return
#' @export
#'
#' @examples
sync_width = function(my_plots){
  stopifnot(class(my_plots) == "list")
  is_ok = sapply(my_plots, function(x){
    "ggplot" %in% class(x) | "grob" %in% class(x)
  })
  stopifnot(all(is_ok))
  my_grobs = lapply(my_plots, function(x){
    if(grid::is.grob(x)){
      x
    }else{
      ggplotGrob(x)
    }
  })

  my_widths = lapply(my_grobs, function(gt){
    gt$widths
  })
  maxWidth = my_widths[[1]]
  if(length(my_widths) > 1){
    for(i in 2:length(my_widths)){
      maxWidth = grid::unit.pmax(maxWidth, my_widths[[i]])
    }
  }
  for(j in 1:length(my_grobs)){
    my_grobs[[j]]$widths = maxWidth
  }
  my_grobs
}

#' Title
#'
#' @param my_plots 
#'
#' @return
#' @export
#'
#' @examples
sync_height = function(my_plots){
  stopifnot(class(my_plots) == "list")
  is_ok = sapply(my_plots, function(x){
    "ggplot" %in% class(x) | "grob" %in% class(x)
  })
  stopifnot(all(is_ok))
  my_grobs = lapply(my_plots, function(x){
    if(grid::is.grob(x)){
      x
    }else{
      ggplotGrob(x)
    }
  })

  my_widths = lapply(my_grobs, function(gt){
    gt$heights
  })
  maxWidth = my_widths[[1]]
  if(length(my_widths) > 1){
    for(i in 2:length(my_widths)){
      maxWidth = grid::unit.pmax(maxWidth, my_widths[[i]])
    }
  }
  for(j in 1:length(my_grobs)){
    my_grobs[[j]]$heights = maxWidth
  }
  my_grobs
}

#' Title
#'
#' @param x 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
sampleCap = function(x, n = 500){
  n = min(n, length(unique(x)))
  out = sample(unique(x), n)
  if(is.factor(out)) out = as.character(out)
  out
}

#' Title
#'
#' @param f 
#'
#' @return
#' @export
#'
#' @examples
get_mapped_reads = function(f){
  stats = Rsamtools::idxstatsBam(f)
  sum(stats[,3])
}


#' guess_feature_file_format
#'
#' @param feature_files 
#'
#' @return
#' @export
#'
#' @examples
guess_feature_file_format = function(feature_files){
  .guess_feature_file_format = function(feature_file){
    file_format = "unknown"
    if(grepl("narrowPeak$", feature_file)){
      file_format = "narrowPeak"
    }else if(grepl("broadPeak$", feature_file)){
      file_format = "broadPeak"
    }else if(grepl("bed$", feature_file)){
      file_format = "bed"
    }else{
      warning("Could not guess file format for feature file: ", feature_file)
    }
    file_format
  }
  
  sapply(feature_files, .guess_feature_file_format)
}

#' Title
#'
#' @param signal_file 
#'
#' @return
#' @export
#'
#' @examples
guess_read_mode = function(signal_file){
  if(signal_file[1] == "null"){
    return("null")
  }
  if(grepl(".bam$", signal_file[1])){
    mode = "bam_SE"
  }else{
    mode = "bigwig"
  }
  message("read_mode has been guessed as ", mode)
  if(mode == "bam_SE"){
    message("Currently ssvQC cannot guess whether a bam file is SE or PE.  Please manually specify bam_PE if appropriate.")
  }
  mode
}

#' Title
#'
#' @param feature_files 
#'
#' @return
#' @export
#'
#' @examples
get_feature_file_load_function = function(feature_files){
  file_types = guess_feature_file_format(feature_files)
  .get_feature_file_load_function = function(file_type){
    switch (file_type,
            narrowPeak = seqsetvis::easyLoad_narrowPeak,
            broadPeak = seqsetvis::easyLoad_broadPeak,
            seacr = seqsetvis::easyLoad_seacr,
            bed = seqsetvis::easyLoad_bed,
            unknown = {
              warning("Treating unknown file type as bed but if you see errors, check file type.")
              seqsetvis::easyLoad_bed
            },
            stop("'", file_type, "' is not a supported file_type")
    )
  }
  
  sapply(file_types, .get_feature_file_load_function)
  
  
}