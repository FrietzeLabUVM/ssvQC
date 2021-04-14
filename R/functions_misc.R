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

sampleCap = function(x, n = 500){
  n = min(n, length(unique(x)))
  out = sample(unique(x), n)
  if(is.factor(out)) out = as.character(out)
  out
}
