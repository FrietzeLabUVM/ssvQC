# library(ssvQC)
# object = readRDS("dev_object.ssvQC.plotSignal.Rds")
# signal_data = object@signal_data$CTCF_features$CTCF_signal
# ###




# stsPlotBinnedHeatmap = function(){
#   
# }
###


#' @import concaveman
stsPlotClusterProfiles = function (profile_dt,
                                   cluster_dt = profile_dt,
                                   cluster_ = "cluster_id",
                                   w = 0.05,
                                   h = 0.05,
                                   x_var = "x",
                                   y_var = "y",
                                   id_var = "id",
                                   wide_var = "name",
                                   tall_var = "tall_none") {
  stopifnot(cluster_ %in% colnames(cluster_dt))
  cent_dt = cluster_dt[, .(tx = median(tx), ty = median(ty)),
                       by = cluster_]
  if(is.null(profile_dt[[tall_var]])){
    profile_dt[[tall_var]] = "none"
  }
  # if(is.null(profile_dt[[tall_var]])){
  #   profile_dt[, `:=`(tid, get(id_var))]
  # }else{
  profile_dt[, `:=`(tid, paste(get(tall_var), get(id_var)))]
  # }
  if(is.null(cluster_dt[["tid"]])){
    if(is.null(cluster_dt[[tall_var]])){
      cluster_dt[, `:=`(tid, get(id_var))]
    }else{
      cluster_dt[, `:=`(tid, paste(get(tall_var), get(id_var)))]
    }
  }
  
  cprof_dt = merge(profile_dt[, setdiff(colnames(profile_dt), cluster_),  with = FALSE],
                   unique(cluster_dt[, c("tid", cluster_), with = FALSE]),
                   by = "tid")
  cprof_dt = cprof_dt[, .(y = mean(y)), c(wide_var, cluster_, x_var)]
  cprof_dt = merge(cprof_dt, cent_dt, by = cluster_)
  cent_dt[, `:=`(xmin, tx - w/2)]
  cent_dt[, `:=`(xmax, tx + w/2)]
  cent_dt[, `:=`(ymin, ty - h/2)]
  cent_dt[, `:=`(ymax, ty + h/2)]
  glyph_dt = as.data.table(GGally::glyphs(cprof_dt,
                                          x_major = "tx",
                                          x_minor = x_var,
                                          y_major = "ty",
                                          y_minor = y_var,
                                          width = w,
                                          height = h))
  my_chull = function(x, y) {
    ch = concaveman::concaveman(cbind(x, y), concavity = 1.4,
                                length_threshold = 0.005)
    list(tx = ch[, 1], ty = ch[, 2])
  }
  ch_dt = cluster_dt[, my_chull(tx, ty), by = c(cluster_)]
  cols = seqsetvis::safeBrew(length(unique(ch_dt[[cluster_]])))
  if (is.numeric(ch_dt[[cluster_]])) {
    ch_dt[[cluster_]] = factor(as.character(ch_dt[[cluster_]]),
                               levels = as.character(sort(unique(ch_dt[[cluster_]]))))
  }
  
  glyph_dt[, group := paste(gid, get(wide_var))]
  
  p_clust_big = ggplot() +
    geom_polygon(data = ch_dt,
                 mapping = aes_string(x = "tx",
                                      y = "ty",
                                      fill = cluster_)) +
    scale_fill_manual(values = cols) +
    geom_rect(data = cent_dt,
              aes(xmin = xmin,
                  xmax = xmax,
                  ymin = ymin,
                  ymax = ymax),
              fill = "white",
              color = "black") +
    geom_path(data = glyph_dt,
              aes_string(x = "gx",
                         y = "gy",
                         group = "group",
                         color = wide_var)) +
    labs(x = "tx", y = "ty") +
    coord_fixed()
  p_clust_big
}
stsPlotClusterProfiles = function (profile_dt,
                                   cluster_dt = profile_dt,
                                   cluster_ = "cluster_id",
                                   w = 0.05,
                                   h = 0.05,
                                   x_var = "x",
                                   y_var = "y",
                                   id_var = "id",
                                   wide_var = "name",
                                   tall_var = "tall_none") {
  stopifnot(cluster_ %in% colnames(cluster_dt))
  cent_dt = cluster_dt[, .(tx = median(tx), ty = median(ty)),
                       by = cluster_]
  if(is.null(profile_dt[[tall_var]])){
    profile_dt[[tall_var]] = "none"
  }
  # if(is.null(profile_dt[[tall_var]])){
  #   profile_dt[, `:=`(tid, get(id_var))]
  # }else{
  profile_dt[, `:=`(tid, paste(get(tall_var), get(id_var)))]
  # }
  if(is.null(cluster_dt[["tid"]])){
    if(is.null(cluster_dt[[tall_var]])){
      cluster_dt[, `:=`(tid, get(id_var))]
    }else{
      cluster_dt[, `:=`(tid, paste(get(tall_var), get(id_var)))]
    }
  }
  
  cprof_dt = merge(profile_dt[, setdiff(colnames(profile_dt), cluster_),  with = FALSE],
                   unique(cluster_dt[, c("tid", cluster_), with = FALSE]),
                   by = "tid")
  cprof_dt = cprof_dt[, .(y = mean(y)), c(wide_var, cluster_, x_var)]
  cprof_dt = merge(cprof_dt, cent_dt, by = cluster_)
  cent_dt[, `:=`(xmin, tx - w/2)]
  cent_dt[, `:=`(xmax, tx + w/2)]
  cent_dt[, `:=`(ymin, ty - h/2)]
  cent_dt[, `:=`(ymax, ty + h/2)]
  glyph_dt = as.data.table(GGally::glyphs(cprof_dt,
                                          x_major = "tx",
                                          x_minor = x_var,
                                          y_major = "ty",
                                          y_minor = y_var,
                                          width = w,
                                          height = h))
  my_chull = function(x, y) {
    ch = concaveman::concaveman(cbind(x, y), concavity = 1.4,
                                length_threshold = 0.005)
    list(tx = ch[, 1], ty = ch[, 2])
  }
  ch_dt = cluster_dt[, my_chull(tx, ty), by = c(cluster_)]
  cols = seqsetvis::safeBrew(length(unique(ch_dt[[cluster_]])))
  if (is.numeric(ch_dt[[cluster_]])) {
    ch_dt[[cluster_]] = factor(as.character(ch_dt[[cluster_]]),
                               levels = as.character(sort(unique(ch_dt[[cluster_]]))))
  }
  
  glyph_dt[, group := paste(gid, get(wide_var))]
  
  p_clust_big = ggplot() +
    geom_polygon(data = ch_dt,
                 mapping = aes_string(x = "tx",
                                      y = "ty",
                                      fill = cluster_)) +
    scale_fill_manual(values = cols) +
    geom_rect(data = cent_dt,
              aes(xmin = xmin,
                  xmax = xmax,
                  ymin = ymin,
                  ymax = ymax),
              fill = "white",
              color = "black") +
    geom_path(data = glyph_dt,
              aes_string(x = "gx",
                         y = "gy",
                         group = "group",
                         color = wide_var)) +
    labs(x = "tx", y = "ty") +
    coord_fixed()
  p_clust_big
}
bin_values = function (x, n_bins, xrng = range(x)) {
  stopifnot(length(xrng) == 2)
  floor(rescale_capped(x, 0:1, xrng) * (n_bins - 1e-05)) + 
    1
}

rescale_capped = function (x, to = c(0, 1), from = range(x, na.rm = TRUE, finite = TRUE)) {
  y = scales::rescale(x, to, from)
  y[y > max(to)] = max(to)
  y[y < min(to)] = min(to)
  y
}


