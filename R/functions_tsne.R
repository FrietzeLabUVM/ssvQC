
#' run_tsne
#'
#' @param profile_dt 
#' @param perplexity 
#' @param n_cores 
#' @param high_topright 
#' @param norm1 
#' @param Y_init 
#' @param verbose 
#' @param force_overwrite 
#' @param wide_var 
#' @param tall_var 
#'
#' @return
#' @import Rtsne
#'
#' @examples
run_tsne = function (profile_dt, 
                     perplexity = 100, 
                     n_cores = getOption("mc.cores", 1), 
                     high_topright = TRUE, 
                     norm1 = TRUE, 
                     Y_init = NULL, 
                     verbose = TRUE, 
                     force_overwrite = FALSE, 
                     x_var = "x", 
                     y_var = "y", 
                     id_var = "id",
                     wide_var = "name", 
                     tall_var = "tall_none") {
  valid_vars = c(ifelse(tall_var == "tall_none", character(), tall_var), wide_var, x_var, y_var, id_var)
  valid_vars = valid_vars[!is.na(valid_vars)]
  stopifnot(valid_vars %in% colnames(profile_dt))
  
  if(is.null(profile_dt[[tall_var]])){
    set(profile_dt, j = tall_var, value = "none")
  }
  tsne_mat = dt2mat(profile_dt, 
                    unique(profile_dt[[wide_var]]), 
                    x_var = x_var,
                    y_var = y_var,
                    id_var = id_var, 
                    wide_var = wide_var,
                    tall_var = tall_var)
  bad_col = apply(tsne_mat, 2, stats::var) == 0
  if (any(bad_col)) {
    warning("zero variance columns detected in tsne matrix input.", 
         "\n", round(sum(bad_col)/length(bad_col) * 100, 
                     2), "% of columns removed.")
    tsne_mat = tsne_mat[, bad_col, drop = FALSE]
  }
  if (is.data.table(Y_init)) {
    if(tall_var != "tall_none"){
      Y_init = as.matrix(Y_init[, .(tx, ty)], rownames.value = paste(Y_init[[id_var]], 
                                                                     Y_init[[tall_var]]))  
    }else{
      Y_init = as.matrix(Y_init[, .(tx, ty)], rownames.value = paste(Y_init[[id_var]], "none"))
    }
    
    Y_init = Y_init[rownames(tsne_mat), ]
    stopifnot(ncol(Y_init) == 2)
  }
  else {
    stopifnot(is.null(Y_init))
  }
  max_perplexity = floor(nrow(tsne_mat)/4)
  if(perplexity > max_perplexity){
    perplexity = max_perplexity
    warning("Reducing perplexity to ", perplexity, " to accommodate data of ", nrow(tsne_mat), " rows.")
  }
  res_tsne = Rtsne::Rtsne(tsne_mat, 
                          Y_init = Y_init, 
                          num_threads = n_cores, 
                          perplexity = perplexity, 
                          check_duplicates = FALSE)
  tsne_dt = as.data.table(res_tsne$Y)
  setnames(tsne_dt, c("tx", "ty"))
  tsne_dt$rn = rownames(tsne_mat)
  
  tsne_dt[, `:=`(c(id_var, tall_var), tstrsplit(rn, " ", keep = seq(2)))]
  
  
  if (norm1) {
    tsne_dt$tx = rescale_capped(tsne_dt$tx) - 0.5
    tsne_dt$ty = rescale_capped(tsne_dt$ty) - 0.5
  }
  if (high_topright) {
    rs = rowSums(tsne_mat)
    tsne_dt$rs = rs[tsne_dt$rn]
    x_cutoff = mean(range(tsne_dt$tx))
    x_flip = sum(tsne_dt[tx > x_cutoff]$rs) < sum(tsne_dt[tx < x_cutoff]$rs)
    if (x_flip) {
      tsne_dt[, `:=`(tx, max(tx) - tx + min(tx))]
    }
    y_cutoff = mean(range(tsne_dt$ty))
    y_flip = sum(tsne_dt[ty > y_cutoff]$rs) < sum(tsne_dt[ty < y_cutoff]$rs)
    if (y_flip) {
      tsne_dt[, `:=`(ty, max(ty) - ty + min(ty))]
    }
    tsne_dt$rs = NULL
  }
  tsne_dt$rn = NULL
  tsne_dt[]
}

rescale_capped = function (x, to = c(0, 1), from = range(x, na.rm = TRUE, finite = TRUE)){
  y = scales::rescale(x, to, from)
  y[y > max(to)] = max(to)
  y[y < min(to)] = min(to)
  y
}

dt2mat = function (prof_dt, 
                   wide_values, 
                   x_var = "x", 
                   y_var = "y", 
                   id_var = "id",
                   wide_var = "name", 
                   tall_var = "tall_none"){
  if(is.null(prof_dt[[tall_var]])){
    if(tall_var == "tall_none"){
      prof_dt[[tall_var]] = "none"
    }
  }
  stopifnot(c(id_var, wide_var, tall_var, x_var, y_var) %in% colnames(prof_dt))
  for (m in wide_values) {
    if (m == wide_values[1]) {
      dt = dcast(prof_dt[get(wide_var) == m], paste0(id_var, "+", tall_var, "~", x_var), value.var = y_var)
      wide_mat = as.matrix(dt[, -seq_len(2)])
      rn = paste(dt[[id_var]], dt[[tall_var]])
    }
    else {
      dt = dcast(prof_dt[get(wide_var) == m], paste(id_var, "+", tall_var, "~", x_var), value.var = y_var)
      stopifnot(all(paste(dt[[id_var]], dt[[tall_var]]) == rn))
      wide_mat = cbind(wide_mat, as.matrix(dt[, -seq_len(2)]))
    }
  }
  rownames(wide_mat) = rn
  wide_mat
}

stsPlotSummaryProfiles = function (profile_dt, 
                                   position_dt, 
                                   x_points, 
                                   y_points = x_points, 
                                   x_var = "x", 
                                   y_var = "y", 
                                   id_var = "id",
                                   wide_var = "name",
                                   tall_var = "tall_none",
                                   q_tall_values = NULL, 
                                   q_wide_values = NULL, 
                                   xrng = range(position_dt$tx), 
                                   yrng = range(position_dt$ty), 
                                   plot_type = c("glyph", "raster")[1], 
                                   rname = NULL, 
                                   odir = file.path(tempdir(), rname), 
                                   force_rewrite = FALSE, 
                                   n_cores = getOption("mc.cores", 1), 
                                   apply_norm = TRUE, 
                                   ylim = c(0, 1), 
                                   ma_size = 2, 
                                   n_splines = 10, 
                                   p = NULL, 
                                   facet_byCell = FALSE, 
                                   line_color_mapping = NULL, 
                                   vertical_facet_mapping = NULL, 
                                   N_floor = 0, 
                                   N_ceiling = NULL, 
                                   min_size = 0.3, 
                                   return_data = FALSE) {
  if(is.null(profile_dt[[tall_var]])){
    profile_dt[[tall_var]] = "none"
  }
  if (is.null(q_tall_values)) {
    if (is.factor(profile_dt[[tall_var]])) {
      q_tall_values = levels(profile_dt[[tall_var]])
    }
    else {
      q_tall_values = sort(unique(profile_dt[[tall_var]]))
    }
  }
  if (is.null(q_wide_values)) {
    if (is.factor(profile_dt[[wide_var]])) {
      q_wide_values = levels(profile_dt[[wide_var]])
    }
    else {
      q_wide_values = sort(unique(profile_dt[[wide_var]]))
    }
  }
  if (is.null(line_color_mapping)) {
    line_color_mapping = seqsetvis::safeBrew(length(unique(profile_dt[[wide_var]])))
    names(line_color_mapping) = unique(profile_dt[[wide_var]])
  }
  line_color_mapping = line_color_mapping[names(line_color_mapping) %in% 
                                            q_wide_values]
  if (is.factor(profile_dt[[tall_var]])) {
    stopifnot(q_tall_values %in% levels(profile_dt[[tall_var]]))
  }
  else {
    stopifnot(q_tall_values %in% unique(profile_dt[[tall_var]]))
  }
  if (is.factor(profile_dt[[wide_var]])) {
    stopifnot(q_wide_values %in% levels(profile_dt[[wide_var]]))
  }
  else {
    stopifnot(q_wide_values %in% unique(profile_dt[[wide_var]]))
  }
  if (!plot_type %in% c("glyph", "raster")) {
    stop("plot_type (\"", plot_type, "\") must be one of \"glyph\" or \"raster\".")
  }
  prof_dt = copy(profile_dt[get(tall_var) %in% q_tall_values & get(wide_var) %in% 
                              q_wide_values])
  pos_dt = copy(position_dt[get(tall_var) %in% q_tall_values])
  if (is.null(rname)) {
    rname = digest::digest(list(prof_dt, pos_dt, x_points, 
                                y_points, xrng, yrng, apply_norm, ylim, line_color_mapping, 
                                n_splines, ma_size, facet_byCell))
  }
  prof_dt[[tall_var]] = factor(prof_dt[[tall_var]], levels = q_tall_values)
  prof_dt[[wide_var]] = factor(prof_dt[[wide_var]], levels = q_wide_values)
  pos_dt[[tall_var]] = factor(pos_dt[[tall_var]], levels = q_tall_values)
  if (plot_type == "raster") {
    if (!facet_byCell) {
      summary_dt = prep_summary(profile_dt = prof_dt, 
                                position_dt = pos_dt, 
                                x_points = x_points, 
                                y_points = y_points, 
                                xrng = xrng, 
                                yrng = yrng, 
                                facet_by = NULL, 
                                x_var = x_var, 
                                y_var = y_var, 
                                id_var = id_var, 
                                wide_var = wide_var, 
                                tall_var = tall_var)
      img_res = prep_images(summary_dt = summary_dt, 
                            xrng = xrng, 
                            yrng = yrng, 
                            x_points = x_points, 
                            y_points = y_points, 
                            rname = rname, 
                            odir = odir, 
                            force_rewrite = force_rewrite, 
                            apply_norm = apply_norm, 
                            ylim = ylim, 
                            ma_size = ma_size, 
                            n_splines = n_splines, 
                            n_cores = n_cores, 
                            line_color_mapping = line_color_mapping, 
                            vertical_facet_mapping = vertical_facet_mapping, 
                            wide_var = wide_var, 
                            x_var = x_var, 
                            y_var = y_var)
      plot_summary_raster(image_dt = img_res$image_dt, 
                          xrng = xrng, 
                          yrng = yrng, 
                          x_points = x_points, 
                          y_points = y_points, 
                          p = p, 
                          line_color_mapping = img_res$line_color_mapping, 
                          N_floor = N_floor, 
                          N_ceiling = N_ceiling, 
                          min_size = min_size, 
                          return_data = return_data, 
                          wide_var = wide_var)
    }
    else {
      summary_dt = prep_summary(profile_dt = prof_dt, 
                                position_dt = pos_dt, 
                                x_points = x_points, 
                                y_points = y_points, 
                                xrng = xrng, 
                                yrng = yrng, 
                                facet_by = tall_var, 
                                x_var = x_var, 
                                y_var = y_var, 
                                id_var = id_var, 
                                wide_var = wide_var)
      img_res = prep_images(summary_dt = summary_dt, 
                            xrng = xrng, 
                            yrng = yrng, 
                            x_points = x_points, 
                            y_points = y_points, 
                            rname = rname, 
                            odir = odir, 
                            force_rewrite = force_rewrite, 
                            apply_norm = apply_norm, 
                            ylim = ylim, 
                            ma_size = ma_size, 
                            n_splines = n_splines, 
                            n_cores = n_cores, 
                            line_color_mapping = line_color_mapping, 
                            vertical_facet_mapping = vertical_facet_mapping, 
                            wide_var = wide_var, 
                            x_var = x_var, 
                            y_var = y_var)
      plot_summary_raster_byCell(image_dt = img_res$image_dt, 
                                 xrng = xrng, 
                                 yrng = yrng, 
                                 x_points = x_points, 
                                 y_points = y_points, 
                                 p = p, 
                                 line_color_mapping = img_res$line_color_mapping, 
                                 N_floor = N_floor, 
                                 N_ceiling = N_ceiling, 
                                 min_size = min_size, 
                                 return_data = return_data)
    }
  }
  else if (plot_type == "glyph") {
    if (!facet_byCell) {
      summary_dt = prep_summary(profile_dt = prof_dt, 
                                position_dt = pos_dt, 
                                x_points = x_points, 
                                y_points = y_points, 
                                xrng = xrng, 
                                yrng = yrng, 
                                facet_by = NULL, 
                                x_var = x_var, 
                                y_var = y_var, 
                                id_var = id_var, 
                                wide_var = wide_var, 
                                tall_var = tall_var)
      plot_summary_glyph(summary_dt = summary_dt, 
                         xrng = xrng, 
                         yrng = yrng, 
                         x_points = x_points, 
                         y_points = y_points, 
                         p = p,
                         ylim = ylim, 
                         N_floor = N_floor, 
                         N_ceiling = N_ceiling, 
                         min_size = min_size, 
                         color_mapping = line_color_mapping, 
                         return_data = return_data,
                         x_var = x_var,
                         y_var = y_var,
                         wide_var = wide_var)
    } else {
      summary_dt_l = lapply(q_tall_values, function(cl) {
        prep_summary(prof_dt[get(tall_var) == cl], 
                     position_dt = pos_dt, 
                     x_points = x_points, 
                     y_points = y_points, 
                     xrng = xrng, 
                     yrng = yrng, 
                     facet_by = NULL, 
                     x_var = x_var, 
                     y_var = y_var, 
                     id_var = id_var, 
                     wide_var = wide_var, 
                     tall_var = tall_var)
      })
      names(summary_dt_l) = q_tall_values
      summary_dt = rbindlist(summary_dt_l, 
                             use.names = TRUE, 
                             idcol = tall_var)
      if (return_data) {
        plot_summary_glyph(summary_dt = summary_dt, 
                           xrng = xrng, 
                           yrng = yrng, 
                           x_points = x_points, 
                           y_points = y_points, 
                           ylim = ylim, 
                           N_floor = N_floor, 
                           N_ceiling = N_ceiling, 
                           min_size = min_size, 
                           color_mapping = line_color_mapping, 
                           return_data = return_data,
                           x_var = x_var,
                           y_var = y_var,
                           wide_var = wide_var)
      }
      else {
        plot_summary_glyph(summary_dt, 
                           x_points = x_points, 
                           y_points = y_points, 
                           xrng = xrng, 
                           yrng = yrng, 
                           ylim = ylim, 
                           N_floor = N_floor, 
                           N_ceiling = N_ceiling, 
                           min_size = min_size, 
                           color_mapping = line_color_mapping,
                           x_var = x_var,
                           y_var = y_var,
                           wide_var = wide_var) + 
          facet_wrap(paste0("~", tall_var))
      }
    }
  }
}

prep_summary = function (profile_dt, 
                         position_dt, 
                         x_points, 
                         y_points = x_points, 
                         xrng = range(position_dt$tx), 
                         yrng = range(position_dt$ty), 
                         facet_by = NULL,
                         x_var = "x", 
                         y_var = "y", 
                         id_var = "id",
                         wide_var = "name", 
                         tall_var = "tall_none"){
  position_dt = copy(position_dt[tx >= min(xrng) & tx <= max(xrng) & 
                                   ty >= min(yrng) & ty <= max(yrng)])
  position_dt = position_dt[get(id_var) %in% unique(profile_dt[[id_var]])]
  if (is.null(position_dt$bx)) 
    position_dt[, `:=`(bx, bin_values(tx, x_points, 
                                      xrng = xrng))]
  if (is.null(position_dt$by)) 
    position_dt[, `:=`(by, bin_values(ty, y_points, 
                                      xrng = yrng))]
  summary_dt = merge(profile_dt, position_dt[, c("bx", "by", tall_var, id_var), with = FALSE], 
                     allow.cartesian = TRUE, 
                     by = intersect(colnames(profile_dt), c(tall_var, id_var)))
  if (is.null(summary_dt[[wide_var]])) 
    summary_dt[[wide_var]] = "signal"
  if (is.null(facet_by)) {
    summary_dt = summary_dt[, list(y_tmp_ = mean(get(y_var))), c("bx", "by", x_var, wide_var)]
  }
  else {
    summary_dt = summary_dt[, list(y_tmp_ = mean(get(y_var))), c("bx", "by", x_var, wide_var, facet_by)]
  }
  setnames(summary_dt, "y_tmp_", y_var)
  N_dt = position_dt[, .(.N), by = .(bx, by)]
  summary_dt = merge(summary_dt, N_dt, by = c("bx", "by"))
  summary_dt[, `:=`(plot_id, paste(bx, by, sep = "_"))]
  summary_dt[]
}

bin_values = function (x, n_bins, xrng = range(x)) 
{
  stopifnot(length(xrng) == 2)
  floor(rescale_capped(x, 0:1, xrng) * (n_bins - 1e-05)) + 
    1
}

bin_values_centers = function (n_bins, rng) 
{
  if (length(rng) != 2) 
    rng = range(rng)
  stopifnot(length(rng) == 2)
  xspc = diff(rng)/n_bins/2
  xs = seq(min(rng) + xspc, max(rng) - xspc, diff(rng)/(n_bins))
  xs
}

nn_clust = function (tsne_res, nsamp = Inf, nn = 100, tall_var = "tall_none", id_var = "id") 
{
  valid_vars = c(ifelse(tall_var == "tall_none", character(), tall_var), id_var)
  valid_vars = valid_vars[!is.na(valid_vars)]
  stopifnot(valid_vars %in% colnames(tsne_res))
  
  if(is.null(tsne_res[[tall_var]])){
    tsne_res[[tall_var]] = "none"
  }
  
  if(nn > .2 * nrow(tsne_res)){
    
    nn = floor(.2 * nrow(tsne_res))
    message("Decreasing nearest-neighbors to ", nn, ".  Original value was too high for dataset.")
  }
  tsne_res[, `:=`(tid, paste(get(tall_var), get(id_var)))]  
  
  mat = t(as.matrix(tsne_res[, .(tx, ty)]))
  colnames(mat) = tsne_res$tid
  mat = mat[, sampleCap(seq(ncol(mat)), nsamp)]
  knn.info <- RANN::nn2(t(mat), k = nn)
  knn <- knn.info$nn.idx
  colnames(knn) = c("tid", paste0("V", seq(nn - 1)))
  knn = as.data.table(knn)
  mknn = melt(knn, id.vars = "tid")
  ADJ = Matrix::Matrix(0, ncol(mat), ncol(mat))
  ADJ[cbind(mknn$tid, mknn$value)] = 1
  rownames(ADJ) = colnames(mat)
  colnames(ADJ) = colnames(mat)
  g <- igraph::graph.adjacency(ADJ, mode = "undirected")
  g <- igraph::simplify(g)
  km <- igraph::cluster_walktrap(g)
  com <- km$membership
  names(com) <- km$names
  com_dt = data.table(tid = names(com), cluster_id = com)
  p_dt = merge(tsne_res, com_dt, by = "tid")
  p = ggplot(p_dt, aes(x = tx, y = ty, color = as.character(cluster_id))) + 
    labs(color = "cluster_id") + geom_point(size = 0.5)
  return(list(data = p_dt, plot = p))
}
