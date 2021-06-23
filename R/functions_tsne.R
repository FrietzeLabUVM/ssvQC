
#' Title
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
                     wide_var = "name", 
                     tall_var = "_none_") {
  valid_vars = c("id", ifelse(tall_var == "_none_", character(), tall_var), wide_var, "x", "y")
  valid_vars = valid_vars[!is.na(valid_vars)]
  stopifnot(valid_vars %in% colnames(profile_dt))
  tsne_mat = dt2mat(profile_dt, unique(profile_dt[[wide_var]]))
  bad_col = apply(tsne_mat, 2, stats::var) == 0
  if (any(bad_col)) {
    stop("zero variance columns detected in tsne matrix input.", 
         "\n", round(sum(bad_col)/length(bad_col) * 100, 
                     2), "% of columns affected.")
  }
  if (is.data.table(Y_init)) {
    if(tall_var != "_none_"){
      Y_init = as.matrix(Y_init[, .(tx, ty)], rownames.value = paste(Y_init$id, 
                                                                     Y_init$tall_var))  
    }else{
      Y_init = as.matrix(Y_init[, .(tx, ty)], rownames.value = paste(Y_init$id))
    }
    
    Y_init = Y_init[rownames(tsne_mat), ]
    stopifnot(ncol(Y_init) == 2)
  }
  else {
    stopifnot(is.null(Y_init))
  }
  res_tsne = Rtsne::Rtsne(tsne_mat, 
                          Y_init = Y_init, 
                          num_threads = n_cores, 
                          perplexity = perplexity, 
                          check_duplicates = FALSE)
  tsne_dt = as.data.table(res_tsne$Y)
  colnames(tsne_dt) = c("tx", "ty")
  tsne_dt$rn = rownames(tsne_mat)
  
  if(tall_var != "_none_"){
    tsne_dt[, `:=`(c("id", "tall_var"), tstrsplit(rn, " ", keep = seq(2)))]
  }else{
    tsne_dt[, `:=`(c("id"), rn)]
  }
  
  
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

dt2mat = function (prof_dt, wide_values, wide_var = "name", tall_var = "_none_"){
  
  if(tall_var != "_none_"){
    for (m in wide_values) {
      if (m == wide_values[1]) {
        dt = dcast(prof_dt[get(wide_var) == m], id + tall_var ~ 
                     x, value.var = "y")
        wide_mat = as.matrix(dt[, -seq_len(2)])
        rn = paste(dt$id, dt$tall_var)
      }
      else {
        dt = dcast(prof_dt[get(wide_var) == m], id + tall_var ~ 
                     x, value.var = "y")
        stopifnot(all(paste(dt$id, dt$tall_var) == rn))
        wide_mat = cbind(wide_mat, as.matrix(dt[, -seq_len(2)]))
      }
    }
  }else{
    for (m in wide_values) {
      if (m == wide_values[1]) {
        dt = dcast(prof_dt[get(wide_var) == m], id ~ 
                     x, value.var = "y")
        wide_mat = as.matrix(dt[, -seq_len(2)])
        rn = dt$id
      }
      else {
        dt = dcast(prof_dt[get(wide_var) == m], id ~ 
                     x, value.var = "y")
        stopifnot(all(dt$id == rn))
        wide_mat = cbind(wide_mat, as.matrix(dt[, -seq_len(2)]))
      }
    }
  }
  
  rownames(wide_mat) = rn
  wide_mat
}

plot_summary_glyph = function (tsne_dt, qtall_vars, id_to_plot = NULL, p = NULL, xrng = c(-0.5, 
                                                                                          0.5), yrng = c(-0.5, 0.5), bg_color = "gray", line_color_mapping = "black", 
                               fill_color_mapping = "gray", label_type = c("text", 
                                                                           "label", "none")[3], bg_points = 5000, arrow_FUN = NULL) 
{
  stopifnot(qtall_vars %in% unique(tsne_dt$tall_var))
  if (is.numeric(label_type)) {
    label_type = c("text", "label", "none")[label_type]
  }
  if (is.null(id_to_plot)) {
    id_to_plot = unique(tsne_dt$id)
  }
  stopifnot(id_to_plot %in% tsne_dt$id)
  lines_dt = tsne_dt[tall_var %in% qtall_vars & id %in% id_to_plot]
  lines_dt$tall_var = factor(lines_dt$tall_var, levels = qtall_vars)
  lines_dt = lines_dt[order(tall_var)][order(id)]
  lo = (seq(id_to_plot)%%length(line_color_mapping)) + 1
  line_color_mapping = line_color_mapping[lo]
  names(line_color_mapping) = id_to_plot
  fo = (seq(id_to_plot)%%length(fill_color_mapping)) + 1
  fill_color_mapping = fill_color_mapping[fo]
  names(fill_color_mapping) = id_to_plot
  if (bg_points < 0) {
    id_tp = character()
  }
  else if (bg_points == 0) {
    id_tp = id_to_plot
  }
  else {
    id_tp = sampleCap(unique(tsne_dt$id), bg_points)
    id_tp = union(id_tp, id_to_plot)
  }
  if (is.null(p)) {
    p = ggplot() + labs(title = paste(qtall_vars, collapse = ", ")) + 
      theme_classic() + coord_cartesian(xlim = xrng, ylim = yrng)
  }
  p = p + annotate("point", x = tsne_dt[id %in% id_tp, 
  ]$tx, y = tsne_dt[id %in% id_tp, ]$ty, color = bg_color)
  ch_dt = lines_dt[, .(ch_i = grDevices::chull(tx, ty)), .(id)]
  lines_dt[, `:=`(ch_i, seq(.N)), by = .(id)]
  ch_res = lines_dt[, .(ch_i = grDevices::chull(tx, ty)), by = .(id)]
  ch_res$o = seq(nrow(ch_res))
  poly_dt = merge(lines_dt, ch_res)
  poly_dt = poly_dt[order(o)]
  for (tid in unique(poly_dt$id)) {
    p = p + annotate("polygon", x = poly_dt[id == tid]$tx, 
                     y = poly_dt[id == tid]$ty, color = line_color_mapping[tid], 
                     fill = fill_color_mapping[tid])
  }
  lab_dt = lines_dt[, .(tx = mean(tx), ty = mean(ty)), by = .(id)]
  switch(label_type, text = {
    p = p + ggrepel::geom_text_repel(data = lab_dt, aes(x = tx, 
                                                        y = ty, label = id), color = "black", show.legend = FALSE)
  }, label = {
    p = p + ggrepel::geom_label_repel(data = lab_dt, aes(x = tx, 
                                                         y = ty, label = id), color = "black", fill = "white", 
                                      show.legend = FALSE)
  }, none = {
    p = p
  })
  p
}

plot_summary_raster = function (image_dt, x_points, y_points = x_points, xrng = c(-0.5, 
                                                                                  0.5), yrng = c(-0.5, 0.5), p = NULL, line_color_mapping = NULL, 
                                N_floor = 0, N_ceiling = NULL, min_size = 0.3, return_data = FALSE) 
{
  image_dt = copy(image_dt)
  image_dt = set_image_rects(image_dt, x_points = x_points, 
                             y_points = y_points, xrng = xrng, yrng = yrng, N_floor = N_floor, 
                             N_ceiling = N_ceiling, min_size = min_size)
  if (return_data) {
    return(image_dt)
  }
  if (is.null(p)) 
    p = ggplot()
  if (!is.null(line_color_mapping)) {
    col_dt = image_dt[, .(tx, ty, wide_var = rep(names(line_color_mapping), 
                                                 each = nrow(image_dt)))]
    p = p + geom_point(data = col_dt, aes(x = tx, y = ty, 
                                          color = wide_var)) + scale_color_manual(values = line_color_mapping)
  }
  p = p + geom_image.rect(data = image_dt, aes(xmin = xmin, 
                                               xmax = xmax, ymin = ymin, ymax = ymax, image = png_file)) + 
    geom_rect(data = image_dt, aes(xmin = xmin, xmax = xmax, 
                                   ymin = ymin, ymax = ymax), fill = NA, color = "black") + 
    coord_cartesian(xlim = xrng, ylim = yrng)
  p
}

stsPlotSummaryProfiles = function (profile_dt, 
                                   position_dt, 
                                   x_points, 
                                   y_points = x_points, 
                                   wide_var = "name",
                                   tall_var = "_none_",
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
  browser()
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
  prof_dt = copy(profile_dt[get(tall_var) %in% q_tall_values & wide_var %in% 
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
      summary_dt = prep_summary(prof_dt, pos_dt, x_points, 
                                y_points, xrng, yrng, NULL)
      img_res = prep_images(summary_dt = summary_dt, xrng = xrng, 
                            yrng = yrng, x_points = x_points, y_points = y_points, 
                            rname = rname, odir = odir, force_rewrite = force_rewrite, 
                            apply_norm = apply_norm, ylim = ylim, ma_size = ma_size, 
                            n_splines = n_splines, n_cores = n_cores, line_color_mapping = line_color_mapping, 
                            vertical_facet_mapping = vertical_facet_mapping)
      plot_summary_raster(image_dt = img_res$image_dt, 
                          xrng = xrng, yrng = yrng, x_points = x_points, 
                          y_points = y_points, p = p, line_color_mapping = img_res$line_color_mapping, 
                          N_floor = N_floor, N_ceiling = N_ceiling, min_size = min_size, 
                          return_data = return_data)
    }
    else {
      summary_dt = prep_summary(prof_dt, pos_dt, x_points, 
                                y_points, xrng, yrng, tall_var)
      img_res = prep_images(summary_dt = summary_dt, xrng = xrng, 
                            yrng = yrng, x_points = x_points, y_points = y_points, 
                            rname = rname, odir = odir, force_rewrite = force_rewrite, 
                            apply_norm = apply_norm, ylim = ylim, facet_by = tall_var, 
                            ma_size = ma_size, n_splines = n_splines, n_cores = n_cores, 
                            line_color_mapping = line_color_mapping, vertical_facet_mapping = vertical_facet_mapping)
      plot_summary_raster_byCell(image_dt = img_res$image_dt, 
                                 xrng = xrng, yrng = yrng, x_points = x_points, 
                                 y_points = y_points, p = p, line_color_mapping = img_res$line_color_mapping, 
                                 N_floor = N_floor, N_ceiling = N_ceiling, min_size = min_size, 
                                 return_data = return_data)
    }
  }
  else if (plot_type == "glyph") {
    if (!facet_byCell) {
      summary_dt = prep_summary(prof_dt, pos_dt, x_points, 
                                y_points, xrng, yrng, NULL)
      plot_summary_glyph(summary_dt, x_points = x_points, 
                         y_points = y_points, xrng = xrng, yrng = yrng, 
                         ylim = ylim, N_floor = N_floor, N_ceiling = N_ceiling, 
                         min_size = min_size, color_mapping = line_color_mapping, 
                         return_data = return_data)
    }
    else {
      summary_dt_l = lapply(q_tall_values, function(cl) {
        prep_summary(prof_dt[get(tall_var) == cl], pos_dt, 
                     x_points, y_points, xrng, yrng, NULL)
      })
      names(summary_dt_l) = q_tall_values
      summary_dt = rbindlist(summary_dt_l, use.names = TRUE, 
                             idcol = tall_var)
      if (return_data) {
        plot_summary_glyph(summary_dt, x_points = x_points, 
                           y_points = y_points, xrng = xrng, yrng = yrng, 
                           ylim = ylim, N_floor = N_floor, N_ceiling = N_ceiling, 
                           min_size = min_size, color_mapping = line_color_mapping, 
                           return_data = return_data)
      }
      else {
        plot_summary_glyph(summary_dt, x_points = x_points, 
                           y_points = y_points, xrng = xrng, yrng = yrng, 
                           ylim = ylim, N_floor = N_floor, N_ceiling = N_ceiling, 
                           min_size = min_size, color_mapping = line_color_mapping) + 
          facet_wrap(paste0("~", tall_var))
      }
    }
  }
}

