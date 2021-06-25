
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
                     x_var = "x", 
                     y_var = "y", 
                     id_var = "id",
                     wide_var = "name", 
                     tall_var = "tall_none") {
  valid_vars = c(ifelse(tall_var == "tall_none", character(), tall_var), wide_var, x_var, y_var, id_var)
  valid_vars = valid_vars[!is.na(valid_vars)]
  stopifnot(valid_vars %in% colnames(profile_dt))
  
  if(is.null(profile_dt[[tall_var]])){
    profile_dt[[tall_var]] = "none"
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
    stop("zero variance columns detected in tsne matrix input.", 
         "\n", round(sum(bad_col)/length(bad_col) * 100, 
                     2), "% of columns affected.")
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
  res_tsne = Rtsne::Rtsne(tsne_mat, 
                          Y_init = Y_init, 
                          num_threads = n_cores, 
                          perplexity = perplexity, 
                          check_duplicates = FALSE)
  tsne_dt = as.data.table(res_tsne$Y)
  colnames(tsne_dt) = c("tx", "ty")
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

plot_summary_glyph.2 = function (tsne_dt, 
                                 qtall_vars, 
                                 id_to_plot = NULL, 
                                 p = NULL, 
                                 xrng = c(-0.5, 0.5), 
                                 yrng = c(-0.5, 0.5), 
                                 bg_color = "gray", 
                                 line_color_mapping = "black", 
                                 fill_color_mapping = "gray", 
                                 label_type = c("text", "label", "none")[3], 
                                 bg_points = 5000, 
                                 arrow_FUN = NULL) 
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

plot_summary_raster = function (image_dt, 
                                x_points, 
                                y_points = x_points, 
                                xrng = c(-0.5, 0.5), 
                                yrng = c(-0.5, 0.5), 
                                p = NULL, 
                                line_color_mapping = NULL, 
                                N_floor = 0, 
                                N_ceiling = NULL, 
                                min_size = 0.3, 
                                return_data = FALSE, 
                                wide_var = "name") 
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
    col_dt = image_dt[, .(tx, ty, wide_var_ = rep(names(line_color_mapping), 
                                                  each = nrow(image_dt)))]
    setnames(col_dt, "wide_var_", wide_var)
    p = p + 
      geom_point(data = col_dt, 
                 aes_string(x = "tx", y = "ty", color = wide_var)) + 
      scale_color_manual(values = line_color_mapping)
  }
  p = p + 
    geom_image.rect(data = image_dt, 
                    aes(xmin = xmin, 
                        xmax = xmax, 
                        ymin = ymin,
                        ymax = ymax, 
                        image = png_file)) + 
    geom_rect(data = image_dt, 
              aes(xmin = xmin, 
                  xmax = xmax, 
                  ymin = ymin, 
                  ymax = ymax), 
              fill = NA, 
              color = "black") + 
    coord_cartesian(xlim = xrng, ylim = yrng)
  p
}

plot_summary_raster_byCell = function (image_dt, 
                                       x_points, 
                                       y_points = x_points, 
                                       p = NULL, 
                                       line_color_mapping = NULL, 
                                       xrng = c(-0.5, 0.5), 
                                       yrng = c(-0.5, 0.5), 
                                       N_floor = 0, 
                                       N_ceiling = NULL, 
                                       min_size = 0.3, 
                                       return_data = FALSE){
  tall_var = bx = by = NULL
  image_dt = set_image_rects(image_dt, x_points = x_points, 
                             y_points = y_points, xrng = xrng, yrng = yrng, N_floor = N_floor, 
                             N_ceiling = N_ceiling, min_size = min_size)
  if (return_data) 
    return(image_dt)
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
    coord_cartesian(xlim = xrng, ylim = yrng) + facet_wrap("tall_var", 
                                                           drop = FALSE)
  p
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

plot_summary_glyph = function (summary_dt, 
                               x_points, 
                               y_points = x_points, 
                               xrng = c(-0.5, 0.5), 
                               yrng = c(-0.5, 0.5), 
                               ylim = NULL, 
                               p = NULL, 
                               N_floor = 0, 
                               N_ceiling = NULL, 
                               min_size = 0.3, 
                               return_data = FALSE, 
                               color_mapping = NULL, 
                               x_var = "x",
                               y_var = "y",
                               wide_var = "name") 
{
  group_size = gx = gy = gid = NULL
  summary_dt = set_size(summary_dt, N_floor, N_ceiling, size.name = "group_size")
  summary_dt = summary_dt[group_size >= min_size]
  if (is.null(ylim)) {
    ylim = range(summary_dt[[y_var]])
  }
  ylim = range(ylim)
  set(summary_dt, i = which(summary_dt[[y_var]] < min(ylim)), j = y_var, value = min(ylim))
  set(summary_dt, i = which(summary_dt[[y_var]] > max(ylim)), j = y_var, value = max(ylim))
  set(summary_dt, j = x_var, value = summary_dt[[x_var]] * summary_dt$group_size)
  set(summary_dt, j = y_var, value = summary_dt[[y_var]] * summary_dt$group_size)
  xs = bin_values_centers(x_points, rng = xrng)
  ys = bin_values_centers(y_points, rng = yrng)
  summary_dt[, `:=`(tx, xs[bx])]
  summary_dt[, `:=`(ty, ys[by])]
  down_scale = max(summary_dt$group_size)
  glyph_dt = as.data.table(GGally::glyphs(summary_dt, x_major = "tx", 
                                          x_minor = x_var, y_major = "ty", y_minor = y_var, 
                                          width = diff(xrng)/x_points * 0.95 * down_scale, height = diff(yrng)/y_points * 
                                            0.95 * down_scale))
  if (return_data) {
    return(glyph_dt)
  }
  if (is.null(color_mapping)) {
    if (is.factor(summary_dt[[wide_var]])) {
      uwide_vars = levels(summary_dt[[wide_var]])
    }
    else if (is.character(summary_dt[[wide_var]])) {
      uwide_vars = unique(summary_dt[[wide_var]])
    }
    color_mapping = seqsetvis::safeBrew(length(uwide_vars))
  }
  if (is.null(p)) {
    p = ggplot()
  }
  glyph_dt[, glyph_group := paste(gid, get(wide_var))]
  p + geom_path(data = glyph_dt, 
                aes_string("gx", 
                           "gy", 
                           group = "glyph_group", 
                           color = wide_var)) + 
    labs(x = "tx", y = "ty") + 
    scale_color_manual(values = color_mapping) + 
    coord_cartesian(xlim = xrng, ylim = yrng)
}

set_size = function (dt, N_floor, N_ceiling, size.name = "img_size") 
{
  tmp_var = NULL
  stopifnot("N" %in% colnames(dt))
  if (is.null(N_ceiling)) {
    N_ceiling = max(dt$N)
  }
  dt[, `:=`(tmp_var, N)]
  dt[tmp_var > N_ceiling, `:=`(tmp_var, N_ceiling)]
  dt[tmp_var < N_floor, `:=`(tmp_var, N_floor)]
  dt[, `:=`(tmp_var, tmp_var - N_floor)]
  dt[, `:=`(tmp_var, tmp_var/N_ceiling)]
  dt[[size.name]] = dt$tmp_var
  dt$tmp_var = NULL
  dt
}

prep_images = function (summary_dt, 
                        x_points, 
                        y_points = x_points, 
                        xrng = c(-0.5, 0.5), 
                        yrng = c(-0.5, 0.5), 
                        rname = NULL, 
                        odir = file.path(tempdir(), rname), 
                        force_rewrite = FALSE, 
                        apply_norm = TRUE, 
                        ylim = c(0, 1), 
                        facet_by = NULL, 
                        ma_size = 2, 
                        n_splines = 10, 
                        n_cores = getOption("mc.cores", 1), 
                        line_color_mapping = NULL, 
                        vertical_facet_mapping = NULL, 
                        wide_var = "name",
                        x_var = "x",
                        y_var = "y") 
{
  if (is.null(line_color_mapping)) {
    line_color_mapping = seqsetvis::safeBrew(length(unique(summary_dt[[wide_var]])))
    names(line_color_mapping) = unique(summary_dt[[wide_var]])
  }
  if (!all(unique(summary_dt[[wide_var]]) %in% names(line_color_mapping))) {
    missing_colors = setdiff(unique(summary_dt[[wide_var]]), 
                             names(line_color_mapping))
    stop("line_color_mapping is missing assignments for: ", 
         paste(missing_colors, collapse = ", "))
  }
  if (is.null(rname)) {
    rname = digest::digest(list(summary_dt, x_points, y_points, 
                                apply_norm, ylim, line_color_mapping, n_splines, 
                                ma_size, facet_by))
  }
  dir.create(odir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(facet_by)) {
    img_dt = unique(summary_dt[, list(bx, by, plot_id)])
    img_dt[, `:=`(png_file, file.path(odir, paste0(plot_id, 
                                                   ".png")))]
  }
  else {
    img_dt = unique(summary_dt[, list(bx, by, plot_id, get(facet_by))])
    colnames(img_dt)[4] = facet_by
    img_dt[, `:=`(png_file, file.path(odir, paste0(get(facet_by), 
                                                   "_", plot_id, ".png")))]
  }
  xs = bin_values_centers(x_points, rng = xrng)
  ys = bin_values_centers(y_points, rng = yrng)
  img_dt[, `:=`(tx, xs[bx])]
  img_dt[, `:=`(ty, ys[by])]
  if (apply_norm) {
    summary_dt[, `:=`(ynorm, get(y_var)/stats::quantile(get(y_var), 0.95)), 
               by = list(wide_var)]
    summary_dt[ynorm > 1, `:=`(ynorm, 1)]
  }
  else {
    summary_dt[, `:=`(ynorm, get(y_var))]
  }
  if (force_rewrite) {
    file.remove(img_dt$png_file[file.exists(img_dt$png_file)])
  }
  if (is.null(facet_by)) {
    count_dt = unique(summary_dt[, .(bx, by, N)])
    img_dt = merge(img_dt, count_dt, by = c("bx", "by"))
  }
  else {
    count_dt = unique(summary_dt[, .(bx, by, get(facet_by), 
                                     N)])
    colnames(count_dt)[colnames(count_dt) == "V3"] = facet_by
    img_dt = merge(img_dt, count_dt, by = c("bx", "by", 
                                            facet_by))
  }
  if (is.null(vertical_facet_mapping)) {
    summary_dt$group = 1
  }
  else {
    if (!all(unique(summary_dt[[wide_var]]) %in% names(vertical_facet_mapping))) {
      missing_groups = setdiff(unique(summary_dt[[wide_var]]), 
                               names(vertical_facet_mapping))
      stop("vertical_facet_mapping is missing assignments for: ", 
           paste(missing_groups, collapse = ", "))
    }
    summary_dt$group = factor(vertical_facet_mapping[summary_dt[[wide_var]]], 
                              levels = unique(vertical_facet_mapping))
  }
  if (any(!file.exists(img_dt$png_file))) {
    plot_info = lapply(which(!file.exists(img_dt$png_file)), 
                       function(i) {
                         fpath = img_dt$png_file[i]
                         p_id = img_dt$plot_id[i]
                         if (is.null(facet_by)) {
                           pdt = summary_dt[plot_id == p_id]
                         }
                         else {
                           pdt = summary_dt[plot_id == p_id & get(facet_by) == 
                                              img_dt[[facet_by]][i]]
                         }
                         pdt[, `:=`(ysm, movingAverage(ynorm, n = ma_size)), 
                             by = c(wide_var)]
                         pdt = seqsetvis::applySpline(pdt, n = n_splines, 
                                                      by_ = wide_var, y_ = "ysm")
                         list(pdt, fpath)
                       })
    hidden = parallel::mclapply(plot_info, function(p_inf) {
      fpath = p_inf[[2]]
      pdt = p_inf[[1]]
      p = ggplot(pdt, aes_string(x = x_var, y = "ysm", ymin = 0, ymax = "ysm", 
                                 color = wide_var, fill = wide_var)) + geom_ribbon(alpha = 0.3) + 
        geom_path(size = 0.6, alpha = 1) + scale_color_manual(values = line_color_mapping) + 
        scale_fill_manual(values = line_color_mapping) + 
        theme_void() + guides(color = "none", fill = "none") + 
        facet_wrap("group", ncol = 1) + theme(strip.background = element_blank(), 
                                              strip.text = element_blank()) + coord_cartesian(ylim = ylim, 
                                                                                              xlim = c(-0.5, 0.5), expand = FALSE)
      ggsave(fpath, p, width = 2, height = 2, units = "cm")
      NULL
    }, mc.cores = n_cores)
  }
  img_dt[, `:=`(png_file, normalizePath(png_file))]
  return(list(image_dt = img_dt[], summary_profile_dt = summary_dt[], 
              x_points = x_points, y_points = y_points, xrng = xrng, 
              yrng = yrng, line_color_mapping = line_color_mapping))
}