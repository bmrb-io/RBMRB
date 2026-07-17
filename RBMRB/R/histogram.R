# =============================================================================
# RBMRB -- Histogram module  (mirrors pybmrb.Histogram)
# =============================================================================

#' Plot a chemical shift distribution histogram
#'
#' Fetches database-wide chemical shift data and renders an interactive
#' distribution plot. Mirrors `pybmrb.Histogram.hist()`.
#'
#' @param residue  Residue name (IUPAC 3-letter). Use "*" for any.
#' @param atom     Atom name (e.g. "CA", "HB*", "*").
#' @param list_of_atoms  Optional character vector of "RES-ATOM" pairs
#'   (e.g. `c("ALA-CB","CYS-N")`). Overrides residue/atom when supplied.
#' @param filtered   Remove outliers. Default TRUE.
#' @param sd_limit   SD factor for outlier removal. Default 10.
#' @param ambiguity  Ambiguity code filter ("*" = no filter).
#' @param ph_min,ph_max  pH range filter.
#' @param t_min,t_max    Temperature range filter (Kelvin).
#' @param standard_amino_acids  Restrict to standard amino acids. Default TRUE.
#' @param histnorm   NULL (count), "percent", "probability", or
#'   "probability density".
#' @param plot_type  "histogram" (default), "box", or "violin".
#' @param output_format  "html" (default), "png", "jpg", "pdf", "svg", "webp".
#' @param output_file  File path to save the plot (optional).
#' @param output_image_width,output_image_height  Image size in pixels.
#' @param show_visualization  Display the figure. Default TRUE.
#' @return Invisibly, a list with elements `values` and `labels`.
#' @export
#' @examples
#' \dontrun{
#'   cs_hist("ALA", "CA")
#'   cs_hist(list_of_atoms = c("ALA-CB","GLY-CA","SER-CB"))
#'   cs_hist("ALA", "CA", histnorm = "probability density",
#'           plot_type = "violin")
#' }
cs_hist <- function(residue = NULL,
                    atom    = NULL,
                    list_of_atoms = NULL,
                    filtered = TRUE, sd_limit = 10,
                    ambiguity = "*",
                    ph_min = NULL, ph_max = NULL,
                    t_min  = NULL, t_max  = NULL,
                    standard_amino_acids = TRUE,
                    histnorm = NULL,
                    plot_type = "histogram",
                    output_format = "html",
                    output_file   = NULL,
                    output_image_width  = 800,
                    output_image_height = 600,
                    show_visualization  = TRUE) {

  plot_type <- match.arg(plot_type, c("histogram","box","violin"))
  specs     <- .resolve_atom_specs(residue, atom, list_of_atoms)

  all_vals   <- numeric(0)
  all_labels <- character(0)

  for (spec in specs) {
    df <- get_cs_data(residue=spec$residue, atom=spec$atom,
                       filtered=filtered, sd_limit=sd_limit,
                       ambiguity=ambiguity,
                       ph_min=ph_min, ph_max=ph_max,
                       t_min=t_min, t_max=t_max,
                       standard_amino_acids=standard_amino_acids)
    if (nrow(df) == 0L) { cli::cli_alert_warning("No data for {spec$label}"); next }
    v <- df$Val[!is.na(df$Val)]
    all_vals   <- c(all_vals,   v)
    all_labels <- c(all_labels, rep(spec$label, length(v)))
  }

  if (length(all_vals) == 0L) {
    cli::cli_alert_danger("No data retrieved.")
    return(invisible(list(values=numeric(0), labels=character(0))))
  }

  fig <- .build_hist_fig(all_vals, all_labels,
                          plot_type=plot_type, histnorm=histnorm,
                          title=.hist_title(specs),
                          x_label="Chemical Shift (ppm)")
  .output_figure(fig, output_format, output_file,
                 output_image_width, output_image_height, show_visualization)
  invisible(list(values=all_vals, labels=all_labels))
}

#' Plot a 2D chemical shift correlation map
#'
#' Mirrors `pybmrb.Histogram.hist2d()`.
#'
#' @param residue   Residue name.
#' @param atom1     First atom (X axis).
#' @param atom2     Second atom (Y axis).
#' @param filtered,sd_limit  See [cs_hist()].
#' @param ambiguity1,ambiguity2  Per-atom ambiguity filters.
#' @param ph_min,ph_max,t_min,t_max  Condition filters.
#' @param histnorm  Normalisation (see [cs_hist()]).
#' @param plot_type "heatmap" (default) or "contour".
#' @param output_format,output_file,output_image_width,output_image_height,show_visualization
#'   See [cs_hist()].
#' @return Invisibly, a list with elements `atom1` and `atom2`.
#' @export
#' @examples
#' \dontrun{
#'   cs_hist2d("CYS", "CB", "N")
#'   cs_hist2d("ALA", "CA", "CB", plot_type = "contour")
#' }
cs_hist2d <- function(residue, atom1, atom2,
                       filtered=TRUE, sd_limit=10,
                       ambiguity1="*", ambiguity2="*",
                       ph_min=NULL, ph_max=NULL,
                       t_min=NULL,  t_max=NULL,
                       histnorm=NULL,
                       plot_type="heatmap",
                       output_format="html",
                       output_file=NULL,
                       output_image_width=800,
                       output_image_height=600,
                       show_visualization=TRUE) {

  plot_type <- match.arg(plot_type, c("heatmap","contour"))
  paired    <- get_2d_cs_data(residue=residue, atom1=atom1, atom2=atom2,
                               filtered=filtered, sd_limit=sd_limit,
                               ambiguity1=ambiguity1, ambiguity2=ambiguity2,
                               ph_min=ph_min, ph_max=ph_max,
                               t_min=t_min, t_max=t_max)

  if (length(paired$atom1) == 0L) {
    cli::cli_alert_danger("No paired data for {residue} {atom1}/{atom2}.")
    return(invisible(paired))
  }

  fig <- .build_hist2d_fig(
    x=paired$atom1, y=paired$atom2,
    plot_type=plot_type, histnorm=histnorm,
    title=sprintf("%s %s vs %s (%d points)",
                   residue, atom1, atom2, length(paired$atom1)),
    x_label=paste(atom1, "Chemical Shift (ppm)"),
    y_label=paste(atom2, "Chemical Shift (ppm)"))

  .output_figure(fig, output_format, output_file,
                 output_image_width, output_image_height, show_visualization)
  invisible(paired)
}

#' Plot a conditional chemical shift histogram
#'
#' Shows the CS distribution of `atom` in `residue` before and after
#' filtering by known chemical shifts of other atoms in the same residue.
#' Mirrors `pybmrb.Histogram.conditional_hist()`.
#'
#' @param residue  Residue name.
#' @param atom     Atom whose distribution to plot.
#' @param filtering_rules  Named list of atom = value pairs used to filter,
#'   e.g. `list(CA = 64.5, H = 7.8)`.
#' @param h_tolerance  Tolerance (ppm) for H filters. Default 0.1.
#' @param c_tolerance  Tolerance (ppm) for C filters. Default 2.0.
#' @param n_tolerance  Tolerance (ppm) for N filters. Default 2.0.
#' @param ph_min,ph_max,t_min,t_max  Condition filters.
#' @param standard_amino_acids  Restrict to standard amino acids. Default TRUE.
#' @param histnorm  Normalisation (see [cs_hist()]).
#' @param plot_type "histogram" (default), "box", or "violin".
#' @param output_format,output_file,output_image_width,output_image_height,show_visualization
#'   See [cs_hist()].
#' @return Invisibly, a list with `all_values`, `filtered_values`, `labels`.
#' @export
#' @examples
#' \dontrun{
#'   conditional_cs_hist("SER", "N",
#'     filtering_rules = list(CA = 55.5, HA = 4.3))
#' }
conditional_cs_hist <- function(residue, atom,
                                 filtering_rules,
                                 ph_min=NULL, ph_max=NULL,
                                 t_min=NULL,  t_max=NULL,
                                 h_tolerance=0.1,
                                 c_tolerance=2.0,
                                 n_tolerance=2.0,
                                 standard_amino_acids=TRUE,
                                 histnorm=NULL,
                                 plot_type="histogram",
                                 output_format="html",
                                 output_file=NULL,
                                 output_image_width=800,
                                 output_image_height=600,
                                 show_visualization=TRUE) {

  plot_type <- match.arg(plot_type, c("histogram","box","violin"))
  if (!is.list(filtering_rules) || length(filtering_rules) == 0L)
    cli::cli_abort("`filtering_rules` must be a non-empty named list.")

  df_all <- get_cs_from_bmrb_db(
    residue=residue, atom=atom,
    ph_min=ph_min, ph_max=ph_max, t_min=t_min, t_max=t_max,
    standard_amino_acids=standard_amino_acids)

  if (nrow(df_all) == 0L) {
    cli::cli_alert_danger("No data for {residue} {atom}.")
    return(invisible(NULL))
  }

  key_cols <- intersect(c("Entry_ID","Comp_index_ID","Entity_assembly_ID"),
                         colnames(df_all))
  df_all[["._key"]] <- do.call(paste, c(df_all[key_cols], sep="_"))
  valid_keys <- df_all[["._key"]]

  for (fa in names(filtering_rules)) {
    ref_val  <- as.numeric(filtering_rules[[fa]])
    tol      <- switch(substr(fa, 1L, 1L), H=h_tolerance, C=c_tolerance,
                       N=n_tolerance, 1.0)
    df_f <- get_cs_from_bmrb_db(residue=residue, atom=fa,
                                  ph_min=ph_min, ph_max=ph_max,
                                  t_min=t_min, t_max=t_max,
                                  standard_amino_acids=standard_amino_acids)
    if (nrow(df_f) == 0L) { cli::cli_alert_warning("No data for {fa}"); next }
    df_f <- df_f[!is.na(df_f$Val) & abs(df_f$Val - ref_val) <= tol, ]
    df_f[["._key"]] <- do.call(paste, c(df_f[key_cols], sep="_"))
    valid_keys <- intersect(valid_keys, df_f[["._key"]])
  }

  all_vals      <- df_all$Val[!is.na(df_all$Val)]
  filtered_vals <- df_all$Val[!is.na(df_all$Val) &
                                df_all[["._key"]] %in% valid_keys]

  cli::cli_alert_info(
    "{length(filtered_vals)} / {length(all_vals)} records match all filters")

  fig <- .build_conditional_hist_fig(
    all_vals=all_vals, filtered_vals=filtered_vals,
    filter_desc=.format_filter_desc(filtering_rules),
    plot_type=plot_type, histnorm=histnorm,
    title=sprintf("%s %s -- Conditional Distribution", residue, atom),
    x_label=paste(atom, "Chemical Shift (ppm)"))

  .output_figure(fig, output_format, output_file,
                 output_image_width, output_image_height, show_visualization)
  invisible(list(all_values=all_vals, filtered_values=filtered_vals,
                 labels=c(rep("All", length(all_vals)),
                          rep("Filtered", length(filtered_vals)))))
}

# =============================================================================
# Figure builders
# =============================================================================

.build_hist_fig <- function(vals, labels, plot_type, histnorm, title, x_label) {
  groups  <- unique(labels)
  palette <- grDevices::colorRampPalette(
    c("#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD",
      "#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")
  )(max(length(groups), 1L))

  fig <- plotly::plot_ly()
  for (i in seq_along(groups)) {
    grp <- groups[[i]]; v <- vals[labels == grp]; col <- palette[[i]]
    if (plot_type == "histogram") {
      fig <- plotly::add_trace(fig, type="histogram", x=v, name=grp,
        histnorm=if (is.null(histnorm)) "" else histnorm,
        marker=list(color=.hex_alpha(col,0.7),
                    line=list(color=col, width=0.5)),
        nbinsx=100, opacity=0.7)
    } else if (plot_type == "violin") {
      fig <- plotly::add_trace(fig, type="violin", y=v, name=grp,
        box=list(visible=TRUE), meanline=list(visible=TRUE),
        fillcolor=.hex_alpha(col,0.5), line=list(color=col))
    } else {
      fig <- plotly::add_trace(fig, type="box", y=v, name=grp,
        marker=list(color=col), boxmean="sd")
    }
  }
  ylab <- if (!is.null(histnorm)) histnorm else "Count"
  plotly::layout(fig,
    title=list(text=title, font=list(size=15)),
    xaxis=list(title=x_label),
    yaxis=list(title=if (plot_type=="histogram") ylab
                     else "Chemical Shift (ppm)"),
    barmode="overlay",
    plot_bgcolor="white", paper_bgcolor="white")
}

.build_hist2d_fig <- function(x, y, plot_type, histnorm, title,
                               x_label, y_label) {
  fig <- plotly::plot_ly()
  if (plot_type == "heatmap") {
    fig <- plotly::add_trace(fig, type="histogram2d", x=x, y=y,
      colorscale="Viridis",
      histnorm=if (is.null(histnorm)) "" else histnorm,
      nbinsx=80, nbinsy=80, showscale=TRUE)
  } else {
    fig <- plotly::add_trace(fig, type="histogram2dcontour", x=x, y=y,
      colorscale="Viridis",
      histnorm=if (is.null(histnorm)) "" else histnorm,
      nbinsx=60, nbinsy=60, showscale=TRUE,
      contours=list(coloring="heatmap"))
    fig <- plotly::add_trace(fig, type="scatter", mode="markers",
      x=x, y=y,
      marker=list(color="rgba(0,0,0,0.08)", size=3),
      showlegend=FALSE, hoverinfo="skip")
  }
  plotly::layout(fig,
    title=list(text=title, font=list(size=15)),
    xaxis=list(title=x_label),
    yaxis=list(title=y_label),
    plot_bgcolor="white", paper_bgcolor="white")
}

.build_conditional_hist_fig <- function(all_vals, filtered_vals, filter_desc,
                                         plot_type, histnorm, title, x_label) {
  fig  <- plotly::plot_ly()
  hn   <- if (is.null(histnorm)) "" else histnorm
  col1 <- "#1F77B4"; col2 <- "#D62728"

  if (plot_type == "histogram") {
    fig <- plotly::add_trace(fig, type="histogram", x=all_vals,
      name="All BMRB", histnorm=hn, nbinsx=100, opacity=0.7,
      marker=list(color=.hex_alpha(col1,0.6),
                  line=list(color=col1, width=0.5)))
    fig <- plotly::add_trace(fig, type="histogram", x=filtered_vals,
      name=paste0("Filtered: ", filter_desc), histnorm=hn, nbinsx=100,
      opacity=0.7,
      marker=list(color=.hex_alpha(col2,0.7),
                  line=list(color=col2, width=0.5)))
  } else if (plot_type == "violin") {
    fig <- plotly::add_trace(fig, type="violin", y=all_vals,
      name="All BMRB", box=list(visible=TRUE), meanline=list(visible=TRUE),
      fillcolor=.hex_alpha(col1,0.5), line=list(color=col1))
    fig <- plotly::add_trace(fig, type="violin", y=filtered_vals,
      name=paste0("Filtered: ", filter_desc),
      box=list(visible=TRUE), meanline=list(visible=TRUE),
      fillcolor=.hex_alpha(col2,0.5), line=list(color=col2))
  } else {
    fig <- plotly::add_trace(fig, type="box", y=all_vals,
      name="All BMRB", marker=list(color=col1), boxmean="sd")
    fig <- plotly::add_trace(fig, type="box", y=filtered_vals,
      name=paste0("Filtered: ", filter_desc),
      marker=list(color=col2), boxmean="sd")
  }

  plotly::layout(fig,
    title=list(text=title, font=list(size=15)),
    barmode="overlay",
    xaxis=list(title=x_label),
    yaxis=list(title=if (plot_type=="histogram") (if (is.null(histnorm)) "Count" else histnorm)
                     else "Chemical Shift (ppm)"),
    plot_bgcolor="white", paper_bgcolor="white")
}

# =============================================================================
# Internal helpers
# =============================================================================

.resolve_atom_specs <- function(residue, atom, list_of_atoms) {
  if (!is.null(list_of_atoms)) {
    return(lapply(list_of_atoms, function(pair) {
      parts <- strsplit(pair, "-", fixed=TRUE)[[1L]]
      if (length(parts) != 2L)
        cli::cli_abort(
          "list_of_atoms entries must be 'RESIDUE-ATOM', got: {pair}")
      list(residue=parts[[1L]], atom=parts[[2L]], label=pair)
    }))
  }
  if (is.null(residue) && is.null(atom))
    cli::cli_abort("Supply at least one of `residue`, `atom`, or `list_of_atoms`.")
  rv <- if (is.null(residue)) "*" else residue
  av <- if (is.null(atom))    "*" else atom
  lb <- if (rv == "*") av else paste0(rv, "-", av)
  list(list(residue=rv, atom=av, label=lb))
}

.hist_title <- function(specs) {
  labs <- vapply(specs, `[[`, character(1L), "label")
  if (length(labs) == 1L) paste("Chemical Shift Distribution:", labs[[1L]])
  else sprintf("Chemical Shift Distribution (%d groups)", length(labs))
}

.format_filter_desc <- function(rules) {
  parts <- mapply(function(n, v) sprintf("%s=%.2f", n, v),
                  names(rules), unlist(rules))
  paste(parts, collapse=", ")
}

.hex_alpha <- function(hex, alpha=0.7) {
  rgb <- grDevices::col2rgb(hex) / 255
  sprintf("rgba(%d,%d,%d,%.2f)",
          round(rgb[1L]*255), round(rgb[2L]*255),
          round(rgb[3L]*255), alpha)
}
