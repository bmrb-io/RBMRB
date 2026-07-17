# =============================================================================
# RBMRB — Spectra module
# Mirrors PyBMRB's pybmrb.Spectra
# All spectrum-simulation functions share the same input/output contract.
# =============================================================================

# ---------------------------------------------------------------------------
# Public spectrum functions
# ---------------------------------------------------------------------------

#' Simulate a \ifelse{html}{\out{<sup>1</sup>H-<sup>15</sup>N HSQC}}{\eqn{^{1}H{-}^{15}N} HSQC} spectrum
#'
#' Plots an interactive ¹H–¹⁵N HSQC spectrum for one or more BMRB entries
#' and/or NMR-STAR files. Overlay multiple datasets as scatter plots; optionally
#' draw trace lines connecting peaks from the same residue position across
#' datasets. Equivalent to `pybmrb.Spectra.n15hsqc()`.
#'
#' @param bmrb_ids         BMRB entry ID(s) — integer, character, or vector.
#' @param input_file_names Path(s) to local NMR-STAR `.str` files.
#' @param auth_tag         Use `Auth_seq_ID` for sequence numbering.
#' @param legend           `NULL` (default), `"residue"`, or `"dataset"`.
#' @param draw_trace       Draw lines connecting same-residue peaks across
#'                         datasets.
#' @param include_sidechain Include sidechain NH₂ peaks (ASN/GLN/ARG/TRP/HIS).
#' @param peak_list        Path to a 2-column CSV file of unassigned peaks
#'                         (`H_ppm`, `N_ppm`) to overlay.
#' @param output_format    `"html"` (default), `"png"`, `"jpg"`, `"pdf"`,
#'                         `"svg"`, or `"webp"`.
#' @param output_file      File path to save the visualization (optional).
#' @param output_image_width,output_image_height Image dimensions in pixels.
#' @param show_visualization Open in browser / display inline. Default `TRUE`.
#'
#' @return Invisibly, a list with elements `x` (H ppm), `y` (N ppm),
#'   `dataset`, `residue`, `info`, and (when `draw_trace = TRUE`) `cs_track`.
#' @export
#' @examples
#' \dontrun{
#' # Single entry
#' n15hsqc(bmrb_ids = 15060, legend = "residue")
#'
#' # Compare two entries with trace lines
#' n15hsqc(bmrb_ids = c(17076, 17077), legend = "dataset", draw_trace = TRUE)
#'
#' # Overlay user's own NMR-STAR file with BMRB entry
#' n15hsqc(bmrb_ids = 15060, input_file_names = "MyData.str",
#'         legend = "dataset")
#' }
n15hsqc <- function(bmrb_ids          = NULL,
                    input_file_names   = NULL,
                    auth_tag           = FALSE,
                    legend             = NULL,
                    draw_trace         = FALSE,
                    include_sidechain  = TRUE,
                    peak_list          = NULL,
                    output_format      = "html",
                    output_file        = NULL,
                    output_image_width  = 800,
                    output_image_height = 600,
                    show_visualization  = TRUE) {

  pl <- create_n15hsqc_peaklist(bmrb_ids = bmrb_ids,
                                 input_file_names = input_file_names,
                                 auth_tag = auth_tag,
                                 draw_trace = draw_trace,
                                 include_sidechain = include_sidechain)
  if (is.null(pl)) { cli::cli_alert_danger("No data to plot."); return(invisible(NULL)) }

  fig <- .plot_2d_spectrum(pl,
    x_label = "\u00b9H Chemical Shift (ppm)",
    y_label = "\u00b9\u2075N Chemical Shift (ppm)",
    title   = "\u00b9H\u2013\u00b9\u2075N HSQC",
    legend  = legend,
    draw_trace = draw_trace,
    extra_peak_list = peak_list)

  .output_figure(fig, output_format, output_file,
                 output_image_width, output_image_height, show_visualization)
  invisible(pl)
}

#' Simulate a \ifelse{html}{\out{<sup>1</sup>H-<sup>13</sup>C HSQC}}{\eqn{^1H{-}^{13}C} HSQC} spectrum
#'
#' Equivalent to `pybmrb.Spectra.c13hsqc()`. See [n15hsqc()] for parameter
#' descriptions.
#'
#' @inheritParams n15hsqc
#' @return Invisibly, a list with the simulated peak list.
#' @export
#' @examples
#' \dontrun{
#' c13hsqc(bmrb_ids = 15060, legend = "residue")
#' }
c13hsqc <- function(bmrb_ids           = NULL,
                    input_file_names    = NULL,
                    auth_tag            = FALSE,
                    legend              = NULL,
                    draw_trace          = FALSE,
                    peak_list           = NULL,
                    output_format       = "html",
                    output_file         = NULL,
                    output_image_width  = 800,
                    output_image_height = 600,
                    show_visualization  = TRUE) {

  pl <- create_c13hsqc_peaklist(bmrb_ids = bmrb_ids,
                                 input_file_names = input_file_names,
                                 auth_tag = auth_tag,
                                 draw_trace = draw_trace)
  if (is.null(pl)) { cli::cli_alert_danger("No data to plot."); return(invisible(NULL)) }

  fig <- .plot_2d_spectrum(pl,
    x_label = "\u00b9H Chemical Shift (ppm)",
    y_label = "\u00b9\u00b3C Chemical Shift (ppm)",
    title   = "\u00b9H\u2013\u00b9\u00b3C HSQC",
    legend  = legend,
    draw_trace = draw_trace,
    extra_peak_list = peak_list)

  .output_figure(fig, output_format, output_file,
                 output_image_width, output_image_height, show_visualization)
  invisible(pl)
}

#' Simulate a \ifelse{html}{\out{<sup>1</sup>H-<sup>1</sup>H TOCSY}}{\eqn{^1H{-}^1H} TOCSY} spectrum
#'
#' Equivalent to `pybmrb.Spectra.tocsy()`. See [n15hsqc()] for shared
#' parameter descriptions.
#'
#' @inheritParams n15hsqc
#' @return Invisibly, a list with the simulated peak list.
#' @export
#' @examples
#' \dontrun{
#' tocsy(bmrb_ids = 15060, legend = "residue")
#' tocsy(bmrb_ids = c(17074, 17076, 17077), legend = "dataset",
#'       draw_trace = TRUE)
#' }
tocsy <- function(bmrb_ids           = NULL,
                  input_file_names    = NULL,
                  auth_tag            = FALSE,
                  legend              = NULL,
                  draw_trace          = FALSE,
                  peak_list           = NULL,
                  output_format       = "html",
                  output_file         = NULL,
                  output_image_width  = 800,
                  output_image_height = 600,
                  show_visualization  = TRUE) {

  pl <- create_tocsy_peaklist(bmrb_ids = bmrb_ids,
                               input_file_names = input_file_names,
                               auth_tag = auth_tag,
                               draw_trace = draw_trace)
  if (is.null(pl)) { cli::cli_alert_danger("No data to plot."); return(invisible(NULL)) }

  fig <- .plot_2d_spectrum(pl,
    x_label = "\u00b9H (ppm)",
    y_label  = "\u00b9H (ppm)",
    title    = "\u00b9H\u2013\u00b9H TOCSY",
    legend   = legend,
    draw_trace = draw_trace,
    extra_peak_list = peak_list)

  .output_figure(fig, output_format, output_file,
                 output_image_width, output_image_height, show_visualization)
  invisible(pl)
}

#' Plot a generic 2D NMR correlation spectrum
#'
#' Generates a 2D spectrum for any two atom types. Supports sequence-walk
#' tracing, preceding/next-residue peaks, and all the same overlay options as
#' the other spectrum functions.
#' Equivalent to `pybmrb.Spectra.generic_2d()`.
#'
#' @param atom_x Atom name for the X axis (IUPAC, e.g. `"N"`, `"CA"`).
#' @param atom_y Atom name for the Y axis.
#' @param include_preceding Include chemical shifts from the preceding (i−1)
#'   residue on the Y axis. Default `FALSE`.
#' @param include_next      Include chemical shifts from the following (i+1)
#'   residue on the Y axis. Default `FALSE`.
#' @param seq_walk Draw a trace connecting sequential i→i±1 pairs for
#'   continuous sequence segments. Default `FALSE`.
#' @param full_walk Like `seq_walk` but ignores missing residues. Default `FALSE`.
#' @inheritParams n15hsqc
#'
#' @return Invisibly, the peak list.
#' @export
#' @examples
#' \dontrun{
#' generic_2d("N", "CA", bmrb_ids = 15000, legend = "residue")
#' generic_2d("N", "CA", bmrb_ids = 15000, include_next = TRUE,
#'            seq_walk = TRUE)
#' }
generic_2d <- function(atom_x,
                       atom_y,
                       bmrb_ids            = NULL,
                       input_file_names    = NULL,
                       auth_tag            = FALSE,
                       legend              = NULL,
                       draw_trace          = FALSE,
                       peak_list           = NULL,
                       include_preceding   = FALSE,
                       include_next        = FALSE,
                       seq_walk            = FALSE,
                       full_walk           = FALSE,
                       output_format       = "html",
                       output_file         = NULL,
                       output_image_width  = 800,
                       output_image_height = 600,
                       show_visualization  = TRUE) {

  pl <- create_2d_peaklist(atom_x = atom_x, atom_y = atom_y,
                            bmrb_ids = bmrb_ids,
                            input_file_names = input_file_names,
                            auth_tag = auth_tag,
                            draw_trace = draw_trace,
                            include_preceding = include_preceding,
                            include_next = include_next)
  if (is.null(pl)) { cli::cli_alert_danger("No data to plot."); return(invisible(NULL)) }

  # Build walk traces if requested
  walk_segments <- NULL
  if (seq_walk || full_walk) {
    walk_segments <- .build_walk(pl, full_walk = full_walk)
  }

  fig <- .plot_2d_spectrum(pl,
    x_label = .axis_label(atom_x),
    y_label  = .axis_label(atom_y),
    title    = paste0(atom_x, " vs ", atom_y, " 2D Spectrum"),
    legend   = legend,
    draw_trace = draw_trace,
    extra_peak_list = peak_list,
    walk_segments = walk_segments)

  .output_figure(fig, output_format, output_file,
                 output_image_width, output_image_height, show_visualization)
  invisible(pl)
}

# ---------------------------------------------------------------------------
# Peaklist builders (no plotting)
# ---------------------------------------------------------------------------

#' Build a ¹H–¹⁵N HSQC peak list
#'
#' Extracts backbone amide (and optionally sidechain NH) chemical shifts from
#' BMRB entries / NMR-STAR files and returns the 2D peak coordinates.
#' Equivalent to `pybmrb.Spectra.create_n15hsqc_peaklist()`.
#'
#' @inheritParams n15hsqc
#' @return A list with elements: `x` (H ppm), `y` (N ppm), `dataset`,
#'   `residue`, `info`, `seq_id`, `cs_track` (if `draw_trace = TRUE`).
#' @export
create_n15hsqc_peaklist <- function(bmrb_ids         = NULL,
                                    input_file_names  = NULL,
                                    auth_tag          = FALSE,
                                    draw_trace        = FALSE,
                                    include_sidechain = TRUE) {
  dfs <- .collect_cs_dfs(bmrb_ids, input_file_names, auth_tag)
  if (length(dfs) == 0L) return(NULL)

  x_out <- y_out <- character(0)
  ds_out <- res_out <- info_out <- sid_out <- character(0)

  for (ds_name in names(dfs)) {
    df <- dfs[[ds_name]]
    seq_col   <- .seq_col(df, auth_tag)
    chain_col <- .chain_col(df, auth_tag)

    for (chain in unique(df[[chain_col]])) {
      cdf <- df[!is.na(df[[chain_col]]) & df[[chain_col]] == chain, ]
      for (sid in sort(unique(cdf[[seq_col]]))) {
        rdf <- cdf[!is.na(cdf[[seq_col]]) & cdf[[seq_col]] == sid, ]
        res <- unique(rdf$Comp_ID)[[1L]]

        # Backbone amide
        H_cs <- .get_atom_cs(rdf, "H")
        N_cs <- .get_atom_cs(rdf, "N")
        if (!is.na(H_cs) && !is.na(N_cs)) {
          x_out   <- c(x_out, H_cs)
          y_out   <- c(y_out, N_cs)
          ds_out  <- c(ds_out, ds_name)
          res_out <- c(res_out, res)
          info_out <- c(info_out, sprintf("%s%s-%s H/N", res, sid, chain))
          sid_out <- c(sid_out, as.character(sid))
        }

        # Sidechain NH groups
        if (include_sidechain && res %in% names(.SIDECHAIN_NH)) {
          for (pair in .SIDECHAIN_NH[[res]]) {
            H2 <- .get_atom_cs(rdf, pair$H)
            N2 <- .get_atom_cs(rdf, pair$N)
            if (!is.na(H2) && !is.na(N2)) {
              x_out    <- c(x_out, H2)
              y_out    <- c(y_out, N2)
              ds_out   <- c(ds_out, ds_name)
              res_out  <- c(res_out, res)
              info_out <- c(info_out, sprintf("%s%s-%s %s/%s", res, sid, chain, pair$H, pair$N))
              sid_out  <- c(sid_out, as.character(sid))
            }
          }
        }
      }
    }
  }

  if (length(x_out) == 0L) return(NULL)
  cs_track <- if (draw_trace) .build_cs_track(x_out, y_out, ds_out, sid_out) else list()

  list(x = as.numeric(x_out), y = as.numeric(y_out),
       dataset  = ds_out, residue = res_out,
       info     = info_out, seq_id = sid_out,
       cs_track = cs_track,
       x_atom   = "H", y_atom = "N")
}

#' Build a ¹H–¹³C HSQC peak list
#'
#' Equivalent to `pybmrb.Spectra.create_c13hsqc_peaklist()`.
#' @inheritParams n15hsqc
#' @return Same structure as [create_n15hsqc_peaklist()].
#' @export
create_c13hsqc_peaklist <- function(bmrb_ids        = NULL,
                                    input_file_names = NULL,
                                    auth_tag         = FALSE,
                                    draw_trace       = FALSE) {
  dfs <- .collect_cs_dfs(bmrb_ids, input_file_names, auth_tag)
  if (length(dfs) == 0L) return(NULL)

  x_out <- y_out <- ds_out <- res_out <- info_out <- sid_out <- character(0)

  for (ds_name in names(dfs)) {
    df <- dfs[[ds_name]]
    seq_col   <- .seq_col(df, auth_tag)
    chain_col <- .chain_col(df, auth_tag)

    # Work on only H and C atoms
    h_df <- df[!is.na(df$Atom_type) & df$Atom_type == "H" & !is.na(df$Val), ]
    c_df <- df[!is.na(df$Atom_type) & df$Atom_type == "C" & !is.na(df$Val), ]

    for (chain in unique(h_df[[chain_col]])) {
      h_chain <- h_df[!is.na(h_df[[chain_col]]) & h_df[[chain_col]] == chain, ]
      c_chain <- c_df[!is.na(c_df[[chain_col]]) & c_df[[chain_col]] == chain, ]

      for (sid in sort(unique(h_chain[[seq_col]]))) {
        h_res <- h_chain[!is.na(h_chain[[seq_col]]) & h_chain[[seq_col]] == sid, ]
        c_res <- c_chain[!is.na(c_chain[[seq_col]]) & c_chain[[seq_col]] == sid, ]
        res   <- unique(h_res$Comp_ID)[[1L]]

        for (k in seq_len(nrow(h_res))) {
          h_atom <- h_res$Atom_ID[[k]]
          h_cs   <- h_res$Val[[k]]
          c_atom <- .h_to_c(h_atom)
          if (is.na(c_atom)) next

          c_cs_row <- c_res[!is.na(c_res$Atom_ID) & c_res$Atom_ID == c_atom, ]
          if (nrow(c_cs_row) == 0L) next
          c_cs <- c_cs_row$Val[[1L]]
          if (is.na(h_cs) || is.na(c_cs)) next

          x_out    <- c(x_out, h_cs)
          y_out    <- c(y_out, c_cs)
          ds_out   <- c(ds_out, ds_name)
          res_out  <- c(res_out, res)
          info_out <- c(info_out, sprintf("%s%s-%s %s/%s", res, sid, chain, h_atom, c_atom))
          sid_out  <- c(sid_out, as.character(sid))
        }
      }
    }
  }

  if (length(x_out) == 0L) return(NULL)
  cs_track <- if (draw_trace) .build_cs_track(x_out, y_out, ds_out, sid_out) else list()

  list(x = as.numeric(x_out), y = as.numeric(y_out),
       dataset  = ds_out, residue = res_out,
       info     = info_out, seq_id = sid_out,
       cs_track = cs_track,
       x_atom   = "H", y_atom = "C")
}

#' Build a ¹H–¹H TOCSY peak list
#'
#' Pairs all protons within each residue (full spin-system). Both diagonal and
#' off-diagonal peaks are included.
#' Equivalent to `pybmrb.Spectra.create_tocsy_peaklist()`.
#' @inheritParams n15hsqc
#' @return Same structure as [create_n15hsqc_peaklist()].
#' @export
create_tocsy_peaklist <- function(bmrb_ids        = NULL,
                                  input_file_names = NULL,
                                  auth_tag         = FALSE,
                                  draw_trace       = FALSE) {
  dfs <- .collect_cs_dfs(bmrb_ids, input_file_names, auth_tag)
  if (length(dfs) == 0L) return(NULL)

  x_out <- y_out <- ds_out <- res_out <- info_out <- sid_out <- character(0)

  for (ds_name in names(dfs)) {
    df <- dfs[[ds_name]]
    seq_col   <- .seq_col(df, auth_tag)
    chain_col <- .chain_col(df, auth_tag)

    h_df <- df[!is.na(df$Atom_type) & df$Atom_type == "H" & !is.na(df$Val), ]

    for (chain in unique(h_df[[chain_col]])) {
      h_chain <- h_df[!is.na(h_df[[chain_col]]) & h_df[[chain_col]] == chain, ]

      for (sid in sort(unique(h_chain[[seq_col]]))) {
        h_res <- h_chain[!is.na(h_chain[[seq_col]]) & h_chain[[seq_col]] == sid, ]
        res   <- unique(h_res$Comp_ID)[[1L]]
        nH    <- nrow(h_res)
        if (nH < 1L) next

        # All H-H pairs (including self-correlations on diagonal)
        for (a in seq_len(nH)) {
          for (b in seq_len(nH)) {
            if (a > b) next  # upper triangle only → symmetric spectrum
            hA <- h_res[a, ]
            hB <- h_res[b, ]
            x_out    <- c(x_out, hA$Val, hB$Val)  # symmetric pair
            y_out    <- c(y_out, hB$Val, hA$Val)
            ds_out   <- c(ds_out, ds_name, ds_name)
            res_out  <- c(res_out, res, res)
            info_out <- c(info_out,
              sprintf("%s%s %s-%s", res, sid, hA$Atom_ID, hB$Atom_ID),
              sprintf("%s%s %s-%s", res, sid, hB$Atom_ID, hA$Atom_ID))
            sid_out  <- c(sid_out, as.character(sid), as.character(sid))
          }
        }
      }
    }
  }

  if (length(x_out) == 0L) return(NULL)
  cs_track <- if (draw_trace) .build_cs_track(x_out, y_out, ds_out, sid_out) else list()

  list(x = as.numeric(x_out), y = as.numeric(y_out),
       dataset  = ds_out, residue = res_out,
       info     = info_out, seq_id = sid_out,
       cs_track = cs_track,
       x_atom   = "H", y_atom = "H")
}

#' Build a generic 2D peak list for any two atoms
#'
#' Equivalent to `pybmrb.Spectra.create_2d_peaklist()`.
#' @param atom_x Atom name for the X axis.
#' @param atom_y Atom name for the Y axis.
#' @param include_preceding Include preceding (i−1) residue Y-axis peaks.
#' @param include_next      Include following  (i+1) residue Y-axis peaks.
#' @inheritParams n15hsqc
#' @return Same structure as [create_n15hsqc_peaklist()].
#' @export
create_2d_peaklist <- function(atom_x,
                                atom_y,
                                bmrb_ids         = NULL,
                                input_file_names  = NULL,
                                auth_tag          = FALSE,
                                draw_trace        = FALSE,
                                include_preceding = FALSE,
                                include_next      = FALSE) {

  dfs <- .collect_cs_dfs(bmrb_ids, input_file_names, auth_tag)
  if (length(dfs) == 0L) return(NULL)

  x_out <- y_out <- ds_out <- res_out <- info_out <- sid_out <- character(0)

  for (ds_name in names(dfs)) {
    df <- dfs[[ds_name]]
    seq_col   <- .seq_col(df, auth_tag)
    chain_col <- .chain_col(df, auth_tag)

    for (chain in unique(df[[chain_col]])) {
      cdf      <- df[!is.na(df[[chain_col]]) & df[[chain_col]] == chain & !is.na(df$Val), ]
      seq_ids  <- sort(unique(cdf[[seq_col]]))

      for (i_idx in seq_along(seq_ids)) {
        sid  <- seq_ids[[i_idx]]
        rdf  <- cdf[!is.na(cdf[[seq_col]]) & cdf[[seq_col]] == sid, ]
        res  <- unique(rdf$Comp_ID)[[1L]]
        x_cs <- .get_atom_cs(rdf, atom_x)
        y_cs <- .get_atom_cs(rdf, atom_y)
        if (is.na(x_cs) || is.na(y_cs)) next

        x_out    <- c(x_out, x_cs)
        y_out    <- c(y_out, y_cs)
        ds_out   <- c(ds_out, ds_name)
        res_out  <- c(res_out, res)
        info_out <- c(info_out, sprintf("%s%s-%s %s/%s", res, sid, chain, atom_x, atom_y))
        sid_out  <- c(sid_out, as.character(sid))

        # Preceding residue y-shift (smaller marker to distinguish)
        if (include_preceding && i_idx > 1L) {
          prev_sid <- seq_ids[[i_idx - 1L]]
          p_rdf    <- cdf[!is.na(cdf[[seq_col]]) & cdf[[seq_col]] == prev_sid, ]
          y_prev   <- .get_atom_cs(p_rdf, atom_y)
          if (!is.na(x_cs) && !is.na(y_prev)) {
            x_out    <- c(x_out, x_cs)
            y_out    <- c(y_out, y_prev)
            ds_out   <- c(ds_out, paste0(ds_name, "_i-1"))
            res_out  <- c(res_out, res)
            info_out <- c(info_out, sprintf("%s%s i-1 %s/%s", res, sid, atom_x, atom_y))
            sid_out  <- c(sid_out, as.character(sid))
          }
        }

        # Next residue y-shift
        if (include_next && i_idx < length(seq_ids)) {
          next_sid <- seq_ids[[i_idx + 1L]]
          n_rdf    <- cdf[!is.na(cdf[[seq_col]]) & cdf[[seq_col]] == next_sid, ]
          y_next   <- .get_atom_cs(n_rdf, atom_y)
          if (!is.na(x_cs) && !is.na(y_next)) {
            x_out    <- c(x_out, x_cs)
            y_out    <- c(y_out, y_next)
            ds_out   <- c(ds_out, paste0(ds_name, "_i+1"))
            res_out  <- c(res_out, res)
            info_out <- c(info_out, sprintf("%s%s i+1 %s/%s", res, sid, atom_x, atom_y))
            sid_out  <- c(sid_out, as.character(sid))
          }
        }
      }
    }
  }

  if (length(x_out) == 0L) return(NULL)
  cs_track <- if (draw_trace) .build_cs_track(x_out, y_out, ds_out, sid_out) else list()

  list(x = as.numeric(x_out), y = as.numeric(y_out),
       dataset  = ds_out, residue = res_out,
       info     = info_out, seq_id = sid_out,
       cs_track = cs_track,
       x_atom   = atom_x, y_atom = atom_y)
}

#' Export a peak list to CSV or Sparky format
#'
#' Converts the output of any `create_*_peaklist()` function to a file or
#' a data.frame. Equivalent to `pybmrb.Spectra.export_peak_list()`.
#'
#' @param peak_list       List returned by any `create_*_peaklist()` function.
#' @param output_file_name Path to write the output file. `NULL` returns only
#'   the data.frame without writing to disk.
#' @param output_format   `"csv"` (default) or `"sparky"`.
#' @param include_side_chain Include sidechain peaks. Default `TRUE`.
#'
#' @return Invisibly, a data.frame with the exported peak list.
#' @export
#' @examples
#' \dontrun{
#' pl <- create_n15hsqc_peaklist(bmrb_ids = 15060)
#' export_peak_list(pl, "peaks.csv", output_format = "csv")
#' export_peak_list(pl, "peaks.sparky", output_format = "sparky")
#' }
export_peak_list <- function(peak_list,
                              output_file_name = NULL,
                              output_format    = "csv",
                              include_side_chain = TRUE) {
  if (is.null(peak_list)) cli::cli_abort("`peak_list` is NULL.")
  output_format <- match.arg(output_format, c("csv", "sparky"))

  df <- data.frame(
    Assignment = peak_list$info,
    x_ppm      = peak_list$x,
    y_ppm      = peak_list$y,
    Dataset    = peak_list$dataset,
    Residue    = peak_list$residue,
    Seq_ID     = peak_list$seq_id,
    stringsAsFactors = FALSE
  )
  # Name x/y columns by actual atom type
  colnames(df)[2L] <- paste0(peak_list$x_atom, "_ppm")
  colnames(df)[3L] <- paste0(peak_list$y_atom, "_ppm")

  if (!is.null(output_file_name)) {
    if (output_format == "sparky") {
      sparky_lines <- c(
        "Assignment\tw1\tw2",
        sprintf("%s\t%.3f\t%.3f", df$Assignment, df[[2L]], df[[3L]])
      )
      writeLines(sparky_lines, output_file_name)
    } else {
      utils::write.csv(df, output_file_name, row.names = FALSE)
    }
    cli::cli_alert_info("Peak list written to {output_file_name}")
  }
  invisible(df)
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Collect CS data.frames from BMRB IDs and/or NMR-STAR files
#' @noRd
.collect_cs_dfs <- function(bmrb_ids, input_file_names, auth_tag) {
  all_dfs <- list()

  if (!is.null(bmrb_ids)) {
    ids <- as.character(.ensure_vector(bmrb_ids))
    for (id in ids) {
      cli::cli_progress_message("Fetching BMRB entry {id} ...")
      entry_json <- tryCatch(
        .fetch_entry_json(id),  # uses ?format=zlib like PyBMRB/pynmrstar
        error = function(e) {
          cli::cli_alert_warning(
            "Could not fetch entry {id}: {conditionMessage(e)}")
          NULL
        }
      )
      if (is.null(entry_json)) next

      dfs <- .entry_json_to_cs_dfs(entry_json)
      if (length(dfs) == 0L) {
        cli::cli_alert_warning(
          "Entry {id}: no Atom_chem_shift loop found. ",
          "Run bmrb_debug_entry({id}) to inspect the API response.")
        next
      }
      df <- do.call(rbind, dfs)
      if (nrow(df) == 0L || !"Val" %in% colnames(df)) {
        cli::cli_alert_warning("Entry {id}: loop found but no Val column.")
        next
      }
      cli::cli_alert_info("Entry {id}: {nrow(df)} chemical shifts loaded.")
      all_dfs[[id]] <- df
    }
  }

  if (!is.null(input_file_names)) {
    fps <- as.character(.ensure_vector(input_file_names))
    raw <- .parse_nmrstar_files(fps, auth_tag = auth_tag)
    for (i in seq_along(raw)) {
      if (!is.null(raw[[i]])) {
        nm <- basename(fps[[i]])
        df <- raw[[i]]
        df[["_dataset_id"]] <- NULL
        all_dfs[[nm]] <- df
      }
    }
  }
  all_dfs
}

#' Get chemical shift value for a specific atom from a residue data.frame
#' @noRd
.get_atom_cs <- function(rdf, atom_id) {
  row <- rdf[!is.na(rdf$Atom_ID) & rdf$Atom_ID == atom_id, ]
  if (nrow(row) == 0L || is.na(row$Val[[1L]])) return(NA_real_)
  row$Val[[1L]]
}

#' Build chemical shift tracking dict (used for draw_trace)
#' @noRd
.build_cs_track <- function(x, y, dataset, seq_id) {
  datasets <- unique(dataset)
  if (length(datasets) < 2L) return(list())

  # For each seq_id, collect x,y across datasets
  sids <- unique(seq_id)
  track <- list()
  for (sid in sids) {
    mask <- seq_id == sid
    dsub_x <- x[mask]
    dsub_y <- y[mask]
    dsub_d <- dataset[mask]
    if (length(unique(dsub_d)) > 1L) {
      key <- paste0("seq_", sid)
      track[[key]] <- list(x = dsub_x, y = dsub_y, dataset = dsub_d)
    }
  }
  track
}

#' Build sequence-walk line segments
#' @noRd
.build_walk <- function(pl, full_walk = FALSE) {
  if (length(pl$x) == 0L) return(NULL)
  ds <- unique(pl$dataset)[[1L]]
  mask <- pl$dataset == ds
  sids <- as.numeric(pl$seq_id[mask])
  xs   <- pl$x[mask]
  ys   <- pl$y[mask]
  ord  <- order(sids)
  sids <- sids[ord]; xs <- xs[ord]; ys <- ys[ord]
  segs <- list()
  for (i in seq_along(sids)[-1L]) {
    gap <- sids[[i]] - sids[[i - 1L]]
    if (gap == 1L || full_walk) {
      segs[[length(segs) + 1L]] <- list(
        x0 = xs[[i - 1L]], y0 = ys[[i - 1L]],
        x1 = xs[[i]],      y1 = ys[[i]]
      )
    }
  }
  segs
}

# ---------------------------------------------------------------------------
# Plotting engine
# ---------------------------------------------------------------------------

#' Build an interactive plotly figure for a 2D peak list
#' @noRd
.plot_2d_spectrum <- function(pl, x_label, y_label, title,
                               legend, draw_trace, extra_peak_list = NULL,
                               walk_segments = NULL) {

  datasets <- unique(pl$dataset)
  residues <- unique(pl$residue)

  # Color logic
  if (!is.null(legend) && legend == "residue") {
    pt_colors <- vapply(pl$residue, .res_color, character(1L))
  } else {
    ds_palette <- grDevices::colorRampPalette(
      c("#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD","#8C564B"))(length(datasets))
    ds_color_map <- stats::setNames(ds_palette, datasets)
    pt_colors <- ds_color_map[pl$dataset]
  }

  # Marker symbols by dataset
  sym_map <- stats::setNames(
    .MARKER_SYMBOLS[((seq_along(datasets) - 1L) %% length(.MARKER_SYMBOLS)) + 1L],
    datasets
  )
  pt_syms <- sym_map[pl$dataset]

  # Build hover text
  hover <- sprintf("Residue: %s<br>Seq: %s<br>Dataset: %s<br>x=%.3f<br>y=%.3f",
                   pl$residue, pl$seq_id, pl$dataset, pl$x, pl$y)

  fig <- plotly::plot_ly()

  # One trace per dataset for clean legends
  for (ds in datasets) {
    mask <- pl$dataset == ds
    fig <- plotly::add_trace(fig,
      type   = "scatter", mode = "markers",
      x      = pl$x[mask],
      y      = pl$y[mask],
      name   = ds,
      marker = list(
        color  = pt_colors[mask],
        symbol = sym_map[[ds]],
        size   = 7,
        line   = list(width = 0.5, color = "rgba(0,0,0,0.3)")
      ),
      text      = hover[mask],
      hovertemplate = "%{text}<extra></extra>",
      showlegend = !is.null(legend)
    )
  }

  # Trace lines between matching residues across datasets
  if (draw_trace && length(pl$cs_track) > 0L) {
    for (track in pl$cs_track) {
      if (length(unique(track$dataset)) > 1L) {
        fig <- plotly::add_trace(fig,
          type = "scatter", mode = "lines",
          x    = track$x, y = track$y,
          line = list(color = "rgba(100,100,100,0.4)", width = 1),
          showlegend = FALSE, hoverinfo = "skip"
        )
      }
    }
  }

  # Sequence walk
  if (!is.null(walk_segments) && length(walk_segments) > 0L) {
    for (seg in walk_segments) {
      fig <- plotly::add_trace(fig,
        type = "scatter", mode = "lines",
        x    = c(seg$x0, seg$x1, NA),
        y    = c(seg$y0, seg$y1, NA),
        line = list(color = "rgba(50,50,200,0.5)", width = 1.5, dash = "dot"),
        showlegend = FALSE, hoverinfo = "skip"
      )
    }
  }

  # Overlay external peak list (CSV)
  if (!is.null(extra_peak_list) && file.exists(extra_peak_list)) {
    ext <- tryCatch(utils::read.csv(extra_peak_list, header = FALSE,
                                     col.names = c("x", "y")),
                    error = function(e) NULL)
    if (!is.null(ext) && nrow(ext) > 0L) {
      fig <- plotly::add_trace(fig,
        type = "scatter", mode = "markers",
        x = ext$x, y = ext$y,
        name = basename(extra_peak_list),
        marker = list(color = "black", symbol = "cross", size = 8),
        showlegend = !is.null(legend)
      )
    }
  }

  # Layout: NMR convention — both axes reversed
  fig <- plotly::layout(fig,
    title  = list(text = title, font = list(size = 16)),
    xaxis  = list(title = x_label, autorange = "reversed"),
    yaxis  = list(title = y_label, autorange = "reversed"),
    legend = list(title = list(text = if (!is.null(legend)) legend else "")),
    plot_bgcolor  = "white",
    paper_bgcolor = "white",
    hovermode = "closest"
  )
  fig
}

#' Save or display a plotly figure
#' @noRd
.output_figure <- function(fig, output_format, output_file,
                            width, height, show_visualization) {
  output_format <- match.arg(output_format,
    c("html", "png", "jpg", "jpeg", "pdf", "svg", "webp"))

  if (!is.null(output_file)) {
    if (output_format == "html") {
      if (!requireNamespace("htmlwidgets", quietly = TRUE))
        cli::cli_abort(
          "Package 'htmlwidgets' is needed to save HTML. ",
          "Install it with: install.packages('htmlwidgets')")
      htmlwidgets::saveWidget(fig, output_file, selfcontained = TRUE)
    } else {
      plotly::save_image(fig, output_file, width = width, height = height)
    }
    cli::cli_alert_info("Saved: {output_file}")
  }

  if (show_visualization) print(fig)
  invisible(fig)
}
