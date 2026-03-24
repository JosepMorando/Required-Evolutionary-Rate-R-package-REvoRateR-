#' Plot REvoRate or genetic offset site trajectories
#'
#' @description
#' Produces a faceted ggplot2 figure (one panel per site) showing genetic
#' offset trajectories (line + points) overlaid on REvoRate heat rectangles,
#' replicating the publication figure from the original analysis script.
#'
#' Sites are ordered by their genetic offset at the last available period
#' (highest offset first).
#'
#' @param revoRate_obj A `REvoRate` object returned by [REvoRate()], or
#'   alternatively a long-format data frame as returned in `$site_long`.
#' @param ncol Integer. Number of columns in the facet grid (default `6`).
#' @param palette Character. Viridis palette for the REvoRate fill scale
#'   (default `"cividis"`).
#' @param save_path Character or `NULL`. If not `NULL`, save the plot to this
#'   file path. File format is inferred from the extension (`.pdf`, `.png`,
#'   `.tiff`, etc.).
#' @param width,height Plot dimensions in inches when saving (defaults
#'   `12` × `9`).
#' @param dpi Resolution in DPI for raster formats (default `600`).
#'
#' @return A `ggplot` object (invisibly when `save_path` is supplied).
#'
#' @examples
#' \dontrun{
#' p <- plot_revoRate(result)
#' print(p)
#'
#' # save directly
#' plot_revoRate(result, save_path = "Fig3_REvoRate.pdf")
#' }
#'
#' @export
plot_revoRate <- function(revoRate_obj,
                          ncol       = 6,
                          palette    = "cividis",
                          save_path  = NULL,
                          width      = 12,
                          height     = 9,
                          dpi        = 600) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plot_revoRate(). ",
         "Install it with: install.packages('ggplot2')")

  # Accept REvoRate object or bare long data frame
  if (inherits(revoRate_obj, "REvoRate")) {
    df <- revoRate_obj$site_long
    if (is.null(df))
      stop("No site coordinates were supplied to REvoRate(); ",
           "cannot plot without `$site_long`.")
  } else if (is.data.frame(revoRate_obj)) {
    df <- revoRate_obj
  } else {
    stop("`revoRate_obj` must be a REvoRate object or a long-format data frame.")
  }

  # ---- numeric period midpoints (for x-axis) ----
  period_mid <- c(
    "2021-2040" = 2030, "2021\u20132040" = 2030,
    "2041-2060" = 2050, "2041\u20132060" = 2050,
    "2061-2080" = 2070, "2061\u20132080" = 2070,
    "2081-2100" = 2090, "2081\u20132100" = 2090
  )

  df$YearMid <- as.numeric(period_mid[as.character(df$Period)])
  df <- df[!is.na(df$YearMid), ]

  # ---- split into GO and RR subsets ----
  go_dat <- df[df$Metric == "GenOffset", c("Site", "YearMid", "Value")]
  names(go_dat)[3] <- "GO"

  rr_rect <- df[df$Metric == "RevoRate", ]
  rr_rect$xmin <- rr_rect$YearMid - 10
  rr_rect$xmax <- rr_rect$YearMid + 10
  rr_rect$RR   <- rr_rect$Value

  # ---- order sites by last GO ----
  site_levels <- go_dat[go_dat$YearMid == max(go_dat$YearMid, na.rm = TRUE), ]
  site_levels <- site_levels[order(-site_levels$GO), "Site"]

  go_dat$Site  <- factor(go_dat$Site,  levels = site_levels)
  rr_rect$Site <- factor(rr_rect$Site, levels = site_levels)

  # ---- build plot ----
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = rr_rect,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                   ymin = -Inf,        ymax = Inf,
                   fill = .data$RR),
      alpha = 0.85, inherit.aes = FALSE
    ) +
    ggplot2::geom_line(
      data = go_dat,
      ggplot2::aes(.data$YearMid, .data$GO),
      color = "white", linewidth = 1.0, lineend = "round"
    ) +
    ggplot2::geom_point(
      data = go_dat,
      ggplot2::aes(.data$YearMid, .data$GO),
      color = "black", size = 1.8, stroke = 0.2
    ) +
    ggplot2::facet_wrap(~ Site, ncol = ncol) +
    ggplot2::scale_x_continuous(
      breaks = c(2030, 2050, 2070, 2090),
      limits = c(2020, 2100),
      expand = ggplot2::expansion(add = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(name = "Genetic Offset") +
    ggplot2::scale_fill_viridis_c(
      option = palette, name = "REvoRate",
      guide  = ggplot2::guide_colorbar(barheight = ggplot2::unit(60, "pt"))
    ) +
    ggplot2::labs(x = "Year (period midpoint)", y = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      strip.text         = ggplot2::element_text(
        face = "bold", margin = ggplot2::margin(b = 2)
      ),
      axis.title.y       = ggplot2::element_text(
        margin = ggplot2::margin(r = 6)
      ),
      legend.position    = "right",
      legend.title       = ggplot2::element_text(size = 9),
      legend.text        = ggplot2::element_text(size = 8)
    )

  if (!is.null(save_path)) {
    ggplot2::ggsave(
      filename = save_path,
      plot     = p,
      width    = width,
      height   = height,
      units    = "in",
      dpi      = dpi
    )
    message("Plot saved to: ", save_path)
    return(invisible(p))
  }

  p
}
