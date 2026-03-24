#' Extract raster metrics at site coordinates
#'
#' @description
#' Given the per-period rasters produced by [REvoRate()], extracts all metric
#' values at a set of sample site coordinates and returns both a wide-format
#' and a long-format (tidy) data frame—matching the output structure of the
#' original analysis script.
#'
#' @param coords A `data.frame` with columns `Long` and `Lat` (decimal
#'   degrees). Row names are used as site IDs; alternatively supply a `Site`
#'   column.
#' @param go_list Named list of genetic offset `RasterLayer` objects (output
#'   of repeated [compute_genetic_offset()] calls).
#' @param rr_list Named list of REvoRate `RasterLayer` objects (output of
#'   [compute_revoRate()]).
#' @param cm_plain Named list of plain climate mismatch `RasterLayer` objects,
#'   or `NULL`.
#' @param cm_gf Named list of GF-weighted climate mismatch `RasterLayer`
#'   objects, or `NULL`.
#' @param period_midpoints Named numeric vector mapping period names to
#'   mid-point years (e.g. `c(p1 = 2030, p2 = 2050, p3 = 2070, p4 = 2090)`).
#' @param raster_crs A CRS object used to create the `SpatialPointsDataFrame`.
#'   Typically `raster::crs(rr_list[[1]])`.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{`wide`}{One row per site; columns for `Site`, `Long`, `Lat`,
#'     then one column per metric × period combination
#'     (e.g. `RevoRate_p1`, `GenOffset_p3`, `ClimMismatch_p2`, ...).}
#'   \item{`long`}{Tidy data frame with columns `Site`, `Long`, `Lat`,
#'     `Metric`, `Period`, `Value`.}
#' }
#'
#' @examples
#' \dontrun{
#' extracted <- extract_site_metrics(
#'   coords           = coord,
#'   go_list          = result$genetic_offset,
#'   rr_list          = result$revoRate,
#'   cm_plain         = result$climate_mismatch,
#'   cm_gf            = result$climate_mismatch_gf,
#'   period_midpoints = c(p1 = 2030, p2 = 2050, p3 = 2070, p4 = 2090),
#'   raster_crs       = raster::crs(result$revoRate[[1]])
#' )
#' head(extracted$wide)
#' head(extracted$long)
#' }
#'
#' @export
extract_site_metrics <- function(coords,
                                 go_list,
                                 rr_list,
                                 cm_plain         = NULL,
                                 cm_gf            = NULL,
                                 period_midpoints,
                                 raster_crs) {

  # ---- site IDs ----
  if ("Site" %in% colnames(coords)) {
    site_ids <- coords$Site
  } else {
    site_ids <- rownames(coords)
    if (is.null(site_ids) || all(site_ids == seq_len(nrow(coords))))
      site_ids <- paste0("Site", seq_len(nrow(coords)))
  }

  # ---- build SpatialPointsDataFrame ----
  pts <- sp::SpatialPointsDataFrame(
    coords      = coords[, c("Long", "Lat")],
    data        = data.frame(Site = site_ids, stringsAsFactors = FALSE),
    proj4string = raster_crs
  )

  # ---- helper: extract a named list of rasters at pts ----
  .extract_named <- function(lst, prefix) {
    if (is.null(lst)) return(NULL)
    mat <- sapply(lst, function(r) raster::extract(r, pts))
    mat <- as.matrix(mat)
    colnames(mat) <- paste0(prefix, "_", names(lst))
    mat
  }

  rr_mat <- .extract_named(rr_list,   "RevoRate")
  go_mat <- .extract_named(go_list,   "GenOffset")
  cm_mat <- .extract_named(cm_plain,  "ClimMismatch")
  cg_mat <- .extract_named(cm_gf,     "ClimMismatchGF")

  # ---- assemble wide table ----
  wide <- data.frame(
    Site = site_ids,
    Long = coords$Long,
    Lat  = coords$Lat,
    stringsAsFactors = FALSE
  )

  for (m in list(rr_mat, go_mat, cm_mat, cg_mat)) {
    if (!is.null(m)) wide <- cbind(wide, m)
  }

  # ---- reshape to long ----
  metric_cols <- setdiff(names(wide), c("Site", "Long", "Lat"))

  long <- tidyr::pivot_longer(
    wide,
    cols      = dplyr::all_of(metric_cols),
    names_to  = c("Metric", "Period"),
    names_pattern = "^(RevoRate|GenOffset|ClimMismatch|ClimMismatchGF)_(.+)$",
    values_to = "Value"
  )

  # ---- decode period labels ----
  if (!is.null(period_midpoints)) {
    period_labels <- .make_period_labels(period_midpoints)
    long$Period   <- ifelse(long$Period %in% names(period_labels),
                            period_labels[long$Period], long$Period)
  }

  list(wide = wide, long = as.data.frame(long))
}


# ---- internal: build "2021-2040" style labels from midpoints ----
#' @keywords internal
.make_period_labels <- function(midpoints) {
  # Assumes 20-year periods centred on midpoint (±10 years)
  # Fallback: just use the midpoint year as string
  labels <- vapply(midpoints, function(m) {
    start <- m - 9   # e.g. 2030 -> 2021
    end   <- m + 10  # e.g. 2030 -> 2040
    # detect the conventional WorldClim 20-yr windows
    if (m == 2030) return("2021-2040")
    if (m == 2050) return("2041-2060")
    if (m == 2070) return("2061-2080")
    if (m == 2090) return("2081-2100")
    paste0(start, "-", end)
  }, character(1))
  setNames(labels, names(midpoints))
}
