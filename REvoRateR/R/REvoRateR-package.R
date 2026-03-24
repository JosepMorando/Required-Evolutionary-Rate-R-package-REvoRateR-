#' REvoRateR: Required Evolutionary Rate from Gradient Forest
#'
#' @description
#' Computes genetic offset and the Required Evolutionary Rate (REvoRate) from a
#' fitted Gradient Forest model and a set of climate raster stacks representing
#' successive future time periods.
#'
#' ## Main workflow
#'
#' The core function is [REvoRate()], which wraps the full computation:
#'
#' 1. Align climate stacks to GF predictor names.
#' 2. Compute **genetic offset** per period via [compute_genetic_offset()].
#' 3. Compute **REvoRate** via [compute_revoRate()].
#' 4. Optionally compute plain and GF-weighted **climate mismatch** via
#'    [compute_climate_mismatch_plain()] and [compute_climate_mismatch_gf()].
#' 5. Optionally extract all raster metrics at site coordinates via
#'    [extract_site_metrics()], returning tidy wide and long data frames.
#' 6. Plot results with [plot_revoRate()].
#'
#' ## Citation
#' If you use this package, please cite the underlying paper:
#'
#' Morando-Milà J, et al. (2025) Re-interpreting genetic offset: quantifying
#' the least required evolutionary rate under climate change at the
#' Mediterranean range margin of European beech. *Evolutionary Applications*.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom raster stack raster nlayers ncell extract crop mask crs calc
#'   cellStats resample extent compareRaster
#' @importFrom sp SpatialPointsDataFrame
#' @importFrom dplyr all_of mutate
#' @importFrom tidyr pivot_longer
#' @importFrom methods is
#' @importFrom stats na.omit
## usethis namespace: end
NULL
