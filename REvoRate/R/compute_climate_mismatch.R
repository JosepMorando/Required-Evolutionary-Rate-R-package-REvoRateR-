#' Compute plain (unweighted) climate mismatch
#'
#' @description
#' Computes the pixel-wise Euclidean distance between future and current
#' climate in raw BIO space:
#' \deqn{CM = \sqrt{\sum_i \Delta BIO_i^2}}
#' Optionally standardizes each layer difference by the SD of the current
#' layer to produce a unitless index.
#'
#' @param current_stack `RasterStack`/`RasterBrick` of present-day climate.
#' @param future_stack  `RasterStack`/`RasterBrick` of future climate for
#'   a single period.
#' @param standardize Logical (default `TRUE`). Divide each layer difference
#'   by the SD of the corresponding current layer before computing the norm.
#'
#' @return A `RasterLayer` of climate mismatch values.
#'
#' @examples
#' \dontrun{
#' cm <- compute_climate_mismatch_plain(current_stack, future_p1, standardize = TRUE)
#' }
#'
#' @export
compute_climate_mismatch_plain <- function(current_stack, future_stack,
                                           standardize = TRUE) {

  stopifnot(raster::nlayers(current_stack) == raster::nlayers(future_stack))

  # ensure consistent layer order
  future_stack <- future_stack[[names(current_stack)]]
  diff_stack   <- future_stack - current_stack

  if (standardize) {
    sds <- sapply(seq_len(raster::nlayers(current_stack)), function(i)
      raster::cellStats(current_stack[[i]], stat = "sd", na.rm = TRUE))
    sds[sds == 0 | !is.finite(sds)] <- 1
    for (i in seq_along(sds))
      diff_stack[[i]] <- diff_stack[[i]] / sds[i]
  }

  raster::calc(diff_stack, fun = function(v) sqrt(sum(v^2, na.rm = TRUE)))
}


#' Compute GF-importance-weighted climate mismatch
#'
#' @description
#' As [compute_climate_mismatch_plain()], but each layer's squared difference
#' is weighted by the relative variable importance from the GF model before
#' summing:
#' \deqn{CM_{GF} = \sqrt{\sum_i w_i \cdot \Delta BIO_i^2}}
#' where \eqn{w_i = imp_i / \sum imp_i}.
#'
#' @param current_stack `RasterStack`/`RasterBrick` of present-day climate.
#' @param future_stack  `RasterStack`/`RasterBrick` of future climate for
#'   a single period.
#' @param gf_model A fitted `gradientForest` model. Variable importances are
#'   extracted via `importance(gf_model)`.
#' @param standardize Logical (default `TRUE`). Standardize each layer
#'   difference by the SD of the current layer before weighting.
#'
#' @return A `RasterLayer` of GF-weighted climate mismatch values.
#'
#' @examples
#' \dontrun{
#' cm_gf <- compute_climate_mismatch_gf(current_stack, future_p1, gfRef)
#' }
#'
#' @export
compute_climate_mismatch_gf <- function(current_stack, future_stack,
                                        gf_model, standardize = TRUE) {

  stopifnot(inherits(gf_model, "gradientForest"))
  stopifnot(raster::nlayers(current_stack) == raster::nlayers(future_stack))

  # normalize importances to sum = 1
  imp <- gradientForest::importance(gf_model)
  imp <- imp / sum(imp)

  # restrict to layers present in both stack and GF
  common_vars  <- intersect(names(current_stack), names(imp))
  imp          <- imp[common_vars]
  current_tmp  <- current_stack[[common_vars]]
  future_tmp   <- future_stack[[common_vars]]

  diff_stack   <- future_tmp - current_tmp

  if (standardize) {
    sds <- sapply(seq_len(raster::nlayers(current_tmp)), function(i)
      raster::cellStats(current_tmp[[i]], stat = "sd", na.rm = TRUE))
    sds[sds == 0 | !is.finite(sds)] <- 1
    for (i in seq_along(sds))
      diff_stack[[i]] <- diff_stack[[i]] / sds[i]
  }

  weights <- as.numeric(imp)
  raster::calc(diff_stack, fun = function(v)
    sqrt(sum((v^2) * weights, na.rm = TRUE)))
}
