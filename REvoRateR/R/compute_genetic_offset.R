#' Compute genetic offset between current and future climate
#'
#' @description
#' Computes the pixel-wise genetic offset as the Euclidean distance in
#' GF-transformed environmental space between present-day and future climate.
#'
#' \deqn{GO = \sqrt{\sum_j \left(GF(x_j^{future}) - GF(x_j^{current})\right)^2}}
#'
#' @param gf_model A fitted `gradientForest` model.
#' @param current_stack A `RasterStack`/`RasterBrick` of present-day climate
#'   variables. Layer names must match predictors in `gf_model`.
#' @param future_stack A `RasterStack`/`RasterBrick` of future climate
#'   variables for a single period.
#'
#' @return A `RasterLayer` of genetic offset values.
#'
#' @examples
#' \dontrun{
#' go_p1 <- compute_genetic_offset(gfRef, current_stack, future_p1)
#' plot(go_p1, main = "Genetic Offset 2021-2040")
#' }
#'
#' @export
compute_genetic_offset <- function(gf_model, current_stack, future_stack) {

  stopifnot(inherits(gf_model, "gradientForest"))
  stopifnot(inherits(current_stack, c("RasterStack", "RasterBrick")))
  stopifnot(inherits(future_stack,  c("RasterStack", "RasterBrick")))

  # --- harmonize NA cells: any NA in future layer 1 propagates to current ---
  clim_tmp  <- current_stack
  na_idx    <- which(is.na(future_stack[[1]][]))
  if (length(na_idx) > 0) {
    for (i in seq_len(raster::nlayers(clim_tmp)))
      clim_tmp[[i]][na_idx] <- NA
  }

  # --- extract env tables ---
  env_cur <- .extract_env(clim_tmp)
  env_fut <- .extract_env(future_stack)

  # --- ensure same cells ---
  env_cur <- .harmonize_cells(env_cur, env_fut)$a
  env_fut <- .harmonize_cells(env_cur, env_fut)$b

  # --- GF predictions ---
  pred_cur <- predict(gf_model, env_cur[, -1, drop = FALSE])
  pred_fut <- predict(gf_model, env_fut[, -1, drop = FALSE])

  # --- Euclidean distance in GF space ---
  go_vals <- sqrt(rowSums((pred_fut - pred_cur)^2))

  # --- map back to raster ---
  out_rast        <- raster::raster(future_stack[[1]])
  out_rast[]      <- NA
  out_rast[env_fut$cell] <- go_vals
  out_rast
}


# ---- helpers shared across functions ----

#' @keywords internal
.extract_env <- function(stack) {
  vals <- raster::extract(stack, seq_len(raster::ncell(stack[[1]])))
  df   <- data.frame(cell = seq_len(raster::ncell(stack[[1]])), vals,
                     stringsAsFactors = FALSE)
  df   <- stats::na.omit(df)
  df
}

#' @keywords internal
.harmonize_cells <- function(a, b) {
  if (identical(a$cell, b$cell)) return(list(a = a, b = b))
  common <- intersect(a$cell, b$cell)
  list(
    a = a[match(common, a$cell), , drop = FALSE],
    b = b[match(common, b$cell), , drop = FALSE]
  )
}
