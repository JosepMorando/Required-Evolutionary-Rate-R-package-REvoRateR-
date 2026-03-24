#' Compute Genetic Offset and Required Evolutionary Rate (REvoRate)
#'
#' @description
#' Given a fitted Gradient Forest model, a raster stack of current climate,
#' a list of future climate raster stacks (one per period), and the midpoint
#' years for each period, `REvoRate()` computes:
#'
#' * **Genetic offset** per pixel and per period (Euclidean distance in
#'   GF-transformed environmental space between present and future).
#' * **REvoRate** per pixel and per period (genetic offset divided by the
#'   number of years elapsed since the previous period).
#' * Optionally, plain and GF-weighted **climate mismatch** (Euclidean
#'   distance in raw / GF-importance-weighted BIO space).
#'
#' If a data frame of site coordinates is supplied, all metrics are also
#' extracted at those points and returned as tidy wide and long data frames.
#'
#' @param gf_model A fitted `gradientForest` model object (from the
#'   `gradientForest` package).
#' @param current_stack A `RasterStack` or `RasterBrick` of present-day
#'   climate variables. Layer names must match the predictor names used to
#'   fit `gf_model`.
#' @param future_stacks A **named list** of `RasterStack` / `RasterBrick`
#'   objects, one per future time period (e.g.
#'   `list(p1 = stack1, p2 = stack2, ...)`). Each stack must have the same
#'   layer names and spatial extent as `current_stack`. The list is processed
#'   in order; names are carried through to the output.
#' @param period_midpoints A **named numeric vector** of mid-point years
#'   corresponding to each entry in `future_stacks` (e.g.
#'   `c(p1 = 2030, p2 = 2050, p3 = 2070, p4 = 2090)`). Names must match
#'   those of `future_stacks`.
#' @param current_midpoint Numeric. The mid-point year of the present-day
#'   climate period (default `1985`, i.e. 1970–2000).  Used only to compute
#'   the time-delta for the **first** future period.
#' @param coords An optional `data.frame` with columns `Long` and `Lat`
#'   (decimal degrees) and row names used as site IDs. When supplied, all
#'   raster metrics are extracted at these locations and wide/long summary
#'   tables are returned. Default `NULL`.
#' @param compute_climate_mismatch Logical. If `TRUE` (default), also compute
#'   plain and GF-weighted climate mismatch rasters.
#' @param standardize_mismatch Logical. If `TRUE` (default), standardize each
#'   BIO layer difference by the SD of the current-climate layer before
#'   computing the Euclidean norm (produces a unitless index).
#' @param verbose Logical. If `TRUE` (default), print progress messages.
#'
#' @return A named list with the following elements:
#'
#' \describe{
#'   \item{`genetic_offset`}{Named list of `RasterLayer` objects—one per
#'     period—containing pixel-wise genetic offset.}
#'   \item{`revoRate`}{Named list of `RasterLayer` objects—one per period—
#'     containing pixel-wise REvoRate (genetic offset / Δt).}
#'   \item{`climate_mismatch`}{Named list of plain (unweighted) climate
#'     mismatch `RasterLayer` objects. `NULL` if
#'     `compute_climate_mismatch = FALSE`.}
#'   \item{`climate_mismatch_gf`}{Named list of GF-importance-weighted
#'     climate mismatch `RasterLayer` objects. `NULL` if
#'     `compute_climate_mismatch = FALSE`.}
#'   \item{`delta_t`}{Named numeric vector of year-deltas used for each
#'     period.}
#'   \item{`site_wide`}{Wide-format `data.frame` with one row per site and
#'     one column per metric × period combination. `NULL` if `coords` is
#'     not supplied.}
#'   \item{`site_long`}{Long-format (tidy) `data.frame` with columns
#'     `Site`, `Long`, `Lat`, `Period`, `Metric`, `Value`. `NULL` if
#'     `coords` is not supplied.}
#' }
#'
#' @details
#' ## Genetic offset
#'
#' Computed as the Euclidean distance in the GF-transformed variable space
#' between the present-day and future-climate predictions at each pixel
#' (Fitzpatrick & Keller 2015; Ellis et al. 2012):
#'
#' \deqn{GO = \sqrt{\sum_j \left(GF(x_j^{future}) - GF(x_j^{current})\right)^2}}
#'
#' ## REvoRate — formal derivation (Morando-Milà et al. 2025, Eqs. 2–5)
#'
#' Let \eqn{r_G = dG/dE} be the rate of genomic change per unit of
#' environmental change, and \eqn{v_E = dE/dt} the velocity of environmental
#' change over time.  The Required Evolutionary Rate is their product:
#'
#' \deqn{REvoRate = \frac{dG}{dt} = \frac{dG}{dE} \cdot \frac{dE}{dt}}
#'
#' This function computes the **discrete baseline** version using
#' bidecadal climate-period midpoints.  For the first period the present-day
#' offset is zero, so:
#'
#' \deqn{RR_{t_1} = \frac{GO_{t_1}}{\Delta t_1}, \quad
#'       \Delta t_1 = midpoint_{t_1} - midpoint_{current}}
#'
#' For subsequent periods, only the incremental change is attributed to that
#' interval (minimum baseline rate):
#'
#' \deqn{RR_{t_k} = \frac{GO_{t_k} - GO_{t_{k-1}}}{\Delta t_k}, \quad k \geq 2}
#'
#' ## Joint interpretation of offset and REvoRate
#'
#' The two metrics define four qualitative risk regimes
#' (Morando-Milà et al. 2025):
#'
#' * **High offset + High REvoRate:** worst case — high adaptive stress; large
#'   mismatch AND rapid change required.
#' * **High offset + Low REvoRate:** chronic lag — large accumulated mismatch,
#'   slow climate velocity.
#' * **Low offset + High REvoRate:** transient stress — rapid environmental
#'   change but small current mismatch.
#' * **Low offset + Low REvoRate:** near equilibrium — population tracking the
#'   optimum.
#'
#' Note that in period 1 offset and REvoRate are linearly related (both
#' measured from the same zero baseline), so differences in ranking only emerge
#' from period 2 onwards.
#'
#' @references
#' Morando-Milà J, Grau O, Ulaszewski B, Vilà-Cabrera A, Peñuelas J, Jump A,
#' Scotti I (2025). Re-interpreting genetic offset: quantifying the least
#' required evolutionary rate under climate change at the Mediterranean range
#' margin of European beech. *Evolutionary Applications*.
#'
#' Fitzpatrick MC, Keller SR (2015). Ecological genomics meets community-level
#' modelling of biodiversity: mapping the genomic landscape of current and
#' future environmental adaptation. *Ecology Letters* 18(1):1–16.
#'
#' Ellis N, Smith SJ, Pitcher CR (2012). Gradient forests: calculating
#' importance gradients on physical predictors. *Ecology* 93(1):156–168.
#'
#' @seealso
#' [compute_genetic_offset()], [compute_revoRate()],
#' [compute_climate_mismatch_plain()], [compute_climate_mismatch_gf()],
#' [extract_site_metrics()], [plot_revoRate()]
#'
#' @examples
#' \dontrun{
#' library(REvoRate)
#' library(gradientForest)
#' library(raster)
#'
#' # Assuming you have:
#' #   gfRef         - a fitted gradientForest model
#' #   current_stack - RasterStack of present-day BIOs (named to match gfRef)
#' #   p1..p4        - RasterStacks of future BIOs for four 20-yr periods
#' #   coord         - data.frame with Long/Lat columns, rownames = site IDs
#'
#' result <- REvoRate(
#'   gf_model          = gfRef,
#'   current_stack     = current_stack,
#'   future_stacks     = list(p1 = p1, p2 = p2, p3 = p3, p4 = p4),
#'   period_midpoints  = c(p1 = 2030, p2 = 2050, p3 = 2070, p4 = 2090),
#'   current_midpoint  = 1985,
#'   coords            = coord,
#'   compute_climate_mismatch = TRUE,
#'   standardize_mismatch     = TRUE
#' )
#'
#' # Wide table ready for export
#' write.csv(result$site_wide, "REvoRate_results_wide.csv", row.names = FALSE)
#' write.csv(result$site_long, "REvoRate_results_long.csv", row.names = FALSE)
#'
#' # Access individual rasters
#' plot(result$revoRate$p1)
#' plot(result$genetic_offset$p3)
#' }
#'
#' @export
REvoRate <- function(gf_model,
                     current_stack,
                     future_stacks,
                     period_midpoints,
                     current_midpoint      = 1985,
                     coords                = NULL,
                     compute_climate_mismatch = TRUE,
                     standardize_mismatch  = TRUE,
                     verbose               = TRUE) {

  # ---- Input validation ----
  .validate_inputs(gf_model, current_stack, future_stacks,
                   period_midpoints, current_midpoint, coords)

  periods <- names(future_stacks)

  # ---- 1. Align current stack layers to GF predictors ----
  gf_vars       <- colnames(gf_model$X)   # predictor names used by GF
  current_stack <- .align_stack(current_stack, gf_vars, "current_stack")

  future_stacks <- lapply(future_stacks, .align_stack, gf_vars = gf_vars,
                          stack_name = "future_stack")

  # ---- 2. Compute genetic offset per period ----
  if (verbose) message("Computing genetic offset...")
  go_list <- vector("list", length(periods))
  names(go_list) <- periods

  for (nm in periods) {
    if (verbose) message("  Period: ", nm)
    go_list[[nm]] <- compute_genetic_offset(
      gf_model      = gf_model,
      current_stack = current_stack,
      future_stack  = future_stacks[[nm]]
    )
  }

  # ---- 3. Compute delta_t per period ----
  midpoints_full <- c(current = current_midpoint, period_midpoints)
  delta_t        <- .compute_delta_t(periods, midpoints_full)

  # ---- 4. Compute REvoRate per period ----
  if (verbose) message("Computing REvoRate...")
  rr_list <- compute_revoRate(go_list, delta_t)

  # ---- 5. Climate mismatch (optional) ----
  cm_plain <- NULL
  cm_gf    <- NULL

  if (compute_climate_mismatch) {
    if (verbose) message("Computing climate mismatch...")
    cm_plain <- lapply(future_stacks, compute_climate_mismatch_plain,
                       current_stack = current_stack,
                       standardize   = standardize_mismatch)
    names(cm_plain) <- periods

    cm_gf <- lapply(future_stacks, compute_climate_mismatch_gf,
                    current_stack = current_stack,
                    gf_model      = gf_model,
                    standardize   = standardize_mismatch)
    names(cm_gf) <- periods
  }

  # ---- 6. Site extraction (optional) ----
  site_wide <- NULL
  site_long <- NULL

  if (!is.null(coords)) {
    if (verbose) message("Extracting metrics at site coordinates...")
    extracted <- extract_site_metrics(
      coords          = coords,
      go_list         = go_list,
      rr_list         = rr_list,
      cm_plain        = cm_plain,
      cm_gf           = cm_gf,
      period_midpoints = period_midpoints,
      raster_crs      = raster::crs(rr_list[[1]])
    )
    site_wide <- extracted$wide
    site_long <- extracted$long
  }

  # ---- Return ----
  structure(
    list(
      genetic_offset       = go_list,
      revoRate             = rr_list,
      climate_mismatch     = cm_plain,
      climate_mismatch_gf  = cm_gf,
      delta_t              = delta_t,
      site_wide            = site_wide,
      site_long            = site_long
    ),
    class = "REvoRate"
  )
}


# ---- print method ----

#' @export
print.REvoRate <- function(x, ...) {
  periods <- names(x$revoRate)
  cat("REvoRate object\n")
  cat("  Periods     :", paste(periods, collapse = ", "), "\n")
  cat("  Delta-t     :", paste(names(x$delta_t), x$delta_t, sep = "=",
                               collapse = ", "), "years\n")
  if (!is.null(x$site_wide)) {
    cat("  Sites       :", nrow(x$site_wide), "\n")
    cat("  Columns     :", paste(names(x$site_wide), collapse = ", "), "\n")
  } else {
    cat("  Sites       : (no coord table supplied)\n")
  }
  invisible(x)
}


# ---- Internal helpers ----

#' @keywords internal
.validate_inputs <- function(gf_model, current_stack, future_stacks,
                              period_midpoints, current_midpoint, coords) {
  if (!inherits(gf_model, "gradientForest"))
    stop("`gf_model` must be a fitted gradientForest object.")
  if (!inherits(current_stack, c("RasterStack", "RasterBrick")))
    stop("`current_stack` must be a RasterStack or RasterBrick.")
  if (!is.list(future_stacks) || is.null(names(future_stacks)))
    stop("`future_stacks` must be a *named* list of RasterStack/RasterBrick objects.")
  if (!is.numeric(period_midpoints) || is.null(names(period_midpoints)))
    stop("`period_midpoints` must be a named numeric vector.")
  if (!all(names(future_stacks) %in% names(period_midpoints)))
    stop("All names in `future_stacks` must appear in `period_midpoints`.")
  if (!is.numeric(current_midpoint) || length(current_midpoint) != 1)
    stop("`current_midpoint` must be a single numeric value.")
  if (!is.null(coords)) {
    if (!is.data.frame(coords))
      stop("`coords` must be a data.frame.")
    if (!all(c("Long", "Lat") %in% colnames(coords)))
      stop("`coords` must have columns named 'Long' and 'Lat'.")
  }
}

#' @keywords internal
.align_stack <- function(stack, gf_vars, stack_name = "stack") {
  available <- names(stack)
  keep       <- intersect(available, gf_vars)
  if (length(keep) == 0)
    stop("No layer names in `", stack_name, "` match the GF predictor names. ",
         "GF predictors: ", paste(gf_vars, collapse = ", "))
  missing <- setdiff(gf_vars, available)
  if (length(missing) > 0)
    warning("These GF predictors are absent from `", stack_name, "` and will ",
            "be ignored: ", paste(missing, collapse = ", "))
  stack[[keep]]
}

#' @keywords internal
.compute_delta_t <- function(periods, midpoints_full) {
  # midpoints_full has names c("current", periods...)
  dt <- numeric(length(periods))
  names(dt) <- periods
  for (i in seq_along(periods)) {
    nm_prev <- if (i == 1) "current" else periods[i - 1]
    dt[periods[i]] <- midpoints_full[periods[i]] - midpoints_full[nm_prev]
  }
  if (any(dt <= 0))
    stop("All period midpoints must be strictly increasing. ",
         "Check `period_midpoints` and `current_midpoint`.")
  dt
}
