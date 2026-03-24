#' Compute Required Evolutionary Rate (REvoRate) from genetic offset rasters
#'
#' @description
#' Converts a list of per-period genetic offset rasters into REvoRate rasters.
#' REvoRate was introduced in Morando-Milà et al. (2025) as a temporal
#' extension of genetic offset that captures the *pace* at which populations
#' must evolve to keep track of shifting climate optima.
#'
#' ## Formal definition (Morando-Milà et al. 2025, Eqs. 2–5)
#'
#' Let \eqn{\Delta G_{G_1 \to G_2}} denote genomic offset between genomic
#' states \eqn{G_1} and \eqn{G_2}, and \eqn{\Delta E_{E_1 \to E_2}} the
#' corresponding environmental change.  A rate of genomic offset is defined as:
#'
#' \deqn{r_G = \frac{\Delta G_{G_1 \to G_2}}{\Delta E_{E_1 \to E_2}} \quad (2)}
#'
#' Moving to infinitesimals:
#'
#' \deqn{r_G = \frac{dG}{dE} \quad (3)}
#'
#' The velocity of environmental change over time is:
#'
#' \deqn{v_E = \frac{dE}{dt} \quad (4)}
#'
#' Finally, the **Required Evolutionary Rate** is the product:
#'
#' \deqn{REvoRate = \frac{dG}{dt} = \frac{dG}{dE} \cdot \frac{dE}{dt} \quad (5)}
#'
#' In practice this function computes the **discrete, baseline** version of
#' Eq. (5).  For the **first** period, the reference genomic state is the
#' present day (offset = 0), so REvoRate collapses to:
#'
#' \deqn{RR_{t_1} = GO_{t_1} \,/\, \Delta t_1}
#'
#' For **subsequent** periods, only the incremental offset change relative to
#' the previous period is attributed to that interval:
#'
#' \deqn{RR_{t_k} = \frac{GO_{t_k} - GO_{t_{k-1}}}{\Delta t_k}, \quad k \geq 2}
#'
#' This baseline rate is the *minimum* rate required: any deviation from the
#' predicted period-by-period optimum would necessitate a higher rate over
#' the same time window.
#'
#' ## Interpreting offset and REvoRate jointly
#'
#' The two metrics define four qualitative regimes
#' (Morando-Milà et al. 2025, Introduction):
#'
#' | Offset | REvoRate | Interpretation |
#' |--------|----------|----------------|
#' | **High** | **High** | Worst case — high adaptive stress; large mismatch AND rapid change |
#' | **High** | Low      | Chronic lag — large accumulated mismatch, slow climate velocity |
#' | Low      | **High** | Transient stress — rapid environmental change, small current mismatch |
#' | Low      | Low      | Near equilibrium — population tracking the optimum |
#'
#' @param go_list Named list of `RasterLayer` objects returned by repeated
#'   calls to [compute_genetic_offset()]. Names must match those of `delta_t`.
#' @param delta_t Named numeric vector of year-deltas for each period.  The
#'   first element is the interval between `current_midpoint` and the first
#'   future period midpoint; subsequent elements are inter-period intervals.
#'   Typically produced internally by [REvoRate()].
#'
#' @return Named list of `RasterLayer` objects — one per period — containing
#'   pixel-wise REvoRate values in units of *genetic-offset units per year*.
#'
#' @references
#' Morando-Milà J, Grau O, Ulaszewski B, Vilà-Cabrera A, Peñuelas J, Jump A,
#' Scotti I (2025). Re-interpreting genetic offset: quantifying the least
#' required evolutionary rate under climate change at the Mediterranean range
#' margin of European beech. *Evolutionary Applications*.
#'
#' @seealso [REvoRate()], [compute_genetic_offset()]
#'
#' @examples
#' \dontrun{
#' go_list <- list(
#'   p1 = compute_genetic_offset(gfRef, current_stack, future_p1),
#'   p2 = compute_genetic_offset(gfRef, current_stack, future_p2),
#'   p3 = compute_genetic_offset(gfRef, current_stack, future_p3),
#'   p4 = compute_genetic_offset(gfRef, current_stack, future_p4)
#' )
#'
#' # delta_t: 45 yrs to p1 midpoint (1985 -> 2030), then 20 yrs between periods
#' rr <- compute_revoRate(go_list,
#'                        delta_t = c(p1 = 45, p2 = 20, p3 = 20, p4 = 20))
#'
#' par(mfrow = c(2, 2))
#' for (nm in names(rr)) plot(rr[[nm]], main = paste("REvoRate", nm))
#' }
#'
#' @export
compute_revoRate <- function(go_list, delta_t) {

  if (!is.list(go_list) || is.null(names(go_list)))
    stop("`go_list` must be a named list.")
  if (!is.numeric(delta_t) || is.null(names(delta_t)))
    stop("`delta_t` must be a named numeric vector.")
  if (!all(names(go_list) %in% names(delta_t)))
    stop("All names in `go_list` must be present in `delta_t`.")

  periods <- names(go_list)
  rr_list <- vector("list", length(periods))
  names(rr_list) <- periods

  for (i in seq_along(periods)) {
    nm <- periods[i]
    dt <- delta_t[nm]

    if (i == 1) {
      # First period: present-day offset = 0, so RR = GO_t1 / delta_t1
      # Discrete version of Eq. 5: dG/dt with dG = GO_t1, dt = delta_t1
      diff_offset <- go_list[[nm]]
    } else {
      # Subsequent periods: incremental change in offset (discrete dG)
      # divided by the inter-period interval (discrete dt) — Eq. 5
      nm_prev     <- periods[i - 1]
      diff_offset <- go_list[[nm]] - go_list[[nm_prev]]
    }

    rr_list[[nm]] <- diff_offset / dt
  }

  rr_list
}
