# Required Evolutionary Rate R-package (REvoRateR)
This R package streamlines the computation of the novel index to complement the genetic offset metric.
# REvoRateR <img src="man/figures/logo.png" align="right" height="120" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/your-github-username/REvoRateR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/your-github-username/REvoRateR/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**REvoRateR** computes the *Required Evolutionary Rate* — the minimum rate of
evolutionary change a population must sustain to remain locally adapted under
projected climate change.  It is derived from **genetic offset** estimated by
a [Gradient Forest](https://doi.org/10.1890/11-0252.1) model over successive
future climate periods.

Developed alongside the paper:

> Morando-Milà J, et al. (2026) Re-interpreting genetic offset: quantifying
> the least required evolutionary rate under climate change at the
> Mediterranean range margin of European beech. *Evolutionary Applications*.

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("your-github-username/REvoRateR")
```

> **Note:** `gradientForest` is not on CRAN. Install it first:
> ```r
> install.packages("gradientForest",
>   repos = "http://R-Forge.R-project.org")
> # or from devtools:
> devtools::install_github("pgreenwell/gradientForest")
> ```

---

## Quick start

```r
library(REvoRateR)

result <- REvoRate(
  gf_model         = gfRef,              # fitted gradientForest model
  current_stack    = current_stack,      # RasterStack, present-day BIOs
  future_stacks    = list(               # one RasterStack per period
    p1 = raster_2021_2040,
    p2 = raster_2041_2060,
    p3 = raster_2061_2080,
    p4 = raster_2081_2100
  ),
  period_midpoints = c(p1=2030, p2=2050, p3=2070, p4=2090),
  current_midpoint = 1985,              # midpoint of reference period
  coords           = coord_df,          # data.frame with Long, Lat columns
  compute_climate_mismatch = TRUE
)

# Tidy tables
write.csv(result$site_wide, "REvoRate_wide.csv", row.names = FALSE)
write.csv(result$site_long, "REvoRate_long.csv", row.names = FALSE)

# Publication figure
plot_revoRate(result, save_path = "Fig_REvoRate.pdf")
```

---

## Output

`REvoRate()` returns a list containing:

| Element | Description |
|---|---|
| `$genetic_offset` | Named list of `RasterLayer` — genetic offset per period |
| `$revoRate` | Named list of `RasterLayer` — REvoRate per period |
| `$climate_mismatch` | Plain (unweighted) climate mismatch rasters |
| `$climate_mismatch_gf` | GF-importance-weighted climate mismatch rasters |
| `$delta_t` | Year-deltas used for each period |
| `$site_wide` | Wide-format table (one row per site, one col per metric × period) |
| `$site_long` | Long-format (tidy) table |

---

## Methods

### Genetic offset
Euclidean distance in GF-transformed environmental space between present-day
and future climate predictions at each pixel:

$$GO = \sqrt{\sum_j \left( GF(x_j^{future}) - GF(x_j^{current}) \right)^2}$$

### REvoRate
For the **first period**, REvoRate is the offset divided by years since the
present-day midpoint:

$$RR_{t_1} = GO_{t_1} \,/\, \Delta t_1$$

For **subsequent periods**, only the incremental increase in offset matters:

$$RR_{t_k} = (GO_{t_k} - GO_{t_{k-1}}) \,/\, \Delta t_k, \quad k \geq 2$$

### Interpreting offset and REvoRate jointly

The two metrics define four qualitative risk regimes:

| Offset | REvoRate | Interpretation |
|--------|----------|----------------|
| **High** | **High** | Worst case — large mismatch AND rapid change required |
| **High** | Low | Chronic lag — large mismatch, slow climate velocity |
| Low | **High** | Transient stress — rapid change, small current mismatch |
| Low | Low | Near equilibrium — population tracking the optimum |

Note that in period 1, offset and REvoRate are linearly related (both measured from the same zero baseline). Differences in ranking between the two metrics only emerge from period 2 onwards.



| Function | Description |
|---|---|
| `REvoRate()` | Main wrapper — runs the full pipeline |
| `compute_genetic_offset()` | Single-period genetic offset raster |
| `compute_revoRate()` | REvoRate from a list of offset rasters |
| `compute_climate_mismatch_plain()` | Unweighted Euclidean BIO mismatch |
| `compute_climate_mismatch_gf()` | GF-importance-weighted BIO mismatch |
| `extract_site_metrics()` | Extract raster values at site coordinates |
| `plot_revoRate()` | Faceted trajectory plot |

---

## Citation

If you use this package, please cite the original paper:

> Morando-Milà J, Grau O, Ulaszewski B, Vilà-Cabrera A, Peñuelas J, Jump A,
> Scotti I (2026). Re-interpreting genetic offset: quantifying the least
> required evolutionary rate under climate change at the Mediterranean range
> margin of European beech. *Evolutionary Applications*.

BibTeX:
```bibtex
@article{MorandoMila2025,
  author  = {Morando-Mil{\`a}, Josep and Grau, Oriol and Ulaszewski, Bartosz
             and Vil{\`a}-Cabrera, Albert and Pe{\~n}uelas, Josep
             and Jump, Alistair and Scotti, Ivan},
  title   = {Re-interpreting genetic offset: quantifying the least required
             evolutionary rate under climate change at the {Mediterranean}
             range margin of {European} beech},
  journal = {Evolutionary Applications},
  year    = {2026}
}
```

---

## License

MIT © Josep Morando-Milà
