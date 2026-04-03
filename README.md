# Evaporative cooling and environmental justice

This repo contains the analysis code for analyses examining who relies on evaporative cooling in Tucson [or beyond], where those households are, and what that means for public health.

## Why this matters

Evaporative coolers are cheap to run but fail when humidity rises. Swamp coolers also pull unfiltered outdoor air indoors (dust, wildfire smoke, ozone).

## Structure

```
data/           Raw data (not tracked; see below)
results/        Tables (.csv) and figures (.pdf, .png)
scripts/
  1_equity_environmental_justice.R
```

## Data
The evaporative cooling dataset is not public yet. The script includes a simulated dataset that mirrors the real data structure so the full pipeline runs end to end. Replace the simulation block (Section 1) with real data reads when available.

## Requirements

R ≥ 4.2. Packages: `here`, `sf`, `ggplot2`, `spdep`, `patchwork`. Install with:

```r
install.packages(c("here", "sf", "ggplot2", "spdep", "patchwork"))
```