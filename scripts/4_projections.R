library(here)
library(sf)
library(terra)
library(ncdf4)
library(ggplot2)
library(patchwork)

bg_sf <- st_read(here("data", "Q1 Data Shapefile", "pima_Q1_data.shp"))

#Dropping empty polygons
bg_sf <- bg_sf[bg_sf$med_inc != 0 & !is.na(bg_sf$med_inc), ]
bg_sf <- bg_sf[!is.na(bg_sf$evp_prp), ]

## Keep only those inside urban areas
UrbanAreas <- st_read(here("data", "2020_Arizona_Census_Urban_Areas", "2020_Arizona_Census_Urban_Areas.shp"))
UrbanAreas <- st_transform(UrbanAreas, st_crs(bg_sf))
UrbanAreas <- st_make_valid(UrbanAreas)
bg_centroids <- st_centroid(bg_sf)
inside <- st_intersects(bg_centroids, st_union(UrbanAreas), sparse = FALSE)[, 1]
bg_sf <- bg_sf[inside, ]

bg <- data.frame(bg_sf)

## Rename columns to match the simulated dataset
colnames(bg)[colnames(bg) == "geoid20"]   <- "GEOID"
colnames(bg)[colnames(bg) == "evp_prp"]   <- "evap_prop"
colnames(bg)[colnames(bg) == "med_inc"]    <- "med_income"
colnames(bg)[colnames(bg) == "pct_mnr"]    <- "pct_minority"
colnames(bg)[colnames(bg) == "ave_age"]    <- "med_year_built"
colnames(bg)[colnames(bg) == "pct_rnt"]    <- "pct_renter"
colnames(bg)[colnames(bg) == "pct_sfr"]    <- "pct_sfh"
colnames(bg)[colnames(bg) == "covennt"]    <- "covenant"

## Compute centroids for lon/lat
bg_sf <- st_as_sf(bg)
bg_sf <- st_transform(bg_sf, 4326)
coords <- st_coordinates(st_centroid(bg_sf))
bg$lon <- coords[, 1]
bg$lat <- coords[, 2]
bg_sf <- st_make_valid(bg_sf)

# 2. OBSERVED BASELINE (from Script 3)

wx <- read.csv(here("data", "weather", "tucson_hourly.csv"))
wx$failure <- wx$relh > 30 & wx$tmpf > 95

## Daily failure counts
daily_obs <- aggregate(failure ~ date, data = wx, FUN = sum)
colnames(daily_obs)[2] <- "failure_hours"
daily_obs$year <- as.integer(format(as.Date(daily_obs$date), "%Y"))
daily_obs$month <- as.integer(format(as.Date(daily_obs$date), "%m"))
daily_obs$failure_day <- as.integer(daily_obs$failure_hours > 0)

## Observed baseline: mean annual failure days (Jun-Sep)
obs_summer <- daily_obs[daily_obs$month >= 6 & daily_obs$month <= 9, ]
obs_annual <- aggregate(failure_day ~ year, data = obs_summer, FUN = sum)
baseline_failure_days <- mean(obs_annual$failure_day)
baseline_failure_se <- sd(obs_annual$failure_day) / sqrt(nrow(obs_annual))

cat("\n=== OBSERVED BASELINE ===\n")
cat("Mean annual failure days (Jun-Sep):", round(baseline_failure_days, 1),
    "+/-", round(baseline_failure_se, 1), "\n")
cat("Years:", paste(range(obs_annual$year), collapse = "-"), "\n")

# 3. DOWNLOAD NEX-GDDP-CMIP6 PROJECTIONS
# NASA NEX-GDDP-CMIP6: daily downscaled climate projections
# at 0.25 degree resolution. Bias-corrected against observations.

dir.create(here("data", "cmip6"), recursive = TRUE, showWarnings = FALSE)
projections_path <- here("data", "cmip6", "tucson_projections.csv")

if (!file.exists(projections_path)) {
  cat("=== Downloading NEX-GDDP-CMIP6 projections ===\n")
  
  ## ~ Tucson coordinates
  tuc_lat <- 32.25
  tuc_lon <- -111.0
  models <- c("ACCESS-CM2", "GFDL-ESM4", "IPSL-CM6A-LR",
              "MPI-ESM1-2-HR", "UKESM1-0-LL")
  ssps <- c("ssp245", "ssp585")
  periods <- list(
    historical = 2005:2014,
    near       = 2025:2034,
    mid        = 2040:2049,
    far        = 2060:2069
  )
  
  variables <- c("tasmax", "hurs")
  s3_base <- "https://nex-gddp-cmip6.s3.us-west-2.amazonaws.com/NEX-GDDP-CMIP6"
  
  all_results <- list()
  
  for (mod in models) {
    for (ssp in ssps) {
      for (period_name in names(periods)) {
        yrs <- periods[[period_name]]
        experiment <- ifelse(period_name == "historical", "historical", ssp)
        
        for (yr in yrs) {
          for (v in variables) {
            
            fname <- sprintf("%s_day_%s_%s_r1i1p1f1_gn_%d.nc", v, mod, experiment, yr)
            url <- sprintf("%s/%s/%s/r1i1p1f1/%s/%s", s3_base, mod, experiment, v, fname)
            
            tmp <- tempfile(fileext = ".nc")
            result <- try(download.file(url, tmp, mode = "wb", quiet = TRUE), silent = TRUE)
            
            if (!inherits(result, "try-error") && file.size(tmp) > 10000) {
              nc_try <- try({
                r <- rast(tmp)
                
                ## Extract the Tucson pixel for Jun-Sep
                r <- rotate(r)
                cell_id <- cellFromXY(r, cbind(tuc_lon, tuc_lat))
                vals <- r[cell_id]
                vals <- as.numeric(vals)
                
                ## Time: layer count = days in year. Subset to Jun-Sep (DOY 152-273)
                is_leap <- (yr %% 4 == 0 & yr %% 100 != 0) | (yr %% 400 == 0)
                jun1 <- ifelse(is_leap, 153, 152)
                sep30 <- ifelse(is_leap, 274, 273)
                
                if (length(vals) >= sep30) {
                  summer_vals <- vals[jun1:sep30]
                  dates <- seq(as.Date(paste0(yr, "-06-01")),
                               as.Date(paste0(yr, "-09-30")), by = "day")
                  if (length(summer_vals) == length(dates)) {
                    d <- data.frame(
                      model = mod, ssp = ssp, period = period_name,
                      year = yr, date = dates, variable = v,
                      value = summer_vals
                    )
                    all_results[[length(all_results) + 1]] <- d
                  }
                }
              }, silent = TRUE)
              unlink(tmp)
            }
          }
          cat("  ", mod, experiment, yr, "\n")
        }
      }
    }
  }
  proj_raw <- do.call(rbind, all_results)
  write.csv(proj_raw, projections_path, row.names = FALSE)
  cat("Projections saved:", nrow(proj_raw), "rows\n")
}

proj <- read.csv(projections_path, stringsAsFactors = FALSE)
proj$date <- as.Date(proj$date)
cat("Projection rows loaded:", nrow(proj), "\n")

# 4. COMPUTE PROJECTED COOLER FAILURE DAYS

proj_tmax <- proj[proj$variable == "tasmax", c("model", "ssp", "period", "year", "date", "value")]
proj_hurs <- proj[proj$variable == "hurs", c("model", "ssp", "period", "year", "date", "value")]
colnames(proj_tmax)[6] <- "tasmax_K"
colnames(proj_hurs)[6] <- "hurs"

proj_wide <- merge(proj_tmax, proj_hurs,
                   by = c("model", "ssp", "period", "year", "date"))

## Convert tasmax from Kelvin to Fahrenheit
proj_wide$tmax_f <- proj_wide$tasmax_K * 9/5 - 459.67

## Failure threshold
proj_wide$failure <- proj_wide$hurs > 30 & proj_wide$tmax_f > 95

## Annual failure days per model x ssp x period
annual_fail <- aggregate(failure ~ model + ssp + period + year,
                         data = proj_wide, FUN = sum)
colnames(annual_fail)[5] <- "failure_days"

## Mean across years within each period
period_fail <- aggregate(failure_days ~ model + ssp + period,
                         data = annual_fail, FUN = mean)

cat("\n=== PROJECTED FAILURE DAYS (mean per summer) ===\n")
print(period_fail)

# 5. TABLE 1: FAILURE DAYS BY PERIOD AND SCENARIO
# Ensemble summary: median [10th-90th percentile] across models

summary_by_period <- function(df) {
  data.frame(
    median = round(median(df$failure_days), 1),
    p10 = round(quantile(df$failure_days, 0.1), 1),
    p90 = round(quantile(df$failure_days, 0.9), 1),
    n_models = length(unique(df$model))
  )
}

table1_list <- list()
for (s in c("ssp245", "ssp585")) {
  for (p in c("historical", "near", "mid", "far")) {
    sub <- period_fail[period_fail$ssp == s & period_fail$period == p, ]
    if (nrow(sub) == 0 && p == "historical") {
      sub <- period_fail[period_fail$period == "historical", ]
    }
    if (nrow(sub) > 0) {
      row <- summary_by_period(sub)
      row$ssp <- s
      row$period <- p
      table1_list[[length(table1_list) + 1]] <- row
    }
  }
}
table1 <- do.call(rbind, table1_list)
table1$label <- paste0(table1$median, " [", table1$p10, "-", table1$p90, "]")

cat("\n=== TABLE 1: Projected failure days ===\n")
print(table1[, c("ssp", "period", "label", "n_models")])

write.csv(table1, here("results", "P4_Table1_projected_failure_days.csv"), row.names = FALSE)

# 6. DELTA METHOD: PROJECT BLOCK-GROUP-LEVEL EXPOSURE
# Delta = projected failure days - model historical failure days.
# Applied to observed baseline to avoid raw model bias.
# Each block group's projected exposure = (observed_baseline + delta) x evap_prop

## Model historical baseline
hist_fail <- period_fail[period_fail$period == "historical", ]
colnames(hist_fail)[3] <- "hist_period"
colnames(hist_fail)[4] <- "hist_failure"

## Compute deltas for each model x ssp x future period
future_fail <- period_fail[period_fail$period != "historical", ]
future_fail <- merge(future_fail, hist_fail[, c("model", "hist_failure")], by = "model")
future_fail$delta <- future_fail$failure_days - future_fail$hist_failure

## Ensemble median delta per ssp x period
delta_summary <- aggregate(delta ~ ssp + period, data = future_fail,
                           FUN = function(x) c(median = median(x),
                                               p10 = quantile(x, 0.1),
                                               p90 = quantile(x, 0.9)))
delta_df <- data.frame(
  ssp = delta_summary$ssp,
  period = delta_summary$period,
  delta_median = delta_summary$delta[, 1],
  delta_p10 = delta_summary$delta[, 2],
  delta_p90 = delta_summary$delta[, 3]
)

cat("\n=== DELTA (change in failure days from historical) ===\n")
print(delta_df)

## Project to block groups
## Projected failure exposure = (baseline + delta) x evap_prop
for (i in seq_len(nrow(delta_df))) {
  s <- delta_df$ssp[i]
  p <- delta_df$period[i]
  d <- delta_df$delta_median[i]
  col_name <- paste0("exposure_", s, "_", p)
  bg[[col_name]] <- (baseline_failure_days + d) * bg$evap_prop
}


# 7. THRESHOLD CROSSING ANALYSIS
# Define inadequacy threshold: compound exposure > X
# where compound exposure = failure_days x evap_prop
# A block group with 40% evap prevalence and 30 failure days
# has exposure = 12 (i.e., 12 household-failure-days per summer).
#
# Threshold: compound exposure > 10 (roughly equivalent to
# 25 failure days at 40% prevalence, or 50 days at 20%).

threshold <- 10

## Current exposure
bg$exposure_current <- baseline_failure_days * bg$evap_prop
bg$above_now <- bg$exposure_current > threshold

## Check each future period x scenario
crossing_summary <- data.frame()
for (s in c("ssp245", "ssp585")) {
  for (p in c("near", "mid", "far")) {
    col <- paste0("exposure_", s, "_", p)
    above <- bg[[col]] > threshold
    newly_above <- above & !bg$above_now
    crossing_summary <- rbind(crossing_summary, data.frame(
      ssp = s, period = p,
      n_above = sum(above),
      n_newly_above = sum(newly_above),
      pct_above = round(mean(above) * 100, 1),
      pct_newly_above = round(mean(newly_above) * 100, 1)
    ))
  }
}

## Add current baseline row
crossing_summary <- rbind(
  data.frame(ssp = "observed", period = "baseline",
             n_above = sum(bg$above_now),
             n_newly_above = 0,
             pct_above = round(mean(bg$above_now) * 100, 1),
             pct_newly_above = 0),
  crossing_summary
)

cat("\n=== TABLE 2: Block groups exceeding threshold ===\n")
print(crossing_summary)

write.csv(crossing_summary, here("results", "P4_Table2_threshold_crossing.csv"), row.names = FALSE)

# 8. TABLE S1: DEMOGRAPHICS OF NEWLY VULNERABLE BLOCK GROUPS

## Under SSP5-8.5 far period (worst case)
col_far <- "exposure_ssp585_far"
newly_vuln <- bg[bg[[col_far]] > threshold & !bg$above_now, ]
already_vuln <- bg[bg$above_now, ]
not_vuln <- bg[!bg$above_now & bg[[col_far]] <= threshold, ]

compare_vars <- c("evap_prop", "med_income", "pct_minority", "pct_renter", "med_year_built")
compare_labels <- c("Evap. prevalence", "Median income ($)", "Minority (%)",
                    "Renter (%)", "Year built")

tableS1 <- data.frame(
  variable = compare_labels,
  newly_mean = sapply(compare_vars, function(v) round(mean(newly_vuln[[v]], na.rm = TRUE), 2)),
  already_mean = sapply(compare_vars, function(v) round(mean(already_vuln[[v]], na.rm = TRUE), 2)),
  not_vuln_mean = sapply(compare_vars, function(v) round(mean(not_vuln[[v]], na.rm = TRUE), 2)),
  p_newly_vs_not = sapply(compare_vars, function(v) {
    if (nrow(newly_vuln) > 2 & nrow(not_vuln) > 2) {
      signif(wilcox.test(newly_vuln[[v]], not_vuln[[v]])$p.value, 3)
    } else { NA }
  })
)

cat("\n=== TABLE S1: Demographics of newly vulnerable block groups ===\n")
print(tableS1)

write.csv(tableS1, here("results", "P4_TableS1_newly_vulnerable_demographics.csv"), row.names = FALSE)

# 9. FIGURES

dir.create(here("results"), recursive = TRUE, showWarnings = FALSE)

## ---- Figure 1: Fan chart of projected failure days ----
## Observed baseline + projected ensemble spread

## Reshape for plotting
period_order <- c("historical", "near", "mid", "far")
period_years <- c(2010, 2030, 2045, 2065)

plot_data <- period_fail
plot_data$period <- factor(plot_data$period, levels = period_order)
plot_data$period_year <- period_years[match(plot_data$period, period_order)]

## Ensemble ribbons
ribbon_data <- aggregate(failure_days ~ ssp + period,
                         data = plot_data,
                         FUN = function(x) c(med = median(x),
                                             lo = quantile(x, 0.1),
                                             hi = quantile(x, 0.9)))
ribbon_df <- data.frame(
  ssp = ribbon_data$ssp,
  period = ribbon_data$period,
  median = ribbon_data$failure_days[, 1],
  lo = ribbon_data$failure_days[, 2],
  hi = ribbon_data$failure_days[, 3]
)
ribbon_df$period_year <- period_years[match(ribbon_df$period, period_order)]

## Observed baseline point
obs_point <- data.frame(period_year = 2020, failure_days = baseline_failure_days)

fig1 <- ggplot() +
  ## SSP2-4.5 ribbon
  geom_ribbon(data = ribbon_df[ribbon_df$ssp == "ssp245", ],
              aes(x = period_year, ymin = lo, ymax = hi),
              fill = "#2c7bb6", alpha = 0.2) +
  geom_line(data = ribbon_df[ribbon_df$ssp == "ssp245", ],
            aes(x = period_year, y = median), color = "#2c7bb6", linewidth = 1) +
  ## SSP5-8.5 ribbon
  geom_ribbon(data = ribbon_df[ribbon_df$ssp == "ssp585", ],
              aes(x = period_year, ymin = lo, ymax = hi),
              fill = "#d73027", alpha = 0.2) +
  geom_line(data = ribbon_df[ribbon_df$ssp == "ssp585", ],
            aes(x = period_year, y = median), color = "#d73027", linewidth = 1) +
  ## Observed baseline
  geom_point(data = obs_point, aes(x = period_year, y = failure_days),
             size = 3, shape = 18) +
  geom_hline(yintercept = baseline_failure_days, linetype = "dashed", color = "grey40") +
  ## Labels
  annotate("text", x = 2068, y = max(ribbon_df$hi[ribbon_df$ssp == "ssp585"]),
           label = "SSP5-8.5", color = "#d73027", size = 3, hjust = 0) +
  annotate("text", x = 2068, y = max(ribbon_df$hi[ribbon_df$ssp == "ssp245"]),
           label = "SSP2-4.5", color = "#2c7bb6", size = 3, hjust = 0) +
  scale_x_continuous(breaks = c(2010, 2020, 2030, 2045, 2065),
                     labels = c("Hist.", "Obs.", "2030", "2045", "2065")) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(x = "", y = "Cooler failure days per summer",
       title = "a")

## ---- Figure 2: Map of projected compound exposure (SSP5-8.5, 2060s) ----
bg_sf2 <- st_as_sf(bg)
if (is.na(st_crs(bg_sf2))) bg_sf2 <- st_set_crs(bg_sf2, 4326)
bg_sf2 <- st_make_valid(bg_sf2)

fig2 <- ggplot(bg_sf2) +
  geom_sf(aes(fill = exposure_ssp585_far), color = "white", size = 0.15) +
  scale_fill_gradientn(colors = c("#2c7bb6", "#abd9e9", "#fee090", "#d73027"),
                       name = "Compound\nexposure\n(2060s, SSP5-8.5)") +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.15, 0.25),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.3, "cm"),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(title = "b")

## ---- Figure 3: Threshold crossing bar chart ----
cross_plot <- crossing_summary[crossing_summary$ssp != "observed", ]
cross_plot$period <- factor(cross_plot$period, levels = c("near", "mid", "far"))
cross_plot$ssp_label <- ifelse(cross_plot$ssp == "ssp245", "SSP2-4.5", "SSP5-8.5")

fig3 <- ggplot(cross_plot, aes(x = period, y = pct_above, fill = ssp_label)) +
  geom_col(position = "dodge", width = 0.6, alpha = 0.8) +
  geom_hline(yintercept = crossing_summary$pct_above[crossing_summary$ssp == "observed"],
             linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("SSP2-4.5" = "#2c7bb6", "SSP5-8.5" = "#d73027"),
                    name = "") +
  scale_x_discrete(labels = c("2030s", "2040s", "2060s")) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.15, 0.85),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(x = "", y = "% block groups above threshold",
       title = "c")

## ---- Figure S1: Map of current vs future exposure side by side ----
figS1a <- ggplot(bg_sf2) +
  geom_sf(aes(fill = exposure_current), color = "white", size = 0.15) +
  scale_fill_gradientn(colors = c("#2c7bb6", "#abd9e9", "#fee090", "#d73027"),
                       name = "Current", limits = range(c(bg$exposure_current, bg$exposure_ssp585_far))) +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(title = "a  Current")

figS1b <- ggplot(bg_sf2) +
  geom_sf(aes(fill = exposure_ssp585_far), color = "white", size = 0.15) +
  scale_fill_gradientn(colors = c("#2c7bb6", "#abd9e9", "#fee090", "#d73027"),
                       name = "2060s SSP5-8.5", limits = range(c(bg$exposure_current, bg$exposure_ssp585_far))) +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(title = "b  SSP5-8.5, 2060s")


# 10. EXPORT
## Figure 1: Fan chart + map + bar chart
fig_main <- (fig1 | fig2) / (fig3 + plot_spacer()) + plot_layout(heights = c(1, 0.7))

pdf(here("results", "P4_Figure1_projections.pdf"), width = 10, height = 8)
print(fig_main)
dev.off()

png(here("results", "P4_Figure1_projections.png"), width = 10, height = 8, units = "in", res = 300)
print(fig_main)
dev.off()

## Figure S1: Current vs future maps
figS1 <- figS1a + figS1b + plot_layout(ncol = 2)

pdf(here("results", "P4_FigureS1_current_vs_future_maps.pdf"), width = 10, height = 5)
print(figS1)
dev.off()

png(here("results", "P4_FigureS1_current_vs_future_maps.png"), width = 10, height = 5, units = "in", res = 300)
print(figS1)
dev.off()
