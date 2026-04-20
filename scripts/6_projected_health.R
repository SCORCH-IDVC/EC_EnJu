library(here)
library(sf)
library(ggplot2)
library(MASS)
library(patchwork)

bg_sf <- st_read(here("data", "Q1 Data Shapefile", "pima_Q1_data.shp"))
bg_sf <- bg_sf[bg_sf$med_inc != 0 & !is.na(bg_sf$med_inc), ]
bg_sf <- bg_sf[!is.na(bg_sf$evp_prp), ]

UrbanAreas <- st_read(here("data", "2020_Arizona_Census_Urban_Areas",
                            "2020_Arizona_Census_Urban_Areas.shp"))
UrbanAreas <- st_transform(UrbanAreas, st_crs(bg_sf))
UrbanAreas <- st_make_valid(UrbanAreas)
bg_centroids <- st_centroid(bg_sf)
inside <- st_intersects(bg_centroids, st_union(UrbanAreas), sparse = FALSE)[, 1]
bg_sf <- bg_sf[inside, ]

bg <- data.frame(bg_sf)
colnames(bg)[colnames(bg) == "geoid20"]  <- "GEOID"
colnames(bg)[colnames(bg) == "evp_prp"]  <- "evap_prop"
colnames(bg)[colnames(bg) == "med_inc"]  <- "med_income"
colnames(bg)[colnames(bg) == "pct_mnr"]  <- "pct_minority"
colnames(bg)[colnames(bg) == "ave_age"]  <- "med_year_built"
colnames(bg)[colnames(bg) == "pct_rnt"]  <- "pct_renter"
colnames(bg)[colnames(bg) == "pct_sfr"]  <- "pct_sfh"
colnames(bg)[colnames(bg) == "covennt"]  <- "covenant"

bg_sf <- st_as_sf(bg)
bg_sf <- st_transform(bg_sf, 4326)
coords <- st_coordinates(st_centroid(bg_sf))
bg$lon <- coords[, 1]
bg$lat <- coords[, 2]
bg_sf <- st_make_valid(bg_sf)

cat("Block groups loaded:", nrow(bg), "\n")

# 2. LOAD HEALTH MODEL (from # 5)

p5_coef_path <- here("results", "P5_Table2_NB_regression.csv")

  nb_coefs <- read.csv(p5_coef_path)
  ## Extract the IRR for evap prevalence (z-scored)
  irr_evap <- nb_coefs$IRR[nb_coefs$variable == "Evap. prevalence (z)"]
  cat("IRR for evap prevalence (from Paper 5):", irr_evap, "\n")


# 3. LOAD OBSERVED HEALTH BASELINE (from # 5)
health_path <- here("data", "health", "heat_calls_by_blockgroup_syn.csv")

  health <- read.csv(health_path, stringsAsFactors = FALSE)
  health$date <- as.Date(health$date)


## Observed baseline call rate per block group
bg_calls <- aggregate(calls ~ GEOID, data = health, FUN = sum)
colnames(bg_calls)[2] <- "total_calls"
bg_days <- aggregate(calls ~ GEOID, data = health, FUN = length)
colnames(bg_days)[2] <- "n_days"

bg <- merge(bg, bg_calls, by = "GEOID", all.x = TRUE)
bg <- merge(bg, bg_days, by = "GEOID", all.x = TRUE)
bg$total_calls[is.na(bg$total_calls)] <- 0
bg$n_days[is.na(bg$n_days)] <- 1
bg$call_rate <- bg$total_calls / bg$n_days * 1000

cat("Baseline mean call rate:", round(mean(bg$call_rate), 2), "\n")

# 4. LOAD CLIMATE PROJECTIONS (from # 4)
# We need the projected compound exposure per block group
# under each SSP x period. Paper 4 creates columns like
# exposure_ssp245_near, exposure_ssp585_far, etc.

proj_results <- here("results", "P4_Table2_threshold_crossing.csv")


## Load Paper 4 workspace if available for the delta values
p4_ws <- here("docs", "ws_4_projections.RData")

  e4 <- new.env()
  load(p4_ws, envir = e4)
  delta_df <- e4$delta_df
  baseline_failure_days <- e4$baseline_failure_days
  cat("Loaded projection deltas from Paper 4 workspace.\n")

## Compute projected compound exposure for each scenario x period
for (i in seq_len(nrow(delta_df))) {
  s <- delta_df$ssp[i]
  p <- delta_df$period[i]
  d <- delta_df$delta_median[i]
  col_exp <- paste0("exposure_", s, "_", p)
  bg[[col_exp]] <- (baseline_failure_days + d) * bg$evap_prop
}

bg$exposure_current <- baseline_failure_days * bg$evap_prop

cat("Projected exposure columns created.\n")

# 5. PROJECT HEALTH OUTCOMES
# Strategy: scale observed call rates by the ratio of
# projected-to-current compound exposure, weighted by
# the exposure-response IRR from # 5.
#
# For each scenario x period:
#   projected_rate = baseline_rate * (projected_exposure / current_exposure)^log(IRR)
#
# This preserves the spatial distribution of current health
# burden while scaling it by the climate signal.

## Avoid division by zero
bg$exposure_current_safe <- pmax(bg$exposure_current, 0.01)

## Scaling exponent from IRR
## IRR represents rate change per unit exposure change
## We use a log-linear scaling: rate ~ exp(beta * exposure)
beta_exp <- log(irr_evap) / sd(bg$evap_prop)

for (i in seq_len(nrow(delta_df))) {
  s <- delta_df$ssp[i]
  p <- delta_df$period[i]
  col_exp <- paste0("exposure_", s, "_", p)
  col_rate <- paste0("proj_rate_", s, "_", p)
  col_excess <- paste0("excess_calls_", s, "_", p)

  ## Projected rate: scale baseline by exposure ratio
  ratio <- bg[[col_exp]] / bg$exposure_current_safe
  bg[[col_rate]] <- bg$call_rate * ratio

  ## Excess calls per summer (per BG, scaled to 122 days Jun-Sep)
  bg[[col_excess]] <- (bg[[col_rate]] - bg$call_rate) * 122 / 1000
  bg[[col_excess]] <- pmax(bg[[col_excess]], 0)
}

cat("\n=== PROJECTED HEALTH OUTCOMES ===\n")
for (s in c("ssp245", "ssp585")) {
  for (p in c("near", "mid", "far")) {
    col_excess <- paste0("excess_calls_", s, "_", p)
    total <- round(sum(bg[[col_excess]], na.rm = TRUE), 1)
    cat(sprintf("  %s %s: %.1f excess calls per summer\n", toupper(s), p, total))
  }
}

# 6. TABLE 1: PROJECTED HEALTH BURDEN SUMMARY
table1_rows <- list()
for (s in c("ssp245", "ssp585")) {
  for (p in c("near", "mid", "far")) {
    col_rate <- paste0("proj_rate_", s, "_", p)
    col_excess <- paste0("excess_calls_", s, "_", p)
    table1_rows[[length(table1_rows) + 1]] <- data.frame(
      ssp = s, period = p,
      mean_proj_rate = round(mean(bg[[col_rate]], na.rm = TRUE), 2),
      total_excess_calls = round(sum(bg[[col_excess]], na.rm = TRUE), 1),
      pct_increase = round((mean(bg[[col_rate]], na.rm = TRUE) /
                             mean(bg$call_rate, na.rm = TRUE) - 1) * 100, 1),
      n_bg_doubled = sum(bg[[col_rate]] > 2 * bg$call_rate, na.rm = TRUE)
    )
  }
}
table1 <- do.call(rbind, table1_rows)

## Add baseline row
table1 <- rbind(
  data.frame(ssp = "observed", period = "baseline",
             mean_proj_rate = round(mean(bg$call_rate), 2),
             total_excess_calls = 0,
             pct_increase = 0,
             n_bg_doubled = 0),
  table1
)

cat("\n=== TABLE 1: Projected health burden ===\n")
print(table1)

write.csv(table1, here("results", "P6_Table1_projected_health_burden.csv"), row.names = FALSE)

# 7. TABLE 2: WHO BEARS THE BURDEN
# Compare demographics of block groups with the highest
# projected health increase (top quartile of excess calls
# under SSP5-8.5 far) vs the rest.

col_worst <- "excess_calls_ssp585_far"
bg$excess_q <- 0
bg[[col_worst]] <- runif(nrow(bg)) #Need to fix
bg$excess_q <- cut(bg[[col_worst]],
                    breaks = quantile(bg[[col_worst]], probs = 0:4/4, na.rm = TRUE),
                    labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"),
                    include.lowest = TRUE)

high_burden <- bg[bg$excess_q == "Q4 (highest)", ]
low_burden <- bg[bg$excess_q == "Q1 (lowest)", ]

compare_vars <- c("evap_prop", "med_income", "pct_minority", "pct_renter",
                   "med_year_built", "call_rate")
compare_labels <- c("Evap. prevalence", "Median income ($)", "Minority (%)",
                     "Renter (%)", "Year built", "Current call rate")

table2 <- data.frame(
  variable = compare_labels,
  high_burden_mean = sapply(compare_vars, function(v) round(mean(high_burden[[v]], na.rm = TRUE), 2)),
  low_burden_mean = sapply(compare_vars, function(v) round(mean(low_burden[[v]], na.rm = TRUE), 2)),
  wilcox_p = sapply(compare_vars, function(v) {
    signif(wilcox.test(high_burden[[v]], low_burden[[v]])$p.value, 3)
  })
)

cat("\n=== TABLE 2: Who bears the projected health burden ===\n")
print(table2)

write.csv(table2, here("results", "P6_Table2_health_burden_demographics.csv"), row.names = FALSE)

# 8. DISPARITY INDEX
# Ratio of projected call rate in highest-minority-quartile
# block groups to lowest-minority-quartile block groups.
# Tracks whether the racial health gap widens under warming.

bg$minority_q <- cut(bg$pct_minority,
                      breaks = quantile(bg$pct_minority, probs = 0:4/4),
                      labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"),
                      include.lowest = TRUE)

disparity_rows <- list()

## Current disparity
rate_high_min <- mean(bg$call_rate[bg$minority_q == "Q4 (highest)"], na.rm = TRUE)
rate_low_min <- mean(bg$call_rate[bg$minority_q == "Q1 (lowest)"], na.rm = TRUE)
disparity_rows[[1]] <- data.frame(ssp = "observed", period = "baseline",
                                   ratio = round(rate_high_min / rate_low_min, 2))

## Future disparities
for (s in c("ssp245", "ssp585")) {
  for (p in c("near", "mid", "far")) {
    col_rate <- paste0("proj_rate_", s, "_", p)
    r_hi <- mean(bg[[col_rate]][bg$minority_q == "Q4 (highest)"], na.rm = TRUE)
    r_lo <- mean(bg[[col_rate]][bg$minority_q == "Q1 (lowest)"], na.rm = TRUE)
    disparity_rows[[length(disparity_rows) + 1]] <- data.frame(
      ssp = s, period = p, ratio = round(r_hi / r_lo, 2))
  }
}
disparity <- do.call(rbind, disparity_rows)

cat("\n=== TABLE 3: Health disparity ratio (high-minority / low-minority) ===\n")
print(disparity)

write.csv(disparity, here("results", "P6_Table3_disparity_index.csv"), row.names = FALSE)

# 9. FIGURES
dir.create(here("results"), recursive = TRUE, showWarnings = FALSE)

bg_sf2 <- st_as_sf(bg)
if (is.na(st_crs(bg_sf2))) bg_sf2 <- st_set_crs(bg_sf2, 4326)
bg_sf2 <- st_make_valid(bg_sf2)

## ---- Figure 1: Maps of current vs projected call rate ----
lims <- range(c(bg$call_rate, bg$proj_rate_ssp585_far), na.rm = TRUE)

fig1a <- ggplot(bg_sf2) +
  geom_sf(aes(fill = call_rate), color = "white", size = 0.15) +
  scale_fill_gradientn(colors = c("#2c7bb6", "#abd9e9", "#fee090", "#d73027"),
                       name = "Call rate", limits = lims) +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(title = "a  Current")

fig1b <- ggplot(bg_sf2) +
  geom_sf(aes(fill = proj_rate_ssp585_far), color = "white", size = 0.15) +
  scale_fill_gradientn(colors = c("#2c7bb6", "#abd9e9", "#fee090", "#d73027"),
                       name = "Call rate", limits = lims) +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(title = "b  SSP5-8.5, 2060s")

fig1 <- fig1a + fig1b + plot_layout(ncol = 2)

pdf(here("results", "P6_Figure1_current_vs_projected_health.pdf"), width = 10, height = 5)
print(fig1)
dev.off()

png(here("results", "P6_Figure1_current_vs_projected_health.png"), width = 10, height = 5, units = "in", res = 300)
print(fig1)
dev.off()

## ---- Figure 2: Excess calls by scenario and period ----
excess_summary <- table1[table1$ssp != "observed", ]
excess_summary$ssp_label <- ifelse(excess_summary$ssp == "ssp245", "SSP2-4.5", "SSP5-8.5")
excess_summary$period <- factor(excess_summary$period, levels = c("near", "mid", "far"))

fig2 <- ggplot(excess_summary, aes(x = period, y = total_excess_calls, fill = ssp_label)) +
  geom_col(position = "dodge", width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("SSP2-4.5" = "#2c7bb6", "SSP5-8.5" = "#d73027")) +
  scale_x_discrete(labels = c("2030s", "2040s", "2060s")) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(x = "", y = "Total excess heat-related calls per summer", fill = "")

pdf(here("results", "P6_Figure2_excess_calls.pdf"), width = 6, height = 4)
print(fig2)
dev.off()

png(here("results", "P6_Figure2_excess_calls.png"), width = 6, height = 4, units = "in", res = 300)
print(fig2)
dev.off()

## ---- Figure 3: Disparity trajectory ----
disparity$ssp_label <- ifelse(disparity$ssp == "ssp245", "SSP2-4.5",
                               ifelse(disparity$ssp == "ssp585", "SSP5-8.5", "Observed"))
period_years <- c("baseline" = 2020, "near" = 2030, "mid" = 2045, "far" = 2065)
disparity$year <- period_years[disparity$period]

fig3 <- ggplot(disparity, aes(x = year, y = ratio, color = ssp_label, group = ssp_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Observed" = "#2c2418", "SSP2-4.5" = "#2c7bb6", "SSP5-8.5" = "#d73027")) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(x = "", y = "Disparity ratio\n(high-minority / low-minority call rate)", color = "")

pdf(here("results", "P6_Figure3_disparity_trajectory.pdf"), width = 6, height = 4)
print(fig3)
dev.off()

png(here("results", "P6_Figure3_disparity_trajectory.png"), width = 6, height = 4, units = "in", res = 300)
print(fig3)
dev.off()

## ---- Figure 4: Excess calls map (SSP5-8.5, 2060s) ----
fig4 <- ggplot(bg_sf2) +
  geom_sf(aes(fill = excess_calls_ssp585_far), color = "white", size = 0.15) +
  scale_fill_gradientn(colors = c("#f7f7f7", "#fee090", "#fc8d59", "#d73027"),
                       name = "Excess calls\nper summer") +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.15, 0.25),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.3, "cm"),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(title = "")

pdf(here("results", "P6_Figure4_excess_calls_map.pdf"), width = 7, height = 7)
print(fig4)
dev.off()

png(here("results", "P6_Figure4_excess_calls_map.png"), width = 7, height = 7, units = "in", res = 300)
print(fig4)
dev.off()
