library(here)
library(sf)
library(ggplot2)
library(MASS)
library(splines)
library(spdep)
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

health_path <- here("data", "health", "heat_calls_by_blockgroup_syn.csv")
health <- read.csv(health_path, stringsAsFactors = FALSE)
health$GEOID <- rep(bg$GEOID, each = length(unique(health$date)))
health$date <- as.Date(health$date)

cat("Health records:", nrow(health), "\n")
cat("Date range:", as.character(range(health$date)), "\n")
cat("Total calls:", sum(health$calls), "\n")

# 3. LOAD WEATHER DATA (from # 3)
wx_path <- here("data", "weather", "tucson_hourly.csv")
if (!file.exists(wx_path)) {
  stop("Run Paper 3 script first to download weather data.")
}

wx <- read.csv(wx_path, stringsAsFactors = FALSE)
wx$tmpf <- as.numeric(wx$tmpf)
wx$relh <- as.numeric(wx$relh)
wx$date <- as.Date(wx$date)
wx <- wx[!is.na(wx$tmpf) & !is.na(wx$relh), ]

## Daily weather summaries
daily_wx <- aggregate(cbind(tmax = tmpf, rh_max = relh) ~ date,
                       data = wx, FUN = max, na.rm = TRUE)
daily_wx$tmean <- aggregate(tmpf ~ date, data = wx, FUN = mean, na.rm = TRUE)$tmpf

## Cooler failure flag (from Paper 3 threshold)
## Stull (2011) wet-bulb approximation
calc_twet <- function(temp_c, rh_pct) {
  temp_c * atan(0.151977 * (rh_pct + 8.313659)^0.5) +
    atan(temp_c + rh_pct) -
    atan(rh_pct - 1.676331) +
    0.00391838 * rh_pct^1.5 * atan(0.023101 * rh_pct) -
    4.686035
}

## Supply air temperature
## ASHRAE Handbook: HVAC Systems and Equipment
## Chapter: Evaporative Air-Cooling Equipment
## saturation effectiveness
## T_supply = T_db - eta * (T_db - T_wb)
## - T_db (dry-bulb temperature)
## - T_wb (wet-bulb temperature)
## - eta (saturation efficiency)
calc_supply <- function(temp_c, rh_pct, eta = 0.85) {
  twet <- calc_twet(temp_c, rh_pct)
  temp_c - eta * (temp_c - twet)
}

#Is supp temp beyond the ASHRAE comfort threshold (>27C)
wx$failure <- calc_supply(wx$tmpf, wx$relh) > 27  # TRUE = failure
#wx$failure <- wx$relh > 30 & wx$tmpf > 95
daily_fail <- aggregate(failure ~ date, data = wx, FUN = sum)
colnames(daily_fail)[2] <- "failure_hours"
daily_wx <- merge(daily_wx, daily_fail, by = "date", all.x = TRUE)
daily_wx$failure_hours[is.na(daily_wx$failure_hours)] <- 0
daily_wx$failure_day <- as.integer(daily_wx$failure_hours > 0)

cat("Daily weather records:", nrow(daily_wx), "\n")

# 4. CROSS-SECTIONAL ANALYSIS (block group level)
# Question: Do block groups with higher evap prevalence have
# higher total heat-related call rates?

## Aggregate total calls per block group across all dates
bg_calls <- aggregate(calls ~ GEOID, data = health, FUN = sum)
colnames(bg_calls)[2] <- "total_calls"

## Number of days observed per block group
bg_days <- aggregate(calls ~ GEOID, data = health, FUN = length)
colnames(bg_days)[2] <- "n_days"

bg <- merge(bg, bg_calls, by = "GEOID", all.x = TRUE)
bg <- merge(bg, bg_days, by = "GEOID", all.x = TRUE)
bg$total_calls[is.na(bg$total_calls)] <- 0
bg$n_days[is.na(bg$n_days)] <- 1

## Call rate per 1000 days (standardized for observation period)
bg$call_rate <- bg$total_calls / bg$n_days * 1000

cat("\n=== CROSS-SECTIONAL SUMMARY ===\n")
cat("Total calls across all BGs:", sum(bg$total_calls), "\n")
cat("Mean calls per BG:", round(mean(bg$total_calls), 2), "\n")
cat("Call rate range:", round(range(bg$call_rate), 2), "\n")

## Spearman correlations: call rate vs predictors
cor_vars <- c("evap_prop", "med_income", "pct_minority", "pct_renter", "med_year_built")
cor_labels <- c("Evap. prevalence", "Median income", "% Minority",
                "% Renter", "Year built")

cor_health <- data.frame(
  variable = cor_labels,
  rho = sapply(cor_vars, function(v) {
    cor(bg$call_rate, bg[[v]], method = "spearman", use = "complete.obs")
  }),
  p = sapply(cor_vars, function(v) {
    cor.test(bg$call_rate, bg[[v]], method = "spearman")$p.value
  })
)
cor_health$rho <- round(cor_health$rho, 3)
cor_health$p <- signif(cor_health$p, 3)

cat("\n=== TABLE 1: Correlations with heat-related call rate ===\n")
print(cor_health)

write.csv(cor_health, here("results", "P5_Table1_health_correlations.csv"), row.names = FALSE)

# 5. NEGATIVE BINOMIAL REGRESSION (cross-sectional)
# Model: total calls ~ evap_prop + income + minority + renter
# Offset by log(n_days) to account for observation period.
# Use negative binomial to handle overdispersion.

bg$z_evap     <- scale(bg$evap_prop)
bg$z_income   <- scale(bg$med_income)
bg$z_minority <- scale(bg$pct_minority)
bg$z_renter   <- scale(bg$pct_renter)

m_nb <- glm.nb(total_calls ~ z_evap + z_income + z_minority + z_renter +
                 offset(log(n_days)),
               data = bg)

cat("\n=== TABLE 2: Negative binomial regression ===\n")
print(summary(m_nb))

coef_nb <- data.frame(
  variable = c("Intercept", "Evap. prevalence (z)", "Income (z)",
               "Minority (z)", "Renter (z)"),
  estimate = round(coef(m_nb), 3),
  se = round(summary(m_nb)$coefficients[, 2], 3),
  IRR = round(exp(coef(m_nb)), 3),
  p = signif(summary(m_nb)$coefficients[, 4], 3)
)

cat("\nIncidence rate ratios:\n")
print(coef_nb)

write.csv(coef_nb, here("results", "P5_Table2_NB_regression.csv"), row.names = FALSE)

# 6. DAILY TIME-SERIES ANALYSIS
# Question: Are heat-related calls higher on cooler-failure days?
# And does evap prevalence modify this effect?
#
# Aggregate daily city-wide calls, merge with weather.
# Poisson regression with spline for seasonality.

daily_calls <- aggregate(calls ~ date, data = health, FUN = sum)
daily_calls <- merge(daily_calls, daily_wx, by = "date", all.x = TRUE)
daily_calls <- daily_calls[!is.na(daily_calls$tmax), ]
daily_calls$doy <- as.integer(format(daily_calls$date, "%j"))
daily_calls$year <- as.integer(format(daily_calls$date, "%Y"))

cat("\n=== DAILY TIME-SERIES ===\n")
cat("Days with calls + weather:", nrow(daily_calls), "\n")
cat("Mean daily calls:", round(mean(daily_calls$calls), 2), "\n")
cat("Failure days:", sum(daily_calls$failure_day), "\n")

## Model: calls ~ failure_day + tmax + ns(doy) + factor(year)
m_ts <- glm(calls ~ failure_day + tmax + ns(doy, df = 4) + factor(year),
            family = poisson, data = daily_calls)

cat("\n=== TABLE 3: Daily time-series Poisson regression ===\n")
print(summary(m_ts))

## Check for overdispersion
disp <- sum(residuals(m_ts, type = "pearson")^2) / m_ts$df.residual
cat("\nDispersion parameter:", round(disp, 2), "\n")

if (disp > 1.5) {
  cat("Overdispersion detected. Refitting with negative binomial.\n")
  m_ts <- glm.nb(calls ~ failure_day + tmax + ns(doy, df = 4) + factor(year),
                  data = daily_calls)
  print(summary(m_ts))
}

## Extract key coefficients
ts_coefs <- summary(m_ts)$coefficients
ts_table <- data.frame(
  variable = c("Failure day", "Tmax (F)"),
  estimate = round(ts_coefs[c("failure_day", "tmax"), 1], 3),
  se = round(ts_coefs[c("failure_day", "tmax"), 2], 3),
  IRR = round(exp(ts_coefs[c("failure_day", "tmax"), 1]), 3),
  p = signif(ts_coefs[c("failure_day", "tmax"), 4], 3)
)
cat("\nKey time-series results (IRR):\n")
print(ts_table)

write.csv(ts_table, here("results", "P5_Table3_timeseries_results.csv"), row.names = FALSE)

# 7. INTERACTION: FAILURE DAY x EVAP PREVALENCE
# Question: Is the health impact of a cooler-failure day
# amplified in neighborhoods with more evap coolers?
#
# Panel model: calls_ij ~ failure_day_j x evap_prop_i + tmax_j
#              + income_i + ns(doy) + offset(log(1))

## Merge health data with weather and block group attributes
panel <- merge(health, daily_wx[, c("date", "tmax", "failure_day", "failure_hours")],
               by = "date", all.x = TRUE)
panel <- merge(panel, bg[, c("GEOID", "evap_prop", "med_income", "pct_minority",
                              "pct_renter")],
               by = "GEOID", all.x = TRUE)
panel <- panel[!is.na(panel$tmax) & !is.na(panel$evap_prop), ]
panel$doy <- as.integer(format(panel$date, "%j"))
panel$year <- as.integer(format(panel$date, "%Y"))

## Standardize
panel$z_evap <- scale(panel$evap_prop)
panel$z_income <- scale(panel$med_income)

cat("\n=== INTERACTION MODEL ===\n")
cat("Panel rows:", nrow(panel), "\n")

## Fit interaction model
m_interact <- glm.nb(calls ~ failure_day * z_evap + tmax + z_income +
                        ns(doy, df = 4) + factor(year),
                      data = panel)

cat("\n=== TABLE 4: Interaction model ===\n")
print(summary(m_interact))

int_coefs <- summary(m_interact)$coefficients
key_rows <- c("failure_day", "z_evap", "tmax", "z_income", "failure_day:z_evap")
key_rows <- key_rows[key_rows %in% rownames(int_coefs)]

int_table <- data.frame(
  variable = key_rows,
  estimate = round(int_coefs[key_rows, 1], 3),
  se = round(int_coefs[key_rows, 2], 3),
  IRR = round(exp(int_coefs[key_rows, 1]), 3),
  p = signif(int_coefs[key_rows, 4], 3)
)

cat("\nKey interaction results:\n")
print(int_table)

write.csv(int_table, here("results", "P5_Table4_interaction_model.csv"), row.names = FALSE)

# 8. SPATIAL AUTOCORRELATION CHECK
coords_sp <- cbind(bg$lon, bg$lat)
nb <- knn2nb(knearneigh(coords_sp, k = 5))
lw <- nb2listw(nb, style = "W")

moran_calls <- moran.test(bg$call_rate, lw)
cat("\n=== SPATIAL AUTOCORRELATION ===\n")
cat("Moran's I on call rate:", round(moran_calls$estimate[1], 3), "\n")
cat("p-value:", signif(moran_calls$p.value, 3), "\n")

# 9. FIGURES
dir.create(here("results"), recursive = TRUE, showWarnings = FALSE)

bg_sf2 <- st_as_sf(bg)
if (is.na(st_crs(bg_sf2))) bg_sf2 <- st_set_crs(bg_sf2, 4326)
bg_sf2 <- st_make_valid(bg_sf2)

## ---- Figure 1a: Map of heat-related call rate ----
fig1a <- ggplot(bg_sf2) +
  geom_sf(aes(fill = call_rate), color = "white", size = 0.15) +
  scale_fill_gradientn(colors = c("#2c7bb6", "#abd9e9", "#fee090", "#d73027"),
                       name = "Call rate\n(per 1000 days)") +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.15, 0.25),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.3, "cm"),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(title = "a  Heat-related calls")

## ---- Figure 1b: Map of evap prevalence (for visual comparison) ----
fig1b <- ggplot(bg_sf2) +
  geom_sf(aes(fill = evap_prop), color = "white", size = 0.15) +
  scale_fill_gradientn(colors = c("#2c7bb6", "#abd9e9", "#fee090", "#d73027"),
                       name = "Evap.\nprevalence") +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.15, 0.25),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.3, "cm"),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(title = "b  Evap. cooler prevalence")

fig1 <- fig1a + fig1b + plot_layout(ncol = 2)

pdf(here("results", "P5_Figure1_call_rate_maps.pdf"), width = 10, height = 5)
print(fig1)
dev.off()

png(here("results", "P5_Figure1_call_rate_maps.png"), width = 10, height = 5, units = "in", res = 300)
print(fig1)
dev.off()

## ---- Figure 2: Scatter of call rate vs evap prevalence ----
rho_ev <- cor(bg$call_rate, bg$evap_prop, method = "spearman")
pval_ev <- cor.test(bg$call_rate, bg$evap_prop, method = "spearman")$p.value
ptext_ev <- ifelse(pval_ev < 0.001, "p < 0.001", paste0("p = ", signif(pval_ev, 2)))

fig2 <- ggplot(bg, aes(x = evap_prop, y = call_rate)) +
  geom_point(size = 1.2, alpha = 0.5, color = "grey30") +
  geom_smooth(method = "lm", se = TRUE, color = "#d73027",
              fill = "#fc8d59", alpha = 0.2, linewidth = 0.7) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("rho == ", round(rho_ev, 2)), parse = TRUE,
           hjust = 1.1, vjust = 1.5, size = 3) +
  annotate("text", x = Inf, y = Inf, label = ptext_ev,
           hjust = 1.1, vjust = 3, size = 2.5, color = "grey40") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(x = "Evap. cooler prevalence",
       y = "Heat-related call rate (per 1000 days)",
       title = "")

pdf(here("results", "P5_Figure2_scatter_calls_evap.pdf"), width = 5, height = 5)
print(fig2)
dev.off()

png(here("results", "P5_Figure2_scatter_calls_evap.png"), width = 5, height = 5, units = "in", res = 300)
print(fig2)
dev.off()

## ---- Figure 3: Daily calls on failure vs non-failure days ----
daily_calls$failure_label <- ifelse(daily_calls$failure_day == 1,
                                     "Failure day", "Non-failure day")
daily_calls$failure_label <- factor(daily_calls$failure_label,
                                     levels = c("Non-failure day", "Failure day"))

fig3 <- ggplot(daily_calls, aes(x = failure_label, y = calls, fill = failure_label)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8, outlier.alpha = 0.4) +
  scale_fill_manual(values = c("Non-failure day" = "#2c7bb6",
                                "Failure day" = "#d73027"), guide = "none") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(x = "", y = "Daily heat-related calls (city-wide)", title = "")

pdf(here("results", "P5_Figure3_failure_day_calls.pdf"), width = 5, height = 4)
print(fig3)
dev.off()

png(here("results", "P5_Figure3_failure_day_calls.png"), width = 5, height = 4, units = "in", res = 300)
print(fig3)
dev.off()

## ---- Figure 4: Call rate by evap quartile ----
bg$evap_q <- cut(bg$evap_prop,
                  breaks = quantile(bg$evap_prop, probs = 0:4/4),
                  labels = c("Q1\n(lowest)", "Q2", "Q3", "Q4\n(highest)"),
                  include.lowest = TRUE)

fig4 <- ggplot(bg, aes(x = evap_q, y = call_rate, fill = evap_q)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8, outlier.alpha = 0.4) +
  scale_fill_manual(values = c("#abd9e9", "#fee090", "#fc8d59", "#d73027"),
                    guide = "none") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(x = "Evap. cooler prevalence quartile",
       y = "Heat-related call rate (per 1000 days)", title = "")

pdf(here("results", "P5_Figure4_calls_by_quartile.pdf"), width = 6, height = 4)
print(fig4)
dev.off()

png(here("results", "P5_Figure4_calls_by_quartile.png"), width = 6, height = 4, units = "in", res = 300)
print(fig4)
dev.off()

