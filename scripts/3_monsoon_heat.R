library(here)
library(sf)
library(ggplot2)
library(spdep)
library(patchwork)
library(splines)

# ============================================================
# 1. LOAD BLOCK GROUP DATA
# ============================================================

tfiles <- list.files(here("data", "Q1 Data Shapefile"))

if (length(tfiles) == 0) {
  stop("Run Paper 1 script first, or place shapefile in data/Q1 Data Shapefile/")
} else {
  bg_sf <- st_read(here("data", "Q1 Data Shapefile", "pima_Q1_data.shp"))
  bg_sf <- bg_sf[bg_sf$med_inc != 0 & !is.na(bg_sf$med_inc), ]
  bg_sf <- bg_sf[!is.na(bg_sf$evp_prp), ]
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
}

cat("Block groups loaded:", nrow(bg), "\n")

# ============================================================
# 2. DOWNLOAD HOURLY WEATHER DATA (AZMET + NWS Tucson)
# ============================================================
# We can also download data from the AZMET raw data archive:
# https://cals.arizona.edu/azmet/
#
# For now, we pull NWS Tucson (TUS) via the Iowa Environmental
# Mesonet (IEM) ASOS archive.

dir.create(here("data", "weather"), recursive = TRUE, showWarnings = FALSE)
wx_path <- here("data", "weather", "tucson_hourly.csv")

if (!file.exists(wx_path)) {
  
  cat("=== Downloading hourly weather data ===\n")
  
  ## NWS Tucson International Airport (KTUS)
  
  years <- 2018:2026
  wx_list <- list()
  
  for (yr in years) {
    url <- sprintf(
      "https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?station=TUS&data=tmpf&data=relh&tz=America/Phoenix&format=onlycomma&latlon=no&elev=no&missing=M&trace=T&direct=no&report_type=3&year1=%d&month1=5&day1=1&year2=%d&month2=10&day2=1",
      yr, yr
    )
    tmp <- tempfile(fileext = ".csv")
    result <- try(download.file(url, tmp, mode = "w", quiet = TRUE), silent = TRUE)
    
    if (!inherits(result, "try-error") && file.size(tmp) > 500) {
      d <- read.csv(tmp, stringsAsFactors = FALSE)
      if (nrow(d) > 0 && "tmpf" %in% colnames(d)) {
        wx_list[[length(wx_list) + 1]] <- d
        cat("  IEM", yr, "downloaded:", nrow(d), "rows\n")
      }
    }
  }
  
  if (length(wx_list) > 0) {
    wx <- do.call(rbind, wx_list)
    
    ## Clean columns
    wx$tmpf <- as.numeric(wx$tmpf)
    wx$relh <- as.numeric(wx$relh)
    wx$datetime <- as.POSIXct(wx$valid, format = "%Y-%m-%d %H:%M", tz = "America/Phoenix")
    wx$date <- as.Date(wx$datetime, tz = "America/Phoenix")
    wx$hour <- as.integer(format(wx$datetime, "%H"))
    wx$year <- as.integer(format(wx$datetime, "%Y"))
    wx$month <- as.integer(format(wx$datetime, "%m"))
    wx$doy <- as.integer(format(wx$datetime, "%j"))
    
    ## Remove rows with missing temp or RH
    wx <- wx[!is.na(wx$tmpf) & !is.na(wx$relh), ]
    
    ## Convert temp to Celsius for consistency
    wx$temp_c <- (wx$tmpf - 32) * 5 / 9
    
    write.csv(wx, wx_path, row.names = FALSE)
    cat("Weather data saved:", nrow(wx), "hourly observations\n")
  } else {
    stop("Weather download failed. Check network connection.")
  }
}

wx <- read.csv(wx_path, stringsAsFactors = FALSE)
wx$datetime <- as.POSIXct(wx$datetime, tz = "America/Phoenix")
wx$date <- as.Date(wx$date)
cat("Weather observations loaded:", nrow(wx), "\n")

# ============================================================
# 3. COOLER FAILURE DAY IDENTIFICATION
# ============================================================
# We flag hours where RH > 30% AND temperature > 95°F
# co-occur. These are the hours when a swamp cooler
# cannot meaningfully cool a home.

wx$failure <- wx$relh > 30 & wx$tmpf > 95

## Aggregate to daily level
daily <- aggregate(
  cbind(failure_hours = failure, max_temp = tmpf, max_rh = relh) ~ date + year + month + doy,
  data = wx,
  FUN = function(x) c(sum(x), max(x), max(x))[1]
)

daily_fail <- aggregate(failure ~ date + year + month + doy, data = wx, FUN = sum)
colnames(daily_fail)[5] <- "failure_hours"

daily_tmax <- aggregate(tmpf ~ date, data = wx, FUN = max, na.rm = TRUE)
colnames(daily_tmax)[2] <- "tmax_f"

daily_rhmax <- aggregate(relh ~ date, data = wx, FUN = max, na.rm = TRUE)
colnames(daily_rhmax)[2] <- "rh_max"

daily <- merge(daily_fail, daily_tmax, by = "date")
daily <- merge(daily, daily_rhmax, by = "date")

## Binary: is this a failure day? (at least 1 hour of co-occurrence)
daily$failure_day <- as.integer(daily$failure_hours > 0)

cat("\n=== COOLER FAILURE SUMMARY ===\n")
cat("Total days analyzed:", nrow(daily), "\n")
cat("Total failure days:", sum(daily$failure_day), "\n")
cat("Failure rate:", round(mean(daily$failure_day) * 100, 1), "%\n")

# ============================================================
# 4. TABLE 1: COOLER FAILURE CLIMATOLOGY
# ============================================================

## Assign season: pre-monsoon (May-Jun), monsoon (Jul-Sep), post-monsoon (Oct)
daily$season <- "Pre-monsoon"
daily$season[daily$month >= 7 & daily$month <= 9] <- "Monsoon"
daily$season[daily$month == 10] <- "Post-monsoon"
daily$season <- factor(daily$season, levels = c("Pre-monsoon", "Monsoon", "Post-monsoon"))

## Annual summary by season
annual_season <- aggregate(
  cbind(failure_days = failure_day, total_failure_hours = failure_hours) ~ year + season,
  data = daily, FUN = sum
)

## Longest consecutive failure streak per year-season
streak_fun <- function(x) {
  if (sum(x) == 0) return(0)
  r <- rle(x)
  max(r$lengths[r$values == 1])
}

streaks <- aggregate(failure_day ~ year + season, data = daily, FUN = streak_fun)
colnames(streaks)[3] <- "longest_streak"

annual_season <- merge(annual_season, streaks, by = c("year", "season"))

## Mean across years
clim <- aggregate(
  cbind(failure_days, total_failure_hours, longest_streak) ~ season,
  data = annual_season, FUN = mean
)
clim_se <- aggregate(
  cbind(failure_days, total_failure_hours, longest_streak) ~ season,
  data = annual_season,
  FUN = function(x) sd(x) / sqrt(length(x))
)
colnames(clim_se)[2:4] <- paste0(colnames(clim_se)[2:4], "_se")

table1 <- merge(clim, clim_se, by = "season")
table1[, 2:7] <- round(table1[, 2:7], 1)

cat("\n=== TABLE 1: Cooler failure climatology ===\n")
print(table1)

write.csv(table1, here("results", "P3_Table1_failure_climatology.csv"), row.names = FALSE)

## Also save full annual breakdown for supplement
write.csv(annual_season, here("results", "P3_TableS1_annual_season_breakdown.csv"), row.names = FALSE)

# ============================================================
# 5. COMPOUND EXPOSURE INDEX
# ============================================================
# Compound exposure = evap cooler prevalence x failure frequency.
# Failure frequency is city-wide (one weather station), so the
# spatial variation comes entirely from cooler prevalence.
# We compute mean annual failure days across all years,
# then multiply by block group evap prevalence.

mean_failure_days <- mean(aggregate(failure_day ~ year, data = daily, FUN = sum)$failure_day)
cat("\nMean annual failure days:", round(mean_failure_days, 1), "\n")

## Compound exposure: proportion of households affected x days of failure
bg$compound_exposure <- bg$evap_prop * mean_failure_days

## Quartiles
bg$exposure_q <- cut(bg$compound_exposure,
                     breaks = quantile(bg$compound_exposure, probs = 0:4/4),
                     labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"),
                     include.lowest = TRUE)

cat("Compound exposure range:", round(range(bg$compound_exposure), 1), "\n")

# ============================================================
# 6. TABLE 2: DEMOGRAPHICS OF HIGH VS LOW EXPOSURE
# ============================================================

high_exp <- bg[bg$exposure_q == "Q4 (highest)", ]
low_exp  <- bg[bg$exposure_q == "Q1 (lowest)", ]

compare_vars <- c("evap_prop", "med_income", "pct_minority", "pct_renter", "med_year_built")
compare_labels <- c("Evap. prevalence", "Median income ($)", "Minority (%)",
                    "Renter (%)", "Year built")

table2 <- data.frame(
  variable = compare_labels,
  high_mean = sapply(compare_vars, function(v) round(mean(high_exp[[v]], na.rm = TRUE), 2)),
  high_se = sapply(compare_vars, function(v) round(sd(high_exp[[v]], na.rm = TRUE) / sqrt(nrow(high_exp)), 2)),
  low_mean = sapply(compare_vars, function(v) round(mean(low_exp[[v]], na.rm = TRUE), 2)),
  low_se = sapply(compare_vars, function(v) round(sd(low_exp[[v]], na.rm = TRUE) / sqrt(nrow(low_exp)), 2)),
  wilcox_p = sapply(compare_vars, function(v) {
    signif(wilcox.test(high_exp[[v]], low_exp[[v]])$p.value, 3)
  })
)

cat("\n=== TABLE 2: High vs. low compound exposure ===\n")
print(table2)

write.csv(table2, here("results", "P3_Table2_demographics_by_exposure.csv"), row.names = FALSE)

# ============================================================
# 7. SPATIAL WEIGHTS AND MORAN'S I ON COMPOUND EXPOSURE
# ============================================================

coords <- cbind(bg$lon, bg$lat)
nb <- knn2nb(knearneigh(coords, k = 5))
lw <- nb2listw(nb, style = "W")

moran_exp <- moran.test(bg$compound_exposure, lw)
cat("\nMoran's I on compound exposure:", round(moran_exp$estimate[1], 3), "\n")
cat("p-value:", signif(moran_exp$p.value, 3), "\n")

# ============================================================
# 8. FIGURES
# ============================================================

dir.create(here("results"), recursive = TRUE, showWarnings = FALSE)

## ---- Figure 1: Heatmap of hourly conditions across summer ----
## x = day of year, y = hour of day, fill = co-occurrence of heat + humidity
## Failure zone shaded

## Aggregate to mean conditions per DOY x hour (across years)
heatmap_data <- aggregate(
  cbind(mean_temp = tmpf, mean_rh = relh, fail_prop = failure) ~ doy + hour,
  data = wx, FUN = mean, na.rm = TRUE
)

## Panel a: temperature
fig1a <- ggplot(heatmap_data, aes(x = doy, y = hour, fill = mean_temp)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("#2c7bb6", "#abd9e9", "#fee090", "#d73027"),
                       name = "Temp (F)") +
  scale_y_continuous(breaks = seq(0, 23, 4)) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.3, "cm"),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(x = "Day of year", y = "Hour", title = "a")

## Panel b: relative humidity
fig1b <- ggplot(heatmap_data, aes(x = doy, y = hour, fill = mean_rh)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("#f7f7f7", "#74add1", "#313695"),
                       name = "RH (%)") +
  scale_y_continuous(breaks = seq(0, 23, 4)) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.3, "cm"),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(x = "Day of year", y = "Hour", title = "b")

## Panel c: failure probability (proportion of years with failure at this DOY x hour)
fig1c <- ggplot(heatmap_data, aes(x = doy, y = hour, fill = fail_prop)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("#f7f7f7", "#fc8d59", "#b30000"),
                       name = "Failure\nprobability",
                       limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 23, 4)) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.3, "cm"),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(x = "Day of year", y = "Hour", title = "c")

fig1 <- fig1a / fig1b / fig1c
pdf(here("results", "P3_Figure1_heatmap.pdf"), width = 8, height = 10)
print(fig1)
dev.off()

png(here("results", "P3_Figure1_heatmap.png"), width = 8, height = 10, units = "in", res = 300)
print(fig1)
dev.off()

## ---- Figure 2: Compound exposure map ----
bg_sf2 <- st_as_sf(bg)
if (is.na(st_crs(bg_sf2))) bg_sf2 <- st_set_crs(bg_sf2, 4326)
bg_sf2 <- st_make_valid(bg_sf2)

fig2 <- ggplot(bg_sf2) +
  geom_sf(aes(fill = compound_exposure), color = "white", size = 0.15) +
  scale_fill_gradientn(colors = c("#2c7bb6", "#abd9e9", "#fee090", "#d73027"),
                       name = "Compound\nexposure") +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.15, 0.25),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.3, "cm"),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(title = "")

pdf(here("results", "P3_Figure2_compound_exposure_map.pdf"), width = 7, height = 7)
print(fig2)
dev.off()

png(here("results", "P3_Figure2_compound_exposure_map.png"), width = 7, height = 7, units = "in", res = 300)
print(fig2)
dev.off()

## ---- Figure S1: Failure days time series ----
annual_total <- aggregate(failure_day ~ year, data = daily, FUN = sum)
colnames(annual_total)[2] <- "failure_days"

figS1 <- ggplot(annual_total, aes(x = year, y = failure_days)) +
  geom_col(fill = "#d73027", alpha = 0.7, width = 0.6) +
  geom_hline(yintercept = mean(annual_total$failure_days),
             linetype = "dashed", color = "grey40") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(x = "Year", y = "Total failure days (May-Sep)", title = "")

pdf(here("results", "P3_FigureS1_annual_failure_days.pdf"), width = 6, height = 4)
print(figS1)
dev.off()

png(here("results", "P3_FigureS1_annual_failure_days.png"), width = 6, height = 4, units = "in", res = 300)
print(figS1)
dev.off()

## ---- Figure S2: Scatter of evap prevalence vs compound exposure ----
figS2 <- ggplot(bg, aes(x = evap_prop, y = compound_exposure)) +
  geom_point(size = 1.2, alpha = 0.5, color = "grey30") +
  geom_smooth(method = "lm", se = TRUE, color = "#d73027",
              fill = "#fc8d59", alpha = 0.2, linewidth = 0.7) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(x = "Evap. cooler prevalence",
       y = "Compound exposure (prevalence x failure days)",
       title = "")

pdf(here("results", "P3_FigureS2_scatter_compound.pdf"), width = 5, height = 5)
print(figS2)
dev.off()

png(here("results", "P3_FigureS2_scatter_compound.png"), width = 5, height = 5, units = "in", res = 300)
print(figS2)
dev.off()
