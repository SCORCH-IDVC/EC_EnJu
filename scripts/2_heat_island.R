library(here)
library(sf)
library(terra)
library(exactextractr)
library(ggplot2)
library(spdep)
library(patchwork)

###############################################################
#Paper 2: Urban Heat Island Interaction                      ##
#Are evaporative-cooler households concentrated in Tucson's   ##
#hottest microclimates, creating a double exposure?           ##
###############################################################

# ============================================================
# 1. LOAD BLOCK GROUP DATA (from Paper 1)
# ============================================================

tfiles <- list.files(here("data", "Q1 Data Shapefile"))

if (length(tfiles) == 0) {
  stop("Run Paper 1 script first, or place shapefile in data/Q1 Data Shapefile/")
} else {
  bg_sf <- st_read(here("data", "Q1 Data Shapefile", "pima_Q1_data.shp"))
  bg_sf <- bg_sf[bg_sf$med_inc != 0 & !is.na(bg_sf$med_inc), ]
  bg_sf <- bg_sf[!is.na(bg_sf$evp_prp), ]
  bg <- data.frame(bg_sf)
  
  ## Rename columns
  colnames(bg)[colnames(bg) == "geoid20"]  <- "GEOID"
  colnames(bg)[colnames(bg) == "evp_prp"]  <- "evap_prop"
  colnames(bg)[colnames(bg) == "med_inc"]  <- "med_income"
  colnames(bg)[colnames(bg) == "pct_mnr"]  <- "pct_minority"
  colnames(bg)[colnames(bg) == "ave_age"]  <- "med_year_built"
  colnames(bg)[colnames(bg) == "pct_rnt"]  <- "pct_renter"
  colnames(bg)[colnames(bg) == "pct_sfr"]  <- "pct_sfh"
  colnames(bg)[colnames(bg) == "covennt"]  <- "covenant"
  
  ## Centroids
  bg_sf <- st_as_sf(bg)
  bg_sf <- st_transform(bg_sf, 4326)
  coords <- st_coordinates(st_centroid(bg_sf))
  bg$lon <- coords[, 1]
  bg$lat <- coords[, 2]
  bg_sf <- st_make_valid(bg_sf)
}

cat("Block groups loaded:", nrow(bg), "\n")

# ============================================================
# 2. DOWNLOAD LST AND NDVI (no account needed)
# ============================================================
# LST: PRISM tmax (4km, monthly) via the prism R package.
#      Downloads directly from Oregon State. No login.
# NDVI: MODIS MOD13Q1 (250m, 16-day) via MODISTools.
#       Downloads from ORNL DAAC API. No login.


dir.create(here("data", "rasters"), recursive = TRUE, showWarnings = FALSE)

lst_path  <- here("data", "rasters", "tucson_lst_summer.tif")
ndvi_path <- here("data", "rasters", "tucson_ndvi_summer.tif")

## ---- LST from PRISM ----
if (!file.exists(lst_path)) {
  
  cat("=== Downloading PRISM tmax (Jun-Sep 2023) ===\n")
  library(prism)
  
  prism_dir <- here("data", "rasters", "prism_raw")
  dir.create(prism_dir, recursive = TRUE, showWarnings = FALSE)
  prism_set_dl_dir(prism_dir)
  
  ## Download monthly tmax for Jun, Jul, Aug, Sep 2023
  get_prism_monthlys(type = "tmax", years = 2023, mon = 6:9, keepZip = FALSE)
  
  ## List downloaded files and read as rasters
  pd <- prism_archive_ls()
  cat("PRISM files downloaded:", length(pd), "\n")
  
  ## Get file paths and stack
  bil_files <- pd_to_file(pd)
  lst_stack <- rast(bil_files)
  
  ## Crop to study area and compute summer mean tmax
  bg_4326 <- st_transform(bg_sf, 4326)
  roi <- ext(as.numeric(st_bbox(bg_4326)))
  lst_stack <- crop(lst_stack, roi)
  lst_stack <- mask(lst_stack, roi)
  
  lst_comp <- app(lst_stack, fun = mean, na.rm = TRUE)
  
  writeRaster(lst_comp, lst_path, overwrite = TRUE)
  cat("LST composite saved (PRISM tmax, Jun-Sep 2023, 4km)\n")
}

## ---- NDVI from MODIS via MODISTools ----
if (!file.exists(ndvi_path)) {
  
  cat("=== Downloading MODIS NDVI (Jun-Sep 2023) ===\n")
  library(MODISTools)
  
  ## Study area center and extent
  bg_4326 <- st_transform(bg_sf, 4326)
  bb <- as.numeric(st_bbox(bg_4326))
  center_lat <- mean(bb[2], bb[4])
  center_lon <- mean(bb[1], bb[3])
  
  ## km extent from center to edges (approximate)
  km_lr <- round((bb[3] - bb[1]) * 111 * cos(center_lat * pi / 180) / 2) + 2
  km_ab <- round((bb[4] - bb[2]) * 111 / 2) + 2
  
  ## Download MOD13Q1 NDVI (250m, 16-day composite)
  ndvi_raw <- mt_subset(
    product   = "MOD13Q1",
    band      = "250m_16_days_NDVI",
    lat       = center_lat,
    lon       = center_lon,
    km_lr     = km_lr,
    km_ab     = km_ab,
    start     = "2023-06-01",
    end       = "2023-09-30",
    site_name = "tucson",
    internal  = TRUE,
    progress  = TRUE
  )
  
  
  cat("NDVI pixels downloaded:", nrow(ndvi_raw), "\n")
  
  ## Scale factor: MODIS NDVI is stored as integer * 10000
  ndvi_raw$value <- ndvi_raw$value * 0.0001
  
  ## Filter out fill values (NDVI should be -0.2 to 1.0)
  ndvi_raw <- ndvi_raw[ndvi_raw$value >= -0.2 & ndvi_raw$value <= 1.0, ]
  
  ## Compute median NDVI per pixel across the summer
  ndvi_median <- aggregate(value ~ pixel + latitude + longitude,
                           data = ndvi_raw, FUN = median, na.rm = TRUE)
  
  ## Convert to spatial points, then rasterize
  ndvi_pts <- st_as_sf(ndvi_median, coords = c("longitude", "latitude"), crs = 4326)
  roi <- ext(bb[1], bb[3], bb[2], bb[4])
  r_template <- rast(roi, res = 0.002, crs = "EPSG:4326")  # ~250m grid
  ndvi_rast <- rasterize(vect(ndvi_pts), r_template, field = "value", fun = mean)
  
  writeRaster(ndvi_rast, ndvi_path, overwrite = TRUE)
  cat("NDVI raster saved (MODIS MOD13Q1, Jun-Sep 2023, 250m)\n")
}

cat("LST raster:", lst_path, "\n")
cat("NDVI raster:", ndvi_path, "\n")

# ============================================================
# 3. ZONAL STATISTICS
# ============================================================

lst_r  <- rast(lst_path)
ndvi_r <- rast(ndvi_path)

## Reproject block groups to match raster CRS
bg_proj <- st_transform(bg_sf, crs(lst_r))

## Extract mean LST and NDVI per block group
bg$mean_lst  <- exact_extract(lst_r, bg_proj, 'mean')
bg$mean_ndvi <- exact_extract(ndvi_r, bg_proj, 'mean')

cat("\n=== ZONAL STATISTICS ===\n")
cat("LST range:", round(range(bg$mean_lst, na.rm = TRUE), 1), "C\n")
cat("NDVI range:", round(range(bg$mean_ndvi, na.rm = TRUE), 3), "\n")

## Drop any block groups with NA raster values
bg <- bg[!is.na(bg$mean_lst) & !is.na(bg$mean_ndvi), ]
cat("Block groups with valid LST + NDVI:", nrow(bg), "\n")

# ============================================================
# 4. DESCRIPTIVE STATISTICS (Table 1)
# ============================================================

## Quartiles based on evap prevalence
bg$evap_q <- cut(bg$evap_prop,
                 breaks = quantile(bg$evap_prop, probs = 0:4/4),
                 labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"),
                 include.lowest = TRUE)

## Summary function
summarize_by_quartile <- function(x, group) {
  means <- tapply(x, group, mean, na.rm = TRUE)
  sds <- tapply(x, group, sd, na.rm = TRUE)
  ns <- tapply(x, group, length)
  ses <- sds / sqrt(ns)
  data.frame(quartile = names(means), mean = round(means, 2),
             sd = round(sds, 2), se = round(ses, 2), n = as.integer(ns))
}

## Build Table 1
vars <- c("mean_lst", "mean_ndvi", "med_income", "pct_renter")
var_labels <- c("Mean LST (C)", "Mean NDVI", "Median income ($)", "Renter (%)")

table1_list <- lapply(vars, function(v) {
  out <- summarize_by_quartile(bg[[v]], bg$evap_q)
  out$variable <- v
  out
})
table1 <- do.call(rbind, table1_list)
table1$var_label <- rep(var_labels, each = 4)

## Kruskal-Wallis tests
kw_tests <- sapply(vars, function(v) {
  kruskal.test(bg[[v]] ~ bg$evap_q)$p.value
})
names(kw_tests) <- var_labels

cat("\n=== TABLE 1: Block group characteristics by evap cooler quartile ===\n")
print(table1[, c("var_label", "quartile", "mean", "se", "n")])
cat("\nKruskal-Wallis p-values:\n")
print(round(kw_tests, 4))

write.csv(table1, here("results", "P2_Table1_quartile_LST_NDVI.csv"), row.names = FALSE)

# ============================================================
# 5. SPEARMAN CORRELATIONS
# ============================================================

cor_vars <- c("mean_lst", "mean_ndvi", "med_income", "pct_renter")
cor_labels <- c("Mean LST", "Mean NDVI", "Median income", "% Renter")

cor_results <- data.frame(
  variable = cor_labels,
  rho = sapply(cor_vars, function(v) cor(bg$evap_prop, bg[[v]], method = "spearman")),
  p = sapply(cor_vars, function(v) cor.test(bg$evap_prop, bg[[v]], method = "spearman")$p.value)
)
cor_results$rho <- round(cor_results$rho, 3)
cor_results$p <- signif(cor_results$p, 3)

cat("\n=== SPEARMAN CORRELATIONS with evap cooler prevalence ===\n")
print(cor_results)

write.csv(cor_results, here("results", "P2_TableS1_correlations.csv"), row.names = FALSE)

# ============================================================
# 6. BIVARIATE LISA: EVAP PREVALENCE x LST
# ============================================================

## Spatial weights
coords <- cbind(bg$lon, bg$lat)
nb <- knn2nb(knearneigh(coords, k = 5))
lw <- nb2listw(nb, style = "W")

## Standardize both variables
z_evap <- scale(bg$evap_prop)
z_lst  <- scale(bg$mean_lst)

## Cross-product for bivariate LISA
## High values = both high or both low; low values = mismatch
bg$cross_evap_lst <- as.numeric(z_evap * lag.listw(lw, z_lst))

## Local Moran's I on evap prevalence with LST as spatial lag
lisa_bi <- localmoran(z_evap, lw)
bg$lisa_bi_I <- lisa_bi[, 1]
bg$lisa_bi_p <- lisa_bi[, 5]

## Classify bivariate LISA clusters
## Using cross-product of local value x spatially lagged LST
lag_lst <- lag.listw(lw, z_lst)

bg$bi_cluster <- "Not significant"
sig <- bg$lisa_bi_p < 0.05
bg$bi_cluster[sig & z_evap > 0 & lag_lst > 0] <- "High evap / High LST"
bg$bi_cluster[sig & z_evap < 0 & lag_lst < 0] <- "Low evap / Low LST"
bg$bi_cluster[sig & z_evap > 0 & lag_lst < 0] <- "High evap / Low LST"
bg$bi_cluster[sig & z_evap < 0 & lag_lst > 0] <- "Low evap / High LST"

cat("\n=== BIVARIATE LISA CLUSTER COUNTS ===\n")
print(table(bg$bi_cluster))

## Double-exposure hotspots
n_double <- sum(bg$bi_cluster == "High evap / High LST")
cat("Double-exposure hotspots (High evap / High LST):", n_double, "\n")

# ============================================================
# 7. GLM: PREDICTORS OF EVAP COOLER PREVALENCE
# ============================================================

## Standardize predictors
bg$z_lst    <- scale(bg$mean_lst)
bg$z_ndvi   <- scale(bg$mean_ndvi)
bg$z_income <- scale(bg$med_income)
bg$z_renter <- scale(bg$pct_renter)

## Quasibinomial GLM
m1 <- glm(evap_prop ~ z_lst + z_ndvi + z_income + z_renter,
          family = quasibinomial, data = bg)

cat("\n=== TABLE 2: GLM RESULTS ===\n")
print(summary(m1))

## Extract coefficients
coef_table <- data.frame(
  variable = c("Intercept", "LST (z)", "NDVI (z)", "Income (z)", "Renter (z)"),
  estimate = round(coef(m1), 3),
  se = round(summary(m1)$coefficients[, 2], 3),
  p = signif(summary(m1)$coefficients[, 4], 3)
)
write.csv(coef_table, here("results", "P2_Table2_GLM_results.csv"), row.names = FALSE)

## Moran's I on residuals
moran_resid <- moran.test(residuals(m1), lw)
cat("\nMoran's I on GLM residuals:", round(moran_resid$estimate[1], 3), "\n")
cat("Moran p-value:", signif(moran_resid$p.value, 3), "\n")

## Spatial error model if needed
if (moran_resid$p.value < 0.05) {
  cat(">> Spatial autocorrelation detected. Fitting spatial error model.\n")
  library(spatialreg)
  m_spatial <- errorsarlm(evap_prop ~ z_lst + z_ndvi + z_income + z_renter,
                          data = bg, listw = lw)
  print(summary(m_spatial))
  
  se_coefs <- summary(m_spatial)$Coef
  coef_table_spatial <- data.frame(
    variable = rownames(se_coefs),
    estimate = round(se_coefs[, 1], 3),
    se = round(se_coefs[, 2], 3),
    p = signif(se_coefs[, 4], 3)
  )
  write.csv(coef_table_spatial, here("results", "P2_Table2b_spatial_error_model.csv"), row.names = FALSE)
}

# ============================================================
# 8. CHARACTERIZE DOUBLE-EXPOSURE HOTSPOTS (Table S2)
# ============================================================

hotspot <- bg[bg$bi_cluster == "High evap / High LST", ]
non_hotspot <- bg[bg$bi_cluster != "High evap / High LST", ]

## Compare demographics
compare_vars <- c("evap_prop", "mean_lst", "mean_ndvi", "med_income",
                  "pct_minority", "pct_renter", "med_year_built")
compare_labels <- c("Evap. prevalence", "Mean LST (C)", "Mean NDVI",
                    "Median income ($)", "Minority (%)", "Renter (%)", "Year built")

hotspot_table <- data.frame(
  variable = compare_labels,
  hotspot_mean = sapply(compare_vars, function(v) round(mean(hotspot[[v]], na.rm = TRUE), 2)),
  hotspot_se = sapply(compare_vars, function(v) round(sd(hotspot[[v]], na.rm = TRUE) / sqrt(nrow(hotspot)), 2)),
  city_mean = sapply(compare_vars, function(v) round(mean(bg[[v]], na.rm = TRUE), 2)),
  city_se = sapply(compare_vars, function(v) round(sd(bg[[v]], na.rm = TRUE) / sqrt(nrow(bg)), 2)),
  wilcox_p = sapply(compare_vars, function(v) {
    signif(wilcox.test(hotspot[[v]], non_hotspot[[v]])$p.value, 3)
  })
)

cat("\n=== TABLE S2: Double-exposure hotspot vs. city-wide ===\n")
print(hotspot_table)

write.csv(hotspot_table, here("results", "P2_TableS2_hotspot_demographics.csv"), row.names = FALSE)

# ============================================================
# 9. FIGURES
# ============================================================

## Merge new variables back to sf for mapping
bg_sf2 <- st_as_sf(bg)
if (is.na(st_crs(bg_sf2))) bg_sf2 <- st_set_crs(bg_sf2, 4326)
bg_sf2 <- st_make_valid(bg_sf2)

## ---- Figure 1: Bivariate choropleth (evap prevalence x LST) ----
## Create bivariate classes: 3x3 grid
bg$evap_cat <- cut(bg$evap_prop,
                   breaks = quantile(bg$evap_prop, probs = c(0, 1/3, 2/3, 1)),
                   labels = c("Low", "Mid", "High"), include.lowest = TRUE)
bg$lst_cat <- cut(bg$mean_lst,
                  breaks = quantile(bg$mean_lst, probs = c(0, 1/3, 2/3, 1)),
                  labels = c("Low", "Mid", "High"), include.lowest = TRUE)
bg$bivar_class <- paste(bg$evap_cat, bg$lst_cat, sep = " / ")

bg_sf2$bivar_class <- bg$bivar_class

## 3x3 bivariate palette (blue-to-red x light-to-dark)
bivar_pal <- c(
  "Low / Low"   = "#e8e8e8", "Low / Mid"   = "#e4acac", "Low / High"   = "#c85a5a",
  "Mid / Low"   = "#b0d5df", "Mid / Mid"   = "#ad9ea5", "Mid / High"   = "#985356",
  "High / Low"  = "#64acbe", "High / Mid"  = "#627f8c", "High / High"  = "#574249"
)

fig1 <- ggplot(bg_sf2) +
  geom_sf(aes(fill = bivar_class), color = "white", size = 0.15) +
  scale_fill_manual(values = bivar_pal, name = "Evap / LST",
                    guide = guide_legend(ncol = 3)) +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(title = "a")

## ---- Figure 1b: Bivariate legend (inset) ----
legend_df <- expand.grid(evap = c("Low", "Mid", "High"),
                         lst = c("Low", "Mid", "High"))
legend_df$fill <- bivar_pal[paste(legend_df$evap, legend_df$lst, sep = " / ")]
legend_df$evap <- factor(legend_df$evap, levels = c("Low", "Mid", "High"))
legend_df$lst  <- factor(legend_df$lst, levels = c("Low", "Mid", "High"))

fig1_legend <- ggplot(legend_df, aes(x = evap, y = lst, fill = fill)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_identity() +
  labs(x = "Evap. prevalence \u2192", y = "LST \u2192") +
  theme_minimal(base_size = 7) +
  theme(axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7, face = "bold"),
        plot.background = element_rect(fill = "white", color = NA))

## ---- Figure 2: LISA cluster map ----
bi_pal <- c("High evap / High LST" = "#d7191c",
            "Low evap / Low LST"    = "#2c7bb6",
            "High evap / Low LST"   = "#fdae61",
            "Low evap / High LST"   = "#abd9e9",
            "Not significant"       = "grey90")

bg_sf2$bi_cluster <- bg$bi_cluster

fig2 <- ggplot(bg_sf2) +
  geom_sf(aes(fill = bi_cluster), color = "white", size = 0.15) +
  scale_fill_manual(values = bi_pal, name = "Bivariate LISA") +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.text = element_text(size = 6.5),
        legend.title = element_text(size = 7),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(title = "b")

## ---- Supplementary Figure S1: Scatterplots ----
make_scatter <- function(xvar, xlab, panel_label) {
  rho <- cor(bg$evap_prop, bg[[xvar]], method = "spearman")
  pval <- cor.test(bg$evap_prop, bg[[xvar]], method = "spearman")$p.value
  ptext <- ifelse(pval < 0.001, "p < 0.001", paste0("p = ", signif(pval, 2)))
  label <- paste0("rho == ", round(rho, 2))
  
  ggplot(bg, aes_string(x = xvar, y = "evap_prop")) +
    geom_point(size = 1.2, alpha = 0.5, color = "grey30") +
    geom_smooth(method = "loess", se = TRUE, color = "#cb181d",
                fill = "#fb6a4a", alpha = 0.2, linewidth = 0.7) +
    annotate("text", x = Inf, y = Inf, label = label, parse = TRUE,
             hjust = 1.1, vjust = 1.5, size = 3) +
    annotate("text", x = Inf, y = Inf, label = ptext,
             hjust = 1.1, vjust = 3, size = 2.5, color = "grey40") +
    theme_minimal(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = 10, face = "bold")) +
    labs(x = xlab, y = "Evap. cooler prevalence", title = panel_label)
}

figS1a <- make_scatter("mean_lst", "Mean land surface temperature (C)", "a")
figS1b <- make_scatter("mean_ndvi", "Mean NDVI", "b")
figS1c <- make_scatter("med_income", "Median household income ($)", "c")
figS1d <- make_scatter("pct_renter", "Renter proportion", "d")

# ============================================================
# 10. EXPORT FIGURES
# ============================================================

dir.create(here("results"), recursive = TRUE, showWarnings = FALSE)

## ---- Figure 1: Bivariate choropleth + legend ----
fig1_final <- fig1 + inset_element(fig1_legend, left = 0.02, bottom = 0.02,
                                   right = 0.25, top = 0.25)

pdf(here("results", "P2_Figure1_bivariate_choropleth.pdf"), width = 7, height = 7)
print(fig1_final)
dev.off()

png(here("results", "P2_Figure1_bivariate_choropleth.png"), width = 7, height = 7, units = "in", res = 300)
print(fig1_final)
dev.off()

## ---- Figure 2: LISA cluster map ----
pdf(here("results", "P2_Figure2_LISA_double_exposure.pdf"), width = 7, height = 7)
print(fig2)
dev.off()

png(here("results", "P2_Figure2_LISA_double_exposure.png"), width = 7, height = 7, units = "in", res = 300)
print(fig2)
dev.off()

## ---- Figure S1: Scatterplots ----
figS1 <- (figS1a + figS1b) / (figS1c + figS1d)

pdf(here("results", "P2_FigureS1_scatterplots.pdf"), width = 8, height = 7)
print(figS1)
dev.off()

png(here("results", "P2_FigureS1_scatterplots.png"), width = 8, height = 7, units = "in", res = 300)
print(figS1)
dev.off()

cat("\n=== DONE ===\n")
cat("Main text outputs:\n")
cat("  P2_Table1_quartile_LST_NDVI.csv\n")
cat("  P2_Table2_GLM_results.csv\n")
cat("  P2_Figure1_bivariate_choropleth.pdf\n")
cat("  P2_Figure2_LISA_double_exposure.pdf\n")
cat("\nSupplement:\n")
cat("  P2_TableS1_correlations.csv\n")
cat("  P2_TableS2_hotspot_demographics.csv\n")
cat("  P2_Table2b_spatial_error_model.csv (if spatial autocorrelation)\n")
cat("  P2_FigureS1_scatterplots.pdf\n")