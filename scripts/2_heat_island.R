library(here)
library(sf)
library(terra)
library(exactextractr)
library(ggplot2)
library(spdep)
library(patchwork)
library(prism)

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

cat("Block groups loaded:", nrow(bg), "\n")

# 2. DOWNLOAD SUMMER TEMPERATURE (PRISM tmax)

dir.create(here("data", "rasters"), recursive = TRUE, showWarnings = FALSE)
lst_path <- here("data", "rasters", "tucson_tmax_summer.tif")

if (!file.exists(lst_path)) {
  
  cat("=== Downloading PRISM tmax (Jun-Sep 2023) ===\n")
  prism_dir <- here("data", "rasters", "prism_raw")
  dir.create(prism_dir, recursive = TRUE, showWarnings = FALSE)
  prism_set_dl_dir(prism_dir)
  
  ## Download monthly tmax for Jun, Jul, Aug, Sep 2023
  get_prism_monthlys(type = "tmax", years = 2023, mon = 6:9, keepZip = FALSE)
  
  ## List downloaded files and stack
  pd <- prism_archive_ls()
  cat("PRISM files downloaded:", length(pd), "\n")
  
  bil_files <- pd_to_file(pd)
  tmax_stack <- rast(bil_files)
  
  ## Crop to study area and compute summer mean tmax
  bg_4326 <- st_transform(bg_sf, 4326)
  roi <- ext(as.numeric(st_bbox(bg_4326)))
  tmax_stack <- crop(tmax_stack, roi)
  tmax_comp <- app(tmax_stack, fun = mean, na.rm = TRUE)
  
  writeRaster(tmax_comp, lst_path, overwrite = TRUE)
  cat("Summer tmax composite saved (PRISM, Jun-Sep 2023, 4km)\n")
}

cat("Temperature raster:", lst_path, "\n")

# 3. ZONAL STATISTICS

tmax_r <- rast(lst_path)

## Reproject block groups to match raster CRS
bg_proj <- st_transform(bg_sf, crs(tmax_r))

## Extract mean summer tmax per block group
bg$mean_tmax <- exact_extract(tmax_r, bg_proj, 'mean')

cat("\n=== ZONAL STATISTICS ===\n")
cat("Tmax range:", round(range(bg$mean_tmax, na.rm = TRUE), 1), "C\n")

## Drop block groups with NA
bg <- bg[!is.na(bg$mean_tmax), ]
cat("Block groups with valid tmax:", nrow(bg), "\n")

# 4. TABLE 1: BLOCK GROUP CHARACTERISTICS BY EVAP QUARTILE

bg$evap_q <- cut(bg$evap_prop,
                 breaks = quantile(bg$evap_prop, probs = 0:4/4),
                 labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"),
                 include.lowest = TRUE)

summarize_by_quartile <- function(x, group) {
  means <- tapply(x, group, mean, na.rm = TRUE)
  sds <- tapply(x, group, sd, na.rm = TRUE)
  ns <- tapply(x, group, length)
  ses <- sds / sqrt(ns)
  data.frame(quartile = names(means), mean = round(means, 2),
             sd = round(sds, 2), se = round(ses, 2), n = as.integer(ns))
}

vars <- c("mean_tmax", "med_income", "pct_minority", "pct_renter")
var_labels <- c("Mean summer tmax (C)", "Median income ($)",
                "Minority (%)", "Renter (%)")

table1_list <- lapply(vars, function(v) {
  out <- summarize_by_quartile(bg[[v]], bg$evap_q)
  out$variable <- v
  out
})
table1 <- do.call(rbind, table1_list)
table1$var_label <- rep(var_labels, each = 4)

kw_tests <- sapply(vars, function(v) {
  kruskal.test(bg[[v]] ~ bg$evap_q)$p.value
})
names(kw_tests) <- var_labels

cat("\n=== TABLE 1 ===\n")
print(table1[, c("var_label", "quartile", "mean", "se", "n")])
cat("\nKruskal-Wallis p-values:\n")
print(round(kw_tests, 4))

write.csv(table1, here("results", "P2_Table1_quartile_tmax.csv"), row.names = FALSE)

# 5. SPEARMAN CORRELATIONS

cor_vars <- c("mean_tmax", "med_income", "pct_minority", "pct_renter")
cor_labels <- c("Mean summer tmax", "Median income", "% Minority", "% Renter")

cor_results <- data.frame(
  variable = cor_labels,
  rho = sapply(cor_vars, function(v) cor(bg$evap_prop, bg[[v]], method = "spearman")),
  p = sapply(cor_vars, function(v) cor.test(bg$evap_prop, bg[[v]], method = "spearman")$p.value)
)
cor_results$rho <- round(cor_results$rho, 3)
cor_results$p <- signif(cor_results$p, 3)

cat("\n=== SPEARMAN CORRELATIONS ===\n")
print(cor_results)

write.csv(cor_results, here("results", "P2_TableS1_correlations.csv"), row.names = FALSE)

# 6. BIVARIATE LISA: EVAP PREVALENCE x TEMPERATURE

## Spatial weights
coords <- cbind(bg$lon, bg$lat)
nb <- knn2nb(knearneigh(coords, k = 5))
lw <- nb2listw(nb, style = "W")

## Standardize
z_evap <- scale(bg$evap_prop)
z_tmax <- scale(bg$mean_tmax)

## Local Moran's I on evap prevalence
lisa_bi <- localmoran(as.numeric(z_evap), lw)
bg$lisa_bi_I <- lisa_bi[, 1]
bg$lisa_bi_p <- lisa_bi[, 5]

## Classify bivariate clusters using local evap and spatially lagged tmax
lag_tmax <- lag.listw(lw, z_tmax)

bg$bi_cluster <- "Not significant"
sig <- bg$lisa_bi_p < 0.05
bg$bi_cluster[sig & z_evap > 0 & lag_tmax > 0] <- "High evap / Hot"
bg$bi_cluster[sig & z_evap < 0 & lag_tmax < 0] <- "Low evap / Cool"
bg$bi_cluster[sig & z_evap > 0 & lag_tmax < 0] <- "High evap / Cool"
bg$bi_cluster[sig & z_evap < 0 & lag_tmax > 0] <- "Low evap / Hot"

cat("\n=== BIVARIATE LISA CLUSTERS ===\n")
print(table(bg$bi_cluster))
cat("Double-exposure hotspots (High evap / Hot):", sum(bg$bi_cluster == "High evap / Hot"), "\n")

# 7. GLM: PREDICTORS OF EVAP COOLER PREVALENCE

bg$z_tmax   <- scale(bg$mean_tmax)
bg$z_income <- scale(bg$med_income)
bg$z_renter <- scale(bg$pct_renter)

## Quasibinomial GLM
m1 <- glm(evap_prop ~ z_tmax + z_income + z_renter,
          family = quasibinomial, data = bg)

cat("\n=== TABLE 2: GLM RESULTS ===\n")
print(summary(m1))

coef_table <- data.frame(
  variable = c("Intercept", "Summer tmax (z)", "Income (z)", "Renter (z)"),
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
  m_spatial <- errorsarlm(evap_prop ~ z_tmax + z_income + z_renter,
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

# 8. CHARACTERIZE DOUBLE-EXPOSURE HOTSPOTS

hotspot <- bg[bg$bi_cluster == "High evap / Hot", ]
non_hotspot <- bg[bg$bi_cluster != "High evap / Hot", ]

compare_vars <- c("evap_prop", "mean_tmax", "med_income",
                  "pct_minority", "pct_renter", "med_year_built")
compare_labels <- c("Evap. prevalence", "Mean summer tmax (C)",
                    "Median income ($)", "Minority (%)",
                    "Renter (%)", "Year built")

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

cat("\n=== TABLE S2: Double-exposure hotspots vs. city-wide ===\n")
print(hotspot_table)

write.csv(hotspot_table, here("results", "P2_TableS2_hotspot_demographics.csv"), row.names = FALSE)

# 9. FIGURES

dir.create(here("results"), recursive = TRUE, showWarnings = FALSE)

bg_sf2 <- st_as_sf(bg)
if (is.na(st_crs(bg_sf2))) bg_sf2 <- st_set_crs(bg_sf2, 4326)
bg_sf2 <- st_make_valid(bg_sf2)

## ---- Figure 1: Bivariate choropleth (evap prevalence x tmax) ----
bg$evap_cat <- cut(bg$evap_prop,
                   breaks = quantile(bg$evap_prop, probs = c(0, 1/3, 2/3, 1)),
                   labels = c("Low", "Mid", "High"), include.lowest = TRUE)
bg$tmax_cat <- cut(bg$mean_tmax,
                   breaks = quantile(bg$mean_tmax, probs = c(0, 1/3, 2/3, 1)),
                   labels = c("Cool", "Mid", "Hot"), include.lowest = TRUE)
bg$bivar_class <- paste(bg$evap_cat, bg$tmax_cat, sep = " / ")

bg_sf2$bivar_class <- bg$bivar_class

bivar_pal <- c(
  "Low / Cool"  = "#e8e8e8", "Low / Mid"  = "#e4acac", "Low / Hot"  = "#c85a5a",
  "Mid / Cool"  = "#b0d5df", "Mid / Mid"  = "#ad9ea5", "Mid / Hot"  = "#985356",
  "High / Cool" = "#64acbe", "High / Mid" = "#627f8c", "High / Hot" = "#574249"
)

fig1 <- ggplot(bg_sf2) +
  geom_sf(aes(fill = bivar_class), color = "white", size = 0.15) +
  scale_fill_manual(values = bivar_pal, name = "Evap / Temp",
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

## Bivariate legend inset
legend_df <- expand.grid(evap = c("Low", "Mid", "High"),
                         tmax = c("Cool", "Mid", "Hot"))
legend_df$fill <- bivar_pal[paste(legend_df$evap, legend_df$tmax, sep = " / ")]
legend_df$evap <- factor(legend_df$evap, levels = c("Low", "Mid", "High"))
legend_df$tmax <- factor(legend_df$tmax, levels = c("Cool", "Mid", "Hot"))

fig1_legend <- ggplot(legend_df, aes(x = evap, y = tmax, fill = fill)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_identity() +
  labs(x = "Evap. prevalence \u2192", y = "Temperature \u2192") +
  theme_minimal(base_size = 7) +
  theme(axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7, face = "bold"),
        plot.background = element_rect(fill = "white", color = NA))

## ---- Figure 2: LISA cluster map ----
bi_pal <- c("High evap / Hot"  = "#d7191c",
            "Low evap / Cool"  = "#2c7bb6",
            "High evap / Cool" = "#fdae61",
            "Low evap / Hot"   = "#abd9e9",
            "Not significant"  = "grey90")

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

## ---- Figure S1: Scatterplots ----
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

figS1a <- make_scatter("mean_tmax", "Mean summer tmax (C)", "a")
figS1b <- make_scatter("med_income", "Median household income ($)", "b")
figS1c <- make_scatter("pct_renter", "Renter proportion", "c")

# 10. EXPORT FIGURES

## Figure 1: Bivariate choropleth + inset legend
fig1_final <- fig1 + inset_element(fig1_legend, left = 0.02, bottom = 0.02,
                                   right = 0.25, top = 0.25)

pdf(here("results", "P2_Figure1_bivariate_choropleth.pdf"), width = 7, height = 7)
print(fig1_final)
dev.off()

png(here("results", "P2_Figure1_bivariate_choropleth.png"), width = 7, height = 7, units = "in", res = 300)
print(fig1_final)
dev.off()

## Figure 2: LISA cluster map
pdf(here("results", "P2_Figure2_LISA_double_exposure.pdf"), width = 7, height = 7)
print(fig2)
dev.off()

png(here("results", "P2_Figure2_LISA_double_exposure.png"), width = 7, height = 7, units = "in", res = 300)
print(fig2)
dev.off()

## Figure S1: Scatterplots (3 panels)
figS1 <- figS1a + figS1b + figS1c + plot_layout(ncol = 3)

pdf(here("results", "P2_FigureS1_scatterplots.pdf"), width = 10, height = 4)
print(figS1)
dev.off()

png(here("results", "P2_FigureS1_scatterplots.png"), width = 10, height = 4, units = "in", res = 300)
print(figS1)
dev.off()
