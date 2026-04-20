library(here)
library(sf)
library(ggplot2)
library(spdep)
library(patchwork)
# 1. LOAD ALL RESULTS
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

# Q1. COOLING ACCESS DESERTS
# Where are the neighborhoods with the highest evaporative
# cooler prevalence AND the lowest incomes? These are the
# places where households cannot afford to transition to
# refrigerated AC. They are stuck.
#
# A cooling access desert is a block group where:
#   - evap prevalence > 75th percentile
#   - median income < 25th percentile
#
# These neighborhoods represent the floor of the
# cooling access hierarchy.

q75_evap <- quantile(bg$evap_prop, 0.75)
q25_income <- quantile(bg$med_income, 0.25)

bg$cooling_desert <- bg$evap_prop >= q75_evap & bg$med_income <= q25_income

cat("\n=== Q1: COOLING ACCESS DESERTS ===\n")
cat("Threshold: evap >=", round(q75_evap, 3), "AND income <=", round(q25_income), "\n")
cat("Cooling deserts:", sum(bg$cooling_desert), "block groups\n")
cat("Population in deserts:",
    round(sum(bg$cooling_desert) / nrow(bg) * 100, 1), "% of block groups\n")

desert <- bg[bg$cooling_desert, ]
non_desert <- bg[!bg$cooling_desert, ]

desert_compare <- data.frame(
  variable = c("Evap. prevalence", "Median income ($)", "Minority (%)",
               "Renter (%)", "Year built"),
  desert_mean = sapply(c("evap_prop", "med_income", "pct_minority", "pct_renter", "med_year_built"),
                        function(v) round(mean(desert[[v]], na.rm = TRUE), 2)),
  city_mean = sapply(c("evap_prop", "med_income", "pct_minority", "pct_renter", "med_year_built"),
                      function(v) round(mean(bg[[v]], na.rm = TRUE), 2)),
  wilcox_p = sapply(c("evap_prop", "med_income", "pct_minority", "pct_renter", "med_year_built"),
                     function(v) signif(wilcox.test(desert[[v]], non_desert[[v]])$p.value, 3))
)

print(desert_compare)
write.csv(desert_compare, here("results", "P7_Table1_cooling_deserts.csv"), row.names = FALSE)

# Q2. TRANSITION COST BURDEN
# What would it cost to replace evaporative coolers with
# refrigerated AC in the most vulnerable block groups?
#
# We estimate the number of evap-cooled households per
# block group and multiply by the average cost of an
# AC unit installation ($3,500-$7,000 depending on
# whether ductwork is needed).
#
# This gives us a dollar figure for the policy conversation.

## Approximate total housing units from ACS (not in our data,
## so we estimate from block group area and Tucson density)
## Conservative estimate: 500 housing units per block group average
est_units_per_bg <- 500

bg$est_evap_units <- round(bg$evap_prop * est_units_per_bg)

## Cost estimates
cost_low <- 3500    # window/wall unit + install
cost_mid <- 5500    # mini-split
cost_high <- 7500   # central AC with ductwork

bg$transition_cost_low <- bg$est_evap_units * cost_low
bg$transition_cost_mid <- bg$est_evap_units * cost_mid
bg$transition_cost_high <- bg$est_evap_units * cost_high

desert$est_evap_units <- round(desert$evap_prop * est_units_per_bg)
desert$transition_cost_low <- desert$est_evap_units * cost_low
desert$transition_cost_mid <- desert$est_evap_units * cost_mid
desert$transition_cost_high <- desert$est_evap_units * cost_high


## City-wide totals
cat("\n=== Q2: TRANSITION COST ESTIMATES ===\n")
cat("Estimated total evap-cooled units:", sum(bg$est_evap_units), "\n")
cat("City-wide transition cost:\n")
cat("  Low ($3,500/unit):  $", formatC(sum(bg$transition_cost_low), format = "d", big.mark = ","), "\n")
cat("  Mid ($5,500/unit):  $", formatC(sum(bg$transition_cost_mid), format = "d", big.mark = ","), "\n")
cat("  High ($7,500/unit): $", formatC(sum(bg$transition_cost_high), format = "d", big.mark = ","), "\n")

## Priority block groups (cooling deserts only)
cat("\nCooling desert transition cost:\n")
cat("  Units:", sum(desert$est_evap_units), "\n")
cat("  Low:  $", formatC(sum(desert$transition_cost_low), format = "d", big.mark = ","), "\n")
cat("  Mid:  $", formatC(sum(desert$transition_cost_mid), format = "d", big.mark = ","), "\n")
cat("  High: $", formatC(sum(desert$transition_cost_high), format = "d", big.mark = ","), "\n")

cost_summary <- data.frame(
  scope = c("City-wide", "Cooling deserts only"),
  units = c(sum(bg$est_evap_units), sum(desert$est_evap_units)),
  cost_low_M = round(c(sum(bg$transition_cost_low), sum(desert$transition_cost_low)) / 1e6, 1),
  cost_mid_M = round(c(sum(bg$transition_cost_mid), sum(desert$transition_cost_mid)) / 1e6, 1),
  cost_high_M = round(c(sum(bg$transition_cost_high), sum(desert$transition_cost_high)) / 1e6, 1)
)
write.csv(cost_summary, here("results", "P7_Table2_transition_costs.csv"), row.names = FALSE)

# Q3. COMPOSITE VULNERABILITY INDEX
# Combine multiple dimensions of vulnerability into a
# single index that identifies the neighborhoods where
# all risk factors converge.
#
# Components (all z-scored, higher = worse):
#   - Evap cooler prevalence (higher = more dependent)
#   - Income (inverted: lower income = higher vulnerability)
#   - Minority proportion (higher = more marginalized)
#   - Renter proportion (higher = less control over housing)
#   - Housing age (inverted: older = worse infrastructure)

bg$z_evap <- as.numeric(scale(bg$evap_prop))
bg$z_income_inv <- as.numeric(scale(-bg$med_income))  # invert: lower = worse
bg$z_minority <- as.numeric(scale(bg$pct_minority))
bg$z_renter <- as.numeric(scale(bg$pct_renter))
bg$z_age_inv <- as.numeric(scale(-bg$med_year_built))  # invert: older = worse

## Equal-weight composite
bg$vulnerability_index <- (bg$z_evap + bg$z_income_inv + bg$z_minority +
                            bg$z_renter + bg$z_age_inv) / 5

## Classify
bg$vuln_q <- cut(bg$vulnerability_index,
                  breaks = quantile(bg$vulnerability_index, probs = 0:4/4),
                  labels = c("Low", "Moderate", "High", "Very high"),
                  include.lowest = TRUE)

cat("\n=== Q3: VULNERABILITY INDEX ===\n")
cat("Range:", round(range(bg$vulnerability_index), 2), "\n")
print(table(bg$vuln_q))

## Demographics by vulnerability quartile
vuln_summary <- data.frame(
  quartile = levels(bg$vuln_q),
  evap = sapply(levels(bg$vuln_q), function(q)
    round(mean(bg$evap_prop[bg$vuln_q == q]) * 100, 1)),
  income = sapply(levels(bg$vuln_q), function(q)
    round(mean(bg$med_income[bg$vuln_q == q]))),
  minority = sapply(levels(bg$vuln_q), function(q)
    round(mean(bg$pct_minority[bg$vuln_q == q]) * 100, 0)),
  renter = sapply(levels(bg$vuln_q), function(q)
    round(mean(bg$pct_renter[bg$vuln_q == q]) * 100, 0))
)
print(vuln_summary)
write.csv(vuln_summary, here("results", "P7_Table3_vulnerability_index.csv"), row.names = FALSE)

# Q4. RENTER-OWNER DIVIDE
# Renters cannot control their cooling infrastructure.
# Landlords choose the cooling system. Renters bear the
# health consequences. This is a structural barrier to
# adaptation. What does the gradient look like?

renter_q <- cut(bg$pct_renter,
                 breaks = quantile(bg$pct_renter, probs = 0:4/4),
                 labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"),
                 include.lowest = TRUE)

renter_evap <- tapply(bg$evap_prop, renter_q, mean, na.rm = TRUE)
renter_income <- tapply(bg$med_income, renter_q, mean, na.rm = TRUE)

cat("\n=== Q4: RENTER-OWNER DIVIDE ===\n")
cat("Mean evap prevalence by renter quartile:\n")
print(round(renter_evap * 100, 1))
cat("Kruskal-Wallis p:", signif(kruskal.test(bg$evap_prop ~ renter_q)$p.value, 3), "\n")

# Q5. COVENANT PERSISTENCE SCORE
# How many dimensions of disadvantage persist in
# historically covenanted areas?

cov <- bg[bg$covenant == 1, ]
non_cov <- bg[bg$covenant == 0, ]

persistence <- data.frame(
  dimension = c("Evap. prevalence", "Income", "Minority %", "Renter %",
                "Year built", "Vulnerability index"),
  covenanted = sapply(c("evap_prop", "med_income", "pct_minority", "pct_renter",
                         "med_year_built", "vulnerability_index"),
                       function(v) round(mean(cov[[v]], na.rm = TRUE), 3)),
  non_covenanted = sapply(c("evap_prop", "med_income", "pct_minority", "pct_renter",
                              "med_year_built", "vulnerability_index"),
                            function(v) round(mean(non_cov[[v]], na.rm = TRUE), 3)),
  wilcox_p = sapply(c("evap_prop", "med_income", "pct_minority", "pct_renter",
                       "med_year_built", "vulnerability_index"),
                     function(v) signif(wilcox.test(cov[[v]], non_cov[[v]])$p.value, 3))
)

cat("\n=== Q5: COVENANT PERSISTENCE ===\n")
print(persistence)
write.csv(persistence, here("results", "P7_Table4_covenant_persistence.csv"), row.names = FALSE)

# Q6. POLICY THRESHOLDS
# At what evaporative cooler prevalence should a neighborhood
# be flagged for intervention? We define thresholds based on
# the health outcome data.
#
# If Paper 5 results exist, identify the evap prevalence
# above which call rates exceed the 90th percentile.

cat("\n=== Q6: POLICY THRESHOLDS ===\n")

## Evap prevalence at various health-risk cutoffs
p90_rate <- quantile(bg$call_rate, 0.90, na.rm = TRUE)
if (!is.na(p90_rate) && p90_rate > 0) {
  high_risk_bg <- bg[bg$call_rate >= p90_rate, ]
  cat("90th percentile call rate:", round(p90_rate, 2), "\n")
  cat("Mean evap prevalence in high-risk BGs:", round(mean(high_risk_bg$evap_prop) * 100, 1), "%\n")
  cat("Suggested intervention threshold: evap prevalence >",
      round(quantile(high_risk_bg$evap_prop, 0.25) * 100, 0), "%\n")
}

## Block groups meeting multiple criteria
bg$priority <- bg$evap_prop >= q75_evap &
  bg$med_income <= q25_income &
  bg$vulnerability_index > quantile(bg$vulnerability_index, 0.75)

cat("Priority neighborhoods (evap Q4 + income Q1 + vulnerability Q4):",
    sum(bg$priority), "\n")

priority_bg <- bg[bg$priority, ]
cat("Mean characteristics of priority neighborhoods:\n")
cat("  Evap:", round(mean(priority_bg$evap_prop) * 100, 1), "%\n")
cat("  Income: $", round(mean(priority_bg$med_income)), "\n")
cat("  Minority:", round(mean(priority_bg$pct_minority) * 100, 0), "%\n")
cat("  Renter:", round(mean(priority_bg$pct_renter) * 100, 0), "%\n")

# 10. FIGURES
bg_sf2 <- st_as_sf(bg)
if (is.na(st_crs(bg_sf2))) bg_sf2 <- st_set_crs(bg_sf2, 4326)
bg_sf2 <- st_make_valid(bg_sf2)

## ---- Figure 1: Cooling access deserts ----
bg_sf2$cooling_desert <- bg$cooling_desert
fig1 <- ggplot(bg_sf2) +
  geom_sf(aes(fill = factor(cooling_desert, levels = c(FALSE, TRUE))),
          color = "white", size = 0.15) +
  scale_fill_manual(values = c("FALSE" = "#e8e8e8", "TRUE" = "#d73027"),
                    labels = c("No", "Yes"), name = "Cooling\ndesert") +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(title = "")

pdf(here("results", "P7_Figure1_cooling_deserts.pdf"), width = 7, height = 7)
print(fig1)
dev.off()

png(here("results", "P7_Figure1_cooling_deserts.png"), width = 7, height = 7, units = "in", res = 300)
print(fig1)
dev.off()

## ---- Figure 2: Vulnerability index map ----
bg_sf2$vulnerability_index <- bg$vulnerability_index

fig2 <- ggplot(bg_sf2) +
  geom_sf(aes(fill = vulnerability_index), color = "white", size = 0.15) +
  scale_fill_gradientn(colors = c("#2c7bb6", "#abd9e9", "#fee090", "#d73027"),
                       name = "Vulnerability\nindex") +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.15, 0.25),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.3, "cm"),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(title = "")

pdf(here("results", "P7_Figure2_vulnerability_index.pdf"), width = 7, height = 7)
print(fig2)
dev.off()

png(here("results", "P7_Figure2_vulnerability_index.png"), width = 7, height = 7, units = "in", res = 300)
print(fig2)
dev.off()

## ---- Figure 3: Priority neighborhoods ----
bg_sf2$priority <- bg$priority

fig3 <- ggplot(bg_sf2) +
  geom_sf(aes(fill = factor(priority, levels = c(FALSE, TRUE))),
          color = "white", size = 0.15) +
  scale_fill_manual(values = c("FALSE" = "#e8e8e8", "TRUE" = "#7b3294"),
                    labels = c("No", "Priority"), name = "") +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(title = "")

pdf(here("results", "P7_Figure3_priority_neighborhoods.pdf"), width = 7, height = 7)
print(fig3)
dev.off()

png(here("results", "P7_Figure3_priority_neighborhoods.png"), width = 7, height = 7, units = "in", res = 300)
print(fig3)
dev.off()
