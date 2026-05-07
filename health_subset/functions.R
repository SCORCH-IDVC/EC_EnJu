
health_script <- function(){
  
  load("WS.RData")
  health_path <- "heat_calls_by_blockgroup_syn.csv"
  health <- read.csv(health_path, stringsAsFactors = FALSE)
  health$date <- as.Date(health$date)

  wx_path <- "tucson_hourly.csv"
  if (!file.exists(wx_path)) {
    stop("Run Paper 3 script first to download weather data.")
  }
  
  wx <- read.csv(wx_path, stringsAsFactors = FALSE)
  wx$tmpf <- as.numeric(wx$tmpf)
  wx$relh <- as.numeric(wx$relh)
  wx$date <- as.Date(wx$date)
  wx <- wx[!is.na(wx$tmpf) & !is.na(wx$relh), ]
  
  daily_wx <- aggregate(cbind(tmax = tmpf, rh_max = relh) ~ date,
                        data = wx, FUN = max, na.rm = TRUE)
  daily_wx$tmean <- aggregate(tmpf ~ date, data = wx, FUN = mean, na.rm = TRUE)$tmpf
  
  calc_twet <- function(temp_c, rh_pct) {
    temp_c * atan(0.151977 * (rh_pct + 8.313659)^0.5) +
      atan(temp_c + rh_pct) -
      atan(rh_pct - 1.676331) +
      0.00391838 * rh_pct^1.5 * atan(0.023101 * rh_pct) -
      4.686035
  }
  
  calc_supply <- function(temp_c, rh_pct, eta = 0.65) {
    twet <- calc_twet(temp_c, rh_pct)
    temp_c - eta * (temp_c - twet)
  }
  
  wx$failure <- calc_supply(wx$temp_c, wx$relh) > 27  # TRUE = failure
  daily_fail <- aggregate(failure ~ date, data = wx, FUN = sum)
  colnames(daily_fail)[2] <- "failure_hours"
  daily_wx <- merge(daily_wx, daily_fail, by = "date", all.x = TRUE)
  daily_wx$failure_hours[is.na(daily_wx$failure_hours)] <- 0
  daily_wx$failure_day <- as.integer(daily_wx$failure_hours > 0)
  
  bg$GEOID <- as.character(bg$GEOID)
  bg$GEOID <- sub("^0", "", bg$GEOID)

  bg_calls <- aggregate(calls ~ GEOID, data = health, FUN = sum)
  colnames(bg_calls)[2] <- "total_calls"
  
  bg_days <- aggregate(calls ~ GEOID, data = health, FUN = length)
  colnames(bg_days)[2] <- "n_days"
  
  bg_calls$GEOID <- as.character(bg_calls$GEOID)
  bg_days$GEOID <- as.character(bg_days$GEOID)
  
  bg <- merge(bg, bg_calls, by = "GEOID", all.x = TRUE)
  bg <- merge(bg, bg_days, by = "GEOID", all.x = TRUE)
  
  names(bg)[grep("total_calls", names(bg))[1]] <- "total_calls"
  names(bg)[grep("n_days", names(bg))[1]] <- "n_days"
  
  bg$total_calls[is.na(bg$total_calls)] <- 0
  bg$n_days[is.na(bg$n_days)] <- 1
  
  bg$call_rate <- bg$total_calls / bg$n_days * 1000
  
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
  
  write.csv(cor_health, "P5_Table1_health_correlations.csv", row.names = FALSE)
  
  bg$z_evap     <- scale(bg$evap_prop)
  bg$z_income   <- scale(bg$med_income)
  bg$z_minority <- scale(bg$pct_minority)
  bg$z_renter   <- scale(bg$pct_renter)
  
  m_nb <- glm.nb(total_calls ~ z_evap + z_income + z_minority + z_renter +
                   offset(log(n_days)),
                 data = bg)
  
  coef_nb <- data.frame(
    variable = c("Intercept", "Evap. prevalence (z)", "Income (z)",
                 "Minority (z)", "Renter (z)"),
    estimate = round(coef(m_nb), 3),
    se = round(summary(m_nb)$coefficients[, 2], 3),
    IRR = round(exp(coef(m_nb)), 3),
    p = signif(summary(m_nb)$coefficients[, 4], 3)
  )
  
  write.csv(coef_nb, "P5_Table2_NB_regression.csv", row.names = FALSE)
  
  daily_calls <- aggregate(calls ~ date, data = health, FUN = sum)
  daily_calls <- merge(daily_calls, daily_wx, by = "date", all.x = TRUE)
  daily_calls <- daily_calls[!is.na(daily_calls$tmax), ]
  daily_calls$doy <- as.integer(format(daily_calls$date, "%j"))
  daily_calls$year <- as.integer(format(daily_calls$date, "%Y"))
  
  m_ts <- glm(calls ~ failure_day + tmax + ns(doy, df = 4) + factor(year),
              family = poisson, data = daily_calls)
  

  disp <- sum(residuals(m_ts, type = "pearson")^2) / m_ts$df.residual

  if (disp > 1.5) {
    m_ts <- glm.nb(calls ~ failure_day + tmax + ns(doy, df = 4) + factor(year),
                   data = daily_calls)
  }
  
  ts_coefs <- summary(m_ts)$coefficients
  ts_table <- data.frame(
    variable = c("Failure day", "Tmax (F)"),
    estimate = round(ts_coefs[c("failure_day", "tmax"), 1], 3),
    se = round(ts_coefs[c("failure_day", "tmax"), 2], 3),
    IRR = round(exp(ts_coefs[c("failure_day", "tmax"), 1]), 3),
    p = signif(ts_coefs[c("failure_day", "tmax"), 4], 3)
  )

  write.csv(ts_table, "P5_Table3_timeseries_results.csv", row.names = FALSE)
  
  panel <- merge(health, daily_wx[, c("date", "tmax", "failure_day", "failure_hours")],
                 by = "date", all.x = TRUE)
  panel <- merge(panel, bg[, c("GEOID", "evap_prop", "med_income", "pct_minority",
                               "pct_renter")],
                 by = "GEOID", all.x = TRUE)
  panel <- panel[!is.na(panel$tmax) & !is.na(panel$evap_prop), ]
  panel$doy <- as.integer(format(panel$date, "%j"))
  panel$year <- as.integer(format(panel$date, "%Y"))
  
  panel$z_evap <- scale(panel$evap_prop)
  panel$z_income <- scale(panel$med_income)

  m_interact <- glm.nb(calls ~ failure_day * z_evap + tmax + z_income +
                         ns(doy, df = 4) + factor(year),
                       data = panel)
  
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
  
  write.csv(int_table, "P5_Table4_interaction_model.csv", row.names = FALSE)
  
  coords_sp <- cbind(bg$lon, bg$lat)
  nb <- knn2nb(knearneigh(coords_sp, k = 5))
  lw <- nb2listw(nb, style = "W")
  
  moran_calls <- moran.test(bg$call_rate, lw)

  bg_sf2 <- st_as_sf(bg)
  if (is.na(st_crs(bg_sf2))) bg_sf2 <- st_set_crs(bg_sf2, 4326)
  bg_sf2 <- st_make_valid(bg_sf2)
  
  ## ---- Figure 1a: Map of heat-related call rate ----
  fig1a <- ggplot(bg_sf2) +
    geom_sf(aes(fill = call_rate), color = "white", size = 0.15) +
    scale_fill_gradientn(colors = c("#2c7bb6", "#abd9e9", "#fee090", "#d73027"),
                         name = "HRI") +
    theme_minimal(base_size = 9) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = c(0.15, 0.25),
          legend.key.height = unit(0.4, "cm"),
          legend.key.width = unit(0.3, "cm"),
          plot.title = element_text(size = 10, face = "bold")) +
    labs(title = "a.")
  
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
    labs(title = "b.")
  
  fig1 <- fig1a + fig1b + plot_layout(ncol = 2)
  
  pdf("P5_Figure1_call_rate_maps.pdf", width = 10, height = 5)
  print(fig1)
  dev.off()
  
  png("P5_Figure1_call_rate_maps.png", width = 10, height = 5, units = "in", res = 300)
  print(fig1)
  dev.off()
  
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
         y = "HRI rate (per 1000 days)",
         title = "")
  
  pdf("P5_Figure2_scatter_calls_evap.pdf", width = 5, height = 5)
  print(fig2)
  dev.off()
  
  png("P5_Figure2_scatter_calls_evap.png", width = 5, height = 5, units = "in", res = 300)
  print(fig2)
  dev.off()
  
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
    labs(x = "", y = "Daily HRI (city-wide)", title = "")
  
  pdf("P5_Figure3_failure_day_calls.pdf", width = 5, height = 4)
  print(fig3)
  dev.off()
  
  png("P5_Figure3_failure_day_calls.png", width = 5, height = 4, units = "in", res = 300)
  print(fig3)
  dev.off()
  
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
         y = "HRI rate (per 1000 days)", title = "")
  
  pdf("P5_Figure4_calls_by_quartile.pdf", width = 6, height = 4)
  print(fig4)
  dev.off()
  
  png("P5_Figure4_calls_by_quartile.png", width = 6, height = 4, units = "in", res = 300)
  print(fig4)
  dev.off()
  
  p5_coef_path <- "P5_Table2_NB_regression.csv"
  health_path <- "heat_calls_by_blockgroup_syn.csv"
  
  nb_coefs <- read.csv(p5_coef_path)
  health <- read.csv(health_path, stringsAsFactors = FALSE)
  
  irr_evap <- nb_coefs$IRR[nb_coefs$variable == "Evap. prevalence (z)"]

  health$date <- as.Date(health$date)
  
  bg_calls <- aggregate(calls ~ GEOID, data = health, FUN = sum)
  
  colnames(bg_calls)[2] <- "total_calls"
  bg_days <- aggregate(calls ~ GEOID, data = health, FUN = length)
  colnames(bg_days)[2] <- "n_days"
  
  bg <- merge(bg, bg_calls, by = "GEOID", all.x = TRUE)
  bg <- merge(bg, bg_days, by = "GEOID", all.x = TRUE)
  
  names(bg)[grep("total_calls", names(bg))[1]] <- "total_calls"
  names(bg)[grep("n_days", names(bg))[1]] <- "n_days"
  
  bg$total_calls[is.na(bg$total_calls)] <- 0
  bg$n_days[is.na(bg$n_days)] <- 1
  bg$call_rate <- bg$total_calls / bg$n_days * 1000
  
  proj_results <- "P4_Table2_threshold_crossing.csv"
  p4_ws <- "ws_4_projections.RData"
  
  e4 <- new.env()
  load(p4_ws, envir = e4)
  delta_df <- e4$delta_df
  baseline_failure_days <- e4$baseline_failure_days

  for (i in seq_len(nrow(delta_df))) {
    s <- delta_df$ssp[i]
    p <- delta_df$period[i]
    d <- delta_df$delta_median[i]
    col_exp <- paste0("exposure_", s, "_", p)
    bg[[col_exp]] <- (baseline_failure_days + d) * bg$evap_prop
  }
  
  bg$exposure_current <- baseline_failure_days * bg$evap_prop
  bg$exposure_current_safe <- pmax(bg$exposure_current, 0.01)
  beta_exp <- log(irr_evap) / sd(bg$evap_prop)
  
  for (i in seq_len(nrow(delta_df))) {
    s <- delta_df$ssp[i]
    p <- delta_df$period[i]
    col_exp <- paste0("exposure_", s, "_", p)
    col_rate <- paste0("proj_rate_", s, "_", p)
    col_excess <- paste0("excess_calls_", s, "_", p)
    
    ratio <- bg[[col_exp]] / bg$exposure_current_safe
    bg[[col_rate]] <- bg$call_rate * ratio
    
    bg[[col_excess]] <- (bg[[col_rate]] - bg$call_rate) * 122 / 1000
    bg[[col_excess]] <- pmax(bg[[col_excess]], 0)
  }
  
  for (s in c("ssp245", "ssp585")) {
    for (p in c("near", "mid", "far")) {
      col_excess <- paste0("excess_calls_", s, "_", p)
      total <- round(sum(bg[[col_excess]], na.rm = TRUE), 1)
    }
  }
  
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
  
  table1 <- rbind(
    data.frame(ssp = "observed", period = "baseline",
               mean_proj_rate = round(mean(bg$call_rate), 2),
               total_excess_calls = 0,
               pct_increase = 0,
               n_bg_doubled = 0),
    table1
  )
  
  write.csv(table1, "P6_Table1_projected_health_burden.csv", row.names = FALSE)
  
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
  
  write.csv(table2, "P6_Table2_health_burden_demographics.csv", row.names = FALSE)
  
  bg$minority_q <- cut(bg$pct_minority,
                       breaks = quantile(bg$pct_minority, probs = 0:4/4),
                       labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"),
                       include.lowest = TRUE)
  
  disparity_rows <- list()
  
  rate_high_min <- mean(bg$call_rate[bg$minority_q == "Q4 (highest)"], na.rm = TRUE)
  rate_low_min <- mean(bg$call_rate[bg$minority_q == "Q1 (lowest)"], na.rm = TRUE)
  disparity_rows[[1]] <- data.frame(ssp = "observed", period = "baseline",
                                    ratio = round(rate_high_min / rate_low_min, 2))
  
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
  
  write.csv(disparity, "P6_Table3_disparity_index.csv", row.names = FALSE)
  
  bg_sf2 <- st_as_sf(bg)
  if (is.na(st_crs(bg_sf2))) bg_sf2 <- st_set_crs(bg_sf2, 4326)
  bg_sf2 <- st_make_valid(bg_sf2)
  
  lims <- range(c(bg$call_rate, bg$proj_rate_ssp585_far), na.rm = TRUE)
  
  fig1a <- ggplot(bg_sf2) +
    geom_sf(aes(fill = call_rate), color = "white", size = 0.15) +
    scale_fill_gradientn(colors = c("#2c7bb6", "#abd9e9", "#fee090", "#d73027"),
                         name = "HRI", limits = lims) +
    theme_minimal(base_size = 9) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 10, face = "bold")) +
    labs(title = "a.")
  
  fig1b <- ggplot(bg_sf2) +
    geom_sf(aes(fill = proj_rate_ssp585_far), color = "white", size = 0.15) +
    scale_fill_gradientn(colors = c("#2c7bb6", "#abd9e9", "#fee090", "#d73027"),
                         name = "HRI", limits = lims) +
    theme_minimal(base_size = 9) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 10, face = "bold")) +
    labs(title = "b.")
  
  fig1 <- fig1a + fig1b + plot_layout(ncol = 2)
  
  pdf("P6_Figure1_current_vs_projected_health.pdf", width = 10, height = 5)
  print(fig1)
  dev.off()
  
  png("P6_Figure1_current_vs_projected_health.png", width = 10, height = 5, units = "in", res = 300)
  print(fig1)
  dev.off()
  
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
    labs(x = "", y = "Total excess HRI per summer", fill = "")
  
  pdf("P6_Figure2_excess_calls.pdf", width = 6, height = 4)
  print(fig2)
  dev.off()
  
  png("P6_Figure2_excess_calls.png", width = 6, height = 4, units = "in", res = 300)
  print(fig2)
  dev.off()
  
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
    labs(x = "", y = "Disparity ratio\n(high-minority / low-minority HRI rate)", color = "")
  
  pdf("P6_Figure3_disparity_trajectory.pdf", width = 6, height = 4)
  print(fig3)
  dev.off()
  
  png("P6_Figure3_disparity_trajectory.png", width = 6, height = 4, units = "in", res = 300)
  print(fig3)
  dev.off()
  
  fig4 <- ggplot(bg_sf2) +
    geom_sf(aes(fill = excess_calls_ssp585_far), color = "white", size = 0.15) +
    scale_fill_gradientn(colors = c("#f7f7f7", "#fee090", "#fc8d59", "#d73027"),
                         name = "Excess HRI\nper summer") +
    theme_minimal(base_size = 9) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = c(0.15, 0.25),
          legend.key.height = unit(0.4, "cm"),
          legend.key.width = unit(0.3, "cm"),
          plot.title = element_text(size = 10, face = "bold")) +
    labs(title = "")
  
  pdf("P6_Figure4_excess_calls_map.pdf", width = 7, height = 7)
  print(fig4)
  dev.off()
  
  png("P6_Figure4_excess_calls_map.png", width = 7, height = 7, units = "in", res = 300)
  print(fig4)
  dev.off()
  
  
}
