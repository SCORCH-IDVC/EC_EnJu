library(here)

dir.create(here("results"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("docs"), recursive = TRUE, showWarnings = FALSE)

## ---- Script 1: Equity & Environmental Justice ----
cat("\n==============================\n")
cat("Running Script 1: Equity\n")
cat("==============================\n")
source(here("scripts", "1_equity_environmental_justice.R"), local = (e1 <- new.env()))
save(list = ls(e1), envir = e1, file = here("docs", "ws_1_equity.RData"))
cat("Workspace saved: ws_1_equity.RData\n")

## ---- Script 2: Urban Heat Island ----
cat("\n==============================\n")
cat("Running Script 2: Heat Island\n")
cat("==============================\n")
source(here("scripts", "2_heat_island.R"), local = (e2 <- new.env()))
save(list = ls(e2), envir = e2, file = here("docs", "ws_2_heat.RData"))
cat("Workspace saved: ws_2_heat.RData\n")

## ---- Script 3: Monsoon Heat Vulnerability ----
cat("\n==============================\n")
cat("Running Script 3: Monsoon\n")
cat("==============================\n")
source(here("scripts", "3_monsoon_heat.R"), local = (e3 <- new.env()))
save(list = ls(e3), envir = e3, file = here("docs", "ws_3_monsoon.RData"))
cat("Workspace saved: ws_3_monsoon.RData\n")

## ---- Script 4: Climate Projections ----
cat("\n==============================\n")
cat("Running Script 4: Projections\n")
cat("==============================\n")
source(here("scripts", "4_projections.R"), local = (e4 <- new.env()))
save(list = ls(e4), envir = e4, file = here("docs", "ws_4_projections.RData"))
cat("Workspace saved: ws_4_projections.RData\n")

## ---- Script 5: Health Outcomes ----
cat("\n==============================\n")
cat("Running Script 5: Health\n")
cat("==============================\n")
source(here("scripts", "5_health_outcomes.R"), local = (e5 <- new.env()))
save(list = ls(e5), envir = e5, file = here("docs", "ws_5_health.RData"))
cat("Workspace saved: ws_5_health.RData\n")

## ---- Script 6: Health Projections ----
cat("\n==============================\n")
cat("Running Script 6: Health Projections\n")
cat("==============================\n")
source(here("scripts", "6_projected_health.R"), local = (e6 <- new.env()))
save(list = ls(e6), envir = e6, file = here("docs", "ws_6_projected_health.RData"))
cat("Workspace saved: ws_6_projected_health.RData\n")

## ---- Script 7: MISC ----
cat("\n==============================\n")
cat("Running Script 7: MISC\n")
cat("==============================\n")
source(here("scripts", "7_synthesis_policy.R"), local = (e7 <- new.env()))
save(list = ls(e7), envir = e7, file = here("docs", "ws_7_synthesis_policy.RData"))
cat("Workspace saved: ws_7_synthesis_policy.RData\n")

## ---- Render website ----
cat("\n==============================\n")
cat("Rendering website\n")
cat("==============================\n")
quarto::quarto_render(here("docs", "results.qmd"),
                      output_format = "html",
                      output_file = "index.html")
cat("Done.\n")
