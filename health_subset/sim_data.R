set.seed(42)

existing <- read.csv("heat_calls_by_blockgroup_old.csv")
geoids <- unique(existing$GEOID)

## Date range
dates <- seq(as.Date("2016-06-01"), as.Date("2024-09-30"), by = "day")
dates <- dates[as.integer(format(dates, "%m")) %in% 6:9]

n_3 <- 5
n_2 <- 26
n_1 <- 2864 - n_3 - n_2

all_combos <- expand.grid(GEOID = geoids, date = dates, stringsAsFactors = FALSE)
idx <- sample(seq_len(nrow(all_combos)), size = 2864)

## Assign call counts
calls_vec <- c(rep(3, n_3), rep(2, n_2), rep(1, n_1))
calls_vec <- sample(calls_vec)

health <- all_combos[idx, ]
health$calls <- calls_vec
health <- health[order(health$date, health$GEOID), ]
rownames(health) <- NULL

write.csv(health, "heat_calls_by_blockgroup_syn.csv", row.names = FALSE)