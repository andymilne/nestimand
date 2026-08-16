for (f in list.files("R", full.names = TRUE)) source(f)
source("dev/demo_data.R")
suppressPackageStartupMessages(library(marginaleffects))
## multi-trial data, so a cells random structure is estimable
dat3 <- do.call(rbind, lapply(1:3, function(i) {
  d <- dat; d$response <- d$response + rnorm(nrow(d), 0, 0.3); d }))
spl <- nesting_spec(dat3, response ~ chord_type * inversion + training +
                    (chord_type * inversion | participant),
                    "inversion %in% chord_type", fit = "lmer")
mc <- nest_fit(spl, random_structure = "chain")
cat("chain RE fit:", class(mc), "\n"); show_code(mc)
e <- estimand(mc, spl, chord_type, policy = "equal", bounds = FALSE, self_check = FALSE)
cat("estimate:", round(as.data.frame(e)$estimate[3], 4), "\n")
cat("\n---- brms dry run ----\n")
spb <- nesting_spec(dat, response ~ chord_type * inversion + training,
                    "inversion %in% chord_type", fit = "brms")
print(nest_fit(spb, dry_run = TRUE, chains = 4, iter = 2000, seed = 1))
cat("\n---- ordinal dry run ----\n")
do2 <- dat; do2$rating <- factor(round(pmin(pmax(do2$response,1),7)), ordered = TRUE)
spo <- nesting_spec(do2, rating ~ chord_type * inversion, "inversion %in% chord_type",
                    fit = "clmm", random = "(1 | participant)")
print(nest_fit(spo, dry_run = TRUE))
