for (f in list.files("R", full.names = TRUE)) source(f)
source("dev/demo_data.R")
suppressPackageStartupMessages(library(brms))
spb <- nesting_spec(dat, response ~ chord_type * inversion + training,
                    "inversion %in% chord_type", fit = "brms")
mb <- nest_fit(spb, chains = 2, iter = 600, warmup = 300, seed = 1, refresh = 0)
cat("draws:", nrow(as.matrix(mb)), " fixef:", nrow(fixef(mb)), "\n")
print(head(rownames(fixef(mb)), 3))
d <- latent_draws(mb, spb, "chord_type", "equal")
cat("posterior of aug - maj: mean", round(mean(d[["aug - maj"]]), 4),
    " sd", round(sd(d[["aug - maj"]]), 4),
    " 95%:", paste(round(quantile(d[["aug - maj"]], c(.025,.975)), 3), collapse = " to "), "\n")
le <- latent_estimand(mb, spb, "chord_type", "equal")
cat("delta-method point/SE:", round(le$estimate[3], 4), round(le$std.error[3], 4), "\n")
