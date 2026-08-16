for (f in list.files("R", full.names = TRUE)) source(f)
source("dev/demo_data.R")
suppressPackageStartupMessages(library(marginaleffects))
sp <- nesting_spec(dat, response ~ chord_type * inversion + training,
                   nests = "inversion %in% chord_type")
m <- nest_fit(sp)
cat("class:", class(m), " rank:", m$rank, "\n"); show_code(m)
cat("\n---- effects mode via emmeans engine ----\n")
me <- nest_fit(sp, engine = "emmeans"); show_code(me)
cat("\n---- mixed, cells RE ----\n")
spl <- nesting_spec(dat, response ~ chord_type * inversion + training +
                    (chord_type * inversion + training | participant),
                    "inversion %in% chord_type", fit = "lmer")
cat("declared:", spl$random_original, "\n")
cat("cells   :", random_terms(spl, "cells"), "\n")
cat("chain   :", random_terms(spl, "chain"), "\n")
cat("\n---- full pipeline code view ----\n")
e <- estimand(m, sp, chord_type, policy = "equal", bounds = FALSE)
show_code(e)
cat("\nestimate:", round(as.data.frame(e)$estimate[3], 4),
    " check:", attr(e,"nestimand")$self_check$status, "\n")
