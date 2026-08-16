for (f in list.files("R", full.names = TRUE)) source(f)
source("dev/demo_data.R")
suppressPackageStartupMessages(library(lme4))
set.seed(9)
dat6 <- do.call(rbind, lapply(1:6, function(i) {
  d <- dat; d$response <- d$response + rnorm(nrow(d), 0, 0.4); d }))
sp <- nesting_spec(dat6, response ~ chord_type * inversion +
                   (chord_type * inversion | participant),
                   "inversion %in% chord_type", fit = "lmer")
m <- suppressWarnings(lme4::lmer(
  stats::as.formula(paste(deparse(cell_formula(sp)), "+", random_terms(sp, "cells"))),
  data = sp$data, control = lme4::lmerControl(check.nobs.vs.nRE = "ignore")))
S <- as.matrix(lme4::VarCorr(m)$participant)
cat("Sigma over realized cells:", nrow(S), "x", ncol(S), "\n")
## the heterogeneity question a reader asks: how much do participants differ in
## the chord-type effect itself? Under cells that is c'Sigma c, with c the SAME
## contrast vector the fixed-effect estimand uses.
pol <- nest_policy(sp, "chord_type", "equal")
cvec <- attr(latent_estimand(m, sp, "chord_type", pol), "nestimand_cvecs")
for (k in c("aug - maj", "dim - maj")) {
  cc <- cvec[[k]]
  cat(sprintf("  participant sd of the %-10s effect: %.3f\n", k,
              sqrt(as.numeric(t(cc) %*% S %*% cc))))
}
cat("\nc for aug - maj (weights over the ten realized cells):\n  ",
    paste(sprintf("%+.2f", cvec[["aug - maj"]]), collapse = " "), "\n")
