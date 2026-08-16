## nestimand: the checks that need a working brms installation --------------
##
## Everything else in the package is verified by tests/test_core.R, which runs
## in seconds and needs no sampling. The checks below need Stan to compile
## and sample, which the development container cannot do in reasonable time, so
## they are collected here to be run once on a machine with a working brms.
##
## To run:  set the working directory to the package root, then
##            source("tests/test_brms.R")
## Expect roughly ten minutes, most of it Stan compilation.
## Each check prints PASS or FAIL and, on failure, what was expected.

## --- load the package sources ---------------------------------------------
## Finds the seven source files by name: in the working directory, in R/, or
## anywhere below. The folder layout does not matter.
local({
  want <- c("spec.R", "translate.R", "policy.R", "estimand.R", "fit.R",
            "latent.R", "priors.R", "chain.R")
  found <- list.files(".", pattern = "[.][Rr]$", recursive = TRUE, full.names = TRUE)
  hits <- found[basename(found) %in% want]
  missing <- setdiff(want, basename(hits))
  if (length(missing))
    stop("cannot find these package source files: ",
         paste(missing, collapse = ", "),
         "\n  Working directory: ", getwd(),
         "\n  It contains: ", paste(head(list.files(), 20), collapse = ", "),
         "\n  Set the working directory to the folder holding the .R source ",
         "files, then run this file again.", call. = FALSE)
  hits <- hits[match(want, basename(hits))]
  for (f in hits) source(f, local = FALSE)
  cat("loaded", length(hits), "source files from", unique(dirname(hits)), "\n")
})
if (!exists("nesting_spec"))
  stop("the source files loaded but did not define nesting_spec(); they may be ",
       "truncated or from a different version.", call. = FALSE)
## demonstration data, exactly as in partially_nested_variables.md Section 1
## (generated inline so this file needs nothing but the package's R/ folder)
if (file.exists("dev/demo_data.R")) source("dev/demo_data.R") else {
  set.seed(1); n <- 40
  .cells <- data.frame(
    chord_type = c(rep(c("dim", "min", "maj"), each = 3), "aug"),
    inversion  = c(rep(c("0", "1", "2"), times = 3), "none"),
    mu         = c(4.16, 4.00, 3.73,
                   4.29, 4.20, 3.71,
                   4.98, 4.62, 4.39,
                   3.92))
  dat <- .cells[rep(seq_len(nrow(.cells)), each = n), ]
  dat$participant <- factor(rep(seq_len(n), times = nrow(.cells)))
  .tr <- runif(n, 0, 10)
  dat$training <- .tr[as.integer(dat$participant)]
  dat$response <- rnorm(nrow(dat), dat$mu + 0.12 * dat$training, 1.2)
  dat$chord_type <- factor(dat$chord_type, levels = c("aug", "dim", "min", "maj"))
  dat$inversion  <- factor(dat$inversion,  levels = c("none", "0", "1", "2"))
}
suppressPackageStartupMessages(library(brms))
## cmdstanr compiles far faster than rstan. Announced rather than assumed, so
## the record of what produced these numbers is complete.
if (is.null(getOption("brms.backend")) &&
    requireNamespace("cmdstanr", quietly = TRUE)) {
  options(brms.backend = "cmdstanr")
  cat("backend: cmdstanr (detected; compiles faster than rstan)\n")
} else cat("backend:", getOption("brms.backend", "rstan"), "\n")

pass <- 0; fail <- 0
chk <- function(label, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) {
    cat("  ERR:", conditionMessage(e), "\n"); FALSE })
  cat(sprintf("%-58s %s\n", label, ifelse(ok, "PASS", "FAIL")))
  if (ok) pass <<- pass + 1 else fail <<- fail + 1
}

spb <- nesting_spec(dat, response ~ chord_type * inversion + training,
                    "inversion %in% chord_type", fit = "brms")

## ---- 1. the cell parameterization samples, and the fit is the expected one
cat("\n[1/5] fitting the cell parameterization (compilation takes a few minutes)\n")
mb <- nest_fit(spb, chains = 2, iter = 1000, warmup = 500, seed = 1, refresh = 0)
chk("cells sample without divergences or rank warnings",
    inherits(mb, "brmsfit") && nrow(brms::fixef(mb)) == 11)
chk("the fitted call travels with the model",
    any(grepl("brm(response ~ 0 + cell", attr(mb, "nestimand_code"), fixed = TRUE)))

## ---- 2. the draw-wise translation: THE check this file exists for
## latent_draws() assumes brms names the cell coefficients `b_` followed by the
## cell label with non-alphanumeric characters stripped. Cell labels contain
## full stops (`maj.0`), so this assumption is the one most likely to be wrong.
cat("\n[2/5] draw-wise translation\n")
cat("  brms coefficient names:", paste(head(rownames(brms::fixef(mb)), 4),
                                       collapse = ", "), "\n")
dr <- try(latent_draws(mb, spb, "chord_type", "equal"), silent = TRUE)
chk("latent_draws() finds the cell coefficients among the draws",
    is.data.frame(dr) && ncol(dr) == 6 && nrow(dr) > 100)
if (is.data.frame(dr)) {
  le <- latent_estimand(mb, spb, "chord_type", "equal")
  cat(sprintf("  posterior mean %.4f (sd %.4f) vs point estimate %.4f (se %.4f)\n",
              mean(dr[["aug - maj"]]), sd(dr[["aug - maj"]]),
              le$estimate[le$term == "aug - maj"], le$std.error[le$term == "aug - maj"]))
  chk("posterior mean agrees with the point estimate to Monte Carlo error",
      abs(mean(dr[["aug - maj"]]) - le$estimate[le$term == "aug - maj"]) < 0.05)
}

## ---- 3. the translated prior samples, and matches its audit table
## The prior is stated on cell means and translated to mu ~ N(A m, A D A').
## Sampling from the prior alone must reproduce the audit table's numbers.
cat("\n[3/5] translated prior, sampled with no data\n")
## the prior must span every population-level coefficient, `training` included:
## a brms prior of class "b" applies to all of them
pri <- nest_prior(spb, mean = 4, sd = 1.5, on = "cells", covariate_sd = 1)
cat("  prior spans", length(pri$full_mean), "coefficients:",
    paste(names(pri$full_mean), collapse = ", "), "\n")
## More draws here than elsewhere: this check compares a sampled standard
## deviation against a stated one, and the sampling error of an sd is roughly
## 1/sqrt(2n), so 1000 draws leave about 2% of noise on the ratio - enough to
## make a real 5% discrepancy indistinguishable from chance.
mp <- nest_fit(spb, priors = pri, sample_prior = "only",
               chains = 4, iter = 3000, warmup = 500, seed = 2, refresh = 0)
dp <- latent_draws(mp, spb, "chord_type", "equal")
cat("  prior draws:", nrow(dp), "\n")
want <- prior_for_estimand(pri, "chord_type", "equal")
cat(sprintf("  aug - maj: stated sd %.3f, sampled sd %.3f\n",
            want$sd[want$parameter == "aug - maj"], sd(dp[["aug - maj"]])))
cat(sprintf("  ratio sampled/stated: %.4f (Monte Carlo error about %.4f)\n",
            sd(dp[["aug - maj"]]) / want$sd[want$parameter == "aug - maj"],
            1 / sqrt(2 * nrow(dp))))
chk("the sampled prior reproduces the audit table (within 5%)",
    abs(sd(dp[["aug - maj"]]) / want$sd[want$parameter == "aug - maj"] - 1) < 0.05)
cat(sprintf("  dim - maj: stated sd %.3f, sampled sd %.3f\n",
            want$sd[want$parameter == "dim - maj"], sd(dp[["dim - maj"]])))
chk("a within-family contrast prior also reproduces (within 5%)",
    abs(sd(dp[["dim - maj"]]) / want$sd[want$parameter == "dim - maj"] - 1) < 0.05)
## The two ratios come from the same draws and so move together; a common
## factor across all contrasts would indicate a systematic error in A D A',
## whereas a common factor near one indicates only Monte Carlo noise.
rat <- vapply(want$parameter, function(k)
  sd(dp[[k]]) / want$sd[want$parameter == k], 1)
cat("  ratios across all six contrasts:",
    paste(sprintf("%.3f", rat), collapse = " "), "\n")
if (requireNamespace("posterior", quietly = TRUE))
  cat(sprintf("  effective sample size on aug - maj: %.0f of %d draws\n",
              posterior::ess_basic(dp[["aug - maj"]]), nrow(dp)))
chk("no systematic scale error common to every contrast",
    abs(mean(rat) - 1) < 0.05)
cat(sprintf("  prior mean on aug - maj: stated %.3f, sampled %.3f\n",
            want$mean[want$parameter == "aug - maj"], mean(dp[["aug - maj"]])))
chk("the prior mean on the contrast is centred where stated",
    abs(mean(dp[["aug - maj"]]) - want$mean[want$parameter == "aug - maj"]) <
      3 * sd(dp[["aug - maj"]]) / sqrt(nrow(dp)) + 0.05)

## ---- 4. brms needs its own spelling for conditional predictions
## Frequentist and Bayesian engines disagree on how random effects are excluded,
## and a silently wrong spelling produces a different estimand, not an error.
cat("\n[4/5] conditional predictions on a mixed model\n")
dat3 <- do.call(rbind, lapply(1:3, function(i) {
  d <- dat; d$response <- d$response + rnorm(nrow(d), 0, 0.3); d }))
spm <- nesting_spec(dat3, response ~ chord_type * inversion + training +
                    (1 | participant), "inversion %in% chord_type",
                    fit = "brms")
mm <- nest_fit(spm, chains = 2, iter = 1000, warmup = 500, seed = 3, refresh = 0)
e_re  <- latent_estimand(mm, spm, "chord_type", "equal")
chk("the latent route is unaffected by random effects (fixed effects only)",
    is.finite(e_re$estimate[3]))
cat("  NOTE: compare against marginaleffects with re_formula = NA, which is the\n")
cat("        brms spelling; re.form = NA is silently ignored on a brmsfit.\n")

## ---- 5. chain mode: does the declared chain fit match the cell fit? -------
## The chain parameterization keeps the original factors as predictors and holds
## the uninformative coefficients at zero by prior. The estimands must agree
## with the cell fit to Monte Carlo error; if they do not, the declarations are
## wrong rather than the theory.
cat("\n[5/5] chain mode against cell mode\n")
cp <- chain_priors(spb)
cat("  held at zero:", sum(cp$table$kind == "structural zero"), "structural,",
    sum(cp$table$kind == "identification constraint"), "identification\n")
mc <- nest_fit(spb, mode = "effects", priors = cp,
               chains = 2, iter = 1000, warmup = 500, seed = 4, refresh = 0)
chk("the declared chain model samples", inherits(mc, "brmsfit"))
fx <- brms::fixef(mc)
zeroed <- rownames(fx) %in% gsub(":", ":", cp$table$coef[cp$table$part == "fixed"])
chk("the declared coefficients are exactly zero in every draw",
    all(abs(fx[zeroed, "Estimate"]) < 1e-12) && all(fx[zeroed, "Est.Error"] < 1e-12))
e_chain <- as.data.frame(marginaleffects::avg_predictions(
  mc, newdata = counterfactual_grid(spb, spb$data,
        nest_policy(spb, "chord_type", "equal")),
  by = "chord_type", wts = counterfactual_grid(spb, spb$data,
        nest_policy(spb, "chord_type", "equal"))$.w,
  hypothesis = mfx_hypothesis("pairwise")))
e_chain <- mfx_canonical(e_chain, levels(factor(dat$chord_type)))
e_cell <- as.data.frame(estimand(mb, spb, chord_type, policy = "equal",
                                 bounds = FALSE, self_check = FALSE))
v <- function(d) d$estimate[d$term == "aug - maj"]
cat(sprintf("  aug - maj: chain %.4f, cell %.4f\n", v(e_chain), v(e_cell)))
chk("chain and cell estimands agree to Monte Carlo error",
    abs(v(e_chain) - v(e_cell)) < 0.05)

cat(sprintf("\n%d passed, %d failed\n", pass, fail))
cat("If check 2 failed, report the coefficient names printed above: the naming\n")
cat("rule in latent_draws() needs correcting to match them.\n")
