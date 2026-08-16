## What does `||` do with factors, on the lme4 you have installed?
## ---------------------------------------------------------------------------
## Standalone: needs only lme4, no other files. Runs in a few seconds.
##   source("doublebar_check.R")
##
## Background. `||` is meant to drop correlations among random effects. It acts
## while the formula is parsed, before any model matrix exists, so historically
## it could not tell a factor from a numeric variable: it split the term and
## left the levels within each part correlated. lme4 2.0-0 (January 2026) added
## a `diag()` term and an option to make `||` behave as expected. This script
## reports which behaviour your installation gives.

suppressPackageStartupMessages(library(lme4))
cat("lme4 version:", as.character(packageVersion("lme4")), "\n")
cat("R version   :", R.version.string, "\n")
has_diag <- exists("getDoublevertDefault", where = asNamespace("lme4"))
cat("getDoublevertDefault() present:", has_diag, "\n")
if (has_diag) cat("current default              :", lme4::getDoublevertDefault(), "\n")
cat("\n")

## ---- balanced data, one three-level factor and one two-level factor -------
set.seed(11)
np <- 40
d <- expand.grid(rep = 1:4, B = factor(c("b1", "b2")),
                 C = factor(c("c1", "c2", "c3")), participant = factor(1:np))
d$x <- rnorm(nrow(d))
d$y <- rnorm(nrow(d), 3 + as.integer(d$C) * 0.2, 1)

npar <- function(f, ...) {
  r <- tryCatch(lFormula(f, data = d,
         control = lmerControl(check.nobs.vs.nRE = "ignore"), ...)$reTrms,
       error = function(e) NULL)
  if (is.null(r)) return("error")
  sprintf("%3d params, blocks %s", length(r$theta),
          paste(lengths(r$cnms), collapse = "+"))
}

cat("--- the reference points ---\n")
cat("continuous, single bar   (1 + x | p)    :", npar(y ~ x + (1 + x | participant)), "\n")
cat("continuous, double bar   (1 + x || p)   :", npar(y ~ x + (1 + x || participant)), "\n")
cat("   expect the double bar to drop 1 parameter: this is the case `||` was built for\n\n")

cat("--- the case at issue: a three-level factor ---\n")
cat("single bar               (1 + C | p)    :", npar(y ~ C + (1 + C | participant)), "\n")
cat("double bar               (1 + C || p)   :", npar(y ~ C + (1 + C || participant)), "\n")
cat("   6 params from the single bar; a genuine diagonal would be 3.\n")
cat("   OLD behaviour ('split')  : 7 params in blocks 1+3 - MORE than the single bar,\n")
cat("                              and the 3-column block is still correlated within.\n")
cat("   NEW behaviour ('diag')   : 3 params, all blocks of width 1.\n\n")

cat("--- two factors ---\n")
cat("single bar               (1 + B + C | p) :", npar(y ~ B * C + (1 + B + C | participant)), "\n")
cat("double bar               (1 + B + C || p):", npar(y ~ B * C + (1 + B + C || participant)), "\n\n")

## ---- does diag() exist, and does it give the intended structure? ----------
cat("--- diag(), if available ---\n")
if (has_diag) {
  cat("diag(1 + C | p)                        :",
      npar(stats::as.formula("y ~ C + diag(1 + C | participant)")), "\n")
  cat("`||` under options(lme4.doublevert.default = 'diag_special'):\n")
  old <- getOption("lme4.doublevert.default")
  options(lme4.doublevert.default = "diag_special")
  cat("  (1 + C || p)                         :", npar(y ~ C + (1 + C || participant)), "\n")
  options(lme4.doublevert.default = old)
} else {
  cat("not available in this version; `||` cannot be made to work on factors here.\n")
}

## ---- the workaround that works on every version ---------------------------
cat("\n--- the version-independent route: one indicator column per level ---\n")
mm <- model.matrix(~ 0 + C, d)
colnames(mm) <- paste0("Cd", seq_len(ncol(mm)))
d <- cbind(d, as.data.frame(mm))
cat("(0 + Cd1 + Cd2 + Cd3 || p)             :",
    npar(y ~ C + (0 + Cd1 + Cd2 + Cd3 || participant)), "\n")
cat("   expect 3 params in blocks 1+1+1 on every version.\n")

## ---- what the fitted summary looks like -----------------------------------
cat("\n--- VarCorr from a single-bar and a double-bar fit ---\n")
m1 <- lmer(y ~ C + (1 + C | participant), data = d,
           control = lmerControl(check.nobs.vs.nRE = "ignore"))
m2 <- lmer(y ~ C + (1 + C || participant), data = d,
           control = lmerControl(check.nobs.vs.nRE = "ignore"))
cat("\nsingle bar:\n"); print(VarCorr(m1), comp = "Std.Dev.")
cat("\ndouble bar:\n"); print(VarCorr(m2), comp = "Std.Dev.")
cat("\nA diagonal structure shows no `Corr` column at all.\n")
