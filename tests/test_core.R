## nestimand: regression suite for the cell-translation core ---------------
## Run from the folder holding the package source files:
##   source("test_core.R")
suppressPackageStartupMessages({library(marginaleffects)})
## --- load the package sources ---------------------------------------------
## Finds the seven source files by name: in the working directory, in R/, or
## anywhere below. The folder layout does not matter.
local({
  want <- c("spec.R", "translate.R", "policy.R", "estimand.R", "fit.R",
            "latent.R", "priors.R", "chain.R", "summary.R")
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

pass <- 0; fail <- 0
chk <- function(label, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) {
    cat("  ERR:", conditionMessage(e), "\n"); FALSE })
  cat(sprintf("%-62s %s\n", label, ifelse(ok, "PASS", "FAIL")))
  if (ok) pass <<- pass + 1 else fail <<- fail + 1
}
err_of <- function(expr) tryCatch({expr; ""}, error = function(e) conditionMessage(e))

## ---- declaration surface (migrated behaviour) ----------------------------
sp <- nesting_spec(dat, response ~ chord_type * inversion + training,
                   nests = "inversion %in% chord_type")
chk("nesting_spec: string nests", inherits(sp, "nesting_spec"))
set.seed(42); nd <- 30
d <- data.frame(arm = rep(c("control","A","A","B","B"), each = nd),
                dose = rep(c("none","low","high","low","high"), each = nd))
d$y <- rnorm(nrow(d), c(1.0,1.8,2.4,1.5,2.9)[match(paste(d$arm, d$dose),
        c("control none","A low","A high","B low","B high"))])
d$arm <- factor(d$arm, levels = c("control","A","B"))
d$dose_1 <- factor(d$dose, levels = c("none","low","high"))
spd <- nesting_spec(d, y ~ arm + arm:dose_1, dose_1 %in% arm)
chk("nesting_spec: bare nests", inherits(spd, "nesting_spec"))
chk("nesting_spec: data name captured ('d')", spd$data_name == "d")
chk("nesting_spec: realized cells counted per family",
    nrow(sp$cells_by_family[[1]]) == 10 && nrow(spd$cells_by_family[[1]]) == 5)
bad <- dat; bad$inversion[bad$chord_type == "aug"] <- NA
chk("nesting_spec: NA refused, apply_sentinel recommended",
    grepl("apply_sentinel", err_of(nesting_spec(bad, response ~ chord_type * inversion,
                                                "inversion %in% chord_type"))))
d3 <- data.frame(v = factor(rep(c("none","p"), 50)), r = c(rep(0,50), runif(50,4,8)))
d3$y <- rnorm(100)
chk("nesting_spec: continuous leaf accepted",
    inherits(nesting_spec(d3, y ~ v * r, "r %in% v"), "nesting_spec"))
chk("nesting_spec: continuous parent refused",
    grepl("cannot nest", err_of(nesting_spec(d3, y ~ v * r, "v %in% r"))))
sp_cl <- dat; sp_cl$cell <- 1
chk("nesting_spec: existing `cell` column refused, cell_name offered",
    grepl("cell_name", err_of(nesting_spec(sp_cl, response ~ chord_type * inversion,
                                           "inversion %in% chord_type"))))
chk("nesting_spec: cell_name honoured",
    identical(nesting_spec(sp_cl, response ~ chord_type * inversion,
              "inversion %in% chord_type", cell_name = "k")$cell_name, "k"))

## ---- cell construction ---------------------------------------------------
chk("cells: 10 realized of 16 crossed", length(sp$cell_levels) == 10)
chk("cells: continuous leaf excluded from the cell factor",
    identical(nesting_spec(d3, y ~ v * r, "r %in% v")$cell_vars, "v"))
chk("cells: cell factor added to the data, no NA",
    is.factor(sp$data$cell) && !anyNA(sp$data$cell) &&
    nlevels(sp$data$cell) == 10)
chk("cell_formula: cells mode",
    identical(paste(deparse(cell_formula(sp)), collapse = " "),
              "response ~ 0 + cell + training"))
chk("cell_formula: effects mode is the identified chain form",
    identical(paste(deparse(cell_formula(sp, "effects")), collapse = " "),
              "response ~ chord_type + chord_type:inversion + training"))

## ---- the two fitting modes ----------------------------------------------
m  <- lm(cell_formula(sp), data = sp$data)
mc <- lm(cell_formula(sp, "effects"), data = sp$data)
chk("fit: cell parameterization is full rank, no aliased coefficients",
    sum(is.na(coef(m))) == 0 && m$rank == length(coef(m)))
chk("fit: chain parameterization aliases 6 coefficients", sum(is.na(coef(mc))) == 6)
chk("fit: the two modes are the same fit", isTRUE(all.equal(logLik(m), logLik(mc))))

## ---- the effect basis ----------------------------------------------------
A <- effect_basis(sp)
chk("effect_basis: square and full rank (10 x 10)",
    nrow(A) == 10 && ncol(A) == 10 && qr(A)$rank == 10)
mu <- predict(m, newdata = transform(cell_grid(sp), training = 0))
b  <- coef(mc)[colnames(A)]
chk("effect_basis: mu = A m reproduces the cell means",
    isTRUE(all.equal(as.numeric(A %*% b), as.numeric(mu), tolerance = 1e-10)))
chk("effect_basis: A^-1 recovers the effect coefficients from cell means",
    isTRUE(all.equal(as.numeric(solve(A) %*% mu), as.numeric(b), tolerance = 1e-10)))

## ---- grid translation ----------------------------------------------------
g <- cell_grid(sp)
chk("cell_grid: one row per realized cell, covariate at its mean",
    nrow(g) == 10 && isTRUE(all.equal(unique(g$training), mean(dat$training))))
full <- expand.grid(chord_type = levels(dat$chord_type),
                    inversion = levels(dat$inversion))
chk("add_cells: unrealized combinations dropped by default",
    nrow(add_cells(sp, full)) == 10 &&
    attr(add_cells(sp, full), "nestimand_dropped") == 6)
chk("add_cells: unrealized combinations refused when asked",
    grepl("not identified by the design",
          err_of(add_cells(sp, full, unrealized = "error"))))
chk("add_cells: a grid missing a nesting variable is refused",
    grepl("cannot be recomputed",
          err_of(add_cells(sp, data.frame(chord_type = "maj")))))
chk("add_cells: recomputed cells agree with the fitted factor",
    identical(as.character(add_cells(sp, dat)$cell), as.character(sp$data$cell)))

## ---- policy weighting ----------------------------------------------------
est <- function(pol, tgt = "chord_type", model = m, spec = sp) {
  gg <- counterfactual_grid(spec, spec$data, pol)
  mfx_canonical(as.data.frame(avg_predictions(model, newdata = gg, by = tgt,
                wts = gg$.w, hypothesis = mfx_hypothesis("pairwise"))),
                levels(factor(spec$data[[tgt]])))
}
val <- function(e, tm) {
  v <- e$estimate[e$term == tm]
  if (!length(v)) {   # label conventions differ between engine versions
    rev_tm <- paste(rev(trimws(strsplit(tm, " - ", fixed = TRUE)[[1]])), collapse = " - ")
    v <- -e$estimate[e$term == rev_tm]
  }
  if (!length(v)) stop("no contrast labelled `", tm, "`; labels present: ",
                       paste(e$term, collapse = ", "))
  v
}
chk("policy equal: aug - maj = -0.6779",
    abs(val(est(nest_policy(sp, "chord_type", "equal")), "aug - maj") + 0.6779) < 1e-4)
chk("policy proportional: -0.6779 (balanced data)",
    abs(val(est(nest_policy(sp, "chord_type", "proportional")), "aug - maj") + 0.6779) < 1e-4)
chk("policy counterfactual: alias of proportional",
    identical(nest_policy(sp, "chord_type", "counterfactual")$p,
              nest_policy(sp, "chord_type", "proportional")$p))
chk("policy hierarchical: -0.6779 at depth one",
    abs(val(est(nest_policy(sp, "chord_type", "hierarchical")), "aug - maj") + 0.6779) < 1e-4)
vert <- vapply(c("0","1","2"), function(iv)
  val(est(nest_policy(sp, "chord_type", "nominated", at = c(inversion = iv))),
      "aug - maj"), 1)
chk("policy nominated: three single-version contrasts",
    all(abs(vert - c(-1.0479, -0.7405, -0.2452)) < 1e-4))
p3 <- nest_policy(sp, "chord_type", c("0" = 0.5, "1" = 0.3, "2" = 0.2))
chk("policy supplied: exactly the convex combination of the vertices",
    abs(val(est(p3), "aug - maj") - sum(vert * c(0.5, 0.3, 0.2))) < 1e-10)
chk("policy: degenerate stratum takes the single realized version",
    identical(unname(p3$p$aug), 1) && identical(names(p3$p$aug), "none"))
chk("policy supplied: incomplete coverage refused, not renormalized",
    grepl("silently", err_of(nest_policy(sp, "chord_type", c("0" = 1)))))
chk("policy nominated: unrealized version refused",
    grepl("must name a version that exists",
          err_of(nest_policy(sp, "chord_type", "nominated", at = c(inversion = "9")))))
chk("policy nominated: `at` required",
    grepl("needs it named", err_of(nest_policy(sp, "chord_type", "nominated"))))
chk("policy standardized: p must be supplied, not aliased",
    grepl("supply it directly", err_of(nest_policy(sp, "chord_type", "standardized"))))
chk("policy within: refused as a policy, redirected",
    grepl("not a policy", err_of(nest_policy(sp, "chord_type", "within"))))
chk("policy: bounds are the range of the vertices",
    abs(min(vert) + 1.0479) < 1e-4 && abs(max(vert) + 0.2452) < 1e-4)

## ---- reorder invariance (belt and braces: impossible by construction) -----
d2 <- dat
set.seed(7)
for (v in c("chord_type", "inversion")) {
  lv <- levels(d2[[v]])
  d2[[v]] <- factor(d2[[v]], levels = c(lv[1], sample(lv[-1])))
}
sp2 <- nesting_spec(d2, response ~ chord_type * inversion + training,
                    nests = "inversion %in% chord_type")
m2 <- lm(cell_formula(sp2), data = sp2$data)
chk("reorder: the permuted fit is the same fit",
    isTRUE(all.equal(logLik(m), logLik(m2))))
cat("marginaleffects", as.character(packageVersion("marginaleffects")),
    "- hypothesis passed as",
    if (mfx_formula_hypothesis()) "~pairwise" else "\"pairwise\"", "\n")
e1 <- est(nest_policy(sp, "chord_type", "equal"))
e2 <- est(nest_policy(sp2, "chord_type", "equal"), model = m2, spec = sp2)
chk("reorder: the estimand is unchanged under level permutation",
    isTRUE(all.equal(sort(round(abs(e1$estimate), 8)),
                     sort(round(abs(e2$estimate), 8)))))

## ---- depth two -----------------------------------------------------------
set.seed(31); nn <- 20
cells2 <- rbind(
  data.frame(chord = c("dim","min","maj"), inv = "0", doub = "octave", mu = c(4.30,4.45,5.10)),
  data.frame(chord = c("dim","min","maj"), inv = "0", doub = "fifth",  mu = c(4.02,4.13,4.86)),
  data.frame(chord = rep(c("dim","min","maj"), 2), inv = rep(c("1","2"), each = 3),
             doub = "none", mu = c(4.00,4.20,4.62, 3.73,3.71,4.39)),
  data.frame(chord = "aug", inv = "none", doub = "none", mu = 3.92))
dd <- cells2[rep(seq_len(nrow(cells2)), each = nn), ]
dd$y <- rnorm(nrow(dd), dd$mu, 1.0)
dd$chord <- factor(dd$chord, levels = c("aug","dim","min","maj"))
dd$inv   <- factor(dd$inv,   levels = c("none","0","1","2"))
dd$doub  <- factor(dd$doub,  levels = c("none","octave","fifth"))
sp2d <- nesting_spec(dd, y ~ chord * inv * doub, c(inv %in% chord, doub %in% inv))
chk("depth 2: 13 realized cells", length(sp2d$cell_levels) == 13)
m2d <- lm(cell_formula(sp2d), data = sp2d$data)
chk("depth 2: cell fit full rank, 13 coefficients",
    sum(is.na(coef(m2d))) == 0 && length(coef(m2d)) == 13)
A2 <- effect_basis(sp2d)
chk("depth 2: effect basis square and full rank (13 x 13)",
    nrow(A2) == 13 && ncol(A2) == 13 && qr(A2)$rank == 13)
e2d <- est(nest_policy(sp2d, "chord", "hierarchical"), tgt = "chord",
           model = m2d, spec = sp2d)
chk("depth 2: hierarchical aug - maj = -0.6283",
    abs(val(e2d, "aug - maj") + 0.6283) < 1e-4)
eq2d <- est(nest_policy(sp2d, "chord", "equal"), tgt = "chord",
            model = m2d, spec = sp2d)
chk("depth 2: equal differs from hierarchical (leaves vs tree)",
    abs(val(eq2d, "aug - maj") - val(e2d, "aug - maj")) > 1e-3)

## ---- estimand(): the function, and the code it kept ----------------------
set.seed(3)
e <- estimand(m, chord_type, spec = sp, policy = "equal")
chk("estimand: aug - maj = -0.6779 in original labels",
    abs(as.data.frame(e)$estimate[as.data.frame(e)$term == "aug - maj"] + 0.6779) < 1e-4)
meta <- attr(e, "nestimand")
chk("estimand: provenance recorded (policy, contrast, build)",
    identical(meta$policy, "equal") && identical(meta$contrast, "pairwise") &&
    identical(meta$build, nestimand_build))
chk("estimand: reorder self-check runs and passes",
    identical(meta$self_check$status, "passed"))
chk("estimand: bounds attached, -1.048 to -0.245 for aug - maj",
    abs(meta$bounds$policy_low[meta$bounds$term == "aug - maj"] + 1.0479) < 1e-4 &&
    abs(meta$bounds$policy_high[meta$bounds$term == "aug - maj"] + 0.2452) < 1e-4)
## fidelity: the code view is the code that ran, so re-running it must agree
env <- new.env(); assign("m", m, env); assign("sp", sp, env)
re <- eval(parse(text = paste(meta$code, collapse = "\n")), envir = env)
chk("show_code: the saved code re-runs and reproduces the estimand",
    isTRUE(all.equal(as.data.frame(re)$estimate, as.data.frame(e)$estimate)))
chk("show_code: prints the lines and returns them invisibly",
    is.character(capture.output(show_code(e))[1]))
## non-core arguments reach the destination function, and appear in the code
e9 <- estimand(m, chord_type, spec = sp, bounds = FALSE, self_check = FALSE,
               conf_level = 0.9)
chk("estimand: `...` passed through to marginaleffects",
    any(grepl("conf_level = 0.9", attr(e9, "nestimand")$code, fixed = TRUE)) &&
    as.data.frame(e9)$conf.low[1] > as.data.frame(e)$conf.low[1])
## within-stratum contrasts
ew <- estimand(m, inversion, spec = sp, contrast = "within", bounds = FALSE)
dw <- as.data.frame(ew)
chk("estimand within: three strata x three contrasts, no sentinel",
    nrow(dw) == 9 && !any(grepl("none", dw$term)) &&
    setequal(unique(dw$stratum), c("dim", "min", "maj")))
chk("estimand within: passes the reorder check",
    identical(attr(ew, "nestimand")$self_check$status, "passed"))
chk("estimand within: refused on a root variable",
    grepl("not nested within anything",
          err_of(estimand(m, chord_type, spec = sp, contrast = "within"))))
chk("estimand: a non-nesting target is refused",
    grepl("need no policy", err_of(estimand(m, training, spec = sp))))
chk("estimand: supplied policy is written into the code verbatim",
    any(grepl('c("0" = 0.5, "1" = 0.3, "2" = 0.2)',
        attr(estimand(m, chord_type, spec = sp, policy = c("0" = .5, "1" = .3, "2" = .2),
             bounds = FALSE, self_check = FALSE), "nestimand")$code, fixed = TRUE)))

## ---- nest_fit(): the fitting side ----------------------------------------
mf <- nest_fit(sp)
chk("nest_fit: fits the cell parameterization, full rank",
    inherits(mf, "lm") && sum(is.na(coef(mf))) == 0)
chk("nest_fit: the call travels with the fit",
    any(grepl("lm(response ~ 0 + cell + training, data = sp$data)",
              attr(mf, "nestimand_code"), fixed = TRUE)))
chk("nest_fit: the parameterization and its reason are stated",
    any(grepl("parameterization: cells", attr(mf, "nestimand_code"))))
chk("nest_fit: emmeans engine selects the effect basis, with its reason",
    identical(attr(nest_fit(sp, engine = "emmeans"), "nestimand_mode"), "effects"))
ef <- estimand(mf, chord_type, spec = sp, policy = "equal", bounds = FALSE, self_check = FALSE)
chk("nest_fit: the code view joins fit and estimand into one pipeline",
    any(grepl("^mf <- lm", attr(ef, "nestimand")$code)) &&
    any(grepl("avg_predictions", attr(ef, "nestimand")$code)))
dr <- nest_fit(sp, dry_run = TRUE, weights = NULL)
chk("nest_fit: dry_run returns code without fitting", inherits(dr, "nestimand_code"))

## random-effects translation
spl <- nesting_spec(dat, response ~ chord_type * inversion + training +
                    (chord_type * inversion + training | participant),
                    "inversion %in% chord_type", fit = "lmer")
chk("random_terms: cells give an unstructured covariance over realized cells",
    identical(random_terms(spl, "cells"), "(0 + cell + training | participant)"))
chk("random_terms: chain gives the progressive grouping submodel",
    identical(random_terms(spl, "chain"),
      paste("(1 + training | participant) + (1 | participant:chord_type) +",
            "(1 | participant:chord_type:inversion)")))
spi <- nesting_spec(dat, response ~ chord_type * inversion + (1 | participant),
                    "inversion %in% chord_type", fit = "lmer")
chk("random_terms: a clean intercept term passes through unchanged",
    identical(random_terms(spi, "cells"), "(1 | participant)"))

## ordinal engines: thresholds already carry an intercept
do2 <- dat; do2$rating <- factor(round(pmin(pmax(do2$response, 1), 7)), ordered = TRUE)
spo <- nesting_spec(do2, rating ~ chord_type * inversion,
                    "inversion %in% chord_type", fit = "clm")
chk("ordinal: cell coding keeps the implicit intercept, no redundant parameter",
    identical(paste(deparse(cell_formula(spo)), collapse = " "), "rating ~ cell"))
mo <- nest_fit(spo)
chk("ordinal: clm fits without the intercept warning, 15 parameters",
    length(coef(mo)) == 15)
chk("ordinal: the scale of the contrast must be stated, not assumed",
    grepl("no single response scale", err_of(estimand(mo, chord_type, spec = spo))))
spbo <- nesting_spec(do2, rating ~ chord_type * inversion,
                     "inversion %in% chord_type", fit = "brms", family = "cumulative()")
chk("ordinal: brms cumulative also takes the threshold-aware coding",
    identical(paste(deparse(cell_formula(spbo)), collapse = " "), "rating ~ cell"))
chk("ordinal: brms family threaded into the emitted call",
    grepl("family = cumulative()", nest_fit(spbo, dry_run = TRUE), fixed = TRUE))
chk("brms: non-core arguments reach the engine and appear in the code",
    grepl("chains = 4, iter = 2000, seed = 1",
          nest_fit(nesting_spec(dat, response ~ chord_type * inversion,
                   "inversion %in% chord_type", fit = "brms"),
                   dry_run = TRUE, chains = 4, iter = 2000, seed = 1), fixed = TRUE))

## ---- latent scale: estimands as linear functionals -----------------------
lp <- latent_estimand(m, "chord_type", "equal", spec = sp)
ap <- as.data.frame(estimand(m, chord_type, spec = sp, bounds = FALSE, self_check = FALSE))
chk("latent: agrees with the prediction route on estimates (1e-10)",
    max(abs(lp$estimate[match(ap$term, lp$term)] - ap$estimate)) < 1e-10)
chk("latent: agrees with the prediction route on standard errors (1e-7)",
    max(abs(lp$std.error[match(ap$term, lp$term)] - ap$std.error)) < 1e-7)
lv <- latent_estimand(m, "chord_type", c("0" = 0.5, "1" = 0.3, "2" = 0.2), spec = sp)
chk("latent: a supplied policy is the convex combination of the vertices",
    abs(lv$estimate[lv$term == "aug - maj"] - sum(vert * c(0.5, 0.3, 0.2))) < 1e-8)
chk("latent: reference and sequential contrasts available",
    nrow(latent_estimand(m, "chord_type", "equal", contrast = "reference", spec = sp)) == 3 &&
    nrow(latent_estimand(m, "chord_type", "equal", contrast = "sequential", spec = sp)) == 3)
lo <- estimand(mo, chord_type, spec = spo, policy = "equal", scale = "latent")
chk("latent: ordinal fit returns one number per contrast, where predictions fail",
    nrow(as.data.frame(lo)) == 6 && is.finite(as.data.frame(lo)$estimate[3]))
chk("latent: ordinal estimand passes the reorder check",
    identical(attr(lo, "nestimand")$self_check$status, "passed"))
chk("latent: bounds computed on the latent scale too",
    !is.null(attr(lo, "nestimand")$bounds) &&
    all(attr(lo, "nestimand")$bounds$policy_low <=
        attr(lo, "nestimand")$bounds$estimate + 1e-10))
chk("latent: the code view names the linear-map route",
    any(grepl("latent_estimand", attr(lo, "nestimand")$code)))
chk("latent: within contrasts have no linear-map form, and say so",
    grepl("no linear-map form",
          err_of(estimand(m, inversion, spec = sp, contrast = "within", scale = "latent"))))
chk("latent: draw-wise translation refuses a frequentist fit",
    grepl("needs a posterior", err_of(latent_draws(m, "chord_type", spec = sp))))
chk("latent: a model not fitted from this spec is refused",
    grepl("refit with nest_fit",
          err_of(latent_estimand(lm(response ~ training, data = dat), "chord_type", spec = sp))))

## ---- prior translation ---------------------------------------------------
pc <- nest_prior(sp, mean = 4, sd = 1.5, on = "cells", covariate_sd = 1)
chk("prior: stated cell sds are recovered exactly after translation",
    isTRUE(all.equal(unname(sqrt(diag(pc$cell_cov))), rep(1.5, 10))))
chk("prior: independent cell priors imply correlated effect priors",
    abs(sqrt(diag(pc$eff_cov))[2] - 1.5 * sqrt(2)) < 1e-10 &&
    abs(pc$eff_cov[2, 3]) > 1e-6)
pe <- nest_prior(sp, mean = c("(Intercept)" = 4, .default = 0),
                 sd = c("(Intercept)" = 2, .default = 0.5), on = "effects",
                 covariate_sd = 1)
chk("prior: `.default` covers the parameters not named",
    length(pe$eff_mean) == 10 && unname(pe$eff_mean[1]) == 4 &&
    all(unname(pe$eff_mean[-1]) == 0))
chk("prior: a silent parameter is refused, not quietly widened",
    grepl("undefined one",
          err_of(nest_prior(sp, mean = c("(Intercept)" = 4), sd = 1,
                            on = "effects", covariate_sd = 1))))
chk("prior: round trip through both spaces is exact",
    isTRUE(all.equal(as.numeric(solve(pc$A) %*% pc$cell_mean),
                     as.numeric(pc$eff_mean))) &&
    isTRUE(all.equal(as.numeric(pc$A %*% pc$eff_cov %*% t(pc$A)),
                     as.numeric(pc$cell_cov))))
au <- prior_audit(pc)
chk("prior_audit: reports both spaces, one row per parameter",
    nrow(au) == 20 && setequal(unique(au$space), c("effects", "cells")))
pfe <- prior_for_estimand(pc, "chord_type", "equal")
chk("prior_for_estimand: aug - maj sd is 1.5 * sqrt(1 + 1/3)",
    abs(pfe$sd[pfe$parameter == "aug - maj"] - 1.5 * sqrt(1 + 1/3)) < 1e-10)
chk("prior_for_estimand: a within-family contrast is tighter",
    abs(pfe$sd[pfe$parameter == "dim - maj"] - 1.5 * sqrt(2/3)) < 1e-10)
spb2 <- nesting_spec(dat, response ~ chord_type * inversion + training,
                     "inversion %in% chord_type", fit = "brms")
prb <- nest_prior(spb2, mean = 4, sd = 1.5, on = "cells", covariate_sd = 1)
chk("prior: brms call carries the statement and the stanvars",
    grepl("prior = prior_statement(prb), stanvars = prior_stanvars(prb)",
          nest_fit(spb2, priors = prb, dry_run = TRUE), fixed = TRUE))
chk("prior: the translation is stated in the emitted code",
    grepl("A D A", nest_fit(spb2, priors = prb, dry_run = TRUE), fixed = TRUE))
chk("prior: refused on an engine that has no prior",
    grepl("no prior to state", err_of(nest_fit(sp, priors = pc, dry_run = TRUE))))
chk("prior: student_t family carried into the statement",
    grepl("multi_student_t(3",
          prior_statement(nest_prior(spb2, mean = 4, sd = 1.5, on = "cells",
                                     covariate_sd = 1,
                                     family = "student_t"))$prior, fixed = TRUE))

## ---- engine compatibility ------------------------------------------------
chk("mfx: formula hypothesis chosen from 0.19.0 onward",
    mfx_formula_hypothesis("0.19.0") && mfx_formula_hypothesis("0.30.2") &&
    !mfx_formula_hypothesis("0.18.0"))
chk("mfx: the emitted code matches the installed version",
    identical(mfx_hypothesis_txt("pairwise", "0.18.0"), "\"pairwise\"") &&
    identical(mfx_hypothesis_txt("pairwise", "0.19.0"), "~pairwise"))
chk("mfx: the emitted estimand code uses the accepted spelling",
    any(grepl(paste0("hypothesis = ", mfx_hypothesis_txt("pairwise")),
              attr(e, "nestimand")$code, fixed = TRUE)))

## ---- engine label conventions --------------------------------------------
## marginaleffects 0.18.x returns `term` = "aug - maj" = -0.6779; 0.32.0 returns
## `hypothesis` = "(maj) - (aug)" = +0.6779. Both must yield the same reported
## estimand, or the same analysis would report opposite signs on two machines.
new_style <- data.frame(
  hypothesis = c("(dim) - (aug)", "(maj) - (aug)", "(maj) - (min)"),
  estimate = c(-0.049, 0.6779, 0.7046), std.error = c(0.21, 0.21, 0.15),
  conf.low = c(-0.47, 0.258, 0.408), conf.high = c(0.37, 1.098, 1.002))
old_style <- data.frame(
  term = c("aug - dim", "aug - maj", "min - maj"),
  estimate = c(0.049, -0.6779, -0.7046), std.error = c(0.21, 0.21, 0.15),
  conf.low = c(-0.37, -1.098, -1.002), conf.high = c(0.47, -0.258, -0.408))
levs <- c("aug", "dim", "min", "maj")
cn <- mfx_canonical(new_style, levs); co <- mfx_canonical(old_style, levs)
chk("labels: the `hypothesis` column is located as well as `term`",
    identical(mfx_term_column(new_style), "hypothesis") &&
    identical(mfx_term_column(old_style), "term"))
chk("labels: parentheses stripped, level order restored",
    identical(cn$term, c("aug - dim", "aug - maj", "min - maj")))
chk("labels: contrast direction fixed by the package, not the engine",
    isTRUE(all.equal(cn$estimate, co$estimate, tolerance = 1e-8)))
chk("labels: confidence bounds swapped and negated with the estimate",
    isTRUE(all.equal(cn$conf.low, co$conf.low, tolerance = 1e-8)) &&
    isTRUE(all.equal(cn$conf.high, co$conf.high, tolerance = 1e-8)))
chk("labels: the number of reversed rows is recorded",
    attr(cn, "nestimand_flipped") == 3 && attr(co, "nestimand_flipped") == 0)
chk("labels: an unrecognized label is left alone rather than mangled",
    identical(mfx_canonical(data.frame(term = "b0", estimate = 1))$term, "b0"))
chk("labels: estimand() reports contrasts in declared level order",
    identical(as.data.frame(e)$term[1:3],
              c("aug - dim", "aug - min", "aug - maj")))

chk("show_code: the normalization is part of the saved code, not applied after",
    any(grepl("mfx_canonical", attr(e, "nestimand")$code)))
chk("show_code: within-contrast code re-runs and reproduces its estimand",
    { ew2 <- estimand(m, inversion, spec = sp, contrast = "within", bounds = FALSE,
                      self_check = FALSE)
      env2 <- new.env(); assign("m", m, env2); assign("sp", sp, env2)
      rr <- eval(parse(text = paste(attr(ew2, "nestimand")$code, collapse = "\n")),
                 envir = env2)
      isTRUE(all.equal(as.data.frame(rr)$estimate, as.data.frame(ew2)$estimate)) })

## ---- priors must span every population-level coefficient ------------------
## brms applies a class = "b" prior to all of them, covariates included. A
## prior stated only on the cells leaves Stan with mismatched dimensions, which
## it reports only as chains finishing unexpectedly.
chk("prior: covariate coefficients identified from the design matrix",
    identical(fitted_covariate_names(spb2), "training") &&
    length(fitted_coef_names(spb2)) == 11)
chk("prior: an incomplete prior is refused before fitting, with the remedy",
    grepl("State `covariate_sd`",
          err_of(nest_prior(spb2, mean = 4, sd = 1.5, on = "cells"))))
chk("prior: the full vector spans the model exactly, in its order",
    identical(names(prb$full_mean), fitted_coef_names(spb2)))
chk("prior: the cell block is untouched, so contrasts are unaffected",
    length(prb$cell_mean) == 10 &&
    abs(prior_for_estimand(prb, "chord_type", "equal")$sd[3] -
        1.5 * sqrt(1 + 1/3)) < 1e-10)
chk("prior: the dimension is checked before sampling starts",
    isTRUE(check_prior_dimension(prb, spb2)))
bad_prior <- prb; bad_prior$full_mean <- bad_prior$full_mean[1:5]
chk("prior: a mismatched prior is caught rather than passed to Stan",
    grepl("must span every population-level",
          err_of(check_prior_dimension(bad_prior, spb2))))
chk("prior: ordinal families refused for now, with the reason",
    grepl("threshold coding",
          err_of(nest_prior(spbo, mean = 4, sd = 1.5, on = "cells"))))

## ---- brms draw naming ----------------------------------------------------
chk("draws: exact `b_` names matched, full stops included",
    identical(draw_names(c("cellaug.none", "training"),
              c("b_cellaug.none", "b_training", "sigma")),
              c("b_cellaug.none", "b_training")))
chk("draws: punctuation-insensitive fallback if the convention changes",
    identical(draw_names("cellaug.none", c("b_cellaugnone")), "b_cellaugnone"))
chk("draws: an unmatched coefficient reports the names available",
    grepl("Draw columns available", err_of(draw_names("cellmaj.9", "b_x"))))

## ---- why cells: the random-effects argument -------------------------------
## A random slope over the nesting structure is rank-deficient by construction
## under the chain form. The engine fits it anyway, so the constraint is
## invisible unless the package says so.
rk <- random_slope_rank(spl)
chk("RE: the chain random slope is rank-deficient (16 columns, rank 10)",
    rk$columns == 16 && rk$rank == 10 && rk$zero_columns == 3 && !rk$identified)
chk("RE: passing it through warns, with the numbers and the remedy",
    { w <- NULL
      withCallingHandlers(random_terms(spl, "as_declared"),
        warning = function(x) { w <<- conditionMessage(x); invokeRestart("muffleWarning") })
      grepl("not identified by any amount of data", w) &&
      grepl("random_structure = \"cells\"", w) })
chk("RE: the cell form puts the whole structure on one identified factor",
    identical(random_terms(spl, "cells"), "(0 + cell + training | participant)"))
chk("RE: a structure with no structural slope is passed through in silence",
    { w <- NULL
      r <- withCallingHandlers(random_terms(spi, "as_declared"),
        warning = function(x) { w <<- conditionMessage(x); invokeRestart("muffleWarning") })
      is.null(w) && identical(r, "(1 | participant)") })

## ---- chain mode: declarations for brms -----------------------------------
spc <- nesting_spec(dat, response ~ chord_type * inversion + training +
                    (chord_type * inversion | participant),
                    "inversion %in% chord_type", fit = "brms")
cp <- chain_priors(spc)
chk("chain: fixed design is 16 informative-plus-empty columns of rank 10",
    cp$fixed$columns == 16 && cp$fixed$rank == 10)
chk("chain: three structural zeros and three identification constraints, fixed",
    sum(cp$table$part == "fixed" & cp$table$kind == "structural zero") == 3 &&
    sum(cp$table$part == "fixed" & cp$table$kind == "identification constraint") == 3)
chk("chain: the random side needs the same six, from the same design",
    sum(cp$table$part == "random" & cp$table$kind == "structural zero") == 3 &&
    sum(cp$table$part == "random" & cp$table$kind == "identification constraint") == 3)
chk("chain: structural zeros are the empty augmented interactions",
    setequal(cp$table$coef[cp$table$part == "fixed" &
                           cp$table$kind == "structural zero"],
             c("chord_typeaug:inversion0", "chord_typeaug:inversion1",
               "chord_typeaug:inversion2")))
chk("chain: the two kinds of zero are distinguished in the emitted code",
    any(grepl("structural zero", cp$code)) &&
    any(grepl("identification constraint", cp$code)) &&
    any(grepl("coding", cp$code)))
chk("chain: the random term matches the fixed side, not the declared crossed form",
    identical(random_terms(spc, "chain_slope"),
              "(0 + chord_type + chord_type:inversion | participant)"))
chk("chain: chain mode selects the chain random term automatically",
    grepl("(0 + chord_type + chord_type:inversion",
          nest_fit(spc, mode = "effects", priors = cp, dry_run = TRUE), fixed = TRUE))
chk("chain: the declarations reach the brms call",
    grepl("prior = chain_prior_object(cp)",
          nest_fit(spc, mode = "effects", priors = cp, dry_run = TRUE), fixed = TRUE))
chk("chain: cells mode is unaffected and needs no declarations",
    grepl("(0 + cell | participant)", nest_fit(spc, dry_run = TRUE), fixed = TRUE))
chk("chain: declarations refused on an engine that cannot express them",
    grepl("brms facility",
          err_of(chain_priors(nesting_spec(dat, response ~ chord_type * inversion,
                              "inversion %in% chord_type", fit = "lmer",
                              random = "(1 | participant)")))))
chk("chain: mismatched mode refused, with the alternative named",
    grepl("nest_prior", err_of(nest_fit(spc, mode = "cells", priors = cp,
                                        dry_run = TRUE))))
chk("chain: regularization is separable from identification",
    is.null(chain_priors(spc, regularize = NULL)$regularize) &&
    !any(grepl("normal", chain_priors(spc, regularize = NULL)$code)))

## ---- the random structure is the user's to simplify -----------------------
## Translation applies only to terms that cross the boundary. A simplified
## structure must come back unchanged: enlarging it would be the package
## overriding a deliberate modelling decision, often one made to get a model
## to converge at all.
resp_re <- function(re) nesting_spec(dat,
  stats::as.formula(paste("response ~ chord_type * inversion + training +", re)),
  "inversion %in% chord_type", fit = "lmer")
chk("RE: a random intercept is left exactly as declared",
    identical(random_terms(resp_re("(1 | participant)"), "cells"),
              "(1 | participant)"))
chk("RE: a slope on the nesting root is identified, so it is left alone",
    identical(random_terms(resp_re("(chord_type | participant)"), "cells"),
              "(chord_type | participant)"))
chk("RE: a hand-written grouping structure passes through untouched",
    identical(random_terms(resp_re(
      "(1 | participant) + (1 | participant:chord_type)"), "cells"),
      "(1 | participant) + (1 | participant:chord_type)"))
chk("RE: the full structure translates, without repeating subsumed terms",
    identical(random_terms(resp_re("(chord_type * inversion | participant)"), "cells"),
              "(0 + cell | participant)"))
chk("RE: covariate slopes survive the translation",
    identical(random_terms(resp_re(
      "(chord_type * inversion + training | participant)"), "cells"),
      "(0 + cell + training | participant)"))
chk("RE: a partial structure is refused, with three named alternatives",
    { e <- err_of(random_terms(resp_re("(inversion | participant)"), "cells"))
      grepl("reduced", e) && grepl("random_structure = \"chain\"", e) &&
      grepl("as_declared", e) })
chk("RE: the chain submodel accepts a partial structure",
    grepl("participant:chord_type",
          random_terms(resp_re("(inversion | participant)"), "chain")))
chk("RE: boundary_vars names the nested variables, not the roots",
    identical(boundary_vars(sp), "inversion"))

## ---- `||` is refused where it misleads ------------------------------------
dat$xnum <- rnorm(nrow(dat))
mk_re <- function(re) nesting_spec(dat,
  stats::as.formula(paste("response ~ chord_type * inversion +", re)),
  "inversion %in% chord_type", fit = "lmer")
chk("double bar: refused on a factor term, with the remedy named",
    { e <- err_of(mk_re("(chord_type * inversion || participant)"))
      grepl("does not give the diagonal covariance", e) && grepl("diag", e) &&
      grepl("expand_re", e) })
chk("double bar: the offending factors are named",
    grepl("`chord_type`", err_of(mk_re("(chord_type || participant)"))))
chk("double bar: allowed on numeric terms, where it behaves as expected",
    inherits(mk_re("(1 + xnum || participant)"), "nesting_spec"))
chk("double bar: a single bar is unaffected",
    inherits(mk_re("(chord_type * inversion | participant)"), "nesting_spec"))

## ---- covariance-structure wrappers survive translation --------------------
## `diag(...)` (lme4 2.0-0 and later) states that the covariance is diagonal.
## Dropping it during translation would turn a diagonal request into an
## unstructured one, silently.
chk("diag(): not double-wrapped in parentheses",
    identical(mk_re("diag(chord_type * inversion | participant)")$random_original,
              "diag(chord_type * inversion | participant)"))
chk("diag(): the wrapper survives translation to cells",
    identical(random_terms(mk_re("diag(chord_type * inversion | participant)"), "cells"),
              "diag(0 + cell | participant)"))
chk("diag(): an untranslated term keeps its wrapper too",
    identical(random_terms(mk_re("diag(1 + chord_type | participant)"), "cells"),
              "diag(1 + chord_type | participant)"))
chk("diag(): a plain bar is unaffected",
    identical(random_terms(mk_re("(chord_type * inversion | participant)"), "cells"),
              "(0 + cell | participant)"))

## ---- apply_sentinel(): `where` is optional ---------------------------------
dna <- dat; dna$inversion[dna$chord_type == "aug"] <- NA
w <- NULL
d_auto <- withCallingHandlers(apply_sentinel(dna, "inversion"),
  warning = function(x) { w <<- conditionMessage(x); invokeRestart("muffleWarning") })
d_expl <- apply_sentinel(dna, "inversion", dna$chord_type == "aug")
chk("apply_sentinel: omitting `where` converts every NA",
    !anyNA(d_auto$inversion) && identical(d_auto$inversion, d_expl$inversion))
chk("apply_sentinel: omitting `where` warns, naming the risk and the remedy",
    grepl("genuine missing data", w) && grepl("`where`", w) && grepl("40", w))
chk("apply_sentinel: supplying `where` converts silently",
    { w2 <- NULL
      withCallingHandlers(apply_sentinel(dna, "inversion", dna$chord_type == "aug"),
        warning = function(x) { w2 <<- conditionMessage(x); invokeRestart("muffleWarning") })
      is.null(w2) })
dmix <- dna; dmix$inversion[5] <- NA
chk("apply_sentinel: with `where`, genuine missingness is still refused",
    grepl("genuine missing data",
          err_of(apply_sentinel(dmix, "inversion", dmix$chord_type == "aug"))))
chk("apply_sentinel: the sentinel is placed first among the levels",
    identical(levels(d_auto$inversion)[1], "none"))
chk("apply_sentinel: numeric variables take 0 by default",
    { dn <- data.frame(v = factor(c("none", "p")), r = c(NA, 5))
      identical(apply_sentinel(dn, "r")$r, c(0, 5)) })

## ---- nest_summary(): the fit in the original parameterization -------------
ns <- nest_summary(mf, spec = sp)
cc <- coef(mc)[!is.na(coef(mc))]
i <- match(names(cc), ns$term)
chk("nest_summary: effect estimates match a directly fitted chain model",
    !anyNA(i) && max(abs(ns$estimate[i] - cc)) < 1e-10)
chk("nest_summary: standard errors match it too",
    max(abs(ns$std.error[i] - sqrt(diag(vcov(mc)))[names(cc)])) < 1e-10)
chk("nest_summary: an additive covariate is a common slope, left as fitted",
    identical(ns$meaning[ns$term == "training"], "common slope") &&
    abs(ns$estimate[ns$term == "training"] - coef(mf)[["training"]]) < 1e-10)
chk("nest_summary: each row states the conditions it equals",
    identical(ns$meaning[ns$term == "chord_typemaj"], "-aug.none + maj.2") &&
    identical(ns$meaning[ns$term == "(Intercept)"], "aug.none"))
nsc <- nest_summary(mf, "cells", spec = sp)
chk("nest_summary: cell space returns the cell means themselves",
    nrow(nsc) == 11 && identical(nsc$meaning[1], "aug.none") &&
    abs(nsc$estimate[1] - coef(mf)[["cellaug.none"]]) < 1e-10)
nso <- nest_summary(mo, spec = spo)
chk("nest_summary: works on an ordinal fit, where the coding differs",
    abs(nso$estimate[nso$term == "chord_typemaj"] - 0.4700) < 1e-3)
chk("nest_summary: the absorbed reference condition is flagged, not left at a bare zero",
    grepl("absorbed into the thresholds", nso$meaning[1]))

## ---- the fit carries its declaration --------------------------------------
chk("nest_fit: the spec travels with the model",
    inherits(attr(mf, "nestimand_spec"), "nesting_spec") &&
    identical(attr(mf, "nestimand_spec_name"), "sp"))
chk("nest_summary: the model alone is enough",
    isTRUE(all.equal(nest_summary(mf)$estimate, nest_summary(mf, spec = sp)$estimate)))
chk("estimand: the model alone is enough; `spec` is only for outside fits",
    isTRUE(all.equal(
      as.data.frame(estimand(mf, chord_type, bounds = FALSE, self_check = FALSE))$estimate,
      as.data.frame(estimand(mf, chord_type, spec = sp, bounds = FALSE,
                             self_check = FALSE))$estimate)))
chk("estimand: the emitted code names the spec as it was named at fitting",
    any(grepl("nest_policy(sp, \"chord_type\"",
        attr(estimand(mf, chord_type, bounds = FALSE, self_check = FALSE),
             "nestimand")$code, fixed = TRUE)))
chk("latent_estimand: the model alone is enough",
    abs(latent_estimand(mf, "chord_type")$estimate[3] + 0.6779) < 1e-4)
chk("a model fitted outside nest_fit is refused, with the reason",
    grepl("does not carry one", err_of(nest_summary(m))))
chk("a model fitted outside nest_fit needs `spec`, and says so",
    grepl("does not carry one",
          err_of(estimand(m, chord_type, bounds = FALSE, self_check = FALSE))))

## ---- covariate slopes translate like the means -----------------------------
## A covariate interacting with the conditions gives one slope per cell, which
## is a vector over the same space as the means and translates the same way: a
## reference-cell slope plus differences from it, exactly as a chain fit reports.
spi <- nesting_spec(dat, response ~ chord_type * inversion * training,
                    "inversion %in% chord_type")
mi <- nest_fit(spi)
si <- nest_summary(mi)
mci <- lm(response ~ chord_type + chord_type:inversion + training +
          chord_type:training + chord_type:inversion:training, data = spi$data)
cci <- coef(mci)[!is.na(coef(mci))]
ii <- match(names(cci), si$term)
chk("slopes: every chain coefficient is reproduced, means and slopes alike",
    !anyNA(ii) && length(cci) == 20 && max(abs(si$estimate[ii] - cci)) < 1e-10)
chk("slopes: the reference-cell slope is named for the covariate alone",
    "training" %in% si$term &&
    abs(si$estimate[si$term == "training"] -
        coef(mi)[[paste0("cell", levels(spi$data$cell)[1], ":training")]]) < 1e-10)
chk("slopes: their meanings mark them as slopes",
    identical(si$meaning[si$term == "chord_typemaj:training"],
              "-aug.none + maj.2 (slope on training)"))
chk("slopes: cell space names them per condition",
    identical(nest_summary(mi, "cells")$term[11],
              "aug.none slope on training"))

## ---- a self-fitted model must correspond to the declaration ---------------
## The estimand functions work from predictions on realized cells, so the
## parameterization does not matter: a crossed or chain fit gives the same
## answer as a cell fit. What does matter is that the model contains the
## declared structure at all.
m_cross <- suppressWarnings(lm(response ~ chord_type * inversion + training, data = sp$data))
m_chain <- lm(cell_formula(sp, "effects"), data = sp$data)
val3 <- function(mm) suppressWarnings(as.data.frame(
  estimand(mm, chord_type, spec = sp, bounds = FALSE, self_check = FALSE))$estimate[3])
chk("self-fitted: a crossed fit gives the same estimand as the cell fit",
    abs(val3(m_cross) + 0.6779) < 1e-4)
chk("self-fitted: a chain fit does too",
    abs(val3(m_chain) + 0.6779) < 1e-4)
chk("self-fitted: a model missing the nested variable is refused",
    grepl("does not contain `inversion`",
          err_of(estimand(lm(response ~ chord_type + training, data = sp$data),
                          chord_type, spec = sp, bounds = FALSE, self_check = FALSE))))
chk("self-fitted: the refusal explains what the answer would have been",
    grepl("weighting the fit already implies",
          err_of(estimand(lm(response ~ chord_type + training, data = sp$data),
                          chord_type, spec = sp, bounds = FALSE, self_check = FALSE))))
chk("nest_summary: refuses a model with no cell coefficients",
    grepl("not fitted in the cell parameterization",
          err_of(nest_summary(m_chain, spec = sp))))

## ---- update() keeps the declaration ---------------------------------------
chk("nest_fit: the fit gains a class, ahead of the engine's own",
    identical(class(mf)[1], "nestimand_fit") && inherits(mf, "lm"))
mu <- update(mf, . ~ . - training)
chk("update: the declaration survives the refit",
    inherits(attr(mu, "nestimand_spec"), "nesting_spec") &&
    identical(class(mu)[1], "nestimand_fit") && length(coef(mu)) == 10)
chk("update: the estimand works from the updated fit, without `spec`",
    abs(as.data.frame(estimand(mu, chord_type, bounds = FALSE,
                               self_check = FALSE))$estimate[3] + 0.6779) < 1e-4)
chk("update: the code view records the new call, not the old one",
    any(grepl("updated from the original call", attr(mu, "nestimand_code"))) &&
    !any(grepl("training", attr(mu, "nestimand_code"))))
chk("update: an update that drops the structure is caught downstream",
    grepl("does not contain",
          err_of(estimand(update(mf, . ~ . - cell), chord_type))))
chk("class: the engine's own methods still work",
    length(predict(mf, newdata = head(sp$data))) == 6 &&
    inherits(summary(mf), "summary.lm") && nrow(vcov(mf)) == 11)

cat(sprintf("\n%d passed, %d failed\n", pass, fail))
