## nestimand: regression suite for the cell-translation core ---------------
## Run from the folder holding the package source files:
##   source("test_core.R")
suppressPackageStartupMessages({library(marginaleffects)})
## --- load the package sources ---------------------------------------------
## Finds the seven source files by name: in the working directory, in R/, or
## anywhere below. The folder layout does not matter.
local({
  want <- c("spec.R", "translate.R", "policy.R", "estimand.R", "fit.R",
            "latent.R", "priors.R", "summary.R", "random.R")
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

## fit_from() takes the arguments a declaration is made from, not a declaration
## object. Most of what follows tests what happens after the fit, and builds one
## declaration to ask several questions of, so this restates that declaration's
## own arguments as a fit_from() call rather than repeating them at every site.
.fit_from_n <- 0L
fit_from <- function(spec, ...) {
  ## the frame the emitted code names has to outlive this call, since a saved
  ## code block is re-run later, in an environment of the test's own making
  .fit_from_n <<- .fit_from_n + 1L
  nm <- paste0(".dat", .fit_from_n)
  assign(nm, spec$data[setdiff(names(spec$data), spec$cell_name)],
         envir = globalenv())
  e <- new.env(parent = parent.frame())
  args <- list(quote(nest_fit), spec$fit, spec$formula_in, as.name(nm))
  nn <- spec_nests(spec)
  if (length(nn)) args$nests <- nn
  if (!is.null(spec$family)) args$family <- spec$family
  if (!is.null(spec$random_original)) args$random <- spec$random_original
  eval(as.call(c(args, as.list(match.call(expand.dots = FALSE)$...))), e)
}

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
## structurally undefined rows are converted now, so what is refused is an NA
## the structure does not account for: one inside a stratum that has values
bad <- dat; bad$inversion[bad$chord_type == "aug"] <- NA
bad$inversion[which(bad$chord_type == "dim")[1]] <- NA
chk("nesting_spec: an NA the structure does not account for is still refused",
    grepl("missing data rather than structure",
          err_of(nesting_spec(bad, response ~ chord_type * inversion,
                              "inversion %in% chord_type")), fixed = TRUE))
d3 <- data.frame(v = factor(rep(c("none","p"), 50)), r = c(rep(0,50), runif(50,4,8)))
d3$y <- rnorm(100)
chk("nesting_spec: continuous leaf accepted",
    inherits(nesting_spec(d3, y ~ v * r, "r %in% v"), "nesting_spec"))
chk("nesting_spec: continuous parent refused",
    grepl("cannot nest", err_of(nesting_spec(d3, y ~ v * r, "v %in% r"))))
sp_cl <- dat; sp_cl$cell <- 1
chk("nesting_spec: a `cell` column already in the data is stepped around",
    identical(suppressMessages(nesting_spec(sp_cl,
      response ~ chord_type * inversion, "inversion %in% chord_type"))$cell_name,
      ".cell"))
chk("nesting_spec: and the user's own column is left as it was",
    identical(suppressMessages(nesting_spec(sp_cl,
      response ~ chord_type * inversion, "inversion %in% chord_type"))$data$cell,
      sp_cl$cell))

## A design with nothing nested is a legitimate thing to hand the package, and
## the answer is that it has nothing to add - but it used to say so by claiming
## every declared nesting variable was continuous, one line after reporting,
## correctly, that nothing was nested. Two arrivals at one branch, one message.
chk("crossed: a design with no nesting says so, and names both ways on",
    { e <- err_of(suppressMessages(nesting_spec(
        within(dat, inversion <- factor(as.character(inversion))),
        response ~ chord_type * inversion)))
      grepl("nothing is nested here", e, fixed = TRUE) &&
        grepl("Fit it directly", e, fixed = TRUE) &&
        grepl("name the design variables in `nests` without `%in%`", e,
              fixed = TRUE) })
chk("crossed: and it names the variables to put there",
    grepl('nests = c("chord_type", "inversion")',
          err_of(suppressMessages(nesting_spec(dat, response ~ chord_type * inversion))),
          fixed = TRUE))
chk("crossed: taking that advice gives a working declaration",
    { sp_x <- suppressMessages(nesting_spec(dat, response ~ chord_type * inversion,
                nests = c("chord_type", "inversion")))
      inherits(sp_x, "nesting_spec") && nrow(sp_x$cells) == 10 })
chk("crossed: the note no longer promises a design the error then refuses",
    { z <- utils::capture.output(
        try(nesting_spec(dat, response ~ chord_type * inversion), silent = TRUE),
        type = "message")
      any(grepl("nothing is nested", z)) &&
        !any(grepl("the design is fully crossed[.]$", z)) })

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

## ---- the two fitting modes ----------------------------------------------
m  <- fit_from(sp)
chk("fit: cell parameterization is full rank, no aliased coefficients",
    sum(is.na(coef(m))) == 0 && m$rank == length(coef(m)))

## ---- the effect basis ----------------------------------------------------
A <- effect_basis(sp)
chk("effect_basis: square and full rank (10 x 10)",
    nrow(A) == 10 && ncol(A) == 10 && qr(A)$rank == 10)
## the basis is the package's own choice of identified effects, so the identity
## is checked against the quantities it produces rather than against whichever
## columns lm's pivot happens to keep
mf0 <- fit_from(sp)
eff <- nest_summary(mf0)
cel <- nest_summary(mf0, "cells")
b  <- eff$estimate[match(colnames(A), eff$term)]
mu <- cel$estimate[match(rownames(A), cel$term)]
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
    abs(val(est(nest_policy(sp, "chord_type", "equal")), "maj - aug") - 0.6779) < 1e-4)
chk("policy proportional: -0.6779 (balanced data)",
    abs(val(est(nest_policy(sp, "chord_type", "proportional")), "maj - aug") - 0.6779) < 1e-4)
chk("policy: `counterfactual` no longer names a policy, and says what to use",
    grepl("route = \"g_computation\"",
          err_of(nest_policy(sp, "chord_type", "counterfactual"))))
chk("policy hierarchical: -0.6779 at depth one",
    abs(val(est(nest_policy(sp, "chord_type", "hierarchical")), "maj - aug") - 0.6779) < 1e-4)
vert <- vapply(c("0","1","2"), function(iv)
  val(est(nest_policy(sp, "chord_type", "nominated", at = c(inversion = iv))),
      "maj - aug"), 1)
chk("policy nominated: three single-version contrasts",
    all(abs(vert - c(1.0479, 0.7405, 0.2452)) < 1e-4))
p3 <- nest_policy(sp, "chord_type", c("0" = 0.5, "1" = 0.3, "2" = 0.2))
chk("policy supplied: exactly the convex combination of the vertices",
    abs(val(est(p3), "maj - aug") - sum(vert * c(0.5, 0.3, 0.2))) < 1e-10)
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
chk("policy: bounds are the range of the vertices",
    abs(max(vert) - 1.0479) < 1e-4 && abs(min(vert) - 0.2452) < 1e-4)

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
## The reorder self-check is gone. It re-ran the estimand with levels permuted
## and warned if the answer moved; the answer cannot move, because the columns a
## pivot drops as redundant leave the column space unchanged, so the fitted
## values on realized conditions are invariant to level order. The property is
## asserted here instead of on every call - which is what the check beneath this
## comment has always done, independently of it.
chk("removed: `self_check` says so rather than being swept into the dots",
    { e <- err_of(suppressMessages(nest_estimand(m, chord_type, self_check = FALSE)))
      grepl("`self_check` has been removed", e, fixed = TRUE) &&
        grepl("Remove the argument", e, fixed = TRUE) })
chk("removed: and so does `spec`",
    grepl("`spec` has been removed",
          err_of(suppressMessages(nest_estimand(m, chord_type, spec = 1))),
          fixed = TRUE))
chk("removed: neither appears in the signature",
    !any(c("self_check", "spec") %in% names(formals(nest_estimand))))

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
    abs(val(e2d, "maj - aug") - 0.6283) < 1e-4)
eq2d <- est(nest_policy(sp2d, "chord", "equal"), tgt = "chord",
            model = m2d, spec = sp2d)
chk("depth 2: equal differs from hierarchical (leaves vs tree)",
    abs(val(eq2d, "maj - aug") - val(e2d, "maj - aug")) > 1e-3)

## ---- nest_estimand(): the function, and the code it kept ----------------------
set.seed(3)
e <- nest_estimand(m, chord_type, policy = "equal")
chk("estimand: aug - maj = -0.6779 in original labels",
    abs(as.data.frame(e)$estimate[as.data.frame(e)$term == "maj - aug"] - 0.6779) < 1e-4)
meta <- attr(e, "nestimand")
chk("estimand: provenance recorded (policy, contrast, build)",
    identical(meta$policy, "equal") && identical(meta$contrast, "pairwise") &&
    identical(meta$build, nestimand_build))
chk("estimand: bounds attached, -1.048 to -0.245 for aug - maj",
    abs(meta$bounds$policy_high[meta$bounds$term == "maj - aug"] - 1.0479) < 1e-4 &&
    abs(meta$bounds$policy_low[meta$bounds$term == "maj - aug"] - 0.2452) < 1e-4)
## fidelity: the code view is the code that ran, so re-running it must agree
env <- new.env(); assign("m", m, env); assign("sp", sp, env)
re <- eval(parse(text = paste(meta$code, collapse = "\n")), envir = env)
chk("show_code: the saved code re-runs and reproduces the estimand",
    isTRUE(all.equal(as.data.frame(re)$estimate, as.data.frame(e)$estimate)))
chk("show_code: prints the lines and returns them invisibly",
    is.character(capture.output(show_code(e))[1]))
## non-core arguments reach the destination function, and appear in the code
e9 <- nest_estimand(m, chord_type, bounds = FALSE,
               conf_level = 0.9)
chk("estimand: `...` passed through to marginaleffects",
    any(grepl("conf_level = 0.9", attr(e9, "nestimand")$code, fixed = TRUE)) &&
    as.data.frame(e9)$conf.low[1] > as.data.frame(e)$conf.low[1])
## contrasts inside each stratum, which is `by =` naming the target's parent
ew <- suppressMessages(nest_estimand(m, inversion, by = "chord_type",
                                bounds = FALSE))
dw <- as.data.frame(ew)
chk("estimand by parent: three strata x three contrasts, no sentinel",
    nrow(dw) == 9 && !any(grepl("none", dw$term)) &&
    setequal(unique(dw$chord_type), c("dim", "min", "maj")))
chk("estimand: a non-nesting target is refused",
    grepl("need no policy", err_of(nest_estimand(m, training))))
chk("estimand: supplied policy is written into the code verbatim",
    any(grepl('c("0" = 0.5, "1" = 0.3, "2" = 0.2)',
        attr(nest_estimand(m, chord_type, policy = c("0" = .5, "1" = .3, "2" = .2),
             bounds = FALSE), "nestimand")$code, fixed = TRUE)))

## ---- declaring and fitting in one call ------------------------------------
## The declaration and the fit are one step unless the user wants them apart,
## and the structure is read off the data when it is not given: a variable is
## structurally undefined exactly where its parent's levels say, and that is
## decidable rather than something to be asked for.
one_d <- local({
  one_rows <- list()
  for (one_ct in c("aug", "dim", "min", "maj")) {
    one_ivs <- if (one_ct == "aug") NA else c("0", "1", "2")
    for (one_iv in one_ivs) for (one_tp in c("t1", "t2"))
      one_rows[[length(one_rows) + 1]] <-
        data.frame(chord_type = one_ct, inversion = one_iv, top = one_tp)
  }
  one_z <- do.call(rbind, one_rows); one_z <- one_z[rep(seq_len(nrow(one_z)), each = 6), ]
  for (one_v in c("chord_type", "inversion", "top")) one_z[[one_v]] <- factor(one_z[[one_v]])
  set.seed(2); one_z$response <- rnorm(nrow(one_z), 4); one_z
})
one_sent <- set_sentinel(one_d, "inversion", where = is.na(one_d$inversion))
chk("infer: the structure is read off the structural NAs",
    identical(infer_nests(one_d, c("chord_type", "inversion", "top"))$nests,
              "inversion %in% chord_type"))
chk("infer: a chain is found, the middle variable being absent too",
    { one_z <- do.call(rbind, lapply(c("aug", "dim"), function(one_ct) {
        one_ivs <- if (one_ct == "aug") NA else c("0", "1", "2")
        do.call(rbind, lapply(one_ivs, function(one_iv) {
          one_zs <- if (is.na(one_iv) || one_iv == "2") NA else c("z1", "z2")
          do.call(rbind, lapply(one_zs, function(one_q)
            data.frame(chord_type = one_ct, inversion = one_iv, Z = one_q))) })) }))
      identical(infer_nests(one_z, c("chord_type", "inversion", "Z"))$nests,
                c("inversion %in% chord_type", "Z %in% inversion")) })
chk("infer: two siblings take the parent, not each other",
    { one_z <- do.call(rbind, lapply(c("aug", "dim", "maj"), function(one_ct) {
        one_a <- if (one_ct == "aug") NA else c("0", "1")
        one_b <- if (one_ct == "aug") NA else c("p", "q")
        expand.grid(chord_type = one_ct, inversion = one_a, X1 = one_b) }))
      identical(infer_nests(one_z, c("chord_type", "inversion", "X1"))$nests,
                c("inversion %in% chord_type", "X1 %in% chord_type")) })
chk("infer: a stray NA is missing data, not structure, and says so",
    { one_z <- one_d; one_z$inversion[which(one_z$chord_type == "dim")[1]] <- NA
      one_r <- infer_nests(one_z, c("chord_type", "inversion", "top"))
      !length(one_r$nests) && grepl("not structural", one_r$notes[1], fixed = TRUE) })
chk("infer: two variables explaining the same gaps is reported, not guessed",
    { one_z <- one_d; one_z$mode <- factor(ifelse(one_z$chord_type == "aug", "x", "y"))
      length(infer_nests(one_z, c("chord_type", "inversion", "top", "mode"))$ambiguous) == 1 })
chk("infer: nesting_spec builds the same declaration either way",
    { one_a <- suppressMessages(nesting_spec(one_d, response ~ chord_type * inversion * top))
      one_b <- suppressMessages(nesting_spec(one_sent, response ~ chord_type * inversion * top,
                                         "inversion %in% chord_type"))
      identical(sort(as.character(one_a$cells$cell)), sort(as.character(one_b$cells$cell))) &&
        identical(sentinel_levels(one_a), sentinel_levels(one_b)) &&
        identical(one_a$parent, one_b$parent) })
chk("infer: it refuses rather than inventing structure from missing data",
    { one_z <- one_d; one_z$inversion[which(one_z$chord_type == "dim")[1]] <- NA
      grepl("not structural",
            err_of(nesting_spec(one_z, response ~ chord_type * inversion * top)),
            fixed = TRUE) })
chk("one call: the structure read from the data fits the model declared by hand",
    { one_a <- suppressMessages(nest_fit("lm",
              response ~ chord_type * (inversion + top), one_d))
      one_b <- suppressMessages(nest_fit("lm",
              response ~ chord_type * (inversion + top), one_sent,
              nests = "inversion %in% chord_type"))
      isTRUE(all.equal(unname(coef(one_a)), unname(coef(one_b)))) && !anyNA(coef(one_a)) })
chk("one call: the declaration is part of the code, which re-runs",
    { one_a <- suppressMessages(nest_fit("lm",
              response ~ chord_type * (inversion + top), one_d))
      one_cd <- attr(one_a, "nestimand_code")
      one_env <- new.env(parent = globalenv()); assign("one_d", one_d, one_env)
      one_r <- suppressMessages(eval(parse(text = paste(c(one_cd, "m"), collapse = "\n")), one_env))
      any(grepl("spec <- nesting_spec", one_cd)) &&
        isTRUE(all.equal(coef(one_a), coef(one_r))) })
chk("one call: the data is named as the user named it",
    { one_cd <- attr(suppressMessages(nest_fit("lm",
              response ~ chord_type * inversion, one_d, dry_run = TRUE)),
              "nestimand_code")
      any(grepl("data = one_d", one_cd, fixed = TRUE)) })
chk("one call: `nests` stays unquoted, as it is written",
    { one_cd <- attr(suppressMessages(nest_fit("lm",
              response ~ chord_type * inversion, one_sent,
              nests = inversion %in% chord_type, dry_run = TRUE)), "nestimand_code")
      any(grepl("nests = inversion %in% chord_type", one_cd, fixed = TRUE)) })
chk("one call: the engine, the formula and the data are all needed",
    grepl("takes an engine, a formula and the data",
          err_of(nest_fit("lm")), fixed = TRUE))
chk("one call: a declaration object is not what nest_fit takes",
    { one_sp <- suppressMessages(nesting_spec(one_sent,
              response ~ chord_type * inversion, "inversion %in% chord_type"))
      nzchar(err_of(suppressMessages(
        nest_fit(one_sp, response ~ chord_type, one_sent)))) })

## ---- nest_fit(): the fitting side ----------------------------------------
mf <- fit_from(sp)
chk("nest_fit: fits the cell parameterization, full rank",
    inherits(mf, "lm") && sum(is.na(coef(mf))) == 0)
chk("nest_fit: the call travels with the fit",
    any(grepl("lm(response ~ 0 + cell + training, data = spec$data)",
              attr(mf, "nestimand_code"), fixed = TRUE)))
chk("nest_fit: the parameterization and its reason are stated",
    any(grepl("parameterization: cells", attr(mf, "nestimand_code"))))
ef <- nest_estimand(mf, chord_type, policy = "equal", bounds = FALSE)
chk("nest_fit: the code view joins fit and estimand into one pipeline",
    any(grepl("^mf <- lm", attr(ef, "nestimand")$code)) &&
    any(grepl("avg_predictions", attr(ef, "nestimand")$code)))
dr <- fit_from(sp, dry_run = TRUE, weights = NULL)
chk("nest_fit: dry_run returns code without fitting", inherits(dr, "nestimand_code"))

## random-effects translation
spl <- nesting_spec(dat, response ~ chord_type * inversion + training +
                    (chord_type * inversion + training | participant),
                    "inversion %in% chord_type", fit = "lmer")
chk("random_terms: a slope over the structure gets the columns the design has",
    identical(random_terms(spl),
      paste0("(1 + ", paste(reduced_columns(spl), collapse = " + "),
             " + training | participant)")))
chk("random_terms: which is the same width as the mean structure",
    length(reduced_columns(spl)) + 1L == nrow(spl$cells))
spi <- nesting_spec(dat, response ~ chord_type * inversion + (1 | participant),
                    "inversion %in% chord_type", fit = "lmer")
chk("random_terms: a clean intercept term passes through unchanged",
    identical(random_terms(spi), "(1 | participant)"))
sp_grp <- nesting_spec(dat, response ~ chord_type * inversion +
                    (1 | participant) + (1 | chord_type:inversion),
                    "inversion %in% chord_type", fit = "lmer")
chk("random_terms: a grouping term is the formula's, and is left as written",
    identical(random_terms(sp_grp),
              "(1 | participant) + (1 | chord_type:inversion)"))

## ordinal engines: thresholds already carry an intercept
do2 <- dat; do2$rating <- factor(round(pmin(pmax(do2$response, 1), 7)), ordered = TRUE)
spo <- nesting_spec(do2, rating ~ chord_type * inversion,
                    "inversion %in% chord_type", fit = "clm")
chk("ordinal: cell coding keeps the implicit intercept, no redundant parameter",
    identical(paste(deparse(cell_formula(spo)), collapse = " "), "rating ~ cell"))
mo <- fit_from(spo)
chk("ordinal: clm fits without the intercept warning, 15 parameters",
    length(coef(mo)) == 15)
chk("ordinal: the response scale is probed, and refused clearly if unavailable",
    { e <- err_of(nest_estimand(mo, chord_type, bounds = FALSE))
      identical(e, "") || grepl("scale = \"latent\"", e) })
spbo <- nesting_spec(do2, rating ~ chord_type * inversion,
                     "inversion %in% chord_type", fit = "brm", family = "cumulative()")
chk("ordinal: brms cumulative also takes the threshold-aware coding",
    identical(paste(deparse(cell_formula(spbo)), collapse = " "), "rating ~ cell"))
chk("ordinal: brms family threaded into the emitted call",
    grepl("family = cumulative()", fit_from(spbo, dry_run = TRUE), fixed = TRUE))
chk("brms: non-core arguments reach the engine and appear in the code",
    grepl("chains = 4, iter = 2000, seed = 1",
          fit_from(nesting_spec(dat, response ~ chord_type * inversion,
                   "inversion %in% chord_type", fit = "brm"),
                   dry_run = TRUE, chains = 4, iter = 2000, seed = 1), fixed = TRUE))

## ---- latent scale: estimands as linear functionals -----------------------
lp <- latent_estimand(m, "chord_type", "equal", spec = sp)
ap <- as.data.frame(nest_estimand(m, chord_type, bounds = FALSE))
chk("latent: agrees with the prediction route on estimates (1e-10)",
    max(abs(lp$estimate[match(ap$term, lp$term)] - ap$estimate)) < 1e-10)
chk("latent: agrees with the prediction route on standard errors (1e-7)",
    max(abs(lp$std.error[match(ap$term, lp$term)] - ap$std.error)) < 1e-7)
lv <- latent_estimand(m, "chord_type", c("0" = 0.5, "1" = 0.3, "2" = 0.2), spec = sp)
chk("latent: a supplied policy is the convex combination of the vertices",
    abs(lv$estimate[lv$term == "maj - aug"] - sum(vert * c(0.5, 0.3, 0.2))) < 1e-8)
chk("latent: reference and sequential contrasts available",
    nrow(latent_estimand(m, "chord_type", "equal", contrast = "reference", spec = sp)) == 3 &&
    nrow(latent_estimand(m, "chord_type", "equal", contrast = "sequential", spec = sp)) == 3)
lo <- nest_estimand(mo, chord_type, policy = "equal")
chk("latent: ordinal fit returns one number per contrast, where predictions fail",
    nrow(as.data.frame(lo)) == 6 && is.finite(as.data.frame(lo)$estimate[3]))
chk("latent: bounds computed on the latent scale too",
    !is.null(attr(lo, "nestimand")$bounds) &&
    all(attr(lo, "nestimand")$bounds$policy_low <=
        attr(lo, "nestimand")$bounds$estimate + 1e-10))
chk("latent: the code view names the linear-map route",
    any(grepl("latent_estimand", attr(lo, "nestimand")$code)))
## A within-stratum contrast crosses no boundary, so no policy weighs anything
## and it is a plain contrast of design rows - a linear map like any other. It
## used to be refused on the latent scale, with a message naming an argument
## `nest_estimand()` does not have.
## `contrast` decides one thing: which comparisons are formed among the target's
## levels. Grouping is `by`, and an interaction is a `:` target. `"within"` was
## the first spelled two ways, `"interaction"` the second.
chk("contrast: it offers the three comparison sets, and forming none",
    identical(eval(formals(nest_estimand)$contrast),
              c("pairwise", "reference", "sequential", "none")))
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
    abs(pfe$sd[pfe$parameter == "maj - aug"] - 1.5 * sqrt(1 + 1/3)) < 1e-10)
chk("prior_for_estimand: a within-family contrast is tighter",
    abs(pfe$sd[pfe$parameter == "maj - dim"] - 1.5 * sqrt(2/3)) < 1e-10)
spb2 <- nesting_spec(dat, response ~ chord_type * inversion + training,
                     "inversion %in% chord_type", fit = "brm")
prb <- nest_prior(spb2, mean = 4, sd = 1.5, on = "cells", covariate_sd = 1)
## A prior on the mean structure is written in brms syntax, in the user's own
## variables, and translated to the coefficients the model is fitted on.
pr_b <- brms::prior_string("normal(4, 1.5)", class = "b")
chk("prior: brms call carries the statement and the stanvars",
    grepl("prior = prior_statement(.prior), stanvars = prior_stanvars(.prior)",
          suppressMessages(fit_from(spb2, prior = pr_b, dry_run = TRUE)),
          fixed = TRUE))
chk("prior: the translation is stated in the emitted code",
    grepl("A D A", suppressMessages(fit_from(spb2, prior = pr_b, dry_run = TRUE)),
          fixed = TRUE))
chk("prior: and is announced, in the variables it was written in",
    any(grepl("translated to the coefficients the model is fitted on",
        utils::capture.output(fit_from(spb2, prior = pr_b, dry_run = TRUE),
                              type = "message"))))
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
              attr(nest_estimand(m, chord_type, bounds = FALSE, p_adjust = "none"),
                   "nestimand")$code, fixed = TRUE)))

## ---- engine label conventions --------------------------------------------
## marginaleffects 0.18.x returns `term` = "maj - aug" = -0.6779; 0.32.0 returns
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
chk("labels: parentheses stripped, declared order restored",
    identical(cn$term, c("dim - aug", "maj - aug", "maj - min")))
chk("labels: contrast direction fixed by the package, not the engine",
    isTRUE(all.equal(cn$estimate, co$estimate, tolerance = 1e-8)))
chk("labels: confidence bounds swapped and negated with the estimate",
    isTRUE(all.equal(cn$conf.low, co$conf.low, tolerance = 1e-8)) &&
    isTRUE(all.equal(cn$conf.high, co$conf.high, tolerance = 1e-8)))
chk("labels: the number of reversed rows is recorded",
    attr(cn, "nestimand_flipped") == 0 && attr(co, "nestimand_flipped") == 3)
chk("labels: an unrecognized label is left alone rather than mangled",
    identical(mfx_canonical(data.frame(term = "b0", estimate = 1))$term, "b0"))
chk("labels: nest_estimand() reports contrasts in declared level order",
    identical(as.data.frame(nest_estimand(m, chord_type, bounds = FALSE))$term[1:3],
              c("dim - aug", "min - aug", "maj - aug")))

chk("show_code: the normalization is part of the saved code, not applied after",
    any(grepl("mfx_canonical",
              attr(nest_estimand(m, chord_type, bounds = FALSE, p_adjust = "none"),
                   "nestimand")$code)))
chk("show_code: grouped-contrast code re-runs and reproduces its estimand",
    { ew2 <- suppressMessages(nest_estimand(m, inversion, by = "chord_type",
                      bounds = FALSE))
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
## ---- ordinal families: the coefficients are contrasts, and translate so ----
po <- nest_prior(spbo, mean = 0, sd = 1, on = "cells")
chk("prior: an ordinal fit takes a translated prior too",
    length(po$full_mean) == 9 &&
    identical(names(po$full_mean), fitted_coef_names(spbo)))
chk("prior: independent cell priors of sd s imply contrast sds of s*sqrt(2)",
    max(abs(sqrt(diag(po$full_cov)) - sqrt(2))) < 1e-10)
chk("prior: neighbouring contrasts share the reference, so they correlate 0.5",
    abs(stats::cov2cor(po$full_cov)[1, 2] - 0.5) < 1e-10)
chk("prior: the dimension matches what the ordinal model will estimate",
    isTRUE(check_prior_dimension(po, spbo)))
chk("prior: the absorbed level is disclosed when printed",
    any(grepl("carried by the .*thresholds",
              capture.output(print(po)))))
chk("prior: a prior on effects translates for an ordinal fit as well",
    length(nest_prior(spbo, mean = 0, sd = 1, on = "effects")$full_mean) == 9)

## ---- brms draw naming ----------------------------------------------------
chk("draws: exact `b_` names matched, full stops included",
    identical(draw_names(c("cellaug.none", "training"),
              c("b_cellaug.none", "b_training", "sigma")),
              c("b_cellaug.none", "b_training")))
chk("draws: punctuation-insensitive fallback if the convention changes",
    identical(draw_names("cellaug.none", c("b_cellaugnone")), "b_cellaugnone"))
chk("draws: an unmatched coefficient reports the names available",
    grepl("Draw columns available", err_of(draw_names("cellmaj.9", "b_x"))))

## ---- the random side follows the formula too -------------------------------
## Translation applies only to terms that reach across the structure. Everything
## else comes back exactly as declared: enlarging it would be the package
## overriding a deliberate modelling decision, often one made to get a model to
## converge at all.
resp_re <- function(re) nesting_spec(dat,
  stats::as.formula(paste("response ~ chord_type * inversion + training +", re)),
  "inversion %in% chord_type", fit = "lmer")
red_bar <- function(spec, extra = "")
  paste0("(1 + ", paste(reduced_columns(spec), collapse = " + "), extra,
         " | participant)")
chk("RE: a random intercept is left exactly as declared",
    identical(random_terms(resp_re("(1 | participant)")), "(1 | participant)"))
chk("RE: a slope on the nesting root is identified, so it is left alone",
    identical(random_terms(resp_re("(chord_type | participant)")),
              "(chord_type | participant)"))
chk("RE: a hand-written grouping structure passes through untouched",
    identical(random_terms(resp_re(
      "(1 | participant) + (1 | participant:chord_type)")),
      "(1 | participant) + (1 | participant:chord_type)"))
chk("RE: a slope over the whole structure gets the columns the design has",
    { sp_re <- resp_re("(chord_type * inversion | participant)")
      identical(random_terms(sp_re), red_bar(sp_re)) })
chk("RE: covariate slopes survive the translation",
    { sp_re <- resp_re("(chord_type * inversion + training | participant)")
      identical(random_terms(sp_re), red_bar(sp_re, " + training")) })
chk("RE: a term the mean structure lacks is refused, and says why",
    { e <- err_of(random_terms(nesting_spec(dat,
        response ~ chord_type + inversion + (chord_type * inversion | participant),
        "inversion %in% chord_type", fit = "lmer")))
      grepl("the model formula does not contain", e, fixed = TRUE) })
chk("RE: boundary_vars names the nested variables, not the roots",
    identical(boundary_vars(sp), "inversion"))
chk("RE: the two sides of the model have the same width",
    { sp_re <- resp_re("(chord_type * inversion | participant)")
      length(reduced_columns(sp_re)) + 1L == nrow(sp_re$cells) })
chk("RE: and the fitted structure is what the formula asked for",
    { sp_re <- suppressMessages(nesting_spec(dat,
        response ~ chord_type + inversion + (chord_type + inversion | participant),
        "inversion %in% chord_type", fit = "lmer"))
      identical(random_terms(sp_re),
                paste0("(1 + ", paste(reduced_columns(sp_re), collapse = " + "),
                       " | participant)")) })

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
chk("diag(): the wrapper survives translation",
    { sp_re <- mk_re("diag(chord_type * inversion | participant)")
      identical(random_terms(sp_re),
        paste0("diag(1 + ", paste(reduced_columns(sp_re), collapse = " + "),
               " | participant)")) })
chk("diag(): an untranslated term keeps its wrapper too",
    identical(random_terms(mk_re("diag(1 + chord_type | participant)")),
              "diag(1 + chord_type | participant)"))
chk("diag(): a plain bar is unaffected",
    { sp_re <- mk_re("(chord_type * inversion | participant)")
      identical(random_terms(sp_re),
        paste0("(1 + ", paste(reduced_columns(sp_re), collapse = " + "),
               " | participant)")) })

## ---- the sentinel is applied by the declaration, not by the user -----------
## `apply_sentinel()` is gone. NA still cannot be left as it is - R deletes those
## rows casewise and silently - but who converts it changed: the declaration says
## which rows should be undefined, or the inference establishes it, and either
## way there is an answer to check the data against. Asking the user to do it by
## hand had no such answer, so `where` could be wrong and the warning for
## omitting it was the whole safeguard.
dna <- dat; dna$inversion[dna$chord_type == "aug"] <- NA
chk("sentinel: a declared structure converts its own undefined rows",
    { sp_na <- suppressMessages(nesting_spec(dna, response ~ chord_type * inversion,
                                             "inversion %in% chord_type"))
      !anyNA(sp_na$data$inversion) &&
        identical(levels(sp_na$data$inversion)[1], "none") &&
        nrow(sp_na$cells) == 10 })
chk("sentinel: and gives the same declaration as data already carrying one",
    { sp_na <- suppressMessages(nesting_spec(dna, response ~ chord_type * inversion,
                                             "inversion %in% chord_type"))
      identical(sort(as.character(sp_na$cells$cell)), sort(as.character(sp$cells$cell))) })
chk("sentinel: an NA the structure does not account for is refused",
    { dmix <- dna; dmix$inversion[5] <- NA
      grepl("missing data rather than structure",
            err_of(nesting_spec(dmix, response ~ chord_type * inversion,
                                "inversion %in% chord_type")), fixed = TRUE) })
chk("sentinel: a value where the structure says there is none is refused",
    { dbad <- dna; dbad$inversion <- as.character(dbad$inversion)
      dbad$inversion[which(dbad$chord_type == "aug")[1]] <- "0"
      dbad$inversion <- factor(dbad$inversion)
      grepl("missing data rather than structure",
            err_of(nesting_spec(dbad, response ~ chord_type * inversion,
                                "inversion %in% chord_type")), fixed = TRUE) })
chk("sentinel: an NA on a variable nested in nothing has no structure to check",
    { droot <- dat; droot$chord_type <- as.character(droot$chord_type)
      droot$chord_type[1] <- NA; droot$chord_type <- factor(droot$chord_type)
      grepl("nested in nothing",
            err_of(nesting_spec(droot, response ~ chord_type * inversion,
                                "inversion %in% chord_type")), fixed = TRUE) })
chk("sentinel: a numeric nested variable takes 0",
    { dn <- data.frame(g = factor(rep(c("a", "b"), each = 2)),
                       r = c(NA, NA, 5, 7), y = rnorm(4))
      sp_n <- suppressMessages(nesting_spec(dn, y ~ g * r, "r %in% g"))
      identical(sp_n$data$r, c(0, 0, 5, 7)) })
chk("sentinel: apply_sentinel is no longer exported",
    !("apply_sentinel" %in% ls(envir = globalenv())) ||
      !is.function(tryCatch(get("apply_sentinel"), error = function(e) NULL)))

## ---- nest_summary(): the fit in the original parameterization -------------
ns <- nest_summary(mf)
chk("nest_summary: an additive covariate is a common slope, left as fitted",
    identical(ns$meaning[ns$term == "training"], "common slope") &&
    abs(ns$estimate[ns$term == "training"] - coef(mf)[["training"]]) < 1e-10)
chk("nest_summary: each row states the conditions it equals",
    identical(ns$meaning[ns$term == "chord_typemaj"],
              paste0("-aug.none + maj.", reference_levels(sp)[["inversion"]])) &&
    identical(ns$meaning[ns$term == "(Intercept)"], "aug.none"))
nsc <- nest_summary(mf, "cells")
chk("nest_summary: cell space returns the cell means themselves",
    nrow(nsc) == 11 && identical(nsc$meaning[1], "aug.none") &&
    abs(nsc$estimate[1] - coef(mf)[["cellaug.none"]]) < 1e-10)
nso <- nest_summary(mo)
chk("nest_summary: works on an ordinal fit, where the coding differs",
    is.finite(nso$estimate[nso$term == "chord_typemaj"]) &&
    identical(nso$meaning[nso$term == "chord_typemaj"],
              paste0("-aug.none + maj.", reference_levels(spo)[["inversion"]])))
chk("nest_summary: the absorbed reference condition is flagged, not left at a bare zero",
    grepl("absorbed into the thresholds", nso$meaning[1]))

## ---- the fit carries its declaration --------------------------------------
chk("nest_fit: the declaration travels with the model",
    inherits(attr(mf, "nestimand_spec"), "nesting_spec") &&
    identical(attr(mf, "nestimand_spec_name"), "spec"))
chk("nest_summary: the model alone is enough",
    is.data.frame(as.data.frame(nest_summary(mf))))
chk("estimand: the model alone is enough",
    is.data.frame(as.data.frame(
      nest_estimand(mf, chord_type, bounds = FALSE))))
chk("estimand: the emitted code names the declaration the fit carries",
    any(grepl("nest_policy(spec, \"chord_type\"",
        attr(nest_estimand(mf, chord_type, bounds = FALSE),
             "nestimand")$code, fixed = TRUE)))
chk("latent_estimand: the model alone is enough",
    abs(latent_estimand(mf, "chord_type")$estimate[3] - 0.6779) < 1e-4)
chk("a model fitted outside nest_fit is refused, with the reason",
    grepl("was not fitted by nest_fit",
          err_of(nest_summary(lm(response ~ training, data = dat)))))
chk("a model fitted outside nest_fit cannot be given a declaration either",
    grepl("was not fitted by nest_fit",
          err_of(nest_estimand(lm(response ~ training, data = dat), chord_type,
                               bounds = FALSE))))

## ---- covariate slopes translate like the means -----------------------------
## A covariate interacting with the conditions gives one slope per cell, which
## is a vector over the same space as the means and translates the same way: a
## reference-cell slope plus differences from it, exactly as a chain fit reports.
spi <- nesting_spec(dat, response ~ chord_type * inversion * training,
                    "inversion %in% chord_type")
mi <- fit_from(spi)
si <- nest_summary(mi)
## the same comparison, on the package's basis rather than lm's pivot
## The slope block is reported as a reference slope plus differences from it,
## the convention a chain fit with a covariate main effect uses. The check that
## matters is therefore the mapping: carrying the effect-space slopes back
## through A must reproduce the per-cell slopes of the cell fit.
Ai <- effect_basis(spi)
sl <- si$estimate[grepl("slope on training", si$meaning)]
cell_sl <- nest_summary(mi, "cells")$estimate[
  grepl("slope on training", nest_summary(mi, "cells")$meaning)]
chk("slopes: 10 effect-space slopes, one per realized condition",
    length(sl) == 10 && length(cell_sl) == 10 && nrow(si) == 20)
chk("slopes: carrying them back through A reproduces the per-cell slopes",
    max(abs(as.numeric(Ai %*% sl) - cell_sl)) < 1e-10)
chk("slopes: and the per-cell slopes are the fitted coefficients themselves",
    max(abs(cell_sl - coef(mi)[grepl(":training$", names(coef(mi)))])) < 1e-10)
chk("slopes: the reference-cell slope is named for the covariate alone",
    "training" %in% si$term &&
    abs(si$estimate[si$term == "training"] -
        coef(mi)[[paste0("cell", levels(spi$data$cell)[1], ":training")]]) < 1e-10)
chk("slopes: their meanings mark them as slopes",
    identical(si$meaning[si$term == "chord_typemaj:training"],
              paste0("-aug.none + maj.", reference_levels(spi)[["inversion"]],
                     " (slope on training)")))
chk("slopes: cell space names them per condition",
    identical(nest_summary(mi, "cells")$term[11],
              "aug.none slope on training"))

## ---- the fit must still correspond to the declaration ---------------------
## A model that does not contain the declared structure cannot respond to it, so
## predictions return whatever weighting the fit already implies rather than the
## policy asked for - and on balanced data the two can coincide, which makes the
## error invisible exactly where it is easiest to make. update() is the way a
## fit gets simplified past that point, so that is where it is checked.
m_thin <- suppressWarnings(update(mf, . ~ . - cell))
chk("correspondence: a refit that drops the declared structure is refused",
    grepl("does not contain `chord_type`",
          err_of(nest_estimand(m_thin, chord_type, bounds = FALSE))))
chk("correspondence: and the refusal says what the answer would have been",
    grepl("weighting the fit already implies",
          err_of(nest_estimand(m_thin, chord_type, bounds = FALSE))))

## ---- update() keeps the declaration ---------------------------------------
chk("nest_fit: the fit gains a class, behind the engine's own",
    identical(class(mf)[1], "lm") && inherits(mf, "nestimand_fit"))
mu <- update(mf, . ~ . - training)
chk("update: the declaration survives the refit",
    inherits(attr(mu, "nestimand_spec"), "nesting_spec") &&
    inherits(mu, "nestimand_fit") && length(coef(mu)) == 10)
chk("update: the estimand works from the updated fit, without `spec`",
    abs(as.data.frame(nest_estimand(mu, chord_type, bounds = FALSE))$estimate[3] - 0.6779) < 1e-4)
chk("update: the code view records the new call, not the old one",
    any(grepl("updated from the original call", attr(mu, "nestimand_code"))) &&
    !any(grepl("training", attr(mu, "nestimand_code"))))
chk("update: an update that drops the structure is caught downstream",
    grepl("does not contain",
          err_of(nest_estimand(update(mf, . ~ . - cell), chord_type))))
chk("class: the engine's own methods still work",
    length(predict(mf, newdata = head(sp$data))) == 6 &&
    inherits(summary(mf), "summary.lm") && nrow(vcov(mf)) == 11)

## ---- marginal contrasts of a nested variable ------------------------------
## `inversion` holds only the sentinel in the augmented stratum, so a contrast
## against that level would compare strata rather than inversions. Such strata
## are excluded, giving the pooled estimand of Section 4.1.
einv <- nest_estimand(mf, inversion, bounds = FALSE)
dinv <- as.data.frame(einv)
chk("nested target: no contrast involves the sentinel level",
    nrow(dinv) == 3 && !any(grepl("none", dinv$term)))
chk("nested target: the pooled contrasts match Section 4.1",
    max(abs(sort(round(dinv$estimate, 4)) -
            sort(c(-0.1632, -0.7337, -0.5706)))) < 1e-3)
chk("nested target: the restriction is stated in the emitted code",
    any(grepl("does not vary in every stratum", attr(einv, "nestimand")$code)) &&
    any(grepl("subset(spec$cells", attr(einv, "nestimand")$code, fixed = TRUE)))
chk("nested target: contrasts inside each parent stratum are unaffected",
    nrow(as.data.frame(suppressMessages(nest_estimand(mf, inversion, by = "chord_type",
                                bounds = FALSE)))) == 9)
chk("root target: no restriction applies, and the estimand is unchanged",
    is.null(degenerate_strata(sp, "chord_type")) &&
    abs(as.data.frame(nest_estimand(mf, chord_type, bounds = FALSE))$estimate[3] - 0.6779) < 1e-4)
chk("versions: a nested variable's versions are the strata it occurs in",
    identical(sort(versions_of(sp, "inversion")[["0"]]), c("dim", "maj", "min")) &&
    identical(versions_of(sp, "inversion")[["none"]], "aug"))
chk("versions: a root variable's versions are unchanged",
    identical(sort(versions_of(sp, "chord_type")[["maj"]]), c("0", "1", "2")))

## ---- bounds: one row per contrast, correctly paired ------------------------
bi <- attr(nest_estimand(mf, inversion), "nestimand")$bounds
chk("bounds: one row per contrast, not recycled against a longer vector",
    nrow(bi) == 3 && !anyDuplicated(bi$term))
chk("bounds: each estimate lies inside its own interval",
    all(bi$estimate >= bi$policy_low - 1e-8) &&
    all(bi$estimate <= bi$policy_high + 1e-8))
chk("bounds: the vertex policies respect the stratum restriction",
    nrow(policy_vertices(sp, "inversion")) == 27)
## the same estimand under permuted levels: the canonical direction follows the
## declared order, so a contrast may be reported the other way round, with its
## interval negated and swapped - the substance is unchanged
dperm <- dat
set.seed(7); dperm$inversion <- factor(dperm$inversion,
                                       levels = c("none", sample(c("0", "1", "2"))))
spp <- nesting_spec(dperm, response ~ chord_type * inversion + training,
                    "inversion %in% chord_type")
bp <- attr(nest_estimand(fit_from(spp), inversion), "nestimand")$bounds
chk("bounds: invariant under level permutation, up to contrast direction",
    isTRUE(all.equal(sort(abs(bi$estimate)), sort(abs(bp$estimate)), tolerance = 1e-8)) &&
    isTRUE(all.equal(sort(abs(c(bi$policy_low, bi$policy_high))),
                     sort(abs(c(bp$policy_low, bp$policy_high))), tolerance = 1e-8)))
chk("labels: the engine's own label column is not left contradicting the estimate",
    { d <- as.data.frame(nest_estimand(mf, chord_type, bounds = FALSE))
      lc <- mfx_term_column(d)
      identical(as.character(d[[lc]]), d$term) })
chk("target: a variable holding the name is read for its value",
    { tv <- "inversion"
      nrow(as.data.frame(nest_estimand(mf, tv, bounds = FALSE))) == 3 })
chk("model: an inline call is named safely in the emitted code",
    any(grepl("^model <- lm",
        attr(nest_estimand(fit_from(sp), chord_type, bounds = FALSE), "nestimand")$code)))

## ---- the declarations must identify the model whatever the level order -----
## The sentinel need not be the reference level. Where it is not, the empty and
## the redundant columns fall differently, and the count of each changes; what
## must not change is that together they leave exactly a full-rank model.
sent_last <- dat
sent_last$inversion <- factor(sent_last$inversion, levels = c("0", "1", "2", "none"))
sp_last <- nesting_spec(sent_last, response ~ chord_type * inversion + training,
                        "inversion %in% chord_type")
m_last <- fit_from(sp_last)
chk("sentinel last: the estimand is unchanged",
    abs(as.data.frame(nest_estimand(m_last, chord_type, bounds = FALSE))$estimate[3] - 0.6779) < 1e-4)
chk("sentinel last: the nested contrast still excludes the sentinel",
    { d <- as.data.frame(nest_estimand(m_last, inversion, bounds = FALSE))
      nrow(d) == 3 && !any(grepl("none", d$term)) })
chk("sentinel last: the effect basis is still square and full rank",
    { A <- effect_basis(sp_last); nrow(A) == 10 && qr(A)$rank == 10 })

## ---- the level order of the sentinel is the package's business -------------
## The chain form is legible only with the sentinel as reference. Rather than
## asking for that, the package finds the sentinel from the design and reorders
## the data the chain form is built from, visibly. Reporting keeps the declared
## order throughout.
chk("sentinel: found from the design, not from its name",
    identical(unlist(sentinel_levels(sp)), c(inversion = "none")) &&
    identical(unlist(sentinel_levels(sp_last)), c(inversion = "none")))
spb_f <- nesting_spec(dat, response ~ chord_type * inversion + training,
                      "inversion %in% chord_type", fit = "brm")
spb_l <- nesting_spec(sent_last, response ~ chord_type * inversion + training,
                      "inversion %in% chord_type", fit = "brm")
chk("sentinel: the cells parameterization needs no reordering",
    !any(grepl("relevel", attr(m_last, "nestimand_code"))))
chk("sentinel: contrasts are still reported in the declared order",
    identical(as.data.frame(nest_estimand(m_last, inversion, bounds = FALSE))$term,
              c("1 - 0", "2 - 0", "2 - 1")))
chk("sentinel: the effect basis is the same either way",
    isTRUE(all.equal(unname(effect_basis(sp)), unname(effect_basis(sp_last)))))

## ---- routes: where the model is evaluated ----------------------------------
## The policy says how versions are weighted; the route says over which rows the
## model is evaluated. On balanced linear data the three agree; they separate
## under covariate imbalance, which is what G-computation exists to handle.
rt3 <- function(mm, spx, rt, pol = "equal")
  as.data.frame(nest_estimand(mm, chord_type, policy = pol, route = rt,
                         bounds = FALSE))$estimate[3]
chk("routes: all three agree on balanced data with a shared covariate",
    max(abs(diff(vapply(c("g_computation", "cells"),
                        function(r) rt3(mf, sp, r), 1)))) < 1e-8)
set.seed(4)
dun <- dat[c(rep(which(dat$inversion == "0"), 2), which(dat$inversion != "0")), ]
dun$training <- rnorm(nrow(dun), ifelse(dun$chord_type == "aug", 7, 3), 1.5)
dun$response <- dun$response + 0.15 * dun$training
spu <- nesting_spec(dun, response ~ chord_type * inversion + training,
                    "inversion %in% chord_type")
mu2 <- fit_from(spu)
chk("routes: cells and g_computation agree in a linear model",
    abs(rt3(mu2, spu, "cells") - rt3(mu2, spu, "g_computation")) < 1e-8)
chk("routes: equal and proportional separate once the design is unbalanced",
    abs(rt3(mu2, spu, "g_computation", "equal") -
        rt3(mu2, spu, "g_computation", "proportional")) > 1e-3)
chk("routes: the route is recorded and printed",
    identical(attr(nest_estimand(mf, chord_type, route = "cells", bounds = FALSE), "nestimand")$route, "cells"))
chk("routes: the emitted code says which rows the model was evaluated over",
    any(grepl("crossed with every realized",
        attr(nest_estimand(mf, chord_type, bounds = FALSE,
                      p_adjust = "none"), "nestimand")$code)) &&
    any(grepl("covariates at their means",
        attr(nest_estimand(mf, chord_type, route = "cells", bounds = FALSE, p_adjust = "none"),
             "nestimand")$code)))

## ---- unit weights: standardizing to another population ---------------------
## The policy weights the versions of a condition; unit weights weight the rows,
## which is how an estimand is standardized to a population other than the
## sample - survey weights, or post-stratification on the covariates.
dw <- dat; dw$wt_hi <- as.numeric(dw$training > 7); dw$wt1 <- 1
spw <- nesting_spec(dw, response ~ chord_type * inversion * training,
                    "inversion %in% chord_type")
mw <- fit_from(spw)
gw <- function(...) as.data.frame(nest_estimand(mw, chord_type, bounds = FALSE, ...))$estimate[3]
chk("weights: uniform weights leave the estimand unchanged",
    abs(gw(weights = "wt1") - gw()) < 1e-10)
chk("weights: reweighting to a subpopulation matches computing on those rows",
    abs(gw(weights = "wt_hi") -
        gw(data = subset(spw$data, training > 7))) < 1e-8)
chk("weights: a numeric vector is accepted as well as a column name",
    abs(gw(weights = spw$data$wt_hi) - gw(weights = "wt_hi")) < 1e-10)
chk("weights: a vector of the wrong length is refused",
    grepl("one weight per row", err_of(gw(weights = c(1, 2)))))
chk("weights: an absent column is refused, and the stale-spec case explained",
    grepl("added after the model was fitted", err_of(gw(weights = "nope"))))
chk("weights: refused on the cells route, which has no units to weight",
    grepl("no units to weight", err_of(gw(route = "cells", weights = "wt_hi"))))
chk("weights: they appear in the emitted code, with their purpose",
    any(grepl("standardized to the population",
        attr(nest_estimand(mw, chord_type, weights = "wt_hi", bounds = FALSE), "nestimand")$code)))

## ---- the reference condition follows the declared level order --------------
## Beyond the structurally empty columns, one column per stratum is redundant.
## Which is dropped is a coding choice that fixes the reference condition, so it
## follows the order declared in the data rather than whichever column a pivot
## reaches last.
ref_of <- function(levs) {
  dd <- dat; dd$inversion <- factor(dd$inversion, levels = levs)
  spr <- nesting_spec(dd, response ~ chord_type * inversion + training,
                      "inversion %in% chord_type")
  ns <- nest_summary(fit_from(spr))
  sub("^-aug.none \\+ dim\\.", "", ns$meaning[ns$term == "chord_typedim"])
}
chk("reference: the first declared non-sentinel level is the reference",
    identical(ref_of(c("0", "1", "2", "none")), "0") &&
    identical(ref_of(c("2", "1", "0", "none")), "2") &&
    identical(ref_of(c("none", "1", "0", "2")), "1"))
chk("reference: reference_levels() reports it",
    identical(reference_levels(sp)[["inversion"]], "0"))
chk("reference: the basis stays square and full rank whichever is chosen",
    { dd <- dat; dd$inversion <- factor(dd$inversion, levels = c("2","1","0","none"))
      A <- effect_basis(nesting_spec(dd, response ~ chord_type * inversion,
                                     "inversion %in% chord_type"))
      nrow(A) == 10 && ncol(A) == 10 && qr(A)$rank == 10 })
chk("reference: the estimand is unaffected by the choice",
    { dd <- dat; dd$inversion <- factor(dd$inversion, levels = c("2","1","0","none"))
      spr <- nesting_spec(dd, response ~ chord_type * inversion + training,
                          "inversion %in% chord_type")
      abs(as.data.frame(nest_estimand(fit_from(spr), chord_type, bounds = FALSE))$estimate[3] - 0.6779) < 1e-4 })

## ---- the reorder check on a mixed model ------------------------------------
## Order instability is a fixed-effects phenomenon, so the check runs on the
## fixed-effects shadow model. Refitting a large random structure is expensive
## and can settle on a different optimum, which would report a failure that has
## nothing to do with level order.

chk("random: a fit with no random effects says so rather than failing",
    grepl("for mixed fits", err_of(random_covariance(mf))))

## ---- every combination of target, route and policy runs --------------------
## Both bugs found here were in combinations never exercised together: the cells
## route with a restricted target, and the observed route with rows in an
## excluded stratum.
combo <- expand.grid(target = c("chord_type", "inversion"),
                     route = c("g_computation", "cells"),
                     policy = c("equal", "proportional"),
                     stringsAsFactors = FALSE)
res <- vapply(seq_len(nrow(combo)), function(i) {
  r <- tryCatch(as.data.frame(nest_estimand(mf, combo$target[i], policy = combo$policy[i],
                  route = combo$route[i], bounds = TRUE))$estimate[1],
                error = function(e) NA_real_)
  r }, 1)
chk("combinations: all twelve target/route/policy pairings run",
    !anyNA(res))
chk("combinations: each target agrees across routes on balanced data",
    max(abs(diff(res[combo$target == "chord_type"]))) < 1e-8 &&
    max(abs(diff(res[combo$target == "inversion"]))) < 1e-8)
chk("combinations: the emitted code runs for each of them",
    all(vapply(seq_len(nrow(combo)), function(i) {
      e <- nest_estimand(mf, combo$target[i], policy = combo$policy[i],
                    route = combo$route[i], bounds = FALSE)
      env <- new.env(); assign("mf", mf, env); assign("sp", sp, env)
      r <- tryCatch(eval(parse(text = paste(attr(e, "nestimand")$code,
                                            collapse = "\n")), envir = env),
                    error = function(x) NULL)
      !is.null(r) }, TRUE)))
chk("policy_weights: a stratum the policy does not cover contributes nothing",
    { pol <- nest_policy(sp, "inversion", "equal",
                         cells = sp$cells[sp$cells$chord_type != "aug", ])
      w <- policy_weights(sp, sp$data, pol)
      all(w[sp$data$chord_type == "aug"] == 0) && any(w > 0) })

## A mixed fit on replicated data, used by the sections that follow.
if (requireNamespace("lme4", quietly = TRUE)) {
  set.seed(9)
  d6 <- do.call(rbind, lapply(1:3, function(i) {
    x <- dat; x$response <- x$response + rnorm(nrow(x), 0, 0.4); x }))
  spm2 <- nesting_spec(d6, response ~ chord_type * inversion + training +
                       (chord_type * inversion | participant),
                       "inversion %in% chord_type", fit = "lmer")
  mm2 <- suppressWarnings(fit_from(spm2))
}

## ---- mixed fits: population-level by default, and grids that carry groups ---
if (requireNamespace("lme4", quietly = TRUE)) {
  rr <- vapply(c("g_computation", "cells"), function(rt)
    tryCatch(suppressWarnings(as.data.frame(
      nest_estimand(mm2, chord_type, route = rt, bounds = FALSE))$estimate[1]), error = function(e) NA_real_), 1)
  chk("mixed: every route runs, the cells grid carrying the grouping column",
      !anyNA(rr) && max(abs(diff(rr))) < 1e-8)
  chk("mixed: cell_grid supplies the grouping factor",
      "participant" %in% names(cell_grid(spm2)))
  chk("mixed: grouping_vars reads them from the declared bars",
      identical(grouping_vars(spm2), "participant"))
  cd <- attr(nest_estimand(mm2, chord_type, route = "cells", bounds = FALSE, p_adjust = "none"), "nestimand")$code
  chk("mixed: the engine's default stands, and is stated",
      any(grepl("default stands", cd)) &&
      !any(grepl("re.form", grep("avg_predictions|hypothesis =", cd, value = TRUE),
                 fixed = TRUE)))
  chk("mixed: a user-supplied re.form is not overridden",
      any(grepl("re.form = NULL",
          attr(nest_estimand(mm2, chord_type, bounds = FALSE,
                        re.form = NULL), "nestimand")$code, fixed = TRUE)))
  chk("mixed: brms would be told in its own spelling",
      { spb_m <- nesting_spec(dat, response ~ chord_type * inversion +
                              (1 | participant), "inversion %in% chord_type",
                              fit = "brm")
        identical(if (identical(spb_m$fit, "brm")) "re_formula" else "re.form",
                  "re_formula") && length(grouping_vars(spb_m)) == 1 })
}
chk("mixed: a fit with no random terms is left alone",
    !any(grepl("re.form", attr(nest_estimand(mf, chord_type, bounds = FALSE, p_adjust = "none"),
                               "nestimand")$code)))

## ---- everything the emitted code calls must be reachable -------------------
## The code runs in the user's environment, so a function that is not exported
## is not there. This walks every combination and checks the whole call set.
exported <- sub("^export\\((.*)\\)$", "\\1",
                grep("^export\\(", readLines("NAMESPACE"), value = TRUE))
combo2 <- expand.grid(target = c("chord_type", "inversion"),
                      route = c("g_computation", "cells"),
                      scale = c("response", "latent"),
                      stringsAsFactors = FALSE)
called <- unlist(lapply(seq_len(nrow(combo2)), function(i) {
  e <- tryCatch(nest_estimand(mf, combo2$target[i], route = combo2$route[i],
                         scale = combo2$scale[i], bounds = TRUE), error = function(x) NULL)
  if (is.null(e)) return(character(0))
  cd <- c(attr(e, "nestimand")$code, attr(mf, "nestimand_code"))
  unlist(regmatches(cd, gregexpr("[A-Za-z._][A-Za-z0-9._]*(?=\\()", cd, perl = TRUE)))
}))
base_or_engine <- c("lm", "glm", "factor", "levels", "c", "subset", "library",
                    "as.character", "paste", "avg_predictions", "function", "if")
missing_exports <- setdiff(setdiff(unique(called), base_or_engine), exported)
missing_exports <- missing_exports[vapply(missing_exports, exists, TRUE)]
chk("emitted code: every nestimand function it calls is exported",
    length(missing_exports) == 0)
chk("routes: a route that averages over each condition's own rows is refused",
    grepl("should be one of", err_of(nest_estimand(mf, chord_type, route = "observed"))))

## ---- printing and multiple targets -----------------------------------------
chk("print: four decimal places by default, more on request",
    { p4 <- capture.output(print(nest_estimand(mf, chord_type, bounds = FALSE)))
      p8 <- capture.output(print(nest_estimand(mf, chord_type, bounds = FALSE), digits = 8))
      any(grepl("0.6779", p4, fixed = TRUE)) &&
      any(grepl("0.67788574", p8, fixed = TRUE)) })
chk("print: the object keeps full precision whatever is shown",
    abs(as.data.frame(nest_estimand(mf, chord_type, bounds = FALSE))$estimate[3] -
        0.677885742) < 1e-8)
chk("print: nest_summary takes digits too",
    any(grepl("4.08435", capture.output(print(nest_summary(mf), digits = 6)),
              fixed = TRUE)))
em <- nest_estimand(mf, c(chord_type, inversion), bounds = FALSE)
chk("targets: several at once, bare or quoted",
    inherits(em, "nestimand_estimands") &&
    identical(names(em), c("chord_type", "inversion")) &&
    identical(names(nest_estimand(mf, c("chord_type", "inversion"), bounds = FALSE)), names(em)))
chk("targets: each element is an ordinary estimand table",
    inherits(em[["chord_type"]], "nestimand_estimand") &&
    nrow(as.data.frame(em[["inversion"]])) == 3)
chk("targets: the nested one is still restricted to the strata it varies in",
    !any(grepl("none", as.data.frame(em[["inversion"]])$term)))
chk("targets: show_code covers every one of them",
    { cd <- capture.output(show_code(em))
      any(grepl("--- chord_type ---", cd)) && any(grepl("--- inversion ---", cd)) })

## ---- interaction contrasts, and the formula forms --------------------------
ei <- nest_estimand(mf, chord_type:inversion)
di <- as.data.frame(ei)
cm <- tapply(dat$response, paste(dat$chord_type, dat$inversion, sep = "."), mean)
chk("interaction: a difference of differences, matching the raw cell means",
    abs(di$estimate[di$term == "(min - dim) x (1 - 0)"] -
        ((cm[["min.1"]] - cm[["min.0"]]) - (cm[["dim.1"]] - cm[["dim.0"]]))) < 1e-4)
chk("interaction: 9 contrasts, none touching the sentinel stratum",
    nrow(di) == 9 && !any(grepl("aug|none", di$term)))
chk("interaction: a policy is dropped rather than refused, and reported",
    { msg <- NULL
      r <- withCallingHandlers(
        nest_estimand(mf, chord_type:inversion, policy = "proportional"),
        message = function(m) { msg <<- conditionMessage(m)
                                invokeRestart("muffleMessage") })
      grepl("does not apply to an interaction", msg) &&
      isTRUE(attr(r, "nestimand")$policy_dropped) &&
      nrow(as.data.frame(r)) == 9 })
chk("interaction: dropping it leaves the marginal contrasts' policy intact",
    { r <- suppressMessages(nest_estimand(mf, chord_type * inversion,
                                     policy = "proportional", bounds = FALSE))
      identical(attr(r[["chord_type"]], "nestimand")$policy, "proportional") &&
      isTRUE(attr(r[["chord_type:inversion"]], "nestimand")$policy_dropped) })
chk("interaction: the printed policy reads as not applicable",
    any(grepl("Policy: not applicable",
        capture.output(print(nest_estimand(mf, chord_type:inversion))))))
chk("interaction: no bounds, since there is no policy to vary",
    is.null(attr(ei, "nestimand")$bounds))
est_star <- nest_estimand(mf, chord_type * inversion, bounds = FALSE)
chk("targets: `a * b` gives a, b, and their interaction",
    identical(names(est_star),
              c("chord_type", "inversion", "chord_type:inversion")) &&
    identical(vapply(est_star, function(z) nrow(as.data.frame(z)), 1L),
              c(chord_type = 6L, inversion = 3L, `chord_type:inversion` = 9L)))
chk("targets: `a:b` gives the interaction alone",
    inherits(ei, "nestimand_estimand") && nrow(di) == 9)
chk("interaction: the emitted code re-runs",
    { env <- new.env(); assign("mf", mf, env); assign("sp", sp, env)
      r <- eval(parse(text = paste(attr(ei, "nestimand")$code, collapse = "\n")),
                envir = env)
      isTRUE(all.equal(as.data.frame(r)$estimate, di$estimate)) })
chk("interaction: a `:` target is how it is asked for",
    { a <- as.data.frame(nest_estimand(mf, chord_type:inversion, bounds = FALSE))
      nrow(a) == nrow(di) && isTRUE(all.equal(a$estimate, di$estimate)) })
chk("interaction: a `*` target gives the marginals and it together",
    { b <- suppressMessages(nest_estimand(mf, chord_type * inversion, bounds = FALSE))
      length(b) == 3 &&
        identical(names(b), c("chord_type", "inversion", "chord_type:inversion")) &&
        isTRUE(all.equal(as.data.frame(b[["chord_type:inversion"]])$estimate,
                         di$estimate)) })

## ---- the bounds follow the route they were asked for -----------------------
## They are the same estimand under other policies, so computing them on a
## different grid from the estimand itself is both inconsistent and, on the
## g-computation grid, far slower than the route requested.
bc <- attr(nest_estimand(mf, inversion, route = "cells"),
           "nestimand")$bounds
bg <- attr(nest_estimand(mf, inversion), "nestimand")$bounds
chk("bounds: the same whichever route computes them, on balanced data",
    isTRUE(all.equal(bc$policy_low, bg$policy_low, tolerance = 1e-8)) &&
    isTRUE(all.equal(bc$policy_high, bg$policy_high, tolerance = 1e-8)))
chk("bounds: the route is carried into the emitted call",
    any(grepl('route = "cells"',
        attr(nest_estimand(mf, inversion, route = "cells",
                      p_adjust = "none"), "nestimand")$code, fixed = TRUE)) &&
    !any(grepl("route =",
        attr(nest_estimand(mf, inversion, p_adjust = "none"),
             "nestimand")$code, fixed = TRUE)))
chk("bounds: the emitted code re-runs on the cells route",
    { e <- nest_estimand(mf, inversion, route = "cells")
      env <- new.env(); assign("mf", mf, env); assign("sp", sp, env)
      r <- eval(parse(text = paste(attr(e, "nestimand")$code, collapse = "\n")),
                envir = env)
      isTRUE(all.equal(as.data.frame(r)$estimate, as.data.frame(e)$estimate)) })

## ---- the estimand does not refit the model ---------------------------------
## The fit belongs in the code view, so that what is shown is a whole analysis,
## but evaluating it there would re-run a fit that has already been done - which
## on a mixed model costs more than everything else together.
counted <- 0
trace_fit <- function() counted <<- counted + 1
sp_c <- nesting_spec(dat, response ~ chord_type * inversion + training,
                     "inversion %in% chord_type")
m_c <- fit_from(sp_c)
local({
  ## a model whose refit would be visible: replace the stored call with one that
  ## increments a counter, then check the counter stays at zero
  attr(m_c, "nestimand_code") <<- c("## nestimand -- fit", "m <- trace_fit()")
})
e_c <- nest_estimand(m_c, chord_type, bounds = FALSE)
chk("estimand: the fit is shown in the code but not re-run",
    counted == 0 && nrow(as.data.frame(e_c)) == 6)
chk("estimand: the code view still contains the fit",
    any(grepl("nestimand -- fit", attr(e_c, "nestimand")$code)))

## ---- brms's own arguments reach the engine ---------------------------------
## `priors` sits after the dots so that R cannot partial-match brms's `prior`
## to it; a brms prior therefore passes straight through.
if (requireNamespace("brms", quietly = TRUE)) {
  cd_b <- suppressMessages(fit_from(spb2, dry_run = TRUE,
                   prior = brms::set_prior("normal(0, 1)", class = "b"),
                   sample_prior = "only", init = 0, file = "somewhere"))
  chk("brms: `prior` is handed to the engine, not captured by `priors`",
      grepl("prior = ", cd_b, fixed = TRUE) &&
      grepl("stanvars = prior_stanvars", cd_b, fixed = TRUE))
  chk("brms: the other engine arguments travel with it",
      grepl('sample_prior = "only"', cd_b, fixed = TRUE) &&
      grepl("init = 0", cd_b, fixed = TRUE) &&
      grepl('file = "somewhere"', cd_b, fixed = TRUE))
}

## ---- brms priors: what belongs where ---------------------------------------
if (requireNamespace("brms", quietly = TRUE)) {
  spb_re <- nesting_spec(dat, response ~ chord_type * inversion + (1 | participant),
                         "inversion %in% chord_type", fit = "brm")
  pb <- nest_prior(spb_re, mean = 4, sd = 1.5, on = "cells")
  chk("brms priors: sd and cor priors pass through in silence",
      { msg <- NULL
        withCallingHandlers(fit_from(spb_re, dry_run = TRUE,
          prior = brms::set_prior("exponential(1)", class = "sd")),
          message = function(m) { msg <<- conditionMessage(m)
                                  invokeRestart("muffleMessage") })
        is.null(msg) })
  chk("brms priors: a class b prior is translated, and the translation reported",
      { msg <- NULL
        cd <- withCallingHandlers(fit_from(spb_re, dry_run = TRUE,
          prior = brms::set_prior("normal(0, 1)", class = "b")),
          message = function(m) { msg <<- conditionMessage(m)
                                  invokeRestart("muffleMessage") })
        grepl("has been translated", msg) &&
        grepl("prior_statement", cd, fixed = TRUE) &&
        grepl("stanvars = prior_stanvars", cd, fixed = TRUE) })
  chk("brms priors: the other classes travel alongside the translated one",
      { cd <- suppressMessages(fit_from(spb_re, dry_run = TRUE,
                prior = c(brms::set_prior("normal(0, 1)", class = "b"),
                          brms::set_prior("exponential(1)", class = "sd"))))
        grepl("prior_statement", cd, fixed = TRUE) &&
        grepl("exponential(1)", cd, fixed = TRUE) })
  chk("brms priors: what the stated prior implies for a contrast is recoverable",
      { tr <- prior_from_brms(spb_re, brms::set_prior("normal(0, 1)", class = "b"))
        abs(prior_for_estimand(tr$translated, "chord_type", "equal")$sd[3] -
            1.1055) < 1e-3 })
  chk("brms priors: a coefficient-specific b prior is refused, with the remedy",
      grepl("whole mean structure",
            err_of(prior_from_brms(spb_re,
              brms::set_prior("normal(0,1)", class = "b", coef = "cellmaj.0")))))
  chk("brms priors: a prior that is not elliptical is refused, with the remedy",
      grepl("translate exactly",
            err_of(prior_from_brms(spb_re,
              brms::set_prior("cauchy(0, 1)", class = "b")))))
}

## ---- sample_prior = "yes" and a translated prior are incompatible ----------
## brms writes one prior draw per class into a scalar, which a multivariate
## prior on the coefficients cannot satisfy: Stan reports an ill-typed
## assignment. Caught here rather than after a compilation.
if (requireNamespace("brms", quietly = TRUE)) {
  chk("sample_prior: \"yes\" with a translated prior is refused before compiling",
      grepl("ill-typed assignment",
            err_of(suppressMessages(fit_from(spb2, dry_run = TRUE,
              prior = brms::set_prior("normal(0, 1)", class = "b"),
              sample_prior = "yes")))))
  chk("sample_prior: the refusal names the way round it",
      { e <- err_of(suppressMessages(fit_from(spb2, dry_run = TRUE,
              prior = brms::set_prior("normal(0, 1)", class = "b"),
              sample_prior = "yes")))
        grepl('sample_prior = "only"', e) })
  chk("sample_prior: \"only\" is accepted, as the prior check uses it",
      inherits(suppressMessages(fit_from(spb2, dry_run = TRUE,
        prior = brms::set_prior("normal(0, 1)", class = "b"),
        sample_prior = "only")), "nestimand_code"))
  chk("sample_prior: \"yes\" is fine when no prior touches the mean structure",
      inherits(suppressMessages(fit_from(spb2, dry_run = TRUE,
        prior = brms::set_prior("exponential(1)", class = "sigma"),
        sample_prior = "yes")), "nestimand_code"))
}

## ---- the multi-target object prints, and interactions pass their check -----
chk("targets: the collection prints even without names",
    { e <- suppressMessages(nest_estimand(mf, chord_type * inversion,
             policy = "proportional", bounds = FALSE))
      length(capture.output(print(unname(e)))) > 0 &&
      length(capture.output(print(e))) > 0 })
chk("targets: show_code covers the collection without names too",
    { e <- suppressMessages(nest_estimand(mf, chord_type * inversion, bounds = FALSE))
      length(capture.output(show_code(unname(e)))) > 0 })
chk("targets: an empty collection says so rather than printing nothing",
    identical(trimws(capture.output(
      print(structure(list(), class = "nestimand_estimands")))), "no estimands"))

## ---- brms objects need reducing before they can be read --------------------
## A brmsformula yields nothing to all.vars(), so a correspondence check written
## for ordinary formulas rejects every brms fit; and brms returns its covariance
## as separate sd and correlation pieces, whose names the matrix product drops.
if (requireNamespace("brms", quietly = TRUE)) {
  bfm <- brms::bf(response ~ cell + cell:training + (0 + cell | participant))
  fm <- bfm
  if (inherits(fm, "bform")) fm <- stats::formula(fm)
  chk("brms: a brmsformula is reduced before its variables are read",
      all(c("cell", "training", "participant") %in% all.vars(fm)) &&
      length(all.vars(bfm)) == 0)
}
if (requireNamespace("lme4", quietly = TRUE)) {
  rcn <- random_covariance(mm2)
  chk("random: the covariance keeps its names through the matrix product",
      !is.null(rownames(rcn[[1]])) &&
      identical(rownames(rcn[[1]]), colnames(rcn[[1]])))
}

## ---- the latent route covers what the prediction route covers --------------
## the linear predictor on an lm is reached by the package's own default,
## since marginaleffects does not accept "link" for that class
lat <- function(...) as.data.frame(nest_estimand(mf, ..., bounds = FALSE))
res <- function(...) as.data.frame(nest_estimand(mf, ..., bounds = FALSE))
chk("latent: an interaction is available, and matches the prediction route",
    { a <- res(chord_type:inversion); b <- lat(chord_type:inversion)
      nrow(b) == 9 && max(abs(a$estimate - b$estimate[match(a$term, b$term)])) < 1e-8 })
chk("latent: a nested target is restricted here too",
    { a <- res(inversion); b <- lat(inversion)
      identical(a$term, b$term) && nrow(b) == 3 &&
      max(abs(a$estimate - b$estimate)) < 1e-8 })
chk("latent: the emitted code names the interaction route",
    any(grepl('contrast = "interaction"',
        attr(nest_estimand(mf, chord_type:inversion, bounds = FALSE), "nestimand")$code, fixed = TRUE)))
chk("latent: an ordinal star form gives all three, correctly sized",
    { eo2 <- suppressMessages(nest_estimand(mo, chord_type * inversion,
               policy = "proportional", bounds = FALSE))
      identical(unname(vapply(eo2, function(z) nrow(as.data.frame(z)), 1L)),
                c(6L, 3L, 9L)) })

## ---- the contrast matrix is shared by both evaluations ---------------------
## Whether a contrast is evaluated at the coefficients or draw by draw, it is
## the same C: only the summary differs.
Cc <- latent_contrast_matrix(mf, sp, "chord_type", "equal")
chk("latent: C has one row per contrast and one column per coefficient",
    nrow(Cc) == 6 && identical(rownames(Cc)[3], "maj - aug") &&
    all(colnames(Cc) %in% names(coef(mf))))
chk("latent: C b reproduces the reported estimates",
    { b <- coef(mf)[colnames(Cc)]
      max(abs(as.numeric(Cc %*% b) -
              latent_estimand(mf, "chord_type", "equal")$estimate)) < 1e-10 })
Ci2 <- latent_contrast_matrix(mf, sp, c("chord_type", "inversion"),
                              contrast = "interaction")
chk("latent: the interaction C has one row per difference of differences",
    nrow(Ci2) == 9 && all(abs(rowSums(Ci2)) < 1e-10))
chk("latent: a frequentist fit still reports a test statistic",
    all(c("statistic", "p.value") %in%
        names(latent_estimand(mf, "chord_type", "equal"))))

## ---- printed alignment ------------------------------------------------------
chk("print: every numeric column is right-aligned, character or not",
    { d <- data.frame(term = "a - b", estimate = 0.5, pd = 0.99)
      x <- structure(d, class = c("nestimand_estimand", "data.frame"),
                     nestimand = list(policy = "equal", route = "g_computation",
                                      contrast = "pairwise", code = "x"))
      ln <- capture.output(print(x))
      ## the value ends where its header ends
      endsw <- function(h, r, nm) {
        i <- regexpr(nm, h, fixed = TRUE)[1] + nchar(nm) - 1
        substr(r, i, i) != " " }
      endsw(ln[1], ln[2], "estimate") && endsw(ln[1], ln[2], "pd") })
chk("print: numeric headers sit over their columns, right-aligned",
    { ln <- capture.output(print(nest_estimand(mf, chord_type, bounds = FALSE)))
      hdr <- ln[1]; row <- ln[2]
      ## the header and the values end at the same column for a numeric field
      regexpr("estimate", hdr, fixed = TRUE)[1] +
        nchar("estimate") - 1 == regexpr("\\s0\\.|\\s-0\\.", row)[1] +
        attr(regexpr("\\s-?0\\.[0-9]+", row), "match.length") - 1 ||
      grepl(" estimate", hdr, fixed = TRUE) })
chk("print: the label column stays left-aligned, header included",
    { ln <- capture.output(print(nest_summary(mf, "cells")))
      any(grepl("^ term ", ln)) })
chk("print: alignment does not disturb the values",
    { d <- as.data.frame(nest_estimand(mf, chord_type, bounds = FALSE))
      abs(d$estimate[3] - 0.677885742) < 1e-8 })

chk("ordinal: the latent scale needs no engine support at all",
    nrow(as.data.frame(nest_estimand(mo, chord_type,
                                bounds = FALSE))) == 6)
chk("ordinal: the probe reports which names it tried",
    { e <- err_of(ordinal_response_type(mo, spo, spo$data))
      identical(e, "") || (grepl("prob", e) && grepl("response", e)) })

## ---- the default scale depends on the family --------------------------------
chk("scale: an ordinal fit defaults to the latent scale",
    identical(attr(nest_estimand(mo, chord_type, bounds = FALSE),
                   "nestimand")$scale, "latent"))
chk("scale: a linear fit still asks for the response scale by default",
    identical(attr(nest_estimand(mf, chord_type, bounds = FALSE),
                   "nestimand")$type, "response"))
chk("scale: the ordinal default gives one number per contrast",
    nrow(as.data.frame(nest_estimand(mo, chord_type, bounds = FALSE))) == 6)
chk("type: response resolves to the name this engine accepts",
    { cd <- suppressMessages(nest_estimand(mo, chord_type, type = "response",
              bounds = FALSE, dry_run = TRUE))
      grepl('type = "prob"', cd, fixed = TRUE) ||
      grepl('type = "response"', cd, fixed = TRUE) })

## ---- the default scale follows the family -----------------------------------
chk("scale: an ordinal fit defaults to the latent scale",
    identical(attr(nest_estimand(mo, chord_type, bounds = FALSE),
                   "nestimand")$scale, "latent"))
chk("scale: everything else defaults to the response scale",
    identical(attr(nest_estimand(mf, chord_type, bounds = FALSE),
                   "nestimand")$type, "response"))
chk("scale: an explicit choice is honoured either way",
    identical(attr(nest_estimand(mf, chord_type, bounds = FALSE), "nestimand")$scale, "latent"))

chk("class: the engine's own class comes first, so dispatch is undisturbed",
    identical(class(mf)[1], "lm") &&
    identical(class(fit_from(sp))[1], "lm"))

## ---- engine arguments on the latent route ----------------------------------
## The latent route never calls the prediction function, so its arguments have
## nowhere to go; an ordinal fit defaults to that route, which makes this easy
## to hit by accident.
chk("latent: conf_level does apply, and widens the interval",
    { a <- as.data.frame(nest_estimand(mf, chord_type,
             conf_level = 0.99, bounds = FALSE))
      b <- as.data.frame(nest_estimand(mf, chord_type,
             conf_level = 0.90, bounds = FALSE))
      all(a$conf.low < b$conf.low) && all(a$conf.high > b$conf.high) })
chk("latent: the response route still takes them",
    nrow(as.data.frame(nest_estimand(mf, chord_type, p_adjust = "holm",
                                bounds = FALSE))) == 6)

## ---- one argument names the quantity ---------------------------------------
meta_of <- function(...) attr(nest_estimand(mf, chord_type, bounds = FALSE, ...), "nestimand")
chk("type: eta is computed from the coefficients, on any class",
    identical(meta_of(type = "eta")$scale, "latent") &&
    identical(attr(nest_estimand(mo, chord_type, type = "eta", bounds = FALSE), "nestimand")$scale, "latent"))
chk("type: an engine name is not reinterpreted as eta",
    { e <- err_of(suppressMessages(nest_estimand(mo, chord_type, type = "latent",
                    bounds = FALSE)))
      identical(e, "") || grepl("type", e) })
chk("type: anything but the linear predictor goes to the prediction function",
    identical(meta_of(type = "response", p_adjust = "none")$scale, "response"))
chk("type: the default is the response scale, and eta for an ordinal fit",
    identical(meta_of()$type, "response") &&
    identical(attr(nest_estimand(mo, chord_type, bounds = FALSE),
                   "nestimand")$type, "eta"))
chk("type: the two agree in a linear model, whichever machinery is used",
    { a <- as.data.frame(nest_estimand(mf, chord_type, bounds = FALSE))
      b <- as.data.frame(nest_estimand(mf, chord_type, type = "response",
                                  bounds = FALSE))
      max(abs(a$estimate - b$estimate[match(a$term, b$term)])) < 1e-8 })
chk("type: it is recorded and printed",
    any(grepl("type: eta",
        capture.output(print(nest_estimand(mf, chord_type, type = "eta",
                                      bounds = FALSE))))))

## ---- the response scale on an ordinal fit is worth a word ------------------
chk("ordinal: the response scale says what it gives and what to do instead",
    { msg <- NULL
      withCallingHandlers(
        tryCatch(nest_estimand(mo, chord_type, type = "response", bounds = FALSE), error = function(e) NULL),
        message = function(m) { msg <<- conditionMessage(m)
                                invokeRestart("muffleMessage") })
      grepl("one per outcome category", msg) && grepl("pairwise", msg) })
chk("ordinal: an engine refusal on that scale is explained, not passed on raw",
    { e <- err_of(suppressMessages(nest_estimand(mo, chord_type, type = "response",
                    bounds = FALSE)))
      identical(e, "") || (grepl("one value per outcome category", e) &&
                           grepl("The engine said", e)) })
chk("ordinal: the default route is unaffected by any of this",
    nrow(as.data.frame(nest_estimand(mo, chord_type, bounds = FALSE))) == 6)

chk("ordinal: the predictive type warns that codes are being averaged",
    { msgs <- character(0)
      withCallingHandlers(
        tryCatch(nest_estimand(mo, chord_type, type = "prediction", bounds = FALSE), error = function(e) NULL),
        message = function(m) { msgs <<- c(msgs, conditionMessage(m))
                                invokeRestart("muffleMessage") })
      any(grepl("category codes", msgs)) && any(grepl("interval one", msgs)) })
chk("prediction: a non-ordinal fit draws no such note",
    { msg <- NULL
      withCallingHandlers(
        tryCatch(nest_estimand(mf, chord_type, type = "prediction", bounds = FALSE), error = function(e) NULL),
        message = function(m) { msg <<- conditionMessage(m)
                                invokeRestart("muffleMessage") })
      is.null(msg) })

chk("ordinal: the size of the comparison set is stated before the work starts",
    { msgs <- character(0)
      withCallingHandlers(
        tryCatch(nest_estimand(mo, chord_type, type = "response", bounds = FALSE), error = function(e) NULL),
        message = function(m) { msgs <<- c(msgs, conditionMessage(m))
                                invokeRestart("muffleMessage") })
      any(grepl("378 contrasts rather than 6", msgs)) })
chk("ordinal: the linear predictor draws no such warning",
    { msgs <- character(0)
      withCallingHandlers(
        nest_estimand(mo, chord_type, bounds = FALSE),
        message = function(m) { msgs <<- c(msgs, conditionMessage(m))
                                invokeRestart("muffleMessage") })
      !any(grepl("contrasts where the linear predictor", msgs)) })

## ---- a hypothesis of the user's own ----------------------------------------
## The package forms one from `contrast`; supplying another would give the
## engine two, so the user's replaces it, and the substitution is announced.
chk("hypothesis: a supplied one replaces the contrast's, without colliding",
    { e <- suppressMessages(nest_estimand(mf, chord_type, hypothesis = "reference",
                                     bounds = FALSE))
      nrow(as.data.frame(e)) == 3 &&
      sum(grepl("hypothesis", attr(e, "nestimand")$code)) == 1 })
chk("hypothesis: the substitution is reported, with what it costs",
    { msg <- NULL
      withCallingHandlers(nest_estimand(mf, chord_type, hypothesis = "reference",
                                   bounds = FALSE),
        message = function(m) { msg <<- conditionMessage(m)
                                invokeRestart("muffleMessage") })
      grepl("instead of contrast", msg) && grepl("subtraction", msg) })
chk("hypothesis: it is refused alongside an interaction, which defines its own",
    grepl("both define which comparisons",
          err_of(nest_estimand(mf, chord_type:inversion, hypothesis = "reference"))))
chk("hypothesis: the emitted code runs",
    { e <- suppressMessages(nest_estimand(mf, chord_type, hypothesis = "reference",
                                     bounds = FALSE))
      env <- new.env(); assign("mf", mf, env); assign("sp", sp, env)
      r <- tryCatch(eval(parse(text = paste(attr(e, "nestimand")$code,
                                            collapse = "\n")), envir = env),
                    error = function(x) NULL)
      !is.null(r) })

chk("messages: one note on the response scale, and none once a hypothesis is given",
    { m1 <- character(0); m2 <- character(0)
      withCallingHandlers(tryCatch(nest_estimand(mo, chord_type, type = "response",
          bounds = FALSE), error = function(e) NULL),
        message = function(x) { m1 <<- c(m1, conditionMessage(x))
                                invokeRestart("muffleMessage") })
      withCallingHandlers(tryCatch(nest_estimand(mo, chord_type, type = "response",
          hypothesis = "reference", bounds = FALSE),
          error = function(e) NULL),
        message = function(x) { m2 <<- c(m2, conditionMessage(x))
                                invokeRestart("muffleMessage") })
      length(m1) == 1 && !any(grepl("comparing them all", m2)) })
chk("messages: the default scale is silent",
    { m <- character(0)
      withCallingHandlers(nest_estimand(mo, chord_type, bounds = FALSE),
        message = function(x) { m <<- c(m, conditionMessage(x))
                                invokeRestart("muffleMessage") })
      length(m) == 0 })

chk("messages: one note on the response scale, and none once a hypothesis is given",
    { m1 <- character(0); m2 <- character(0)
      withCallingHandlers(tryCatch(nest_estimand(mo, chord_type, type = "response",
          bounds = FALSE), error = function(e) NULL),
        message = function(x) { m1 <<- c(m1, conditionMessage(x))
                                invokeRestart("muffleMessage") })
      withCallingHandlers(tryCatch(nest_estimand(mo, chord_type, type = "response",
          hypothesis = "reference", bounds = FALSE),
          error = function(e) NULL),
        message = function(x) { m2 <<- c(m2, conditionMessage(x))
                                invokeRestart("muffleMessage") })
      length(m1) == 1 && !any(grepl("comparing them all", m2)) })
chk("messages: the default scale is silent",
    { m <- character(0)
      withCallingHandlers(nest_estimand(mo, chord_type, bounds = FALSE),
        message = function(x) { m <<- c(m, conditionMessage(x))
                                invokeRestart("muffleMessage") })
      length(m) == 0 })

chk("messages: the substitution note says what changes, in plain terms",
    { msg <- NULL
      withCallingHandlers(nest_estimand(mf, chord_type, hypothesis = "reference",
                                   bounds = FALSE),
        message = function(m) { msg <<- conditionMessage(m)
                                invokeRestart("muffleMessage") })
      grepl("using your `hypothesis`", msg) &&
      grepl("which way round each subtraction goes", msg) })
chk("messages: a small job draws no size warning",
    { msg <- character(0)
      withCallingHandlers(tryCatch(nest_estimand(mo, chord_type, type = "response",
          route = "cells", hypothesis = "reference", bounds = FALSE), error = function(e) NULL),
        message = function(m) { msg <<- c(msg, conditionMessage(m))
                                invokeRestart("muffleMessage") })
      !any(grepl("take a while", msg)) })

## ---- the star form with a hypothesis of the user's own ---------------------
chk("targets: a * b with a hypothesis returns all three, and says how they divide",
    { z <- character(0)
      r <- withCallingHandlers(nest_estimand(mf, chord_type * inversion,
             hypothesis = "reference", bounds = FALSE),
           message = function(m) { z <<- c(z, conditionMessage(m))
                                   invokeRestart("muffleMessage") })
      identical(names(r), c("chord_type", "inversion", "chord_type:inversion")) &&
      any(grepl("expands as a formula does: 3 results", z)) })
chk("targets: without a hypothesis all three parts come back",
    { r <- nest_estimand(mf, chord_type * inversion, bounds = FALSE)
      length(r) == 3 })
chk("messages: a note is said once, not once per part",
    { z <- character(0)
      withCallingHandlers(nest_estimand(mf, chord_type * inversion,
        hypothesis = "reference", bounds = FALSE),
        message = function(m) { z <<- c(z, conditionMessage(m))
                                invokeRestart("muffleMessage") })
      sum(grepl("expands as a formula does", z)) == 1 &&
      !any(grepl("using your `hypothesis`", z)) })
chk("messages: the ordinal note does not fire on a gaussian fit",
    { z <- character(0)
      withCallingHandlers(nest_estimand(mf, chord_type, type = "response",
        bounds = FALSE),
        message = function(m) { z <<- c(z, conditionMessage(m))
                                invokeRestart("muffleMessage") })
      !any(grepl("outcome category", z)) })

## ---- interactions within a grouped output ----------------------------------
## Where the engine returns one row per condition per group - the categories of
## an ordinal fit - a difference of differences is only interpretable within a
## group, so the matrix is block-diagonal: the same columns, once per group.
dg <- data.frame(chord_type = rep(c("dim", "min", "maj"), each = 4),
                 inversion  = rep(c("0", "1"), times = 6),
                 group      = rep(c("1", "1", "2", "2"), times = 3))
Hg <- interaction_matrix(dg, c("chord_type", "inversion"))
chk("interaction: one block of contrasts per group",
    ncol(Hg) == 6 && nrow(Hg) == 12)
chk("interaction: every contrast stays inside one group",
    all(apply(Hg, 2, function(col) length(unique(dg$group[col != 0])) == 1)))
chk("interaction: the group is named in the label",
    all(grepl(" \\| [12]$", colnames(Hg))))
chk("interaction: each is still a difference of differences",
    all(abs(colSums(Hg)) < 1e-12) && all(colSums(abs(Hg)) == 4))
chk("interaction: ungrouped output is unchanged",
    ncol(interaction_matrix(dg[dg$group == "1", c("chord_type", "inversion")],
                            c("chord_type", "inversion"))) == 3)

## ---- the route reaches every part of a call --------------------------------
## The interaction had its own grid construction, which ignored `route`: with
## route = "cells" the marginal contrasts used 9 rows while the interaction
## quietly built the whole G-computation grid.
chk("route: the interaction honours it too",
    { cd <- suppressMessages(nest_estimand(mf, chord_type:inversion, route = "cells",
              dry_run = TRUE, bounds = FALSE, p_adjust = "none"))
      grepl("cell_grid(spec)", cd, fixed = TRUE) &&
      !grepl("counterfactual_grid", cd, fixed = TRUE) })
chk("route: and gives the same answer either way in a linear model",
    { a <- as.data.frame(nest_estimand(mf, chord_type:inversion, route = "cells",
                                  bounds = FALSE))
      b <- as.data.frame(nest_estimand(mf, chord_type:inversion, bounds = FALSE))
      max(abs(a$estimate - b$estimate)) < 1e-8 })
chk("dry_run: several targets give one script, not a list of them",
    { cd <- suppressMessages(nest_estimand(mf, chord_type * inversion,
              dry_run = TRUE, bounds = FALSE))
      inherits(cd, "nestimand_code") &&
      length(capture.output(show_code(cd))) > 10 &&
      grepl("--- chord_type ---", cd, fixed = TRUE) })

## ---- the two routes are different quantities on a nonlinear scale ----------
## On the link scale averaging the design rows and averaging the predictions
## are the same operation, so the routes agree exactly. On the response scale
## the link stands between them and they do not.
do3 <- dat
do3$rating <- factor(round(pmin(pmax(do3$response, 1), 4)), ordered = TRUE)
do3$z <- as.numeric(scale(do3$training))
sp3 <- nesting_spec(do3, rating ~ chord_type * inversion + z,
                    "inversion %in% chord_type", fit = "clm")
m3 <- fit_from(sp3)
lat3 <- function(rt) as.data.frame(nest_estimand(m3, chord_type, route = rt,
                                            bounds = FALSE))
chk("routes: identical on the link scale, to machine precision",
    { a <- lat3("cells"); b <- lat3("g_computation")
      max(abs(a$estimate - b$estimate[match(a$term, b$term)])) < 1e-12 })
resp3 <- function(...) suppressMessages(as.data.frame(
  nest_estimand(m3, chord_type, type = "response", hypothesis = "reference",
           bounds = FALSE, ...)))
chk("routes: they differ on the response scale, as a nonlinear link requires",
    { a <- resp3(route = "cells"); b <- resp3()
      max(abs(a$estimate - b$estimate[match(a$term, b$term)])) > 1e-4 })
chk("subsample: estimates the same quantity as the full grid",
    { set.seed(2)
      max(abs(resp3()$estimate - resp3(subsample = 200)$estimate)) < 0.02 })
chk("subsample: the emitted script draws the same rows again",
    { set.seed(2)
      e1 <- suppressMessages(nest_estimand(m3, chord_type, subsample = 50,
              bounds = FALSE))
      env <- new.env(parent = globalenv())
      assign("m3", m3, env); assign("sp3", sp3, env)
      r <- eval(parse(text = paste(attr(e1, "nestimand")$code, collapse = "\n")),
                envir = env)
      isTRUE(all.equal(as.data.frame(r)$estimate,
                       as.data.frame(e1)$estimate)) })
chk("subsample: it is announced, since it trades exactness for time",
    { msg <- character(0); set.seed(2)
      withCallingHandlers(tryCatch(nest_estimand(m3, chord_type, type = "response",
          subsample = 50, bounds = FALSE),
          error = function(e) NULL),
        message = function(m) { msg <<- c(msg, conditionMessage(m))
                                invokeRestart("muffleMessage") })
      any(grepl("Monte Carlo error", msg)) })
chk("print: the outcome category is shown where the output is grouped",
    { d <- data.frame(group = "3", term = "a - b", estimate = 0.5)
      x <- structure(d, class = c("nestimand_estimand", "data.frame"),
                     nestimand = list(policy = "equal", route = "cells",
                                      contrast = "pairwise", type = "response",
                                      code = "x"))
      any(grepl("group", capture.output(print(x)))) })

## ---- Bayesian summaries on the prediction route too ------------------------
chk("posterior: a frequentist fit keeps its statistic and p-value",
    { d <- as.data.frame(nest_estimand(mf, chord_type, bounds = FALSE))
      all(c("statistic", "p.value") %in% names(d)) && !"pd" %in% names(d) })
chk("posterior: the summariser leaves a non-Bayesian fit alone",
    { e0 <- nest_estimand(mf, chord_type, bounds = FALSE)
      identical(names(as.data.frame(add_posterior_summary(e0, mf))),
                names(as.data.frame(e0))) })
chk("posterior: pd is the larger tail, computed from draws",
    { z <- c(-1, 2, 3, 4)            # 3 of 4 positive
      abs(max(mean(z > 0), mean(z < 0)) - 0.75) < 1e-12 })

## ---- the two levers on a large posterior -----------------------------------
chk("subsample: the counts reported respect it",
    { z <- character(0)
      withCallingHandlers(nest_estimand(mf, chord_type, subsample = 50,
                                   bounds = FALSE),
        message = function(m) { z <<- c(z, conditionMessage(m))
                                invokeRestart("muffleMessage") })
      any(grepl("random 50 of", z)) })
chk("subsample: the code records the draw so the same rows come back",
    { cd <- suppressMessages(nest_estimand(mf, chord_type, subsample = 50,
              dry_run = TRUE, bounds = FALSE))
      grepl("set.seed(", cd, fixed = TRUE) &&
      grepl("sample.int(nrow(spec$data), 50)", cd, fixed = TRUE) &&
      max(nchar(strsplit(as.character(cd), "\n")[[1]])) < 200 })
chk("subsample: it estimates the same quantity",
    { set.seed(11)
      a <- as.data.frame(nest_estimand(mf, chord_type, bounds = FALSE))$estimate
      b <- suppressMessages(as.data.frame(nest_estimand(mf, chord_type,
             subsample = 200, bounds = FALSE))$estimate)
      max(abs(a - b)) < 0.05 })
chk("ndraws: harmless where there is no posterior to thin",
    nrow(as.data.frame(suppressMessages(nest_estimand(mf, chord_type,
      ndraws = 100, bounds = FALSE)))) == 6)

## ---- the reorder check should not cry wolf ---------------------------------

chk("arguments: subsample, ndraws and a hypothesis together",
    { r <- suppressMessages(nest_estimand(mf, chord_type, type = "response",
             subsample = 100, ndraws = 100, hypothesis = "reference",
             bounds = FALSE))
      nrow(as.data.frame(r)) == 3 })
chk("arguments: every combination of the levers runs",
    { ok <- TRUE
      for (sub in list(NULL, 100)) for (dr in list(NULL, 100))
        for (hyp in list(NULL, "reference")) {
          cl <- quote(nest_estimand(mf, "chord_type", bounds = FALSE))
          if (!is.null(sub)) cl$subsample <- sub
          if (!is.null(dr))  cl$ndraws <- dr
          if (!is.null(hyp)) cl$hypothesis <- hyp
          r <- try(suppressMessages(eval(cl)), silent = TRUE)
          if (inherits(r, "try-error")) ok <- FALSE
        }
      ok })

chk("ndraws: the note quotes the number asked for, not a fixed one",
    { body <- paste(deparse(nest_estimand), collapse = " ")
      grepl("1/sqrt", body) && !grepl("sqrt\\(500\\)", body) })
chk("ndraws: the size note counts the draws in use",
    { body <- paste(deparse(nest_estimand), collapse = " ")
      grepl("dots\\$ndraws", body) })
chk("messages: advice already taken is not repeated",
    { body <- paste(deparse(nest_estimand), collapse = " ")
      grepl("is.null\\(subsample\\)", body) &&
      grepl("is.null\\(dots\\$ndraws\\)", body) })

chk("ndraws: the engine's own spelling works on either route",
    { a <- suppressMessages(nest_estimand(mo, chord_type, ndraws = 500, dry_run = TRUE, bounds = FALSE))
      b <- suppressMessages(nest_estimand(mf, chord_type, type = "response",
             ndraws = 500, dry_run = TRUE, bounds = FALSE))
      grepl("ndraws = 500", a, fixed = TRUE) &&
      grepl("ndraws = 500", b, fixed = TRUE) })
chk("hypothesis: refused on the linear predictor, for the reason that applies",
    grepl("no groups to compare within",
          err_of(nest_estimand(mf, chord_type, type = "eta",
                          hypothesis = ~ pairwise | group))))

## ---- a model class the prediction machinery does not handle ----------------
## clmm is one: marginaleffects supports no type for it, so the response scale
## is unavailable however it is asked for. The linear predictor is unaffected,
## being taken from the coefficients.
if (requireNamespace("ordinal", quietly = TRUE)) {
  spmm <- nesting_spec(do2, rating ~ chord_type * inversion + (1 | participant),
                       "inversion %in% chord_type", fit = "clmm")
  mmm <- suppressWarnings(fit_from(spmm))
  chk("clmm: the linear predictor works",
      nrow(as.data.frame(nest_estimand(mmm, chord_type, bounds = FALSE))) == 6)
  chk("clmm: the response scale says the class is unsupported, and what does work",
      { e <- err_of(nest_estimand(mmm, chord_type, type = "response", bounds = FALSE))
        grepl("does not support models of class", e) &&
        grepl('type = "eta" works', e) })
  chk("clmm: the message no longer names the retired `scale` argument",
      { body <- paste(deparse(ordinal_response_type), collapse = " ")
        !grepl('scale = ..latent', body) })
}

chk("messages: a star call says once what its parts would each repeat",
    { z <- character(0)
      withCallingHandlers(nest_estimand(mf, chord_type * inversion,
        policy = "proportional", hypothesis = "reference", bounds = FALSE),
        message = function(m) { z <<- c(z, conditionMessage(m))
                                invokeRestart("muffleMessage") })
      length(z) == 1 && grepl("expands as a formula does: 3 results", z) &&
      grepl("hypothesis", z) && grepl("policy", z) })
chk("messages: a single target still explains itself",
    { z <- character(0)
      withCallingHandlers(nest_estimand(mf, chord_type, hypothesis = "reference",
        bounds = FALSE),
        message = function(m) { z <<- c(z, conditionMessage(m))
                                invokeRestart("muffleMessage") })
      length(z) == 1 && grepl("using your `hypothesis`", z) })
chk("messages: the star flag is cleared afterwards",
    { nest_estimand(mf, chord_type * inversion, bounds = FALSE)
      z <- character(0)
      withCallingHandlers(nest_estimand(mf, chord_type, hypothesis = "reference",
        bounds = FALSE),
        message = function(m) { z <<- c(z, conditionMessage(m))
                                invokeRestart("muffleMessage") })
      length(z) == 1 })

## ---- the reorder check does not inherit the estimand's subsample -----------
## A subsample is taken to make the estimand affordable. The check asks a
## different question - does the estimand move when levels are permuted - and
## must work from the full declaration: a small sample may not contain every
## cell, and a spec rebuilt from it would have fewer of them.

## ---- what a subsample does and does not touch ------------------------------
## The grid crosses every sampled row with every cell, so cell coverage in the
## sample is irrelevant. The policy is another matter: it counts how often each
## version is realized, and must be taken from the whole data.
chk("subsample: cell coverage in the sample does not matter",
    { few <- sp$data[sp$data$chord_type == "aug", ][1:5, ]
      g <- counterfactual_grid(sp, few, nest_policy(sp, "chord_type", "equal"))
      length(unique(as.character(g$cell))) == nrow(sp$cells) })
chk("subsample: a proportional policy is taken from the whole data",
    { set.seed(4)
      a <- as.data.frame(nest_estimand(mf, chord_type, policy = "proportional",
             bounds = FALSE))
      b <- suppressMessages(as.data.frame(nest_estimand(mf, chord_type,
             policy = "proportional", subsample = 40, bounds = FALSE)))
      !anyNA(b$estimate) &&
      max(abs(a$estimate - b$estimate[match(a$term, b$term)])) < 1e-10 })

## ---- every emitted script must run, on every route and scale ---------------
## The latent branch referred to `cells` without emitting the line that makes
## it: the restriction was applied internally but not written down.
chk("emitted code: runs for every target, scale and bounds setting",
    { ok <- TRUE
      for (tgt in c("chord_type", "inversion"))
        for (ty in c("link", "response"))
          for (bd in c(TRUE, FALSE)) {
            cl <- bquote(nest_estimand(mf, .(tgt), type = .(ty), bounds = .(bd)))
            e <- tryCatch(suppressMessages(eval(cl)), error = function(x) NULL)
            if (is.null(e)) next
            env <- new.env(parent = globalenv())
            assign("mf", mf, env); assign("sp", sp, env)
            r <- tryCatch(eval(parse(text = paste(attr(e, "nestimand")$code,
                                                  collapse = "\n")), envir = env),
                          error = function(x) NULL)
            if (is.null(r)) ok <- FALSE
          }
      ok })
chk("emitted code: the latent branch writes the restriction it uses",
    { cd <- nest_estimand(mf, inversion, bounds = TRUE,
                     dry_run = TRUE)
      grepl("cells <- subset(", cd, fixed = TRUE) &&
      grepl("cells = cells", cd, fixed = TRUE) })

## ---- a random term with nothing to translate just prints -------------------
if (requireNamespace("lme4", quietly = TRUE)) {
  sp_sl <- nesting_spec(d6, response ~ chord_type * inversion +
                        (chord_type | participant), "inversion %in% chord_type",
                        fit = "lmer")
  m_sl <- suppressWarnings(fit_from(sp_sl))
  out_sl <- capture.output(print(random_covariance(m_sl, space = "effects")))
  chk("random: a term that is not a cell covariance prints without commentary",
      !any(grepl("no effect-space counterpart|shown as fitted|translation",
                 out_sl)))
  chk("random: its header claims no parameterization it does not have",
      grepl("^Group `participant`  \\(4 dimensions\\)$", out_sl[1]))
  chk("random: a cell covariance says which parameterization it is in",
      { out2 <- capture.output(print(random_covariance(mm2, space = "effects")))
        grepl("effects parameterization", out2[1]) })
  chk("random: the fact is still recorded for code that needs it",
      identical(attr(random_covariance(m_sl)[[1]], "translated"), FALSE))
}

## ---- the summary says which engine produced it -----------------------------
chk("nest_summary: the header names the fitting function",
    { ln <- capture.output(print(nest_summary(mf)))[1]
      grepl("fitted with lm\\(\\)", ln) })
chk("nest_summary: and the formula it was fitted with",
    { ln <- capture.output(print(nest_summary(mf)))[2]
      grepl("response ~", ln) && grepl("cell", ln) })
chk("nest_summary: the engine is recorded, not only printed",
    identical(attr(nest_summary(mf), "nestimand_fit"), "lm"))
if (requireNamespace("ordinal", quietly = TRUE)) {
  chk("nest_summary: an ordinal fit names its own engine",
      grepl("ordinal::clm\\(\\)",
            capture.output(print(nest_summary(mo)))[1]))
}
chk("nest_summary: the engine names map to real functions",
    identical(engine_call("brm"), "brms::brm()") &&
    identical(engine_call("glmer"), "lme4::glmer()") &&
    identical(engine_call("something"), "something()"))

## ---- an identity link makes the two scales one quantity --------------------
## The contrast can then be taken from the coefficients, which averages the
## design matrix rather than one prediction per row. Gated on the link, not the
## family: gaussian(link = "log") is not identity.
chk("shortcut: an identity link routes the response scale through the coefficients",
    identical(attr(nest_estimand(mf, chord_type, type = "response", bounds = FALSE), "nestimand")$scale, "latent"))
chk("shortcut: the type asked for is still what is reported",
    identical(attr(nest_estimand(mf, chord_type, type = "response", bounds = FALSE), "nestimand")$type, "response"))
chk("shortcut: it gives the same numbers as the prediction route",
    { a <- as.data.frame(nest_estimand(mf, chord_type, type = "response",
                                  bounds = FALSE))
      b <- as.data.frame(nest_estimand(mf, chord_type, bounds = FALSE))
      max(abs(a$estimate - b$estimate[match(a$term, b$term)])) < 1e-10 })
chk("shortcut: the emitted code says why it was taken",
    any(grepl("the link is the identity",
              attr(nest_estimand(mf, chord_type, type = "response", bounds = FALSE), "nestimand")$code)))
chk("shortcut: a non-identity link does not take it",
    { d3 <- dat; d3$y <- exp(dat$response / 3)
      sp3b <- nesting_spec(d3, y ~ chord_type * inversion,
                           "inversion %in% chord_type", fit = "glm",
                           family = gaussian(link = "log"))
      m3b <- fit_from(sp3b)
      identical(attr(nest_estimand(m3b, chord_type, type = "response", bounds = FALSE), "nestimand")$scale,
                "response") })
chk("shortcut: it stands aside for a hypothesis of the user's own",
    identical(attr(suppressMessages(nest_estimand(mf, chord_type, type = "response",
                     hypothesis = "reference", bounds = FALSE)), "nestimand")$scale, "response"))
chk("model_link: reads the link, defaulting to identity where there is none",
    identical(model_link(mf), "identity"))

chk("shortcut: unit weights are folded into the design-row average",
    identical(attr(nest_estimand(mw, chord_type, weights = "wt_hi", bounds = FALSE), "nestimand")$scale, "latent"))
chk("weights: the two routes agree once they are",
    { a <- as.data.frame(nest_estimand(mw, chord_type, weights = "wt_hi",
             bounds = FALSE))
      b <- as.data.frame(nest_estimand(mw, chord_type, weights = "wt_hi",
             bounds = FALSE, p_adjust = "none"))
      max(abs(a$estimate - b$estimate[match(a$term, b$term)])) < 1e-8 })
chk("weights: a wrong length is refused on this route too",
    grepl("one weight per row",
          err_of(latent_estimand(mw, "chord_type", spec = spw,
                                 weights = c(1, 2)))))
chk("shortcut: an ordinal fit is not mistaken for an identity link",
    !identical(model_link(mo, spo), "identity"))
chk("shortcut: a self-fitted model without cell coefficients cannot take it",
    !linear_map_available(lm(response ~ chord_type * inversion, data = sp$data),
                          sp, sp$data))
chk("latent: `data` reaches the contrast, so a subset restricts it",
    { a <- as.data.frame(nest_estimand(mw, chord_type, weights = "wt_hi",
             bounds = FALSE))$estimate
      b <- as.data.frame(nest_estimand(mw, chord_type,
             data = subset(mw_data_hi <- spw$data, training > 7),
             bounds = FALSE))$estimate
      max(abs(a - b)) < 1e-8 })

## ---- the equality does not depend on balance -------------------------------
set.seed(7)
cells_u <- data.frame(chord_type = c(rep(c("dim", "min", "maj"), each = 3), "aug"),
                      inversion = c(rep(c("0", "1", "2"), 3), "none"))
n_u <- c(120, 40, 25, 90, 30, 15, 200, 60, 35, 50)
du <- cells_u[rep(seq_len(10), n_u), ]
du$chord_type <- factor(du$chord_type, levels = c("aug", "dim", "min", "maj"))
du$inversion  <- factor(du$inversion,  levels = c("none", "0", "1", "2"))
du$x <- rnorm(nrow(du), ifelse(du$chord_type == "aug", 6, 3), 2)
du$y <- rnorm(nrow(du), 4 + 0.5 * (du$chord_type == "maj") + 0.3 * du$x, 1)
spu2 <- nesting_spec(du, y ~ chord_type * inversion * x, "inversion %in% chord_type")
mu3 <- fit_from(spu2)
chk("shortcut: agrees with the prediction route on badly unbalanced data",
    { ok <- TRUE
      for (pol in c("equal", "proportional")) {
        a <- as.data.frame(nest_estimand(mu3, chord_type, policy = pol, bounds = FALSE, p_adjust = "none"))
        b <- as.data.frame(nest_estimand(mu3, chord_type, policy = pol, bounds = FALSE))
        if (max(abs(a$estimate - b$estimate[match(a$term, b$term)])) > 1e-10)
          ok <- FALSE
      }
      ok })

## ---- a type the engine rejects is the engine's to refuse -------------------
## The package recognises "link" and its synonyms, but not in defiance of the
## engine: a type the user typed that marginaleffects will not accept for this
## model is passed through, so the error comes from the package that owns the
## vocabulary. A type the package chooses for itself is not checked, since some
## classes the engine cannot predict from at all.
chk("type: an engine name it rejects reaches it",
    { e <- err_of(suppressMessages(nest_estimand(mf, chord_type, type = "link",
                    bounds = FALSE)))
      grepl("Assertion on 'type'", e) || grepl("Must be element of set", e) })
chk("type: an engine name it accepts is used",
    { d_b <- dat; d_b$bin <- as.numeric(dat$response > 4)
      sp_b <- nesting_spec(d_b, bin ~ chord_type * inversion,
                           "inversion %in% chord_type", fit = "glm",
                           family = binomial())
      nrow(as.data.frame(suppressMessages(nest_estimand(fit_from(sp_b), chord_type,
             type = "link", bounds = FALSE)))) == 6 })
chk("type: the package's own default is not put to that test",
    nrow(as.data.frame(nest_estimand(mf, chord_type, bounds = FALSE))) == 6 &&
    nrow(as.data.frame(nest_estimand(mo, chord_type, bounds = FALSE))) == 6)
chk("engine_accepts: reports what this class takes",
    engine_accepts(mf, sp, "response") && !engine_accepts(mf, sp, "link"))
chk("type: eta is the package's own name, and no engine claims it",
    !("eta" %in% marginaleffects:::type_dictionary$type))
chk("type: a slower engine equivalent draws a note pointing at eta",
    { z <- character(0)
      d_g <- dat
      sp_g <- nesting_spec(d_g, response ~ chord_type * inversion,
                           "inversion %in% chord_type", fit = "glm",
                           family = gaussian())
      withCallingHandlers(nest_estimand(fit_from(sp_g), chord_type, type = "link",
                                   bounds = FALSE),
        message = function(m) { z <<- c(z, conditionMessage(m))
                                invokeRestart("muffleMessage") })
      any(grepl('type = "eta" gives these same numbers', z)) })

chk("type: no equivalence is offered for a type about to be refused",
    { z <- character(0)
      withCallingHandlers(
        tryCatch(nest_estimand(mf, chord_type, type = "link", bounds = FALSE), error = function(e) NULL),
        message = function(m) { z <<- c(z, conditionMessage(m))
                                invokeRestart("muffleMessage") })
      length(z) == 0 })
chk("type: the equivalence is offered where the engine will produce it",
    { sp_g2 <- nesting_spec(dat, response ~ chord_type * inversion,
                            "inversion %in% chord_type", fit = "glm",
                            family = gaussian())
      z <- character(0)
      withCallingHandlers(
        nest_estimand(fit_from(sp_g2), chord_type, type = "link", bounds = FALSE),
        message = function(m) { z <<- c(z, conditionMessage(m))
                                invokeRestart("muffleMessage") })
      any(grepl('type = "eta" gives these same numbers', z)) })

## ---- engines that return a Matrix rather than a matrix ---------------------
## lme4 returns a dpoMatrix, and arithmetic on it dispatches to methods whose
## result is not always a base array. Coercing once keeps everything downstream
## ordinary; the failure it caused was in the standard errors, not the estimates.
if (requireNamespace("lme4", quietly = TRUE)) {
  chk("vcov: a Matrix is coerced to a base matrix",
      { V <- vcov_beta(mm2, names(coef_vector(mm2)))
        is.matrix(V) && !isS4(V) })
  chk("vcov: standard errors come through finite on a mixed fit",
      { d <- as.data.frame(nest_estimand(mm2, chord_type, bounds = FALSE))
        all(is.finite(d$std.error)) && all(d$std.error > 0) })
  chk("vcov: the bounds path works on one too",
      { b <- attr(nest_estimand(mm2, chord_type, bounds = TRUE),
                  "nestimand")$bounds
        nrow(b) == 6 && all(is.finite(b$policy_low)) })
}
chk("quad_form_diag: the row-wise form equals the diagonal of the product",
    { set.seed(2); C <- matrix(rnorm(12), 3, 4); V <- crossprod(matrix(rnorm(16), 4))
      max(abs(quad_form_diag(C, V) - diag(C %*% V %*% t(C)))) < 1e-12 })
chk("vcov: a covariance whose names do not match is refused, not multiplied",
    grepl("appear in the model's covariance",
          err_of(vcov_beta(mf, c("nope", "also_nope")))))

## ---- what each engine accepts, read from the installed version -------------
et <- engine_types()
chk("engine_types: reports a row per class, with eta first",
    identical(as.data.frame(et)$class[1], "any") &&
    grepl("eta", as.data.frame(et)$types[1]) &&
    all(c("lm", "glm", "clm") %in% as.data.frame(et)$class))
chk("engine_types: matches what the engine actually does",
    { row_lm <- as.data.frame(et)$types[as.data.frame(et)$class == "lm"]
      grepl("response", row_lm) && !grepl("\\blink\\b", row_lm) &&
      engine_accepts(mf, sp, "response") && !engine_accepts(mf, sp, "link") })
chk("engine_types: a single fit narrows it to that class",
    { e1 <- as.data.frame(engine_types(mf))
      nrow(e1) == 2 && identical(e1$class[2], "lm") })
chk("engine_types: an unsupported class says so rather than showing nothing",
    { e2 <- as.data.frame(engine_types(structure(list(), class = "nosuchclass")))
      grepl("no entry", e2$types[nrow(e2)]) })
chk("engine_types: it names the version it read",
    grepl(as.character(utils::packageVersion("marginaleffects")),
          capture.output(print(et))[1], fixed = TRUE))

## ---- a transformation written into the formula ----------------------------
## all.vars() sees through I(x^2) to x, so such a term was silently reduced to
## the untransformed column - and its interaction with the conditions lost.
chk("declaration: an in-place transformation is refused, not reduced",
    { e <- err_of(nesting_spec(dat, response ~ chord_type * inversion *
                                 I(training^2), "inversion %in% chord_type"))
      grepl("transforms a variable in place", e) &&
      grepl("as a column of the data", e) })
chk("declaration: the same model as a column is accepted",
    { d_t <- dat; d_t$t2 <- d_t$training^2
      sp_t <- nesting_spec(d_t, response ~ chord_type * inversion * t2,
                           "inversion %in% chord_type")
      grepl("cell:t2", paste(deparse(cell_formula(sp_t)), collapse = ""),
            fixed = TRUE) })
chk("declaration: random-effects terms are not mistaken for transformations",
    { sp_r <- nesting_spec(dat, response ~ chord_type * inversion +
                             (1 | participant), "inversion %in% chord_type",
                           fit = "lmer")
      inherits(sp_r, "nesting_spec") })

## ---- the eta route honours the route it was given --------------------------
chk("eta: route reaches the contrast matrix",
    { pol_e <- nest_policy(sp, "chord_type", "equal")
      Mc <- policy_contrast_matrix(sp, "chord_type", pol_e, sp$data, mf,
                                   route = "cells")
      Mg <- policy_contrast_matrix(sp, "chord_type", pol_e, sp$data, mf)
      nrow(Mc) == nrow(Mg) && max(abs(Mc - Mg)) < 1e-8 })
chk("eta: and the emitted code records it",
    grepl('route = "cells"',
          nest_estimand(mf, chord_type, type = "eta", route = "cells",
                   dry_run = TRUE, bounds = FALSE), fixed = TRUE))

## ---- the engines are named for their fitting functions ---------------------
chk("declaration: the brms engine is `brm`, as the others are `lm`, `clm`",
    identical(nesting_spec(dat, response ~ chord_type * inversion,
                           "inversion %in% chord_type", fit = "brm")$fit, "brm"))
chk("declaration: `brms`, the package's name, is accepted for it",
    identical(nesting_spec(dat, response ~ chord_type * inversion,
                           "inversion %in% chord_type", fit = "brms")$fit, "brm"))
chk("declaration: the emitted call is brms::brm either way",
    { cd <- fit_from(nesting_spec(dat, response ~ chord_type * inversion,
                     "inversion %in% chord_type", fit = "brms"), dry_run = TRUE)
      grepl("brms::brm(", cd, fixed = TRUE) })
chk("declaration: every engine name maps to a fitting function",
    all(vapply(c("lm", "glm", "lmer", "glmer", "clm", "clmm", "brm"),
               function(k) nzchar(engine_call(k)), TRUE)))

## ---- random effects: the engine's default, stated not overridden -----------
## Whether the group-level effects enter is the engine's decision. The package
## says which default is in force, since it settles whether the estimand
## describes a typical group or the sampled ones.
if (requireNamespace("lme4", quietly = TRUE)) {
  chk("re.form: nothing is injected into the call",
      { cd <- attr(suppressMessages(suppressWarnings(nest_estimand(mm2, chord_type,
                     type = "response", p_adjust = "none", bounds = FALSE))), "nestimand")$code
        call_lines <- grep("avg_predictions|hypothesis =", cd, value = TRUE)
        !any(grepl("re.form", call_lines, fixed = TRUE)) })
  chk("re.form: the default in force is announced",
      { z <- character(0)
        withCallingHandlers(suppressWarnings(nest_estimand(mm2, chord_type,
          type = "response", p_adjust = "none", bounds = FALSE)),
          message = function(m) { z <<- c(z, conditionMessage(m))
                                  invokeRestart("muffleMessage") })
        any(grepl("default stands", z)) && any(grepl("re.form = NA", z)) })
  chk("re.form: supplying one silences the note and reaches the engine",
      { z <- character(0)
        cd <- attr(withCallingHandlers(suppressWarnings(nest_estimand(mm2, chord_type,
                type = "response", re.form = NA, p_adjust = "none",
                bounds = FALSE)),
                message = function(m) { z <<- c(z, conditionMessage(m))
                                        invokeRestart("muffleMessage") }),
              "nestimand")$code
        !any(grepl("default stands", z)) &&
        any(grepl("re.form = NA",
                  grep("avg_predictions|hypothesis =", cd, value = TRUE),
                  fixed = TRUE)) })
  chk("re.form: the question does not arise on the linear predictor",
      { z <- character(0)
        withCallingHandlers(suppressWarnings(nest_estimand(mm2, chord_type,
          type = "eta", bounds = FALSE)),
          message = function(m) { z <<- c(z, conditionMessage(m))
                                  invokeRestart("muffleMessage") })
        !any(grepl("default stands", z)) })
}

## ---- what eta does about the random effects --------------------------------
## It takes the fixed effects, which is the typical group - the equivalent of
## re.form = NA. On a linear scale the sampled-group answer is the same number,
## since the group deviations average to zero.
if (requireNamespace("lme4", quietly = TRUE)) {
  g_re <- function(...) as.data.frame(suppressMessages(suppressWarnings(
    nest_estimand(mm2, chord_type, bounds = FALSE, ...))))$estimate
  chk("eta: equals the prediction route with the deviations held at zero",
      max(abs(g_re(type = "eta") -
              g_re(type = "response", p_adjust = "none", re.form = NA))) < 1e-8)
  chk("eta: and, on a linear scale, with them included too",
      max(abs(g_re(type = "eta") -
              g_re(type = "response", p_adjust = "none", re.form = NULL))) < 1e-8)
  chk("eta: because the group deviations average to zero",
      max(abs(colMeans(lme4::ranef(mm2)[[1]]))) < 1e-6)
  chk("re.form: the note says eta is the typical group, not the sampled ones",
      { z <- character(0)
        withCallingHandlers(suppressWarnings(nest_estimand(mm2, chord_type,
          type = "response", p_adjust = "none", bounds = FALSE)),
          message = function(m) { z <<- c(z, conditionMessage(m))
                                  invokeRestart("muffleMessage") })
        any(grepl("average group", z)) && any(grepl("at zero", z)) })
}

if (requireNamespace("lme4", quietly = TRUE)) {
  se_of <- function(rf) as.data.frame(suppressMessages(suppressWarnings(
    nest_estimand(mm2, chord_type, type = "response", p_adjust = "none",
             re.form = rf, bounds = FALSE))))$std.error
  ## marginaleffects forms the interval by the delta method from the
  ## fixed-effect covariance, so the two settings give the same standard error
  ## up to numerical differentiation - they differ at 1e-9, not at all in kind.
  chk("re.form: a frequentist fit gives the same interval either way",
      max(abs(se_of(NA) - se_of(NULL))) < 1e-6)
}

## ---- eta and the random-effects form ---------------------------------------
if (requireNamespace("lme4", quietly = TRUE)) {
  chk("eta: a random-effects form is refused, since eta is the average group",
      grepl("ask for different things",
            err_of(nest_estimand(mm2, chord_type, type = "eta", re.form = NULL))))
  chk("eta: on a frequentist fit the refusal says nothing is lost",
      grepl("not available on any route",
            err_of(nest_estimand(mm2, chord_type, type = "eta", re.form = NULL))))
}

## ---- eta over the sampled groups, where the draws allow it -----------------
## The group deviations average to zero across draws but not within one, and
## that is where the extra posterior width comes from. Tested on constructed
## draws, since a fitted brms model is not available here.
set.seed(1); nd_g <- 800
Dg <- cbind(b_cellA = rnorm(nd_g, 1), b_cellB = rnorm(nd_g, 2))
for (g in c("g1", "g2", "g3")) for (tm in c("cellA", "cellB"))
  Dg <- cbind(Dg, matrix(rnorm(nd_g, 0, 0.5), nd_g, 1,
              dimnames = list(NULL, sprintf("r_participant[%s,%s]", g, tm))))
Ug <- group_mean_draws(Dg, c("b_cellA", "b_cellB"))
chk("group draws: one column per coefficient, one row per draw",
    nrow(Ug) == nd_g && ncol(Ug) == 2)
chk("group draws: they average to zero across draws",
    max(abs(colMeans(Ug))) < 0.05)
chk("group draws: but vary within a draw, which is the point",
    all(apply(Ug, 2, stats::sd) > 0.1))
Cg <- matrix(c(-1, 1), 1, 2, dimnames = list("B - A", c("b_cellA", "b_cellB")))
s_avg <- Dg[, colnames(Cg)] %*% t(Cg)
s_grp <- (Dg[, colnames(Cg)] + Ug) %*% t(Cg)
chk("group draws: including them widens the posterior without moving it",
    abs(mean(s_avg) - mean(s_grp)) < 0.05 && stats::sd(s_grp) > stats::sd(s_avg))
chk("group draws: a coefficient with no group-level counterpart contributes nothing",
    all(group_mean_draws(Dg, "b_training") == 0))
chk("eta: a non-Bayesian fit still refuses the sampled-group request",
    grepl("ask for different things",
          err_of(nest_estimand(mf, chord_type, type = "eta", re.form = NULL))))

## ---- deeper and branching structures --------------------------------------
## Two things the one-level chain never exercises: a stratum key with more than
## one ancestor, and a parent holding two children. Both are places where a
## position in the family vector was mistaken for ancestry.
deep <- local({
  rows <- list()
  for (ct in c("aug", "dim", "min", "maj")) {
    invs <- if (ct == "aug") "none" else c("0", "1", "2")
    x1s  <- if (ct == "aug") "none" else c("a", "b")
    for (iv in invs) for (x in x1s) {
      zs <- if (iv %in% c("none", "2")) "none" else c("z1", "z2")
      for (z in zs)
        rows[[length(rows) + 1]] <-
          data.frame(chord_type = ct, inversion = iv, X1 = x, Z = z)
    }
  }
  g <- do.call(rbind, rows)
  d <- g[rep(seq_len(nrow(g)), each = 12), ]
  d$chord_type <- factor(d$chord_type, levels = c("aug", "dim", "min", "maj"))
  d$inversion  <- factor(d$inversion,  levels = c("none", "0", "1", "2"))
  d$X1 <- factor(d$X1, levels = c("none", "a", "b"))
  d$Z  <- factor(d$Z,  levels = c("none", "z1", "z2"))
  set.seed(11); d$response <- rnorm(nrow(d), 4)
  d
})
sp_d <- nesting_spec(deep, response ~ chord_type * inversion * X1 * Z,
                     c("inversion %in% chord_type", "X1 %in% chord_type",
                       "Z %in% inversion"))
m_d <- fit_from(sp_d)
chk("branching: a parent may hold two nested variables",
    identical(unname(sp_d$parent[c("inversion", "X1", "Z")]),
              c("chord_type", "chord_type", "inversion")))
chk("branching: c(a, b) %in% parent declares them together",
    identical(nesting_spec(deep, response ~ chord_type * inversion * X1,
                           c(inversion, X1) %in% chord_type)$parent[c("inversion", "X1")],
              c(inversion = "chord_type", X1 = "chord_type")))
chk("branching: a + b %in% parent is the same declaration",
    identical(nesting_spec(deep, response ~ chord_type * inversion * X1,
                           inversion + X1 %in% chord_type)$parent[c("inversion", "X1")],
              c(inversion = "chord_type", X1 = "chord_type")))
chk("branching: the left of %in% may be bracketed for precedence",
    identical(nesting_spec(deep, response ~ chord_type * inversion * X1,
                           (inversion + X1) %in% chord_type)$parent[c("inversion", "X1")],
              c(inversion = "chord_type", X1 = "chord_type")))
chk("branching: nested brackets and c() mix freely on the left",
    identical(nesting_spec(deep, response ~ chord_type * inversion * X1,
                           c((inversion), X1) %in% chord_type)$parent[c("inversion", "X1")],
              c(inversion = "chord_type", X1 = "chord_type")))
chk("branching: an interaction on the left is refused, and named",
    grepl("inversion:X1",
          err_of(nesting_spec(deep, response ~ chord_type * inversion * X1,
                              inversion:X1 %in% chord_type))))
chk("branching: a variable inside two parents is refused",
    grepl("more than one",
          err_of(nesting_spec(deep, response ~ chord_type * inversion * X1,
                              c("X1 %in% chord_type", "X1 %in% inversion")))))
chk("branching: chained %in% is refused with the one-level-per-entry remedy",
    grepl("one level per", err_of(nesting_spec(deep, response ~ chord_type,
                                               "Z %in% inversion %in% chord_type"))))
chk("ancestry: the strata of a sibling are its own parents, not its siblings",
    identical(nest_ancestors(sp_d, "X1"), "chord_type") &&
    identical(nest_ancestors(sp_d, "Z"), c("chord_type", "inversion")))
chk("ancestry: degenerate strata are keyed by every ancestor",
    identical(degenerate_strata(sp_d, "X1")$vars, "chord_type") &&
    identical(degenerate_strata(sp_d, "Z")$vars, c("chord_type", "inversion")))
chk("depth: the effect basis is square and full rank over the realized cells",
    { A <- effect_basis(sp_d)
      nrow(A) == nrow(sp_d$cells) && qr(A)$rank == nrow(A) })
chk("depth: the fit has one coefficient per realized cell, none aliased",
    length(coef(m_d)) == nrow(sp_d$cells) && !anyNA(coef(m_d)))
## The same structure without the sibling, small enough that the bounds are
## enumerated rather than skipped: three levels of nesting is the point here.
sp_z <- nesting_spec(subset(deep, X1 %in% c("none", "a")),
                     response ~ chord_type * inversion * Z,
                     c("inversion %in% chord_type", "Z %in% inversion"))
m_z <- fit_from(sp_z)
chk("depth: the pooled estimand of the deepest variable runs, with bounds",
    { e <- nest_estimand(m_z, Z, policy = "equal")
      b <- attr(e, "nestimand_bounds")
      !is.null(b) && b$policy_low < b$estimate && b$estimate < b$policy_high })
chk("depth: its emitted code runs on its own and gives the same estimate",
    { e <- nest_estimand(m_z, Z, policy = "equal")
      cd <- attr(nest_estimand(m_z, Z, policy = "equal", dry_run = TRUE),
                 "nestimand_code")
      en <- new.env(parent = globalenv())
      assign("sp_z", sp_z, en); assign("m_z", m_z, en)
      out <- eval(parse(text = paste(c(cd, "est"), collapse = "\n")), envir = en)
      isTRUE(all.equal(as.data.frame(out)$estimate, as.data.frame(e)$estimate)) })
chk("branching: contrasts of a sibling group by its own parent, not its sibling",
    { w <- as.data.frame(suppressMessages(nest_estimand(m_d, X1, by = "chord_type",
                                                   bounds = FALSE)))
      setequal(w$chord_type, c("dim", "min", "maj")) })
chk("branching: the declaration round-trips through spec_nests()",
    setequal(spec_nests(sp_d),
             c("inversion %in% chord_type", "X1 %in% chord_type",
               "Z %in% inversion")))
## Hierarchical weighting over a branching family: each variable's split is
## conditional on its own ancestors, so siblings are independent choices whose
## probabilities multiply, and no order of declaration is privileged.
h_a <- nesting_spec(deep, response ~ chord_type * inversion * X1 * Z,
                    c("inversion %in% chord_type", "X1 %in% chord_type",
                      "Z %in% inversion"))
h_b <- nesting_spec(deep, response ~ chord_type * inversion * X1 * Z,
                    c("X1 %in% chord_type", "inversion %in% chord_type",
                      "Z %in% inversion"))
chk("hierarchical: the split is conditional on ancestry, not on the order declared",
    { pa <- nest_policy(h_a, "chord_type", "hierarchical")$p[["maj"]]
      pb <- nest_policy(h_b, "chord_type", "hierarchical")$p[["maj"]]
      isTRUE(all.equal(sort(unname(pa)), sort(unname(pb)))) })
chk("hierarchical: over a set of siblings it agrees with equal",
    { sp_s <- nesting_spec(deep[deep$Z %in% c("none", "z1"), ],
                           response ~ chord_type * inversion * X1,
                           c(inversion, X1) %in% chord_type)
      p <- nest_policy(sp_s, "chord_type", "hierarchical")$p[["maj"]]
      max(abs(p - 1 / length(p))) < 1e-12 })
chk("hierarchical: a ragged support is renormalized over what exists",
    { d <- deep[!(deep$inversion != "0" & deep$X1 == "b"), ]   # X1 = b only at inversion 0
      sp_r <- nesting_spec(d[d$Z %in% c("none", "z1"), ],
                           response ~ chord_type * inversion * X1,
                           c(inversion, X1) %in% chord_type)
      p <- nest_policy(sp_r, "chord_type", "hierarchical")$p[["maj"]]
      abs(sum(p) - 1) < 1e-12 && max(abs(p - 1 / length(p))) < 1e-12 })
chk("hierarchical: it still separates from equal where one variable is inside another",
    { ph <- nest_policy(h_a, "chord_type", "hierarchical")$p[["maj"]]
      pe <- nest_policy(h_a, "chord_type", "equal")$p[["maj"]]
      max(abs(ph - pe)) > 1e-6 })

## A branching family on the random side: the design's own columns again, and
## the same ones whichever order the siblings were declared in.
re_dat <- local({
  d <- deep[deep$Z %in% c("none", "z1"), c("chord_type", "inversion", "X1")]
  d <- d[rep(seq_len(nrow(d)), each = 2), ]
  d$participant <- factor(rep(rep(1:8, each = 2), length.out = nrow(d)))
  set.seed(12); d$response <- rnorm(nrow(d)); d
})
re_f <- response ~ chord_type * inversion * X1 +
  (chord_type * inversion * X1 | participant)
re_a <- nesting_spec(re_dat, re_f, c(inversion, X1) %in% chord_type, fit = "lmer")
re_b <- nesting_spec(re_dat, re_f, c(X1, inversion) %in% chord_type, fit = "lmer")
chk("random branching: the bar gets the columns the design identifies",
    identical(random_terms(re_a),
              paste0("(1 + ", paste(reduced_columns(re_a), collapse = " + "),
                     " | participant)")))
chk("random branching: one column per realized condition, and no more",
    length(reduced_columns(re_a)) + 1L == nrow(re_a$cells))
chk("random branching: the order the siblings were declared in changes nothing",
    identical(sort(reduced_columns(re_a)), sort(reduced_columns(re_b))))

## The declared random terms are read from the parsed expression. A bar whose
## left side is bracketed - the natural way to write a slope over a set of
## variables - hid the bar from the regex that used to match parentheses, and
## the term arrived with no grouping factor.
chk("bars: a bracketed left side is still one term, with its grouping factor",
    { b <- bar_terms_of("(chord_type * (inversion + X1) | participant)")
      length(b) == 1 && b[[1]]$grp == "participant" &&
        b[[1]]$wrapper == "" && grepl("inversion", b[[1]]$lhs) })
chk("bars: several terms are separated, however they are bracketed",
    { b <- bar_terms_of("(1 | participant) + ((a + b) | stimulus)")
      length(b) == 2 &&
        all(vapply(b, `[[`, "", "grp") %in% c("participant", "stimulus")) })
chk("bars: a covariance wrapper is kept, not read as the grouping factor",
    { b <- bar_terms_of("diag(a * b | participant)")
      length(b) == 1 && b[[1]]$wrapper == "diag" && b[[1]]$grp == "participant" })
chk("bars: a nested grouping factor is read whole",
    identical(bar_terms_of("(1 | school:class)")[[1]]$grp, "school:class"))
chk("bars: the translation is the same however the left side is bracketed",
    { mk <- function(bar) suppressMessages(nesting_spec(re_dat,
        stats::as.formula(paste("response ~ chord_type * inversion * X1 +", bar)),
        "inversion %in% chord_type", fit = "lmer"))
      identical(random_terms(mk("(chord_type * (inversion + X1) | participant)")),
                random_terms(mk("(chord_type * (inversion + X1) | participant)"))) &&
      grepl("^\\(1 \\+ dm_",
            random_terms(mk("(chord_type * (inversion + X1) | participant)"))) })
chk("bars: grouping_vars reads them from the parse too",
    identical(grouping_vars(suppressMessages(nesting_spec(re_dat,
        response ~ chord_type * inversion * X1 +
          (chord_type * (inversion + X1) | participant),
        "inversion %in% chord_type", fit = "lmer"))), "participant"))


## A categorical variable of the design that is nested in nothing. It belongs
## in the cell factor all the same - a variable left out of it is a covariate,
## which cannot be the target of an estimand.
cross_dat <- local({
  rows <- list()
  for (ct in c("aug", "dim", "min", "maj")) {
    invs <- if (ct == "aug") "none" else c("0", "1", "2")
    for (iv in invs) for (tp in c("t1", "t2"))
      rows[[length(rows) + 1]] <-
        data.frame(chord_type = ct, inversion = iv, top = tp)
  }
  g <- do.call(rbind, rows); d <- g[rep(seq_len(nrow(g)), each = 12), ]
  d$chord_type <- factor(d$chord_type, levels = c("aug", "dim", "min", "maj"))
  d$inversion  <- factor(d$inversion, levels = c("none", "0", "1", "2"))
  d$top <- factor(d$top)
  set.seed(15); d$response <- rnorm(nrow(d), 4); d
})
sp_x <- nesting_spec(cross_dat, response ~ chord_type * inversion * top,
                     c("inversion %in% chord_type", "top"))
m_x <- fit_from(sp_x)
chk("crossed: a bare name declares a variable nested in nothing",
    identical(sp_x$crossed, "top") &&
      identical(sp_x$cell_vars, c("chord_type", "inversion", "top")))
chk("crossed: the same declaration unquoted",
    identical(nesting_spec(cross_dat, response ~ chord_type * inversion * top,
                           c(inversion %in% chord_type, top))$cell_vars,
              c("chord_type", "inversion", "top")))
chk("crossed: it joins the cell factor, which stays full rank",
    { A <- effect_basis(sp_x)
      nrow(A) == 20 && qr(A)$rank == 20 && !anyNA(coef(m_x)) })
chk("crossed: and is a target like any other",
    { e <- as.data.frame(nest_estimand(m_x, top, policy = "equal", bounds = FALSE))
      identical(e$term, "t2 - t1") })
chk("crossed: a numeric one is refused, since that is a covariate",
    grepl("is a covariate",
          err_of(nesting_spec(within(cross_dat, top <- as.numeric(top)),
                              response ~ chord_type * inversion * top,
                              c("inversion %in% chord_type", "top")))))
chk("crossed: declaring it both ways is refused",
    grepl("one position in the structure",
          err_of(nesting_spec(cross_dat, response ~ chord_type * inversion * top,
                              c("top %in% chord_type", "top")))))
chk("crossed: a factor crossed with the structure needs no declaration",
    { sp <- suppressMessages(nesting_spec(cross_dat,
              response ~ chord_type * inversion * top, "inversion %in% chord_type"))
      identical(sp$cell_vars, c("chord_type", "inversion", "top")) &&
        identical(sp$crossed, "top") })
chk("crossed: and the fold is said aloud",
    grepl("part of the categorical design",
          paste(utils::capture.output(nesting_spec(cross_dat,
            response ~ chord_type * inversion * top, "inversion %in% chord_type"),
            type = "message"), collapse = " ")))
chk("crossed: one entering additively is left as a covariate",
    { sp <- nesting_spec(cross_dat, response ~ chord_type * inversion + top,
                         "inversion %in% chord_type")
      identical(sp$cell_vars, c("chord_type", "inversion")) &&
        identical(sp$covariates, "top") })
sp_undeclared <- nesting_spec(cross_dat, response ~ chord_type * inversion + top,
                              "inversion %in% chord_type")
m_undeclared <- fit_from(sp_undeclared)
chk("target: an additive factor is told how to become a target",
    { e <- err_of(nest_estimand(m_undeclared, top, policy = "equal"))
      grepl("additively", e, fixed = TRUE) })

## `a * b * c` parses as `(a * b) * c`, so the operands have to be gathered
## through the nesting: reading the top call alone left `a * b` as a name.
chk("target: a * b * c names three targets, not two",
    { e <- err_of(nest_estimand(mf, chord_type * inversion * training,
                           policy = "equal", bounds = FALSE))
      grepl("`training` is not among", e, fixed = TRUE) })
chk("target: a * b still names two",
    inherits(nest_estimand(mf, chord_type * inversion, policy = "equal",
                      bounds = FALSE), "nestimand_estimands"))

## A factor covariate crossed with the cells is named `cell<k>:x<level>`, one
## set of columns per non-reference level. Matching the variable name exactly
## found none of them, and every column fell through to the leftover block,
## where it was labelled a threshold.
fac_dat <- local({ d <- dat; set.seed(14)
  ## ordered, so it stays a covariate: its contrasts say it is meant as a
  ## quantity, and it is crossed with the cells rather than folded into them
  d$top <- factor(rep(c("t1", "t2", "t3"), length.out = nrow(d)), ordered = TRUE); d })
sp_fac <- nesting_spec(fac_dat, response ~ chord_type * inversion * top * training,
                       "inversion %in% chord_type")
m_fac <- fit_from(sp_fac)
s_fac <- as.data.frame(nest_summary(m_fac))
chk("summary: a factor covariate crossed with the cells is translated",
    sum(grepl("slope on top", s_fac$meaning)) == 2 * nrow(sp_fac$cells))
chk("summary: and nothing of it is left over as a threshold",
    !any(grepl("threshold", s_fac$meaning)))
chk("summary: a real threshold is still called one",
    any(grepl("threshold", as.data.frame(nest_summary(mo))$meaning)))
chk("summary: every row still translates one coefficient",
    nrow(s_fac) == length(coef(m_fac)))

## A Bayesian fit is summarized from its draws: a p-value computed from a
## normal approximation to the posterior answers a question the model was not
## fitted to ask. The mapping is exercised here; the brms fit itself is in
## tests/test_brms.R.
set.seed(21)
S_draws <- cbind(up = rnorm(4000, 1, 0.3), flat = rnorm(4000, 0, 1),
                 held = rep(0, 4000))
ds <- draws_summary(S_draws, colnames(S_draws), 0.9)
chk("nest_summary posterior: the summary is the posterior's own, not an approximation",
    { isTRUE(all.equal(unname(ds$estimate), unname(colMeans(S_draws)))) &&
      isTRUE(all.equal(unname(ds$std.error), unname(apply(S_draws, 2, sd)))) })
chk("nest_summary posterior: the interval is a quantile interval at the level asked for",
    isTRUE(all.equal(c(ds$conf.low[1], ds$conf.high[1]),
                     unname(quantile(S_draws[, "up"], c(0.05, 0.95))))))
chk("nest_summary posterior: pd is the mass on the side of zero holding more of it",
    { isTRUE(all.equal(ds$pd[1], mean(S_draws[, "up"] > 0))) &&
      ds$pd[2] > 0.5 && ds$pd[2] < 0.6 })
chk("nest_summary posterior: a row held at zero has no direction to report",
    is.na(ds$pd[3]))
chk("nest_summary posterior: no p-value or test statistic is offered",
    !any(c("p.value", "statistic") %in% names(ds)))
chk("nest_summary posterior: the printed table shows pd where a frequentist shows p",
    { fake <- ds; class(fake) <- c("nestimand_summary", class(fake))
      attr(fake, "nestimand_space") <- "effects"
      attr(fake, "nestimand_fit") <- "brm"
      fake$meaning <- c("a", "b", "c")
      txt <- paste(utils::capture.output(print(fake)), collapse = " ")
      grepl(" pd ", txt) && !grepl("p.value", txt, fixed = TRUE) })
chk("nest_summary posterior: a frequentist fit still reports the p-value",
    { d <- as.data.frame(nest_summary(mf))
      "p.value" %in% names(d) && !("pd" %in% names(d)) })

## An interaction contrast over three variables: the corners are the eight
## conditions two levels of each variable define, and the sign is the product
## of the three simple contrasts. The matrix was written for two variables and
## silently matched nothing when given three.
sp_i3 <- suppressMessages(nesting_spec(cross_dat,
  response ~ chord_type * inversion * top, "inversion %in% chord_type"))
m_i3 <- fit_from(sp_i3)
chk("interaction: three variables give a difference of differences of differences",
    { H <- interaction_matrix(sp_i3$cells,
                              c("chord_type", "inversion", "top"))
      ncol(H) == 9 && all(colSums(H) == 0) && all(colSums(abs(H)) == 8) })
chk("interaction: the sentinel stratum is left out, having no two inversions",
    { H <- interaction_matrix(sp_i3$cells, c("chord_type", "inversion", "top"))
      !any(grepl("aug", colnames(H))) })
chk("interaction: `a * b * c` crosses, as a formula does",
    { e <- nest_estimand(m_i3, chord_type * inversion * top, policy = "equal",
                    bounds = FALSE)
      identical(names(e), c("chord_type", "inversion", "top",
                            "chord_type:inversion", "chord_type:top",
                            "inversion:top", "chord_type:inversion:top")) })
chk("interaction: the sentinel never enters, whichever target is named last",
    { e <- as.data.frame(nest_estimand(m_i3, inversion:top, type = "eta",
                                  bounds = FALSE))
      nrow(e) == 3 && !any(grepl("none", e$term)) })
chk("interaction: the restriction is the targets' together, not the last one's",
    { d <- degenerate_strata_multi(sp_i3, c("inversion", "top"))
      !any(grepl("aug", d$keep)) && any(grepl("aug", d$drop)) })
chk("interaction: a stratum keeps its place where every target does vary in it",
    { e <- as.data.frame(nest_estimand(m_i3, chord_type:top, type = "eta",
                                  bounds = FALSE))
      nrow(e) == 6 && any(grepl("aug", e$term)) })
chk("interaction: an interaction over some of the design variables averages over the rest",
    { ## the target names two of the three cell variables, so each combination
      ## covers two cells; the contrast is formed on their average, and taking
      ## one of them would answer at an arbitrary level of the third
      e <- as.data.frame(nest_estimand(m_i3, chord_type:inversion, type = "eta",
                                  bounds = FALSE))
      b <- coef(m_i3)
      av <- function(a, i) mean(c(b[[paste0("cell", a, ".", i, ".t1")]],
                                  b[[paste0("cell", a, ".", i, ".t2")]]))
      hand <- (av("min", "1") - av("min", "0")) - (av("dim", "1") - av("dim", "0"))
      isTRUE(all.equal(e$estimate[e$term == "(min - dim) x (1 - 0)"], hand)) })
chk("interaction: the two routes agree on it, as they must on a linear scale",
    { a <- as.data.frame(nest_estimand(m_i3, chord_type:inversion, type = "eta",
                                  bounds = FALSE))
      b <- as.data.frame(nest_estimand(m_i3, chord_type:inversion,
                                  bounds = FALSE))
      isTRUE(all.equal(a$estimate[order(a$term)], b$estimate[order(b$term)])) })
chk("interaction: a design with no such set of corners says so, and why",
    { ## every variable has two levels, but no eight conditions carry two of
      ## each: the design cannot answer a three-way comparison
      d <- data.frame(chord_type = c("dim", "dim", "min", "min"),
                      inversion  = c("0", "1", "0", "1"),
                      top        = c("t1", "t2", "t2", "t1"))
      grepl("realizes no set of 8 conditions",
            err_of(interaction_matrix(d, c("chord_type", "inversion", "top")))) })
chk("interaction: a variable with one level here is named rather than pivoted on",
    grepl("has one here",
          err_of(interaction_matrix(sp_i3$cells[sp_i3$cells$chord_type == "dim", ],
                                    c("chord_type", "inversion", "top")))))

## A variable given a random slope keeps its fixed effect, and a covariate
## declared crossed with the structure keeps that crossing on both sides.
slope_dat <- local({
  d <- deep[deep$Z %in% c("none", "z1"), ]
  d$participant <- factor(rep(1:8, length.out = nrow(d)))
  set.seed(13); d$GMSI <- rnorm(nrow(d)); d
})
sp_sl <- suppressMessages(nesting_spec(slope_dat,
  response ~ chord_type * inversion * X1 * GMSI +
    (chord_type * inversion * X1 * GMSI | participant),
  "inversion %in% chord_type", fit = "lmer"))
chk("covariates: a variable with a random slope keeps its fixed effect",
    "GMSI" %in% sp_sl$covariates)
chk("covariates: the grouping factor is not one of them",
    !("participant" %in% sp_sl$covariates))
chk("covariates: crossed with the structure on the fixed side",
    grepl("cell:GMSI", paste(deparse(cell_formula(sp_sl)), collapse = " "),
          fixed = TRUE))
chk("random: a covariate crossed with the structure keeps the crossing",
    { cols <- reduced_columns(sp_sl)
      identical(random_terms(sp_sl),
        paste0("(1 + ", paste(cols, collapse = " + "), " + GMSI + ",
               paste(paste0(cols, ":GMSI"), collapse = " + "),
               " | participant)")) })
chk("random: one left additive in the bar stays additive",
    { sp_ad <- suppressMessages(nesting_spec(slope_dat,
        response ~ chord_type * inversion * X1 * GMSI +
          (chord_type * inversion * X1 + GMSI | participant),
        "inversion %in% chord_type", fit = "lmer"))
      identical(random_terms(sp_ad),
        paste0("(1 + ", paste(reduced_columns(sp_ad), collapse = " + "),
               " + GMSI | participant)")) })

## A numeric nested variable is a slope inside the cells, not a cell variable.
## It stays visible in the printed structure, since a variable whose values are
## labels is easily left numeric by accident.
cont_dat <- local({
  d <- deep[deep$Z %in% c("none", "z1"), ]
  d$X1n <- as.numeric(factor(d$X1)); d })
chk("continuous nested: the declaration is said aloud rather than absorbed",
    grepl("slope within the realized cells",
          paste(utils::capture.output(
            nesting_spec(cont_dat, response ~ chord_type * inversion * X1n,
                         c(inversion, X1n) %in% chord_type),
            type = "message"), collapse = " ")))
chk("continuous nested: it is a covariate crossed with the cells, not a cell variable",
    { sp_c <- suppressMessages(
        nesting_spec(cont_dat, response ~ chord_type * inversion * X1n,
                     c(inversion, X1n) %in% chord_type))
      identical(sp_c$cell_vars, c("chord_type", "inversion")) &&
        identical(unname(sp_c$cont_nested), "X1n") &&
        grepl("cell:X1n", paste(deparse(cell_formula(sp_c)), collapse = " ")) })
chk("continuous nested: the printed structure still shows it",
    grepl("X1n (continuous)",
          paste(utils::capture.output(print(suppressMessages(
            nesting_spec(cont_dat, response ~ chord_type * inversion * X1n,
                         c(inversion, X1n) %in% chord_type)))), collapse = " "),
          fixed = TRUE))

## A declaration that asks for less than the saturated structure is fitted as
## written. The cell factor is saturated by construction, so the restriction can
## only be expressed in the effects parameterization, and the fitting mode
## follows the declaration rather than the default.
res_sp <- suppressMessages(nesting_spec(cross_dat,
  response ~ chord_type * (inversion + top), "inversion %in% chord_type"))
res_m <- fit_from(res_sp)
chk("restricted: the declared terms are the formula's own, nothing added",
    identical(declared_terms(res_sp),
              c("chord_type", "inversion", "top",
                "chord_type:inversion", "chord_type:top")))
chk("restricted: a `+` between two variables is not upgraded to a `*`",
    { sp <- suppressMessages(nesting_spec(cross_dat,
        response ~ chord_type + inversion, "inversion %in% chord_type"))
      ## the pooled model: one inversion effect, shared by every chord type
      identical(declared_terms(sp), c("chord_type", "inversion")) &&
        term_span(sp, declared_terms(sp)) == 6 &&
        identical(as.character(fitting_mode(sp)), "reduced") })
chk("restricted: and the count reported is the model that is fitted",
    { sp <- suppressMessages(nesting_spec(cross_dat,
        response ~ chord_type + inversion, "inversion %in% chord_type"))
      m <- suppressMessages(fit_from(sp))
      length(coef(m)) == term_span(sp, declared_terms(sp)) && !anyNA(coef(m)) })
chk("restricted: the sentinel is absence, so the contrasts are among real levels",
    { sp <- suppressMessages(nesting_spec(cross_dat,
        response ~ chord_type + inversion, "inversion %in% chord_type"))
      nm <- unname(attr(reduced_design(sp), "effect_names"))
      identical(nm, c("(Intercept)", "chord_typedim", "chord_typemin",
                      "chord_typemaj", "inversion1", "inversion2")) })
chk("restricted: a row where the variable is undefined contributes nothing",
    { sp <- suppressMessages(nesting_spec(cross_dat,
        response ~ chord_type + inversion, "inversion %in% chord_type"))
      X <- reduced_design(sp)
      aug <- grepl("^aug", rownames(X))
      all(X[aug, grepl("inversion", colnames(X))] == 0) })
chk("restricted: they span less than the cells, and the mode follows",
    { sp_full <- suppressMessages(nesting_spec(cross_dat,
        response ~ chord_type * inversion * top, "inversion %in% chord_type"))
      term_span(res_sp, declared_terms(res_sp)) < nrow(res_sp$cells) &&
        identical(as.character(fitting_mode(res_sp)), "reduced") &&
        identical(as.character(fitting_mode(sp_full)), "cells") })
chk("restricted: the fit is the model the formula names, not the saturated one",
    { hand <- lm(response ~ chord_type + top + chord_type:inversion +
                   chord_type:top, data = cross_dat)
      isTRUE(all.equal(unname(fitted(res_m)), unname(fitted(hand)))) &&
        df.residual(res_m) == df.residual(hand) })
chk("restricted: the aliased columns the chain form carries are dropped, not fitted",
    { b <- coef(res_m)
      sum(!is.na(b)) == term_span(res_sp, declared_terms(res_sp)) })
chk("restricted: estimands run, and the two routes agree",
    { a <- as.data.frame(nest_estimand(res_m, chord_type, policy = "equal",
                                  bounds = FALSE))
      b <- as.data.frame(nest_estimand(res_m, chord_type, type = "eta",
                                  policy = "equal", bounds = FALSE))
      nrow(a) == 6 && isTRUE(all.equal(a$estimate, b$estimate)) })
chk("restricted: nest_summary reports the coefficients rather than refusing",
    { d <- as.data.frame(nest_summary(res_m))
      nrow(d) == sum(!is.na(coef(res_m))) &&
        all(d$meaning == "as fitted" | grepl("held at zero", d$meaning)) })
chk("restricted: the printed spec shows the structure that will be fitted",
    { txt <- paste(utils::capture.output(print(res_sp)), collapse = " ")
      grepl("Structure fitted", txt, fixed = TRUE) &&
        grepl("chord_type + inversion + top", txt, fixed = TRUE) &&
        grepl("14 of the 20 realized cells", txt, fixed = TRUE) &&
        !grepl("0 \\+ cell", txt) })

## brms cannot drop an uninformative column the way lm does, so the effects
## parameterization needs its constant(0) block or the posterior is improper
## along it. That form is now reached only when something asks for it - the
## emmeans engine, or a prior stated coordinate-wise in effect space - since a
## restricted declaration goes to the reduced form, which carries no such
## column. The block is derived and applied rather than left to be discovered.
res_brm <- suppressMessages(nesting_spec(cross_dat,
  response ~ chord_type * (inversion + top), "inversion %in% chord_type",
  fit = "brm"))
## The message and the fit are one number now: term_span() is the width of the
## design fit_from() builds, so the count reported cannot drift from the model
## fitted. It used to, and said "fitted as written" while fitting the saturated
## structure.
chk("saturation: the message counts the realized cells the formula spans",
    { txt <- paste(utils::capture.output(
        nesting_spec(cross_dat, response ~ chord_type * (inversion + top),
                     "inversion %in% chord_type"), type = "message"),
        collapse = " ")
      grepl("spans 14 of the 20 realized conditions", txt, fixed = TRUE) &&
        grepl("one coefficient per effect it names", txt, fixed = TRUE) })
chk("saturation: and the count is the number of coefficients fitted",
    { sp <- suppressMessages(nesting_spec(cross_dat,
        response ~ chord_type * (inversion + top), "inversion %in% chord_type"))
      m <- suppressMessages(fit_from(sp))
      length(coef(m)) == term_span(sp, declared_terms(sp)) && !anyNA(coef(m)) })


## The declaration is carried by the fit, so these functions take the model first
## and the target second. A spec passed positionally lands in `target` and used
## to fail several frames later on `if (target %in% f)`, which mentions neither
## argument. Every route to a policy runs through versions_of(), so the check
## sits there and each entry point inherits it.
chk("target: a spec in the target slot is named as one",
    { e <- err_of(nest_policy(res_sp, res_sp))
      grepl("nesting_spec was passed in its place", e, fixed = TRUE) &&
        grepl("carries it", e, fixed = TRUE) })
chk("target: and reported before the engine is looked at",
    { e <- err_of(latent_draws(res_m, res_sp))
      grepl("nesting_spec was passed in its place", e, fixed = TRUE) })
chk("target: several variables are refused with the function that takes them",
    { e <- err_of(nest_policy(res_sp, c("chord_type", "top")))
      grepl("of length 2", e, fixed = TRUE) && grepl("nest_estimand()", e, fixed = TRUE) })
chk("target: a non-string says what it was",
    grepl("numeric of length 1", err_of(nest_policy(res_sp, 3)), fixed = TRUE))
## `latent_draws()` left the cells unresolved and passed NULL on, where
## `latent_estimand()` dropped the strata in which the target does not vary. One
## function decides it now, so the two routes are over the same cells.
chk("cells: one function decides which cells an estimand is over",
    { a <- estimand_cells(res_sp, "top")
      b <- estimand_cells(res_sp, "inversion")
      nrow(a) == nrow(res_sp$cells) && nrow(b) < nrow(res_sp$cells) &&
        identical(estimand_cells(res_sp, "top", res_sp$cells[1:2, ]),
                  res_sp$cells[1:2, ]) })

## ---- the random side follows the declaration too -------------------------
## The mean structure stopped being silently saturated when the reduced design
## arrived; the random side went on being saturated for a while longer. A
## declaration asking for less than the full crossing asks for less on both
## sides, and the cell factor can express that on neither.
re_red_d <- local({
  rows <- list()
  for (ct in c("aug", "dim", "min", "maj")) {
    invs <- if (ct == "aug") "none" else c("0", "1", "2")
    for (iv in invs) for (tp in c("t1", "t2"))
      rows[[length(rows) + 1]] <-
        data.frame(chord_type = ct, inversion = iv, top = tp)
  }
  g <- do.call(rbind, rows); d <- g[rep(seq_len(nrow(g)), each = 6), ]
  d$participant <- factor(rep(1:6, times = nrow(g)))
  d$chord_type <- factor(d$chord_type, levels = c("aug", "dim", "min", "maj"))
  d$inversion <- factor(d$inversion, levels = c("none", "0", "1", "2"))
  d$top <- factor(d$top)
  set.seed(21); d$response <- rnorm(nrow(d), 4); d$GMSI <- rnorm(nrow(d)); d
})
re_red_f <- response ~ chord_type * (inversion + top) +
  (chord_type * (inversion + top) | participant)
re_red <- suppressMessages(nesting_spec(re_red_d, re_red_f,
  "inversion %in% chord_type", fit = "lmer"))
re_sat <- suppressMessages(nesting_spec(re_red_d,
  response ~ chord_type * inversion * top +
    (chord_type * inversion * top | participant),
  "inversion %in% chord_type", fit = "lmer"))
chk("random reduced: the same columns the mean structure is fitted on",
    { r <- random_terms(re_red)
      cols <- setdiff(colnames(reduced_design(re_red)), "(Intercept)")
      identical(r, sprintf("(1 + %s | participant)",
                           paste(cols, collapse = " + "))) })
chk("random reduced: so the two sides have the same dimension",
    { r <- random_terms(re_red)
      n_re <- length(strsplit(sub("^\\(", "", sub(" \\|.*$", "", r)), " \\+ ")[[1]])
      n_fx <- ncol(reduced_design(re_red))
      n_re == n_fx && n_fx == 14 })
chk("random reduced: which is smaller than one dimension per condition",
    ncol(reduced_design(re_red)) == 14 && nrow(re_red$cells) == 20)
chk("random reduced: and it is what the fit uses",
    { cd <- suppressMessages(attr(fit_from(re_red, dry_run = TRUE),
                                  "nestimand_code"))
      any(grepl("(1 + dm_chord_typedim", cd, fixed = TRUE)) &&
        !any(grepl("0 + cell | participant", cd, fixed = TRUE)) })
chk("random reduced: a saturated declaration gets one column per condition",
    { r <- random_terms(re_sat)
      length(strsplit(sub("^\\(", "", sub(" \\|.*$", "", r)), " \\+ ")[[1]]) ==
        nrow(re_sat$cells) })
chk("random reduced: a term the mean structure lacks is refused, and says why",
    { sp_x <- suppressMessages(nesting_spec(re_red_d,
        response ~ chord_type * (inversion + top) +
          (chord_type * inversion * top | participant),
        "inversion %in% chord_type", fit = "lmer"))
      e <- err_of(random_terms(sp_x))
      grepl("which the model formula does not contain", e, fixed = TRUE) &&
        grepl("Add the term to the model formula", e, fixed = TRUE) })
chk("random reduced: a crossed covariate gets a slope per column and one more",
    { sp_x <- suppressMessages(nesting_spec(re_red_d,
        response ~ chord_type * (inversion + top) * GMSI +
          (chord_type * (inversion + top) * GMSI | participant),
        "inversion %in% chord_type", fit = "lmer"))
      r <- random_terms(sp_x)
      cols <- setdiff(colnames(reduced_design(sp_x)), "(Intercept)")
      grepl(" GMSI + ", r, fixed = TRUE) &&
        all(vapply(paste0(cols, ":GMSI"), function(z)
          grepl(z, r, fixed = TRUE), TRUE)) })
chk("random reduced: the emitted note says what it did and why",
    { cd <- suppressMessages(attr(fit_from(re_red, dry_run = TRUE),
                                  "nestimand_code"))
      txt <- paste(cd, collapse = " ")
      grepl("random structure", txt, fixed = TRUE) &&
        grepl("translated to", txt, fixed = TRUE) &&
        grepl("carries columns the", txt, fixed = TRUE) })
chk("random reduced: the covariance is reported in the effects, not the columns",
    { m_rr <- suppressMessages(fit_from(re_red))
      r <- random_covariance(m_rr, space = "cells")
      nm <- rownames(r$participant)
      !any(grepl("^dm_", nm)) })

## The shadow is fitted here rather than by nest_fit(), so it carries no
## parameterization label unless one is put on it - and an unlabelled fit is
## read as the cell form. For a restricted declaration the shadow is not in the
## cell form, so the check built design rows over conditions whose coefficients
## the shadow does not have and reported the model as one nest_fit() had not
## produced. Reported as `reorder check: error`, on an estimand that was itself
## perfectly sound.

## ---- contrast = "none": the predictions the contrasts are formed from -------
## Each row of the policy matrix is already the policy-weighted average of the
## design rows for one level of the target. The contrasts are differences
## between its rows, so asking for none is asking for that matrix itself - and
## the marginal mean of a partially nested design is policy-dependent in
## exactly the way a contrast is, which is the argument for taking both from
## one call rather than from two that could drift apart.
pred_e <- suppressMessages(nest_estimand(mf, chord_type, policy = "equal",
                                         contrast = "none"))
chk("none: one row per level of the target, named by the level",
    { d <- as.data.frame(pred_e)
      nrow(d) == 4 && identical(d$term, c("aug", "dim", "min", "maj")) })
chk("none: the predictions are what the contrasts are differences of",
    { d <- as.data.frame(pred_e)
      cd <- as.data.frame(suppressMessages(nest_estimand(mf, chord_type,
              policy = "equal", bounds = FALSE)))
      v <- function(z, t) z$estimate[z$term == t]
      abs((v(d, "maj") - v(d, "aug")) - v(cd, "maj - aug")) < 1e-10 })
chk("none: the two routes agree, as they must on an identity link",
    { a <- as.data.frame(suppressMessages(nest_estimand(mf, chord_type,
             contrast = "none", type = "eta", bounds = FALSE)))
      b <- as.data.frame(pred_e)
      max(abs(a$estimate - b$estimate[match(a$term, b$term)])) < 1e-8 })
chk("none: a level with a single version has bounds that collapse to a point",
    { b <- attr(pred_e, "nestimand")$bounds
      abs(b$policy_low[b$term == "aug"] - b$policy_high[b$term == "aug"]) < 1e-10 &&
        b$policy_high[b$term == "maj"] > b$policy_low[b$term == "maj"] })
chk("none: the emitted code re-runs and reproduces it",
    { p2 <- suppressMessages(nest_estimand(mf, chord_type, contrast = "none",
              bounds = FALSE))
      env2 <- new.env(parent = globalenv())
      assign("mf", mf, env2); assign("spec", attr(mf, "nestimand_spec"), env2)
      r <- suppressMessages(eval(parse(text = paste(attr(p2, "nestimand")$code,
                                                    collapse = "\n")), env2))
      isTRUE(all.equal(as.data.frame(r)$estimate, as.data.frame(p2)$estimate)) })
chk("none: a nested target reports the levels that exist, not the sentinel",
    { d <- as.data.frame(suppressMessages(nest_estimand(mf, inversion,
             contrast = "none", bounds = FALSE)))
      nrow(d) == 3 && !any(grepl("none", d$term)) })
chk("none: a target naming two variables gives one row per combination",
    { d <- as.data.frame(suppressMessages(nest_estimand(mf, chord_type:inversion,
             contrast = "none", bounds = FALSE)))
      nrow(d) == 9 && all(grepl(", ", d$term, fixed = TRUE)) })
chk("none: and the combinations exclude the stratum where one does not vary",
    { d <- as.data.frame(suppressMessages(nest_estimand(mf, chord_type:inversion,
             contrast = "none", bounds = FALSE)))
      !any(grepl("aug", d$term)) })
chk("none: the footer says what came back rather than `none`",
    grepl("no contrast: policy-weighted predictions",
          paste(utils::capture.output(print(pred_e)), collapse = " "), fixed = TRUE))

## The engine's own result object is the handle its own accessors work through -
## get_draws() above all, which is how a posterior is taken elsewhere to be
## processed or plotted. Canonicalizing the contrast directions used to coerce
## it to a plain data frame first, so a Bayesian estimand came back with no
## draws reachable. Keeping the object means keeping the draws in step with it:
## a contrast written the other way round has its estimate negated, and draws
## that disagreed with the estimate printed beside them would be worse than none.
mfx_obj <- local({
  d <- data.frame(term = c("aug - maj", "maj - dim"), estimate = c(-1, 2),
                  conf.low = c(-1.5, 1.5), conf.high = c(-0.5, 2.5))
  attr(d, "posterior_draws") <- matrix(c(-1.1, -0.9, 2.1, 1.9), nrow = 2,
                                       byrow = TRUE)
  class(d) <- c("predictions", "data.frame")
  mfx_canonical(d, levs = c("aug", "dim", "min", "maj"))
})
chk("draws: the engine's own result object survives canonicalization",
    inherits(mfx_obj, "predictions"))
chk("draws: a flipped contrast takes its draws with it",
    isTRUE(all.equal(rowMeans(attr(mfx_obj, "posterior_draws")),
                     mfx_obj$estimate)))
chk("draws: and the label is the flipped one",
    identical(mfx_obj$term, c("maj - aug", "maj - dim")))
chk("draws: draws that cannot be aligned give up the object, not the agreement",
    { d <- data.frame(term = c("aug - maj", "maj - dim"), estimate = c(-1, 2))
      attr(d, "posterior_draws") <- matrix(1:3, nrow = 3)   # wrong length
      class(d) <- c("predictions", "data.frame")
      o <- mfx_canonical(d, levs = c("aug", "dim", "min", "maj"))
      identical(class(o), "data.frame") && identical(o$estimate, c(1, 2)) })
chk("draws: a fit whose route runs through the engine keeps the class",
    { d_b <- dat; d_b$bin <- as.integer(d_b$response > stats::median(d_b$response))
      m_b <- suppressMessages(nest_fit("glm", bin ~ chord_type * inversion, d_b,
               family = "binomial", nests = "inversion %in% chord_type"))
      inherits(suppressMessages(nest_estimand(m_b, chord_type, policy = "equal",
               bounds = FALSE)), "predictions") })
## One accessor for the posterior of the estimand, whichever route produced it.
## Reading the engine's own draws through the engine's own accessor would work
## on one version of marginaleffects and not the next: the function was renamed
## and the object grew properties where it had attributes. The draws are
## snapshotted instead, in one orientation, under one name.
draws_obj <- local({
  o <- data.frame(term = c("maj - aug", "maj - dim"), estimate = c(1, 2))
  attr(o, "posterior_draws") <- matrix(c(0.9, 1.1, 1.8, 2.2), nrow = 2,
                                       byrow = TRUE)
  class(o) <- c("predictions", "data.frame")
  attr(o, "nestimand_draws") <-
    name_draw_columns(estimand_draws_from_engine(o), o)
  o
})
chk("draws: the engine's draws are snapshotted as draws x contrasts",
    identical(dim(attr(draws_obj, "nestimand_draws")), c(2L, 2L)) &&
      identical(colnames(attr(draws_obj, "nestimand_draws")),
                c("maj - aug", "maj - dim")))
chk("draws: long is one row per draw per contrast",
    { d <- nest_draws(draws_obj)
      is.data.frame(d) && identical(names(d), c("drawid", "term", "draw")) &&
        nrow(d) == 4 && is.numeric(d$draw) && is.character(d$term) })
chk("draws: on a collection the target is its own column, not glued into term",
    { coll <- structure(list(a = draws_obj, b = draws_obj),
                        class = "nestimand_estimands")
      d <- nest_draws(coll)
      is.data.frame(d) &&
        identical(names(d), c("drawid", "target", "term", "draw")) &&
        setequal(unique(d$target), c("a", "b")) &&
        setequal(unique(d$term), c("maj - aug", "maj - dim")) })
chk("draws: a stratum is a column too, where the estimand has one",
    { S <- attr(draws_obj, "nestimand_draws")
      attr(S, "nestimand_labels") <-
        data.frame(stratum = c("dim", "dim"), term = colnames(S),
                   stringsAsFactors = FALSE)
      z <- draws_obj; attr(z, "nestimand_draws") <- S
      identical(names(nest_draws(z)), c("drawid", "stratum", "term", "draw")) })
chk("draws: DxP is a data frame, one row per draw and one column per contrast",
    { d <- nest_draws(draws_obj, "DxP")
      is.data.frame(d) && nrow(d) == 2 &&
        identical(names(d), c("maj - aug", "maj - dim")) &&
        isTRUE(all.equal(d[["maj - aug"]], c(0.9, 1.1))) })
chk("draws: PxD is its transpose, rows named for the contrasts",
    { d <- nest_draws(draws_obj, "PxD")
      is.data.frame(d) && identical(rownames(d), c("maj - aug", "maj - dim")) &&
        identical(unname(as.matrix(d)),
                  t(unname(as.matrix(nest_draws(draws_obj, "DxP"))))) })
chk("draws: the column means are the estimates they came from",
    isTRUE(all.equal(unname(colMeans(nest_draws(draws_obj, "DxP"))),
                     as.data.frame(draws_obj)$estimate)))
## brms::hypothesis() splits its argument on the operators before parsing it, so
## a column called `maj - aug` is read as arithmetic on parameters `maj` and
## `aug`, and backticks do not rescue it. `names = "syntactic"` is for handing
## the draws to it.
chk("draws: syntactic names carry the direction of the comparison",
    { d <- nest_draws(draws_obj, "DxP", names = "syntactic")
      identical(names(d), c("maj_vs_aug", "maj_vs_dim")) })
chk("draws: a label starting with a digit takes the target's name in front",
    identical(syntactic_terms(c("1 - 0", "2 - 0"), "inversion"),
              c("inversion_1_vs_0", "inversion_2_vs_0")))
chk("draws: and a bare V where there is no target to use",
    identical(syntactic_terms("1 - 0"), "V_1_vs_0"))
chk("draws: every syntactic name parses as one symbol",
    { z <- syntactic_terms(c("dim: maj - aug", "chord_type:inversion", "1 - 0"),
                           "t")
      all(vapply(z, function(k)
        tryCatch(is.name(str2lang(k)), error = function(e) FALSE), TRUE)) })
chk("draws: the labels asked in travel with the renamed draws",
    { d <- nest_draws(draws_obj, "DxP", names = "syntactic")
      identical(unname(attr(d, "nestimand_terms")),
                c("maj - aug", "maj - dim")) })
chk("draws: and the map survives whichever shape is taken",
    all(vapply(c("long", "DxP", "PxD"), function(sh)
      !is.null(attr(nest_draws(draws_obj, sh, names = "syntactic"),
                    "nestimand_terms")), TRUE)))
chk("draws: names default to the terms, unchanged",
    identical(names(nest_draws(draws_obj, "DxP")),
              c("maj - aug", "maj - dim")))
chk("draws: a collection stacks its parts, each named by its target",
    { coll <- structure(list(chord_type = draws_obj, inversion = draws_obj),
                        class = "nestimand_estimands")
      cn <- colnames(nest_draws(coll, "DxP"))
      length(cn) == 4 && all(grepl("^chord_type: |^inversion: ", cn)) })
chk("draws: a frequentist estimand is refused with the reason",
    grepl("exist for a Bayesian fit",
          err_of(nest_draws(structure(data.frame(a = 1),
                                      class = c("nestimand_estimand", "data.frame"))))))
chk("draws: the long form is read when the wide one is not reachable",
    { o <- data.frame(term = c("a - b", "c - b"), estimate = c(1, 2))
      class(o) <- c("predictions", "data.frame")
      lng <- data.frame(rowid = rep(1:2, each = 3), drawid = rep(1:3, 2),
                        draw = c(1, 1, 1, 2, 2, 2))
      S <- with(list(), {
        r <- as.integer(as.factor(lng$rowid)); d <- as.integer(as.factor(lng$drawid))
        M <- matrix(NA_real_, max(d), max(r)); M[cbind(d, r)] <- lng$draw; M })
      identical(dim(S), c(3L, 2L)) && identical(colMeans(S), c(1, 2)) })
chk("draws: the accessor takes whichever spelling the installed version has",
    is.function(mfx_draws) &&
      is.null(mfx_draws(structure(list(), class = "not_a_result"))))

## An argument the coefficient route cannot take is one thing; an argument that
## is nothing's is another. Reporting the second as the first sends the reader
## to the wrong documentation - `method` was described as an argument of
## avg_predictions, which has no such argument.
chk("dots: an argument of the prediction function is named as one",
    { e <- err_of(nest_estimand(m_i3, top, type = "eta", vcov = "HC3"))
      grepl("is an argument of marginaleffects", e, fixed = TRUE) })
chk("dots: an argument of neither is named as neither",
    { e <- err_of(nest_estimand(m_i3, top, type = "eta", method = "within"))
      grepl("is not an argument of nest_estimand() or of marginaleffects", e,
            fixed = TRUE) })
chk("dots: a near miss of a real argument is suggested",
    grepl("did you mean `contrast`",
          err_of(nest_estimand(m_i3, top, type = "eta", contrst = "within")),
          fixed = TRUE))
chk("dots: a name close to nothing draws no guess",
    !grepl("did you mean",
           err_of(nest_estimand(m_i3, top, type = "eta", method = "within")),
           fixed = TRUE))
chk("dots: the message names the vocabulary that does exist",
    { e <- err_of(nest_estimand(m_i3, top, type = "eta", method = "within"))
      grepl("`contrast`", e, fixed = TRUE) && grepl("`by`", e, fixed = TRUE) &&
        grepl("`:` target", e, fixed = TRUE) })

## ---- `by`: the estimand within each level of something --------------------
## As marginaleffects groups an average, `by` groups an estimand: the target's
## contrasts are formed inside each group, the policy weighted over that group's
## conditions alone. It is the only way to group an estimand.
by_e <- suppressMessages(nest_estimand(m_i3, top, by = chord_type, policy = "proportional",
                                  bounds = FALSE))
chk("by: one row per group, labelled by the grouping variable",
    { d <- as.data.frame(by_e)
      "chord_type" %in% names(d) && nrow(d) == 4 &&
        setequal(d$chord_type, c("aug", "dim", "min", "maj")) })
chk("by: each group is the estimand computed over that group's cells",
    { d <- as.data.frame(by_e)
      hand <- vapply(c("aug", "dim", "min", "maj"), function(k) {
        cs <- sp_i3$cells[sp_i3$cells$chord_type == k, , drop = FALSE]
        pol <- nest_policy(sp_i3, "top", "proportional", cells = cs)
        as.data.frame(latent_estimand(m_i3, "top", pol, spec = sp_i3,
                                      cells = cs))$estimate }, 1)
      isTRUE(all.equal(d$estimate[match(names(hand), d$chord_type)],
                       unname(hand))) })
chk("by: the two routes agree on it",
    { a <- as.data.frame(suppressMessages(nest_estimand(m_i3, top, by = chord_type,
             policy = "proportional", bounds = FALSE)))
      b <- as.data.frame(suppressMessages(nest_estimand(m_i3, top, by = chord_type,
             policy = "proportional", type = "eta", bounds = FALSE)))
      isTRUE(all.equal(a$estimate, b$estimate)) })
chk("by: on a nested target it groups by the parent, as `within` used to",
    { a <- as.data.frame(suppressMessages(nest_estimand(m_i3, inversion, by = chord_type,
             bounds = FALSE)))
      nrow(a) == 9 && setequal(a$chord_type, c("dim", "min", "maj")) })
chk("by: several grouping variables give one row per realized combination",
    { d <- as.data.frame(suppressMessages(nest_estimand(m_i3, top,
             by = c(chord_type, inversion), bounds = FALSE)))
      all(c("chord_type", "inversion") %in% names(d)) &&
        nrow(d) == nrow(unique(sp_i3$cells[, c("chord_type", "inversion")])) })
chk("by: a nested target's own restriction still applies inside the grouping",
    { d <- as.data.frame(suppressMessages(nest_estimand(m_i3, inversion,
             by = chord_type, bounds = FALSE)))
      !("aug" %in% as.character(d$chord_type)) })
chk("by: a group in which the target does not vary is left out, and said",
    { ## `top` is t1 only within aug, so aug has no top contrast to form
      d2 <- cross_dat[!(cross_dat$chord_type == "aug" & cross_dat$top == "t2"), ]
      sp2 <- suppressMessages(nesting_spec(d2, response ~ chord_type * inversion * top,
                                           "inversion %in% chord_type"))
      m2 <- fit_from(sp2)
      z <- character(0)
      d <- as.data.frame(withCallingHandlers(
        nest_estimand(m2, top, by = chord_type, bounds = FALSE),
        message = function(m) { z <<- c(z, conditionMessage(m))
                                invokeRestart("muffleMessage") }))
      !("aug" %in% as.character(d$chord_type)) &&
        any(grepl("does not vary", z)) })
chk("by: the bounds are reported per group",
    { e <- suppressMessages(nest_estimand(m_i3, top, by = chord_type,
             policy = "proportional"))
      b <- attr(e, "nestimand")$bounds
      !is.null(b) && "chord_type" %in% names(b) && nrow(b) == 4 })
chk("by: a covariate cannot group an estimand",
    grepl("not one", err_of(nest_estimand(m_i3, top, by = response))))
chk("by: nor can the target group itself",
    grepl("also the target", err_of(nest_estimand(m_i3, top, by = top))))
chk("by: the emitted code is one runnable block per group",
    { cd <- attr(suppressMessages(nest_estimand(m_i3, top, by = chord_type,
             bounds = FALSE, dry_run = TRUE)),
             "nestimand_code")
      sum(grepl("^## group:", cd)) == 4 &&
        !inherits(try(parse(text = paste(cd, collapse = "\n")), silent = TRUE),
                  "try-error") })

## ---- the target is expanded by the formula machinery ----------------------
## A hand-rolled walk over the target expression handled `a * b`, then `a * b *
## c` once it was taught to recurse, then failed on `a * (b + c)`. R already
## knows how to expand a right-hand side; the results must be its term labels,
## whatever the bracketing.
for (.tg in list(quote(chord_type * inversion),
                 quote(chord_type * inversion * top),
                 quote(chord_type * (inversion + top)),
                 quote((chord_type + inversion) * top),
                 quote(chord_type + top))) {
  .labs <- attr(stats::terms(stats::as.formula(
    paste("~", paste(deparse(.tg), collapse = " ")))), "term.labels")
  chk(paste0("target `", paste(deparse(.tg), collapse = " "),
             "` gives the formula's terms"),
      { e <- suppressMessages(eval(bquote(nest_estimand(m_i3, .(.tg), policy = "equal",
               bounds = FALSE))))
        identical(names(e), .labs) })
}
chk("target: a bare interaction is still one result, not a list",
    { e <- suppressMessages(nest_estimand(m_i3, chord_type:inversion, bounds = FALSE))
      !inherits(e, "nestimand_estimands") &&
        identical(attr(e, "nestimand")$contrast, "interaction") })


## The effect basis reparameterizes the cell fit, which is saturated whatever
## the formula says, so the basis is saturated too: built from the declared
## terms it was smaller than the thing it is meant to translate, and
## `effect_basis()` refused a spec that was perfectly well formed.
chk("basis: a restricted formula still gives a square, full-rank basis",
    { sp <- suppressMessages(nesting_spec(cross_dat,
        response ~ chord_type * (inversion + top), "inversion %in% chord_type"))
      A <- effect_basis(sp)
      nrow(A) == nrow(sp$cells) && ncol(A) == nrow(A) && qr(A)$rank == nrow(A) })
chk("basis: and the same one the fully crossed formula gives",
    { r <- suppressMessages(nesting_spec(cross_dat,
        response ~ chord_type * (inversion + top), "inversion %in% chord_type"))
      x <- suppressMessages(nesting_spec(cross_dat,
        response ~ chord_type * inversion * top, "inversion %in% chord_type"))
      identical(chain_terms(r), chain_terms(x)) })
chk("basis: no term names a nested variable without its parent",
    { sp <- suppressMessages(nesting_spec(cross_dat,
        response ~ chord_type * inversion * top, "inversion %in% chord_type"))
      all(vapply(strsplit(chain_terms(sp), ":"), function(vs)
        all(unlist(lapply(vs, function(v) nest_ancestors(sp, v))) %in% vs), TRUE)) })
chk("basis: a chain is still the prefix ladder it always was",
    identical(chain_terms(sp2d), c("chord", "chord:inv", "chord:inv:doub")))

chk("saturation: a formula spanning less than the realized cells says so",
    grepl("one coefficient per effect it names",
          paste(utils::capture.output(
            nesting_spec(dat, response ~ chord_type + inversion,
                         "inversion %in% chord_type"), type = "message"),
            collapse = " ")))
chk("saturation: a formula that does span them is silent",
    !length(utils::capture.output(
      nesting_spec(dat, response ~ chord_type * inversion,
                   "inversion %in% chord_type"), type = "message")))

cat(sprintf("\n%d passed, %d failed\n", pass, fail))
