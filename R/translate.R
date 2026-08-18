## nestimand: the translation module ---------------------------------------
## This file replaces the prototype's `model_formula()`. It owns the maps
## between the original variable space, in which questions, priors, and
## reports are stated, and the realized-cell space, in which estimation is
## unconditionally well posed.

## --- the sentinel ----------------------------------------------------------
## A nested variable takes some level only where it is undefined - the sentinel.
## It is found from the design rather than from its name: the level or levels
## realized only in strata where the variable does not vary.
sentinel_levels <- function(spec) {
  out <- list()
  for (fam in spec$cat_families) {
    if (length(fam) < 2) next
    for (k in seq_along(fam)[-1]) {
      v <- fam[k]; above <- fam[seq_len(k - 1)]
      tab <- spec$cells
      key <- do.call(paste, c(unname(lapply(above, function(a)
        as.character(tab[[a]]))), sep = "."))
      varies <- names(which(tapply(as.character(tab[[v]]), key,
                                   function(z) length(unique(z))) > 1))
      free <- unique(as.character(tab[[v]])[key %in% varies])
      s <- setdiff(unique(as.character(tab[[v]])), free)
      if (length(s)) out[[v]] <- s
    }
  }
  out
}

## The chain parameterization is legible only when the sentinel is the reference
## level: the structurally impossible cells are then exactly the all-zero
## columns named stratum:level. With a real level as reference the same
## information smears across a different set of zero and dependent columns.
## Reporting is unaffected - contrasts are stated in the declared order - so the
## reordering is confined to the design the chain form is built from, and is
## written into the emitted code rather than done silently.
sentinel_first <- function(spec, data = spec$data) {
  s <- sentinel_levels(spec)
  for (v in names(s)) {
    lv <- levels(factor(data[[v]]))
    data[[v]] <- factor(data[[v]], levels = c(s[[v]], setdiff(lv, s[[v]])))
  }
  data
}

sentinel_relevel_code <- function(spec, data_name) {
  s <- sentinel_levels(spec)
  if (!length(s)) return(NULL)
  c("## the chain form is legible only with the sentinel as reference level: the",
    "## impossible cells are then exactly the all-zero columns. Contrasts are",
    "## still reported in the order declared in the data.",
    vapply(names(s), function(v)
      sprintf('%s[["%s"]] <- relevel(factor(%s[["%s"]]), "%s")',
              data_name, v, data_name, v, s[[v]][1]), ""))
}

## --- fitting formulae ------------------------------------------------------

## Ordinal engines carry an intercept-like term already - the thresholds - so
## a zero-intercept cell coding is one parameter too many there. `clm()` warns
## and reinstates the intercept silently; brms would sample a weakly identified
## model. The single cell factor under default contrasts is full rank either
## way, since dummy coding of one factor cannot produce structural zeros, so
## only the intercept convention changes and no estimand is affected.
has_thresholds <- function(spec)
  spec$fit %in% c("clm", "clmm") ||
  (identical(spec$fit, "brm") && !is.null(spec$family) &&
     grepl("cumulative|sratio|cratio|acat", spec$family))

cell_formula <- function(spec, mode = c("cells", "effects"), intercept = NULL) {
  mode <- match.arg(mode)
  if (is.null(intercept)) intercept <- has_thresholds(spec)
  rhs <- if (mode == "cells") {
    cn <- spec$cell_name
    cov_terms <- unlist(lapply(spec$covariates, function(cv)
      if (isTRUE(spec$cov_by_cell[[cv]])) paste0(cn, ":", cv) else cv))
    paste(c(if (intercept) cn else paste0("0 + ", cn), cov_terms), collapse = " + ")
  } else {
    paste(c(chain_terms(spec), spec$covariates[!spec$cov_by_cell],
            unlist(lapply(spec$covariates[spec$cov_by_cell], function(cv)
              paste0(chain_terms(spec), ":", cv)))), collapse = " + ")
  }
  stats::as.formula(paste(spec$outcome, "~", rhs), env = parent.frame())
}

## The identified chain basis, closed under ancestry: retained because the
## effect-basis fitting mode and the emmeans engine both require the original
## factors as predictors (see `fitting_mode()`).
chain_terms <- function(spec) {
  fams <- spec$cat_families
  rank_of <- function(v) { for (f in fams) if (v %in% f) return(match(v, f)); Inf }
  ancestors <- function(v) {
    for (f in fams) if (v %in% f) return(f[seq_len(match(v, f) - 1)])
    character(0)
  }
  canon <- function(vars) { vars <- unique(vars); vars[order(vapply(vars, rank_of, 1), vars)] }
  labs <- spec$term_labels
  labs <- lapply(strsplit(labs, ":"), function(vs) vs[vs %in% spec$cell_vars])
  labs <- labs[lengths(labs) > 0]
  closed <- lapply(labs, function(vs) canon(unique(c(vs, unlist(lapply(vs, ancestors))))))
  used <- unique(unlist(closed))
  for (f in fams) if (any(f %in% used))
    for (k in seq_along(f)) if (any(f[seq_len(k)] %in% used))
      closed <- c(closed, list(f[seq_len(k)]))
  out <- unique(vapply(closed, paste, "", collapse = ":"))
  out[order(lengths(strsplit(out, ":")), out)]
}

## --- grid translation ------------------------------------------------------

## Recompute the cell factor from the crossed original factors. Mandatory for
## every grid: a grid built in original space carries no cell column, and a
## grid inherited from the data carries a stale one.
add_cells <- function(spec, newdata, unrealized = c("drop", "error")) {
  unrealized <- match.arg(unrealized)
  miss <- setdiff(spec$cell_vars, names(newdata))
  if (length(miss))
    stop("the grid does not contain the nesting variable(s) `",
         paste(miss, collapse = "`, `"), "`, so the realized cell cannot be ",
         "recomputed. Grids are constructed in the original variable space and ",
         "translated; supply every declared nesting variable.")
  key <- do.call(paste, c(unname(lapply(spec$cell_vars, function(v)
    as.character(newdata[[v]]))), sep = "."))
  bad <- !key %in% spec$cell_levels
  if (any(bad)) {
    if (unrealized == "error")
      stop(sum(bad), " row(s) of the grid name combinations that do not exist ",
           "in the design (", paste(utils::head(unique(key[bad]), 3), collapse = ", "),
           if (length(unique(key[bad])) > 3) ", ...", "). No model can predict ",
           "them, and any estimand averaging over them is not identified by the ",
           "design. Restrict the grid to realized cells, or declare the nesting ",
           "structure that makes these combinations impossible.")
    newdata <- newdata[!bad, , drop = FALSE]
    key <- key[!bad]
  }
  newdata[[spec$cell_name]] <- factor(key, levels = spec$cell_levels)
  attr(newdata, "nestimand_dropped") <- sum(bad)
  newdata
}

## The realized grid in original space: one row per cell, covariates held at
## their means (the `datagrid` convention), constructed so that no unrealized
## combination is ever formed in the first place.
cell_grid <- function(spec, data = spec$data, covariates = c("mean", "keep")) {
  covariates <- match.arg(covariates)
  g <- spec$cells
  if (covariates == "mean")
    for (cv in spec$covariates)
      g[[cv]] <- if (is.numeric(data[[cv]])) mean(data[[cv]], na.rm = TRUE)
                 else stats::na.omit(data[[cv]])[1]
  ## A grid for a mixed model still needs its grouping columns present: the
  ## prediction machinery reads them even when the random effects are excluded.
  ## Their value is immaterial, since the estimand is population-level.
  for (gv in grouping_vars(spec))
    if (gv %in% names(data) && !gv %in% names(g))
      g[[gv]] <- stats::na.omit(data[[gv]])[1]
  g
}

## the grouping factors named to the right of a bar
grouping_vars <- function(spec) {
  bars <- spec$random_original
  if (is.null(bars)) return(character(0))
  bl <- regmatches(bars, gregexpr("\\|[^)]*", bars))[[1]]
  unique(unlist(lapply(gsub("^\\|\\s*", "", bl), function(z)
    all.vars(stats::as.formula(paste("~", z))))))
}

## --- which columns carry the identification constraint ----------------------
## Beyond the structurally empty columns, one column per stratum is redundant:
## the stratum's own indicator is the sum of its cells. Which one is dropped is
## a coding choice, and it fixes the reference condition, so it follows the
## level order declared in the data - the first level that is not the sentinel -
## rather than whichever column a pivot happens to reach last.
reference_levels <- function(spec, data = spec$data) {
  s <- sentinel_levels(spec)
  out <- list()
  for (fam in spec$cat_families) for (v in fam[-1]) {
    lv <- levels(factor(data[[v]]))
    free <- setdiff(lv, s[[v]])
    if (length(free)) out[[v]] <- free[1]
  }
  out
}

## Which columns to constrain, given those reference choices. One column per
## stratum is redundant, and the choice is made constructively rather than by
## pattern: a candidate is dropped only if the remaining columns still span the
## same space. Deeper terms are offered first, so the constraint falls on the
## most specific coefficient, and a pivot settles anything left over.
identification_columns <- function(spec, X, data = spec$data) {
  refs <- reference_levels(spec, data)
  cols <- colnames(X)
  target_rank <- qr(X)$rank
  cand <- cols[vapply(cols, function(k)
    any(vapply(names(refs), function(v)
      grepl(paste0("(^|:)", v, refs[[v]], "($|:)"), k), TRUE)), TRUE)]
  cand <- cand[order(-lengths(strsplit(cand, ":")), match(cand, cols))]
  keep <- cols
  dropped <- character(0)
  for (k in cand) {
    if (length(keep) == target_rank) break
    trial <- setdiff(keep, k)
    if (qr(X[, trial, drop = FALSE])$rank == target_rank) {
      keep <- trial; dropped <- c(dropped, k)
    }
  }
  dropped
}

## --- the effect basis ------------------------------------------------------

## A: cell means as a linear function of identified chain-basis effects,
## mu = A m, with A square and full rank by construction. The inverse gives
## the effect coefficients implied by any set of cell means, which is what
## the prior and draw translations use.
effect_basis <- function(spec) {
  tab <- spec$cells
  ref <- sentinel_first(spec)
  for (v in spec$cell_vars) tab[[v]] <- factor(tab[[v]], levels = levels(factor(ref[[v]])))
  f <- stats::as.formula(paste("~", paste(chain_terms(spec), collapse = " + ")))
  X <- stats::model.matrix(f, tab)
  ## drop the structurally empty columns, then the declared reference columns;
  ## a pivot settles anything those two leave over
  empty <- colSums(X != 0) == 0
  Xne <- X[, !empty, drop = FALSE]
  drop_ref <- identification_columns(spec, Xne, ref)
  Xr <- Xne[, setdiff(colnames(Xne), drop_ref), drop = FALSE]
  q <- qr(Xr)
  A <- Xr[, q$pivot[seq_len(q$rank)], drop = FALSE]
  rownames(A) <- as.character(tab[[spec$cell_name]])
  if (nrow(A) != ncol(A) || qr(A)$rank != ncol(A))
    stop("the effect basis is not square and full rank (", nrow(A), " cells, ",
         ncol(A), " identified effects). This indicates an undeclared nesting ",
         "family or a cell table that is not the saturated realized structure.")
  structure(A, class = c("nestimand_basis", "matrix", "array"),
            cells = rownames(A), effects = colnames(A))
}

## --- fitting mode ----------------------------------------------------------

## The parameterization is itself a translation: chosen by the layer, never a
## silent difference. Cells by default; the effect basis when the requested
## engine or the declared priors can only be expressed there.
fitting_mode <- function(spec, engine = "marginaleffects", priors = NULL) {
  if (identical(engine, "emmeans"))
    return(structure("effects", reason = paste(
      "engine = \"emmeans\" is formula-driven and requires the original factors",
      "as predictors")))
  if (!is.null(priors) && isTRUE(attr(priors, "requires_effect_basis")))
    return(structure("effects", reason = paste(
      "the declared priors are coordinate-aligned in effect space (independent",
      "non-elliptical, or constrained support), which pins the fitting basis")))
  structure("cells", reason = "default: estimation is unconditionally well posed")
}

## --- the translation object ------------------------------------------------

translation <- function(spec) {
  A <- effect_basis(spec)
  structure(list(spec = spec, A = A, Ainv = solve(A),
                 cells = rownames(A), effects = colnames(A),
                 cell_name = spec$cell_name, cell_vars = spec$cell_vars),
            class = "nestimand_translation")
}

print.nestimand_translation <- function(x, ...) {
  cat("nestimand translation\n")
  cat("  original space :", paste(x$cell_vars, collapse = ", "), "\n")
  cat("  cell space     :", length(x$cells), "realized cells\n")
  cat("  effect basis   :", length(x$effects), "identified effects, full rank\n")
  invisible(x)
}

## --- engine compatibility --------------------------------------------------
## marginaleffects accepted `hypothesis = "pairwise"` up to 0.18.x and requires
## `hypothesis = ~pairwise` from 0.19.0. The package emits whichever the
## installed version accepts, rather than assuming one, so the same script runs
## either side of that change.
mfx_formula_hypothesis <- function(
    version = if (requireNamespace("marginaleffects", quietly = TRUE))
      utils::packageVersion("marginaleffects") else "0.19.0")
  utils::compareVersion(as.character(version), "0.19.0") >= 0

mfx_hypothesis <- function(contrast = "pairwise", version = NULL) {
  use_f <- if (is.null(version)) mfx_formula_hypothesis() else
    mfx_formula_hypothesis(version)
  if (use_f) stats::as.formula(paste("~", contrast)) else contrast
}

mfx_hypothesis_txt <- function(contrast = "pairwise", version = NULL) {
  use_f <- if (is.null(version)) mfx_formula_hypothesis() else
    mfx_formula_hypothesis(version)
  if (use_f) paste0("~", contrast) else sprintf('"%s"', contrast)
}

## Engine label conventions have moved between marginaleffects versions: the
## contrast label lived in `term` in 0.18.x and lives in `hypothesis` in recent
## releases; labels were bare (`aug - maj`) and are now parenthesised
## (`(maj) - (aug)`); and the direction of the difference reversed. The last is
## not cosmetic - inheriting it would mean the same analysis reporting
## opposite-signed effects on two machines - so this package fixes the
## convention itself: contrasts run in declared factor-level order, earlier
## level minus later, whatever the engine returns, and the estimate is negated
## where the engine ran it the other way.
mfx_term_column <- function(d) {
  d <- as.data.frame(d)
  cand <- c("term", "hypothesis", "contrast")
  cand <- c(intersect(cand, names(d)), setdiff(names(d), cand))
  for (nm in cand) {
    v <- d[[nm]]
    if ((is.character(v) || is.factor(v)) && any(grepl(" - ", as.character(v), fixed = TRUE)))
      return(nm)
  }
  NA_character_
}

## strip the engine's decoration: "(maj) - (aug)" -> c("maj", "aug")
mfx_split_label <- function(x) {
  x <- trimws(as.character(x))
  parts <- strsplit(x, " - ", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  sub("^\\((.*)\\)$", "\\1", parts)
}

mfx_canonical <- function(d, levs = NULL) {
  df <- as.data.frame(d)
  nm <- mfx_term_column(df)
  if (is.na(nm)) return(df)
  lab <- lapply(df[[nm]], mfx_split_label)
  ok <- lengths(lab) == 2
  if (!any(ok)) { df$term <- as.character(df[[nm]]); return(df) }
  if (is.null(levs)) levs <- unique(unlist(lab))
  rank_of <- function(z) { i <- match(z, levs); if (is.na(i)) length(levs) + 1L else i }
  ## canonical direction: later declared level minus earlier
  flip <- vapply(seq_along(lab), function(i)
    ok[i] && rank_of(lab[[i]][1]) < rank_of(lab[[i]][2]), TRUE)
  term <- vapply(seq_along(lab), function(i) {
    if (!ok[i]) return(as.character(df[[nm]][i]))
    p <- lab[[i]]; if (flip[i]) p <- rev(p)
    paste(p[1], "-", p[2])
  }, "")
  s <- ifelse(flip, -1, 1)
  for (cl in intersect(c("estimate", "statistic"), names(df))) df[[cl]] <- df[[cl]] * s
  if (all(c("conf.low", "conf.high") %in% names(df))) {
    lo <- df$conf.low; hi <- df$conf.high
    df$conf.low  <- ifelse(flip, -hi, lo)
    df$conf.high <- ifelse(flip, -lo, hi)
  }
  ## the engine's own label column would otherwise still read in its original
  ## direction, contradicting the estimate beside it
  if (!identical(nm, "term")) df[[nm]] <- term
  df$term <- term
  attr(df, "nestimand_flipped") <- sum(flip)
  df
}

## kept for callers that only need the label column normalized
mfx_with_term <- function(d, levs = NULL) mfx_canonical(d, levs)
