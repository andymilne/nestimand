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
    for (v in fam[-1]) {
      ## the strata of `v` are its ancestors, which is not the same as
      ## everything before it in the family once a parent holds two children
      above <- nest_ancestors(spec, v)
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

cell_formula <- function(spec, mode = c("cells", "reduced", "effects"),
                         intercept = NULL) {
  mode <- match.arg(mode)
  if (is.null(intercept)) intercept <- has_thresholds(spec)
  rhs <- if (mode == "reduced") {
    z <- setdiff(colnames(reduced_design(spec)), "(Intercept)")
    ## A covariate crossed with the structure has one slope per column of the
    ## design *and* one for the column the design does not carry - its intercept,
    ## which is the reference condition. Without that term the crossing would be
    ## one slope short of the structure it is crossed with.
    cov_terms <- unlist(lapply(spec$covariates, function(cv)
      if (isTRUE(spec$cov_by_cell[[cv]])) c(cv, paste0(z, ":", cv)) else cv))
    paste(c(z, cov_terms), collapse = " + ")
  } else if (mode == "cells") {
    cn <- spec$cell_name
    cov_terms <- unlist(lapply(spec$covariates, function(cv)
      if (isTRUE(spec$cov_by_cell[[cv]])) paste0(cn, ":", cv) else cv))
    paste(c(if (intercept) cn else paste0("0 + ", cn), cov_terms), collapse = " + ")
  } else {
    ## The effects parameterization reads its coefficients as a chain, which a
    ## term naming a nested variable without its parent has no reading as, so
    ## this is the one place the ancestry-closed form of the declaration is
    ## used. It is a niche path - emmeans, and priors stated coordinate-wise in
    ## effect space - and it is why closed_terms() still exists.
    tm <- closed_terms(spec)
    paste(c(tm, spec$covariates[!spec$cov_by_cell],
            unlist(lapply(spec$covariates[spec$cov_by_cell], function(cv)
              paste0(tm, ":", cv)))), collapse = " + ")
  }
  stats::as.formula(paste(spec$outcome, "~", rhs), env = parent.frame())
}

## The structure the formula declares, closed under ancestry: the terms the user
## wrote, with each nested variable given its parent, since a nested variable
## has no effect outside the strata it varies in - `inversion` alone is not a
## term this design admits, and becomes `chord_type:inversion`. This is what is
## fitted when the declaration asks for less than the saturated structure, and
## it is a subset of `chain_terms()`, never more.
## The terms the formula asks for, as written: each term reduced to the design
## variables it names, and nothing added. This is the whole of what the package
## fits. A `+` between two variables the design admits an interaction between is
## a restriction the user stated, and stating it is the point of writing `+`.
declared_terms <- function(spec, labels = spec$term_labels) {
  labs <- lapply(strsplit(labels, ":"), function(vs) vs[vs %in% spec$cell_vars])
  out <- unique(vapply(Filter(length, labs), paste, "", collapse = ":"))
  out[order(lengths(strsplit(out, ":")), out)]
}

## The ancestry-closed form of the same thing, in which every nested variable is
## named alongside its parent. The reduced design does not need it - the design
## over the realized cells says what is estimable without being told - and it is
## retained only for the effects parameterization, whose coefficients are read
## as a chain and where a term naming a nested variable without its parent has
## no such reading.
closed_terms <- function(spec, labels = spec$term_labels) {
  fams <- spec$cat_families
  rank_of <- function(v) { for (f in fams) if (v %in% f) return(match(v, f)); Inf }
  canon <- function(vs) { vs <- unique(vs); vs[order(vapply(vs, rank_of, 1), vs)] }
  labs <- lapply(strsplit(labels, ":"), function(vs)
    vs[vs %in% spec$cell_vars])
  labs <- labs[lengths(labs) > 0]
  closed <- lapply(labs, function(vs)
    canon(unique(c(vs, unlist(lapply(vs, function(v) nest_ancestors(spec, v)))))))
  out <- unique(vapply(closed, paste, "", collapse = ":"))
  out[order(lengths(strsplit(out, ":")), out)]
}

## How many dimensions of the realized-cell space a set of terms spans, which is
## what decides whether the declaration is the saturated structure or less. It is
## the width of the design those terms give, so that the count reported and the
## model fitted can never disagree.
term_span <- function(spec, terms) {
  if (!length(terms)) return(0L)
  ncol(design_over_cells(spec, terms))
}

## The sentinel is not a value the nested variable takes: it records that the
## variable is undefined there. So it is ordered last rather than first - a real
## level is the reference, and the contrasts are among the levels that exist -
## and the rows where it stands contribute nothing to any column the variable
## helps to build. Written as a level instead, it would make `chord_type +
## inversion` report `inversion 0 versus an augmented chord`, which is not an
## inversion effect, and would duplicate the parent's own term.
sentinel_absent <- function(spec, tab) {
  s <- sentinel_levels(spec)
  for (v in names(s)) {
    real <- setdiff(levels(factor(tab[[v]])), s[[v]])
    tab[[v]] <- factor(as.character(tab[[v]]), levels = c(real, s[[v]]))
  }
  tab
}

## The design a set of terms gives over the realized cells: one row per realized
## condition, the undefined conditions contributing nothing, and the columns the
## design cannot inform removed - identically zero first, then linearly
## redundant. What remains is full rank by construction, so the model fitted is
## the model the formula names and no coefficient has to be held at zero.
design_over_cells <- function(spec, terms, tab = spec$cells) {
  cellid <- as.character(tab[[spec$cell_name]])
  tab <- sentinel_absent(spec, tab)
  s <- sentinel_levels(spec)
  f <- stats::as.formula(paste("~", if (length(terms))
    paste(terms, collapse = " + ") else "1"))
  X <- stats::model.matrix(f, tab)
  tl <- c("(Intercept)", attr(stats::terms(f), "term.labels"))[attr(X, "assign") + 1L]
  for (v in names(s)) {
    involves <- vapply(strsplit(tl, ":", fixed = TRUE), function(z) v %in% z, TRUE)
    X[as.character(tab[[v]]) %in% s[[v]], involves] <- 0
  }
  keep <- colSums(X != 0) > 0                  # conditions that do not exist
  X <- X[, keep, drop = FALSE]; tl <- tl[keep]
  q <- qr(X); piv <- sort(q$pivot[seq_len(q$rank)])
  X <- X[, piv, drop = FALSE]; tl <- tl[piv]   # and the redundant ones
  rownames(X) <- cellid
  attr(X, "term_of") <- stats::setNames(tl, colnames(X))
  X
}

## The identified chain basis, closed under ancestry: retained because the
## effect-basis fitting mode and the emmeans engine both require the original
## factors as predictors (see `fitting_mode()`).
chain_terms <- function(spec) {
  vars <- spec$cell_vars
  if (length(vars) > 12)
    stop(length(vars), " categorical design variables is more structure than ",
         "the effect basis enumerates. The cell parameterization is unaffected; ",
         "it is the reading of the coefficients as effects that has no compact ",
         "form at this size.")
  fams <- spec$cat_families
  rank_of <- function(v) { for (f in fams) if (v %in% f) return(match(v, f)); Inf }
  canon <- function(vs) { vs <- unique(vs); vs[order(vapply(vs, rank_of, 1), vs)] }
  ## The basis is the saturated one over the realized cells, and does not follow
  ## the declared formula. Cells mode fits `~ 0 + cell` whatever the formula
  ## says, so a basis built from the declared terms would be smaller than the
  ## thing it is meant to reparameterize - `chord_type * (inversion + top)`
  ## spans 14 of 20 cells, and the translation would have no basis for the other
  ## six. What the structure does decide is which terms may appear: a variable
  ## is named only alongside its own ancestors, so every term is an effect
  ## within a stratum that exists, and a term naming a nested variable without
  ## its parent - which would be a marginal effect the design cannot support -
  ## never arises.
  closed <- Filter(length, lapply(seq_len(2^length(vars)) - 1L, function(i) {
    vs <- vars[bitwAnd(i, 2L^(seq_along(vars) - 1L)) > 0]
    if (!length(vs)) return(NULL)
    if (!all(unlist(lapply(vs, function(v) nest_ancestors(spec, v))) %in% vs))
      return(NULL)                       # not closed under ancestry
    canon(vs)
  }))
  out <- unique(vapply(closed, paste, "", collapse = ":"))
  out[order(lengths(strsplit(out, ":")), out)]
}

## --- the reduced design ----------------------------------------------------
## A declaration that asks for less than the saturated structure cannot be
## written on the cell factor, but it does not need the chain form either. Its
## design over the realized cells, with the columns the data cannot inform
## removed - identically zero, then redundant - is full rank by construction and
## needs nothing held at zero. It is the same object the cell factor is: one row
## per realized condition, computed once and carried in the data. The cell
## factor is its special case, the saturated structure coded as indicators.
##
## Column names are the effects they stand for, with `:` written `.` so that
## they are ordinary names in a formula and survive every engine unaltered.
reduced_design <- function(spec) {
  X <- design_over_cells(spec, declared_terms(spec))
  eff <- colnames(X)
  term_of <- attr(X, "term_of")
  colnames(X) <- reduced_names(eff)
  ## the effect each column stands for, kept verbatim so that reporting can
  ## name it as the user wrote it rather than as the syntactic column name
  attr(X, "effect_names") <- stats::setNames(eff, colnames(X))
  attr(X, "term_of") <- stats::setNames(unname(term_of), colnames(X))
  X
}

## The columns of the reduced design that carry a given set of terms. The random
## side takes its columns from the same design the fixed side is built from,
## rather than building its own: a term coded in isolation gets different columns
## from the same term coded beside its relatives - marginality decides that - and
## two designs built separately would not be two views of one model.
reduced_columns <- function(spec, terms = NULL) {
  X <- reduced_design(spec)
  cols <- setdiff(colnames(X), "(Intercept)")
  if (is.null(terms)) return(cols)
  tm <- attr(X, "term_of")
  cols[tm[cols] %in% terms]
}

## The readable name of a reduced column, for reporting. Unknown names - the
## covariate crossings the fitting formula forms, and anything the engine adds -
## are returned with the prefix removed and the separator restored, which is
## right for `dm_a.b:x` and harmless otherwise.
reduced_labels <- function(spec, terms) {
  map <- attr(reduced_design(spec), "effect_names")
  vapply(terms, function(z) {
    parts <- strsplit(z, ":", fixed = TRUE)[[1]]
    paste(ifelse(parts %in% names(map), map[parts], parts), collapse = ":")
  }, "", USE.NAMES = FALSE)
}

## The column names the design is carried under. `:` becomes `.` because a colon
## in a variable name would be read as an interaction operator, and the `dm_`
## prefix keeps a generated column from colliding with one already in the data -
## `with_reduced()` refuses rather than overwriting if it would. It stands for
## design matrix, which is all these columns are: nothing is centred or scaled.
reduced_names <- function(x)
  ifelse(x == "(Intercept)", "(Intercept)",
         paste0("dm_", gsub("[^A-Za-z0-9._]", ".", x)))

## The data with the reduced design carried alongside it, one column per
## identified effect, looked up by the cell each row belongs to.
with_reduced <- function(spec, data = spec$data) {
  keep <- setdiff(colnames(reduced_design(spec)), "(Intercept)")
  clash <- intersect(keep, names(data))
  if (length(clash))
    stop("the reduced design would overwrite the column(s) `",
         paste(clash, collapse = "`, `"), "` of the data. Rename them.")
  reduced_augment(spec, data)
}

## The same columns written onto a frame that may already carry them, and may
## carry them wrongly: a counterfactual grid is built by moving rows to other
## cells, and a column copied across with the row would still describe the cell
## the row came from. They are a function of the cell and are recomputed from
## it, every time, wherever a grid is formed.
reduced_augment <- function(spec, data) {
  if (!spec$cell_name %in% names(data))
    stop("the reduced design is looked up by realized cell, and this frame has ",
         "no `", spec$cell_name, "` column. Grids are translated with ",
         "add_cells() before the design is written onto them.")
  X <- reduced_design(spec)
  keep <- setdiff(colnames(X), "(Intercept)")
  data[keep] <- as.data.frame(X[as.character(data[[spec$cell_name]]), keep,
                                drop = FALSE], check.names = FALSE)
  data
}

## Does a declared random structure ask for less than the saturated one? Asked
## of the random side on its own: a formula can cross the structure fully in the
## mean and restrict it by group, and the cell factor can express the second no
## better than the first.
random_restricted <- function(spec) {
  if (is.null(spec$random_original)) return(FALSE)
  bl <- tryCatch(bar_terms_of(spec$random_original), error = function(e) NULL)
  if (is.null(bl)) return(FALSE)
  any(vapply(bl, function(b) {
    lhs <- tryCatch(attr(stats::terms(stats::as.formula(paste("~", b$lhs))),
                         "term.labels"), error = function(e) character(0))
    tm <- declared_terms(spec, lhs)
    length(tm) > 0 && term_span(spec, tm) < nrow(spec$cells)
  }, TRUE))
}

## Does this declaration need the reduced form at all? Asked of the spec rather
## than of a fit, so that a grid can be built the same way wherever it is built.
is_reduced <- function(spec)
  isTRUE(tryCatch(term_span(spec, declared_terms(spec)) < nrow(spec$cells),
                  error = function(e) FALSE))

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
  unique(unlist(lapply(bar_terms_of(bars), function(b)
    all.vars(stats::as.formula(paste("~", b$grp))))))
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
  ## `~ 0 + cell` is saturated by construction, so a declaration that asks for
  ## less cannot be written on it: fitting it there would enlarge the model the
  ## user wrote. The cell form is kept wherever the two agree, which is whenever
  ## the formula crosses the structure fully - there it is the same design in
  ## coordinates that read as cell means.
  sp <- term_span(spec, declared_terms(spec))
  if (sp < nrow(spec$cells))
    return(structure("reduced", reason = paste(
      "the formula spans", sp, "of the", nrow(spec$cells), "realized cells, so",
      "it asks for less than the saturated structure: the cell factor cannot",
      "express that, so the design the formula gives over the realized cells is",
      "carried per cell, with the columns it cannot inform removed")))
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
