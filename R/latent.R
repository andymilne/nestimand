## nestimand: estimands as linear functionals of the cell coefficients -----
## On a linear scale - the latent scale of an ordinal model, the link scale of
## a generalized linear model - a policy contrast is exactly c'b, with c a
## weighted difference of design-matrix rows. Some engines will compute the same
## thing through their own prediction machinery - brms under type = "link",
## whose ordinal linear predictor excludes the thresholds - but clm's does not,
## the accepted type names have moved between releases, and the contrast vector
## c is needed in its own right for the heterogeneity of an effect. This is what the translation
## layer is for: the map from cells to named effects is a matrix, so the
## estimand needs no prediction machinery, is immune to the per-category
## expansion that ordinal fits provoke, and costs one matrix product.

## The contrast vector c for a policy, formed on the counterfactual grid so
## that the covariate distribution is common to every stratum. On a linear
## scale, averaging predictions over that grid and averaging its design-matrix
## rows are the same operation, so this is G-computation, not an approximation.
policy_contrast_matrix <- function(spec, target, policy, data = spec$data,
                                   model = NULL, cells = spec$cells,
                                   weights = NULL, route = "g_computation") {
  g <- policy_grid(spec, policy, data, cells, weights, route)
  X <- design_rows(spec, g, fit_mode(model))
  if (!is.null(model)) X <- align_design(X, model)
  lev <- target_stratum(spec, target, g)
  levs <- levels(factor(spec$data[[target]]))
  levs <- levs[levs %in% lev]
  M <- t(vapply(levs, function(s) {
    i <- lev == s
    ws <- g$.w[i] / sum(g$.w[i])
    as.numeric(crossprod(ws, X[i, , drop = FALSE]))
  }, numeric(ncol(X))))
  colnames(M) <- colnames(X)
  M
}

## The rows an estimand is averaged over, each weighted by the policy. The route
## says which rows, here as on the prediction side: every observed row crossed
## with every condition, or one row per condition with the covariates at their
## means. The two coincide unless a design column is a nonlinear function of a
## covariate.
policy_grid <- function(spec, policy, data = spec$data, cells = spec$cells,
                        weights = NULL, route = "g_computation") {
  g <- if (identical(route, "cells")) {
    d <- cell_grid(spec, data)
    d <- d[as.character(d[[spec$cell_name]]) %in%
           as.character(cells[[spec$cell_name]]), , drop = FALSE]
    d$.w <- policy_weights(spec, d, policy)
    d$.row <- seq_len(nrow(d))
    d
  } else counterfactual_grid(spec, data, policy, cells = cells)
  ## Unit weights standardize to another population. They multiply the row
  ## weights, exactly as they do on the prediction route: the contrast is a
  ## weighted average of design rows either way, and which weights are used is
  ## the only difference.
  if (!is.null(weights)) {
    wv <- if (is.character(weights) && length(weights) == 1L) data[[weights]]
          else weights
    if (length(wv) != nrow(data))
      stop("`weights` has ", length(wv), " values but the data has ",
           nrow(data), " rows: one weight per row is needed.")
    g$.w <- g$.w * wv[g$.row]
  }
  g
}

## --- a numeric target ------------------------------------------------------
## A step along a numeric target is named in marginaleffects' vocabulary, so
## that `contrast` means the same thing whichever route computes it: "sd" and
## "2sd" span one and two standard deviations around the mean, "iqr" the
## quartiles, "minmax" the range, a number a step of that size up from each
## row's own value, two numbers the step between them, and a function of the
## values the two columns it returns. The moments are those of the data the
## model was fitted to, as marginaleffects takes them.
step_values <- function(spec, target, step, x) {
  check_step(step, target)
  ref <- spec$data[[target]]
  if (is.function(step)) {
    v <- as.data.frame(step(x))
    if (ncol(v) != 2L || nrow(v) != length(x))
      stop("a function given as `contrast` takes the values of `", target,
           "` and returns two columns, the low and the high value to compare, ",
           "one row per value.", call. = FALSE)
    return(list(lo = as.numeric(v[[1]]), hi = as.numeric(v[[2]])))
  }
  if (is.character(step)) {
    m <- mean(ref, na.rm = TRUE); s <- stats::sd(ref, na.rm = TRUE)
    lh <- switch(step,
      sd = c(m - s / 2, m + s / 2), "2sd" = c(m - s, m + s),
      iqr = stats::quantile(ref, c(0.25, 0.75), names = FALSE, na.rm = TRUE),
      minmax = range(ref, na.rm = TRUE))
    return(list(lo = rep(lh[1], length(x)), hi = rep(lh[2], length(x))))
  }
  if (length(step) == 1L) return(list(lo = x, hi = x + step))
  list(lo = rep(min(step), length(x)), hi = rep(max(step), length(x)))
}

check_step <- function(step, target) {
  ok <- is.function(step) ||
    (is.character(step) && length(step) == 1L &&
       step %in% c("sd", "2sd", "iqr", "minmax")) ||
    (is.numeric(step) && all(is.finite(step)) &&
       ((length(step) == 1L && step != 0) ||
          (length(step) == 2L && step[1] != step[2])))
  if (!ok)
    stop("`", target, "` is numeric, so `contrast` is the step compared along ",
         "it: \"sd\" or \"2sd\" for one or two standard deviations, \"iqr\" or ",
         "\"minmax\", a number for a step of that size, two different numbers ",
         "for the step between them, or a function of the values returning the ",
         "low and the high value to compare. ", paste(deparse(step), collapse = " "),
         " is none of these.", call. = FALSE)
  invisible(TRUE)
}

## What a step is called in the estimand's `term` column, whichever route
## computed it.
step_label <- function(step) {
  if (is.function(step)) return("custom step")
  if (is.character(step))
    return(switch(step, sd = "+1 SD", "2sd" = "+2 SD", iqr = "Q3 - Q1",
                  minmax = "Max - Min"))
  if (length(step) == 1L) return(sprintf("%+g", step))
  sprintf("%g - %g", max(step), min(step))
}

## The size of a step in the target's own units, for saying what was compared.
## A function's step can differ from row to row, so it has no single size.
step_size <- function(spec, target, step) {
  if (is.function(step)) return(NA_real_)
  v <- step_values(spec, target, step, spec$data[[target]][1])
  v$hi - v$lo
}

## The contrast vector for a step: the change in each design row when the target
## moves from the low value to the high one, averaged over the policy-weighted
## grid - every observed row in every realized condition. On a linear scale
## that is exactly the average change in prediction, so this is G-computation,
## as it is for a categorical target. The two designs are reduced one at a time,
## so the grid's design is never held twice.
## With `by_cell`, one row per realized condition instead: the step's effect in
## each, which is what the bounds range over, from one pass over the grid.
step_contrast_matrix <- function(spec, target, policy, step, data = spec$data,
                                 model = NULL, cells = spec$cells,
                                 weights = NULL, route = "g_computation",
                                 by_cell = FALSE) {
  g <- policy_grid(spec, policy, data, cells, weights, route)
  v <- step_values(spec, target, step, g[[target]])
  key <- as.character(g[[spec$cell_name]])
  ws <- if (by_cell) g$.w / stats::ave(g$.w, key, FUN = sum) else g$.w / sum(g$.w)
  side <- function(val) {
    gs <- g; gs[[target]] <- val
    X <- design_rows(spec, gs, fit_mode(model))
    if (!is.null(model)) X <- align_design(X, model)
    if (by_cell) rowsum(X * ws, key)
    else matrix(as.numeric(crossprod(ws, X)), nrow = 1L,
                dimnames = list(step_label(step), colnames(X)))
  }
  hi <- side(v$hi)
  lo <- side(v$lo)
  hi - lo[rownames(hi), colnames(hi), drop = FALSE]
}

## The names an engine gives the columns of a design. brms rewrites the
## characters it does not allow in a parameter name - `I(x^2)` becomes `IxE2`,
## `poly(x, 2)1` becomes `polyx21` - and its own rule is applied here, so that
## design columns and the fit's coefficients are matched by name rather than by
## position.
engine_names <- function(model, cols) {
  if (!inherits(model, "brmsfit")) return(cols)
  rn <- tryCatch(get("rename", envir = asNamespace("brms")),
                 error = function(e) NULL)
  if (is.null(rn)) cols else rn(cols)
}

coef_vector <- function(model) {
  if (inherits(model, "clm") || inherits(model, "clmm")) {
    b <- model$beta
    if (is.null(b)) b <- coef(model)
    return(b)
  }
  if (inherits(model, "merMod") || methods::is(model, "merMod"))
    return(lme4::fixef(model))
  if (inherits(model, "glmmTMB")) return(unlist(glmmTMB::fixef(model)$cond))
  if (inherits(model, "brmsfit")) {
    b <- brms::fixef(model)[, "Estimate"]
    names(b) <- rownames(brms::fixef(model))
    return(b)
  }
  stats::coef(model)
}

vcov_beta <- function(model, nm) {
  ## Several engines return a Matrix rather than a base matrix - lme4 among
  ## them - and the arithmetic then dispatches to methods whose result is not
  ## always a base array. Coercing once here keeps everything downstream
  ## ordinary, and costs nothing at these sizes.
  V <- tryCatch(as.matrix(stats::vcov(model)),
                error = function(e) stats::vcov(model))
  if (inherits(model, "brmsfit")) {
    dn <- colnames(V)
    V <- V[dn %in% nm, dn %in% nm, drop = FALSE]
    return(V)
  }
  keep <- intersect(nm, colnames(V))
  if (!length(keep))
    stop("none of the ", length(nm), " coefficients the contrast needs appear ",
         "in the model's covariance matrix, which has ",
         if (is.null(colnames(V))) "no column names" else
           paste(ncol(V), "columns"), ". Report the model class (",
         paste(class(model), collapse = "/"), ").")
  as.matrix(V[keep, keep, drop = FALSE])
}

## Pairwise, reference, or sequential differences of the stratum rows. Each runs
## later declared level minus earlier, so that it reads as a departure from the
## reference condition.
contrast_pairs <- function(levs, contrast = "pairwise") {
  switch(contrast,
    pairwise = { ij <- utils::combn(seq_along(levs), 2)
                 lapply(seq_len(ncol(ij)), function(k) rev(ij[, k])) },
    reference = lapply(seq_along(levs)[-1], function(j) c(j, 1L)),
    sequential = lapply(seq_along(levs)[-1], function(j) c(j, j - 1L)),
    ## no comparisons: the levels stand as they are, and what is reported is the
    ## policy-weighted prediction for each rather than a difference between two
    none = list(),
    stop("contrast must be pairwise, reference, sequential or none here."))
}

## The design of a grid, in the parameterization the fit is in. The effects
## form is fitted on data with the sentinel as reference level - `nest_fit()`
## relevels it, and writes the relevel into the emitted code - so a grid coded
## from the data as it stands would carry `inversionnone` columns the fit does
## not have. The columns must be built the way the fit built them, or the map
## between them is not a map at all.
design_rows <- function(spec, g, mode = "cells") {
  ref <- spec$data
  ## The grid may cover only part of the design - one stratum of it, when an
  ## estimand is asked for within each level of something - and a factor with
  ## one level in that slice has no contrasts to code, so `model.matrix()`
  ## refuses. The columns have to be the ones the fit has either way, so every
  ## factor keeps the levels the fit saw rather than the levels this grid
  ## happens to contain.
  for (v in intersect(names(g), names(ref)))
    if (is.factor(ref[[v]]) || is.character(ref[[v]]))
      g[[v]] <- factor(as.character(g[[v]]), levels = levels(factor(ref[[v]])))
  ## Where the fit names columns of the design rather than the original factors
  ## - the reduced form, or slopes crossed with part of the structure - the grid
  ## carries them the same way the fitted data does: looked up by the cell each
  ## row belongs to.
  if (uses_reduced(spec)) g <- reduced_augment(spec, g)
  ## the bases the fit fixed, not ones recomputed from this grid
  X <- fitting_design(spec, g, mode)
  if (colnames(X)[1] == "(Intercept)") X[, -1, drop = FALSE] else X
}

## Which parameterization a fit is in, so that the design rows built here match
## its coefficients. A fit from `nest_fit()` carries it; anything else is read
## as the cell form.
fit_mode <- function(model) {
  md <- attr(model, "nestimand_mode")
  if (is.null(md)) "cells" else md
}

## A design column whose coefficient the engine dropped as aliased is held at
## zero and contributes nothing, so it is removed from the design rather than
## carried into the product as NA. Only a column with no coefficient at all -
## which means the model was not fitted from this declaration - is an error.
align_design <- function(X, model, what = "cell") {
  colnames(X) <- engine_names(model, colnames(X))
  b <- coef_vector(model)
  aliased <- names(b)[is.na(b)]
  missing_cols <- setdiff(colnames(X), names(b))
  if (length(missing_cols))
    stop("the fitted model has no coefficient for ",
         length(missing_cols), " column(s) of the ", what, " design matrix (",
         paste(utils::head(missing_cols, 3), collapse = ", "),
         if (length(missing_cols) > 3) ", ...", "), so the linear map cannot ",
         "be formed. This happens when the model was not fitted from this ",
         "declaration; refit with nest_fit().")
  X[, setdiff(colnames(X), aliased), drop = FALSE]
}

## One design row per realized cell, covariates averaged over the grid, so that
## a contrast among cells is c'b in the same way a contrast among strata is.
cell_design_rows <- function(spec, data, cells, model) {
  pol <- structure(list(kind = "uniform", target = NULL,
    p = NULL, at = NULL), class = "nestimand_policy")
  g <- counterfactual_grid(spec, data, cells = cells)
  X <- design_rows(spec, g, fit_mode(model))
  colnames(X) <- engine_names(model, colnames(X))
  b <- coef_vector(model)
  X <- X[, intersect(colnames(X), names(b)[!is.na(b)]), drop = FALSE]
  key <- as.character(g[[spec$cell_name]])
  lv <- as.character(cells[[spec$cell_name]])
  M <- t(vapply(lv, function(k) colMeans(X[key == k, , drop = FALSE]),
                numeric(ncol(X))))
  rownames(M) <- lv
  M
}

## The contrast matrix C, one row per reported contrast, columns the model's
## coefficients: the whole estimand is C b, whether that is evaluated at the
## coefficients or draw by draw.
## Which realized cells an estimand of `target` is taken over: all of them,
## less any stratum in which the target does not vary, since such a stratum
## offers no version to weight and would enter the average as a constant. One
## function decides it, so that every route to a contrast is over the same
## cells - `latent_draws()` used to leave the choice unmade and pass NULL on.
estimand_cells <- function(spec, target, cells = NULL) {
  if (!is.null(cells)) return(cells)
  dg <- degenerate_strata_multi(spec, target)
  if (is.null(dg) || !length(dg$drop)) spec$cells else
    spec$cells[deg_key(spec$cells, dg$vars) %in% dg$keep, , drop = FALSE]
}

latent_contrast_matrix <- function(model, spec, target, policy = "equal",
                                   at = NULL, contrast = "pairwise",
                                   data = spec$data, cells = NULL,
                                   weights = NULL, route = "g_computation") {
  ## A numeric target is compared at two of its values, and `contrast` is the
  ## step between them; the policy weights the conditions it is compared in.
  if (is_numeric_target(spec, target)) {
    if (is.null(cells)) cells <- spec$cells
    pol <- if (inherits(policy, "nestimand_policy")) policy
           else nest_policy(spec, target, policy, at, data, cells = cells)
    return(step_contrast_matrix(spec, target, pol, contrast, data, model,
                                cells = cells, weights = weights, route = route))
  }
  ## Several variables named at once are one condition, not several targets:
  ## an interaction contrasts differences between their combinations, and no
  ## contrast reports the combinations themselves. Both work over the same
  ## cells, and neither weighs anything, so neither takes a policy.
  combo <- identical(contrast, "interaction") ||
    (identical(contrast, "none") && length(target) > 1L)
  if (!combo) cells <- estimand_cells(spec, target, cells)
  if (combo) {
    cells <- estimand_cells(spec, target)
    M <- cell_design_rows(spec, data, cells, model)
    ## Several cells can share one combination of the targets, whenever the
    ## design holds a variable the interaction does not name. The contrast is
    ## then formed on their average, which is what the prediction route does by
    ## averaging predictions; matching a combination to the first cell that
    ## carries it would answer the question at one arbitrary level of the
    ## others, silently.
    key <- do.call(paste, c(unname(lapply(target, function(v)
      as.character(cells[[v]]))), sep = "\r"))
    if (anyDuplicated(key)) {
      idx <- split(seq_along(key), key)
      M <- do.call(rbind, lapply(idx, function(i) colMeans(M[i, , drop = FALSE])))
      cells <- cells[vapply(idx, `[`, 1L, 1L), target, drop = FALSE]
      rownames(M) <- NULL; rownames(cells) <- NULL
    }
    if (identical(contrast, "none")) {
      rownames(M) <- do.call(paste, c(unname(lapply(target, function(v)
        as.character(cells[[v]]))), sep = ", "))
      return(M)
    }
    H <- interaction_matrix(cells, target)
    C <- t(H) %*% M
    rownames(C) <- colnames(H)
    return(C)
  }
  pol <- if (inherits(policy, "nestimand_policy")) policy
         else nest_policy(spec, target, policy, at, data, cells = cells)
  M <- policy_contrast_matrix(spec, target, pol, data, model, cells = cells,
                              weights = weights, route = route)
  ## Each row of M is already the policy-weighted average of the design rows for
  ## one level of the target - the marginal mean. The contrasts are differences
  ## between its rows, so asking for none is asking for M itself.
  if (identical(contrast, "none")) return(M)
  levs <- rownames(M)
  prs <- contrast_pairs(levs, contrast)
  C <- do.call(rbind, lapply(prs, function(p) M[p[1], ] - M[p[2], ]))
  rownames(C) <- vapply(prs, function(p) paste(levs[p[1]], "-", levs[p[2]]), "")
  colnames(C) <- colnames(M)
  C
}

latent_estimand <- function(model, target, policy = "equal", at = NULL,
                            contrast = "pairwise", data = NULL,
                            conf_level = 0.95, spec = NULL, cells = NULL,
                            ndraws = NULL, weights = NULL,
                            route = "g_computation", re_formula = NA) {
  spec <- resolve_spec(model, spec)
  if (is.null(data)) data <- spec$data
  C <- latent_contrast_matrix(model, spec, target, policy, at, contrast, data,
                              cells, weights, route)
  ## One row per reported contrast. Anything else means the contrast matrix has
  ## been built over grid rows rather than over conditions, and the result would
  ## be meaningless as well as enormous.
  if (nrow(C) > 10000)
    stop("the contrast matrix has ", format(nrow(C), big.mark = ","),
         " rows, where one per comparison was expected. This is a bug in ",
         "nestimand rather than in the model: please report the model class (",
         paste(class(model), collapse = "/"), "), the number of realized cells (",
         nrow(spec$cells), ") and the contrast (", contrast, ").")
  lo <- (1 - conf_level) / 2

  ## A Bayesian fit is summarized from its draws. Applying the delta method to a
  ## posterior covariance would give a normal approximation to the posterior and
  ## a frequentist test statistic computed from it, which has no Bayesian
  ## reading; the draws give the posterior of the contrast itself.
  if (inherits(model, "brmsfit")) {
    D <- as.matrix(brms::as_draws_matrix(model))
    if (!is.null(ndraws) && ndraws < nrow(D))
      D <- D[sort(sample.int(nrow(D), ndraws)), , drop = FALSE]
    nm <- draw_names(colnames(C), colnames(D))
    B <- D[, nm, drop = FALSE]
    ## `re_formula = NULL` asks for the sampled groups rather than the average
    ## one. Draw by draw the group deviations do not average to zero, and that
    ## is where the extra width comes from: the estimate is unchanged, the
    ## posterior is wider. Only the group mean is needed, since the estimand
    ## averages over groups anyway.
    if (is.null(re_formula)) {
      U <- group_mean_draws(D, colnames(C))
      B <- B + U
    }
    S <- B %*% t(C)
    pd <- apply(S, 2, function(z) max(mean(z > 0), mean(z < 0)))
    out <- data.frame(
      term = rownames(C),
      estimate = apply(S, 2, mean), std.error = apply(S, 2, stats::sd),
      conf.low = apply(S, 2, stats::quantile, probs = lo),
      conf.high = apply(S, 2, stats::quantile, probs = 1 - lo),
      pd = pd, row.names = NULL)
    if (!is.null(attr(C, "strata")))
      out <- cbind(stratum = attr(C, "strata"), out)
    attr(out, "nestimand_scale") <- "linear predictor (latent / link), posterior"
    attr(out, "nestimand_cvecs") <- stats::setNames(
      lapply(seq_len(nrow(C)), function(i) stats::setNames(C[i, ], colnames(C))),
      rownames(C))
    ## labelled as the table is, stratum included where there is one
    attr(out, "nestimand_draws") <- name_draw_columns(S, out)
    return(out)
  }

  b <- coef_vector(model)[colnames(C)]
  V <- vcov_beta(model, colnames(C))[colnames(C), colnames(C), drop = FALSE]
  ## C should have one row per contrast and one column per coefficient. If it
  ## does not, the product below would allocate a matrix the size of the grid
  ## squared, and R would report a long-vector error from deep inside it.
  check_contrast_shape(C, V, b)
  est <- as.numeric(C %*% b)
  ## diag(C V C') row by row: forming the full product would allocate one
  ## square per contrast pair, which on a large contrast set exceeds what R
  ## can address, and every off-diagonal entry is discarded anyway.
  se  <- sqrt(quad_form_diag(C, V))
  z <- stats::qnorm(1 - lo)
  out <- data.frame(term = rownames(C), estimate = est, std.error = se,
                    statistic = est / se,
                    p.value = 2 * stats::pnorm(-abs(est / se)),
                    conf.low = est - z * se, conf.high = est + z * se,
                    row.names = NULL)
  if (!is.null(attr(C, "strata")))
    out <- cbind(stratum = attr(C, "strata"), out)
  attr(out, "nestimand_scale") <- "linear predictor (latent / link)"
  attr(out, "nestimand_cvecs") <- stats::setNames(
    lapply(seq_len(nrow(C)), function(i) stats::setNames(C[i, ], colnames(C))),
    rownames(C))
  out
}

## The Bayesian counterpart: the same c, applied draw by draw, so the result is
## a posterior rather than a delta-method interval.
latent_draws <- function(model, target, policy = "equal", at = NULL,
                         contrast = "pairwise", data = NULL, spec = NULL,
                         cells = NULL, weights = NULL,
                         route = "g_computation") {
  spec <- resolve_spec(model, spec)
  ## before the engine check, so that a misplaced argument is reported as one
  ## rather than as whatever the engine happens to be
  check_target(spec, target)
  if (is.null(data)) data <- spec$data
  if (!inherits(model, "brmsfit"))
    stop("draw-wise translation needs a posterior; this is a frequentist fit. ",
         "Use latent_estimand() for the delta-method interval.")
  cells <- estimand_cells(spec, target, cells)
  pol <- if (inherits(policy, "nestimand_policy")) policy
         else nest_policy(spec, target, policy, at, data, cells = cells)
  M <- policy_contrast_matrix(spec, target, pol, data, model, cells = cells,
                              weights = weights, route = route)
  D <- as.matrix(brms::as_draws_matrix(model))
  nm <- draw_names(colnames(M), colnames(D))
  B <- D[, nm, drop = FALSE]
  levs <- rownames(M)
  prs <- contrast_pairs(levs, contrast)
  out <- do.call(cbind, lapply(prs, function(p)
    as.numeric(B %*% (M[p[1], ] - M[p[2], ]))))
  colnames(out) <- vapply(prs, function(p) paste(levs[p[1]], "-", levs[p[2]]), "")
  as.data.frame(out)
}

## brms names population-level coefficients `b_` followed by the design-matrix
## column name, punctuation included (`b_cellaug.none`). Matching is done on
## the exact name first and, only if that fails, on a punctuation-insensitive
## comparison, so a future change to either convention is survivable.
draw_names <- function(cols, draw_cols) {
  exact <- paste0("b_", cols)
  if (all(exact %in% draw_cols)) return(exact)
  flat <- function(z) tolower(gsub("[^A-Za-z0-9]", "", z))
  key <- flat(sub("^b_", "", draw_cols))
  hit <- match(flat(cols), key)
  if (anyNA(hit))
    stop("these cell coefficients were not found among the posterior draws: ",
         paste(cols[is.na(hit)], collapse = ", "),
         ".\n  Draw columns available: ",
         paste(utils::head(grep("^b_", draw_cols, value = TRUE), 6), collapse = ", "),
         "\n  Report these names: the matching rule needs updating.")
  draw_cols[hit]
}


## A shape check with a legible failure: the alternative is an allocation error
## from inside a matrix product, which says nothing about what went wrong.
check_contrast_shape <- function(C, V, b) {
  if (nrow(C) > 1000)
    stop("the contrast matrix has ", format(nrow(C), big.mark = ","), " rows, ",
         "which is not a set of contrasts: it has the shape of the grid ",
         "itself. Something upstream has failed to aggregate. Report this with ",
         "dim(C) = ", paste(dim(C), collapse = " x "), " and ",
         length(b), " coefficients.", call. = FALSE)
  if (ncol(C) != length(b) || ncol(C) != ncol(V))
    stop("the contrast matrix has ", ncol(C), " columns, against ", length(b),
         " coefficients and a ", paste(dim(V), collapse = " x "),
         " covariance: they cannot be multiplied. This usually means the fit ",
         "does not carry the cell design's coefficients - a model fitted by ",
         "hand in crossed or chain form - in which case ask for a quantity ",
         "computed by prediction instead.", call. = FALSE)
  invisible(TRUE)
}


## The diagonal of C V C', computed row by row. Forming the full product would
## allocate one entry per pair of contrasts and discard all but the diagonal;
## on a large contrast set that exceeds what R can address.
quad_form_diag <- function(C, V) {
  C <- as.matrix(C)
  V <- as.matrix(V)
  as.numeric(rowSums((C %*% V) * C))
}


## The draw-wise mean over groups of each group-level coefficient, aligned to
## the fixed-effect columns a contrast uses. A coefficient with no group-level
## counterpart contributes nothing, and one that is common to every condition
## cancels in a contrast, so both are left at zero.
group_mean_draws <- function(D, coefs) {
  U <- matrix(0, nrow(D), length(coefs), dimnames = list(NULL, coefs))
  rcols <- grep("^r_", colnames(D), value = TRUE)
  if (!length(rcols)) return(U)
  ## r_<group>[<level>,<term>]
  term <- sub("^r_[^\\[]+\\[[^,]+,(.*)\\]$", "\\1", rcols)
  for (k in seq_along(coefs)) {
    want <- sub("^b_", "", coefs[k])
    i <- which(punct_free(term) == punct_free(want))
    if (length(i)) U[, k] <- rowMeans(D[, rcols[i], drop = FALSE])
  }
  U
}

punct_free <- function(x) gsub("[^A-Za-z0-9]", "", x)
