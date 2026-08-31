## nestimand: the weighting policy -----------------------------------------
## Every weighting scheme is a distribution p over the versions of a compound
## condition - the levels of the nested variable realized within a stratum -
## and each such p names a stochastic intervention. The estimand is therefore
## a causal contrast between two well-defined policies rather than an
## approximation to an atomic contrast.

policy_aliases <- c("equal", "proportional", "hierarchical",
                    "nominated", "standardized", "within")

## Versions of the compound condition within each stratum of `target`.
## The stratum is a level of `target`; the versions are the realized
## combinations of the variables nested below it in the same family.
versions_of <- function(spec, target, cells = spec$cells) {
  fam <- NULL
  for (f in spec$cat_families) if (target %in% f) fam <- f
  if (is.null(fam))
    stop("`", target, "` is not a categorical variable of any declared nesting ",
         "family (", paste(unlist(spec$cat_families), collapse = ", "), ").")
  ## The versions of a condition are the realized cells it occurs in, described
  ## by the other categorical variables. For a nesting root those are the
  ## variables below it; for a nested variable they are the strata it appears
  ## in, which is the same construction read the other way.
  others <- setdiff(spec$cell_vars, target)
  key <- if (length(others))
    do.call(paste, c(unname(lapply(others, function(v) as.character(cells[[v]]))),
                     sep = ".")) else rep("", nrow(cells))
  split(key, as.character(cells[[target]]))
}

## Strata in which the target does not vary: there the comparison would leave
## the target's own levels and compare strata instead, under the target's label.
degenerate_strata <- function(spec, target) {
  above <- nest_ancestors(spec, target)
  if (!length(above)) return(NULL)
  tab <- spec$cells
  key <- do.call(paste, c(unname(lapply(above, function(v) as.character(tab[[v]]))),
                          sep = "."))
  n <- tapply(as.character(tab[[target]]), key, function(z) length(unique(z)))
  list(vars = above, keep = names(n)[n > 1], drop = names(n)[n == 1])
}

## The restriction for a set of targets. Each nested variable is compared only
## where it varies, and an interaction exists only where every one of its
## variables does, so the restrictions are intersected. With a chain the deepest
## target's restriction implied the others', which is why taking it alone was
## enough until a target could sit outside that chain - a variable crossed with
## the structure, whose own restriction is empty, was being read as no
## restriction at all, and the sentinel level of the nested target survived into
## the comparison.
degenerate_strata_multi <- function(spec, targets, cells = spec$cells) {
  degs <- lapply(targets, function(v) degenerate_strata(spec, v))
  degs <- Filter(function(d) !is.null(d) && length(d$drop), degs)
  if (!length(degs)) return(NULL)
  vars <- unique(unlist(lapply(degs, `[[`, "vars")))
  keep_row <- rep(TRUE, nrow(cells))
  for (d in degs) keep_row <- keep_row & (deg_key(cells, d$vars) %in% d$keep)
  key <- deg_key(cells, vars)
  list(vars = vars, keep = unique(key[keep_row]), drop = unique(key[!keep_row]))
}

## The stratum key `degenerate_strata()` reports is composite - one part per
## ancestor - so a cell table is restricted by rebuilding that key, not by
## matching the first ancestor alone. With one ancestor the two coincide,
## which is why the difference only shows at a nesting depth of three or more.
deg_key <- function(tab, vars)
  do.call(paste, c(unname(lapply(vars, function(v) as.character(tab[[v]]))),
                   sep = "."))

## p as a per-stratum named vector, from an alias or a supplied distribution.
nest_policy <- function(spec, target, policy = "equal", at = NULL, data = spec$data,
                        cells = spec$cells) {
  vs <- versions_of(spec, target, cells)
  if (is.character(policy) && length(policy) == 1L) {
    if (identical(policy, "counterfactual"))
      stop("`counterfactual` named two different things and no longer names a ",
           "policy. For the empirical version frequencies use ",
           "policy = \"proportional\"; for the grid on which the model is ",
           "evaluated use route = \"g_computation\", which is the default.")
    kind <- match.arg(policy, policy_aliases)
    if (kind == "within")
      stop("`within` is not a policy: it emits per-stratum contrasts, which do ",
           "not cross the structural boundary and so require no distribution ",
           "over versions. Use estimand(..., contrast = \"within\").")
    if (kind == "nominated" && is.null(at))
      stop("policy = \"nominated\" places all mass on one version and needs it ",
           "named, e.g. at = c(inversion = \"0\").")
    if (kind == "standardized")
      stop("policy = \"standardized\" takes p from outside the design; supply it ",
           "directly, e.g. policy = c(\"0\" = 0.5, \"1\" = 0.3, \"2\" = 0.2).")
    p <- switch(kind,
      equal = lapply(vs, function(v) stats::setNames(rep(1 / length(v), length(v)), v)),
      proportional = {
        key <- do.call(paste, c(unname(lapply(spec$cell_vars, function(x)
          as.character(data[[x]]))), sep = "."))
        cellcount <- table(factor(key, levels = spec$cell_levels))
        vkey <- stats::setNames(
          do.call(paste, c(unname(lapply(setdiff(spec$cell_vars, target), function(x)
            as.character(cells[[x]]))), sep = ".")),
          as.character(cells[[spec$cell_name]]))
        lapply(names(vs), function(s) {
          rows <- names(vkey)[as.character(cells[[target]]) == s]
          cnt <- as.numeric(cellcount[rows])
          stats::setNames(cnt / sum(cnt), vs[[s]])
        }) |> stats::setNames(names(vs))
      },
      hierarchical = {
        hw <- hierarchical_weights(spec, target, cells)
        stats::setNames(lapply(names(vs), function(s) hw[[s]][vs[[s]]]), names(vs))
      },
      nominated = lapply(vs, function(v) {
        if (length(v) == 1L) return(stats::setNames(1, v))  # degenerate stratum
        pick <- version_key(spec, target, at, v)
        if (!pick %in% v)
          stop("the nominated version `", pick, "` is not realized in a stratum ",
               "that has more than one version (versions present: ",
               paste(v, collapse = ", "), "). A nominated contrast must name a ",
               "version that exists wherever a choice is available.")
        stats::setNames(as.numeric(v == pick), v)
      }))
  } else {
    ## a supplied distribution over version labels, applied in every stratum in
    ## which those versions are realized
    p <- lapply(names(vs), function(s) {
      v <- vs[[s]]
      ## A degenerate stratum has one version and no choice to make: the
      ## policy is vacuous there, and the single realized condition takes all
      ## the mass. This is the sentinel stratum of a partially nested design.
      if (length(v) == 1L) return(stats::setNames(1, v))
      missing_v <- setdiff(v, names(policy))
      if (length(missing_v))
        stop("the supplied policy names no mass for version(s) `",
             paste(missing_v, collapse = "`, `"), "`, which are realized in ",
             "stratum `", s, "`. Restricting p to the named versions and ",
             "renormalizing would be a different intervention, silently ",
             "chosen; state the mass for every realized version, including a ",
             "zero where a version is deliberately excluded.")
      w <- as.numeric(policy[match(v, names(policy))])
      if (sum(w) == 0)
        stop("the supplied policy places zero mass on every realized version ",
             "of stratum `", s, "`.")
      stats::setNames(w / sum(w), v)
    })
    names(p) <- names(vs)
    kind <- "supplied"
  }
  structure(list(kind = kind, target = target, p = p, at = at),
            class = "nestimand_policy")
}

## uniform at each level of the chain, rather than uniform over the leaves:
## the two coincide at depth one and differ at depth two and beyond.
## Hierarchical weighting: mass splits equally at each node of the declared
## structure, not equally over the leaves. The split for a variable is made
## conditional on that variable's own ancestors - never on its position among
## the other variables - so variables sharing a parent are independent choices
## whose probabilities multiply, and the result does not depend on the order the
## declarations were written. On a chain this is the successive split down the
## levels, unchanged. Where the realized combinations are fewer than that
## product, the weights are renormalized over what exists, as everywhere else in
## the package: a version that does not occur cannot carry mass.
hierarchical_weights <- function(spec, target, cells = spec$cells) {
  vars <- setdiff(spec$cell_vars, target)
  lab <- if (length(vars))
    do.call(paste, c(unname(lapply(vars, function(v) as.character(cells[[v]]))),
                     sep = ".")) else rep("", nrow(cells))
  lapply(split(seq_len(nrow(cells)), as.character(cells[[target]])), function(i) {
    w <- rep(1, length(i))
    for (v in vars) {
      anc <- intersect(nest_ancestors(spec, v), c(vars, target))
      key <- if (length(anc))
        do.call(paste, c(unname(lapply(anc, function(a)
          as.character(cells[[a]][i]))), sep = "\r")) else rep("", length(i))
      w <- w / ave_unique(as.character(cells[[v]][i]), key)
    }
    stats::setNames(w / sum(w), lab[i])
  })
}

ave_unique <- function(lev, key)
  as.numeric(stats::ave(seq_along(lev), key,
                        FUN = function(i) length(unique(lev[i]))))

version_key <- function(spec, target, at, versions) {
  below <- setdiff(spec$cell_vars, target)
  miss <- setdiff(below, names(at))
  if (length(miss))
    stop("`at` must name a level for every variable nested below `", target,
         "`: missing `", paste(miss, collapse = "`, `"), "`.")
  paste(as.character(unlist(at[below])), collapse = ".")
}

print.nestimand_policy <- function(x, ...) {
  cat("nestimand policy over versions of `", x$target, "`  [", x$kind, "]\n", sep = "")
  for (s in names(x$p)) {
    pv <- x$p[[s]]
    cat("  ", s, ": ", paste(sprintf("%s=%.3g", names(pv), pv), collapse = "  "),
        "\n", sep = "")
  }
  invisible(x)
}

## --- the counterfactual grid ----------------------------------------------
## Every row of the data, crossed with every realized cell: G-computation, so
## that the version distribution is not confounded with the covariate
## distribution. The observed-rows subset coincides with this only under
## balance and covariate independence.
counterfactual_grid <- function(spec, data = spec$data, policy = NULL,
                                cells = spec$cells) {
  keep_cols <- setdiff(names(data), c(spec$cell_vars, spec$cell_name))
  g <- data[rep(seq_len(nrow(data)), each = nrow(cells)), keep_cols, drop = FALSE]
  for (v in c(spec$cell_vars, spec$cell_name))
    g[[v]] <- rep(cells[[v]], times = nrow(data))
  g[[".row"]] <- rep(seq_len(nrow(data)), each = nrow(cells))
  rownames(g) <- NULL
  if (!is.null(policy)) g[[".w"]] <- policy_weights(spec, g, policy)
  g
}

policy_weights <- function(spec, grid, policy) {
  target <- policy$target
  ## the same construction as versions_of(): a version is named by the other
  ## categorical variables, whichever side of the target they sit on
  others <- setdiff(spec$cell_vars, target)
  vkey <- if (length(others))
    do.call(paste, c(unname(lapply(others, function(v) as.character(grid[[v]]))), sep = "."))
  else rep("", nrow(grid))
  strat <- as.character(grid[[target]])
  w <- numeric(nrow(grid))
  for (s in unique(strat)) {
    i <- strat == s
    pv <- policy$p[[s]]
    ## a stratum the policy does not cover - one excluded because the target
    ## does not vary there - contributes nothing rather than erroring
    w[i] <- if (is.null(pv)) 0 else as.numeric(pv[match(vkey[i], names(pv))])
  }
  w[is.na(w)] <- 0
  w
}

## --- partial identification ------------------------------------------------
## Every mixture estimand is a convex combination of the single-version
## contrasts, so their range is the identification region over all admissible
## policies. It is a companion to any across-boundary contrast, not an
## alternative to it.
policy_vertices <- function(spec, target) {
  vs <- versions_of(spec, target)
  expand.grid(vs, stringsAsFactors = FALSE)
}

## --- routes to an estimand -------------------------------------------------
## A policy says how the versions are weighted. It does not say over which rows
## the model is evaluated, and that is a separate choice.
##
##   g_computation   every observed row crossed with every realized cell, so
##                   each condition is averaged over the same covariate
##                   distribution: the population-averaged effect, and the only
##                   route with a causal reading when covariates differ across
##                   conditions.
##   cells           one row per realized cell, covariates at their means. The
##                   conditional effect at an average covariate value.
##
## On the linear predictor the two coincide whenever every design column is a
## linear function of the variables in the grid, which is the ordinary case:
## averaging the columns over the observed rows and evaluating them at the
## covariate means are then the same operation.

## --- interaction contrasts --------------------------------------------------
## A contrast of contrasts: (a1 b1 - a1 b2) - (a2 b1 - a2 b2). It crosses no
## structural boundary, since every cell it uses exists, and so takes no policy.
## The matrix is built from the table the engine returns, so its rows are in
## whatever order that engine used.
## When the engine returns grouped output - the outcome categories of an
## ordinal fit, or any other `group` column - a difference of differences is
## only interpretable within a group. The matrix is then block-diagonal: the
## same columns, once per group, labelled with it.
interaction_matrix <- function(d, vars, group = NULL) {
  d <- as.data.frame(d)
  if (is.null(group) && "group" %in% names(d) &&
      length(unique(as.character(d$group))) > 1)
    group <- "group"
  if (!is.null(group) && !is.na(group)) {
    g <- as.character(d[[group]])
    blocks <- lapply(unique(g), function(k) {
      idx <- which(g == k)
      Hk <- interaction_matrix(d[idx, , drop = FALSE], vars, group = NA_character_)
      M <- matrix(0, nrow(d), ncol(Hk),
                  dimnames = list(NULL, paste0(colnames(Hk), " | ", k)))
      M[idx, ] <- Hk
      M
    })
    return(do.call(cbind, blocks))
  }
  key <- do.call(paste, c(unname(lapply(vars, function(v) as.character(d[[v]]))),
                          sep = "\r"))
  lv <- lapply(vars, function(v) unique(as.character(d[[v]])))
  names(lv) <- vars
  ## An interaction contrast picks two levels of every variable and takes the
  ## product of the simple contrasts over the corners they define: with two
  ## variables that is a difference of differences, with three a difference of
  ## those, and so on. The corner cells must all exist, which is what rules a
  ## comparison out in a partially nested design.
  single <- vars[lengths(lv) < 2]
  if (length(single))
    stop("an interaction contrast needs two levels of every variable, and `",
         paste(single, collapse = "`, `"), "` has one here.")
  pairs <- lapply(lv, function(l) utils::combn(seq_along(l), 2, simplify = FALSE))
  ## the last variable varies fastest, so a two-variable interaction comes out
  ## in the order it always has
  grid <- expand.grid(rev(lapply(pairs, seq_along)), KEEP.OUT.ATTRS = FALSE)
  grid <- grid[, rev(seq_len(ncol(grid))), drop = FALSE]
  n_corner <- 2^length(vars)
  if (nrow(grid) > 500)
    stop("this interaction would form ", nrow(grid), " contrasts over ",
         length(vars), " variables. Name fewer targets, or `at` a level of one ",
         "of them, or supply your own `hypothesis`.")
  cols <- list(); nm <- character(0)
  for (r in seq_len(nrow(grid))) {
    sel <- lapply(seq_along(vars), function(k)
      lv[[k]][pairs[[k]][[grid[r, k]]]])          # c(earlier, later)
    corners <- expand.grid(sel, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    kc <- do.call(paste, c(unname(as.list(corners)), sep = "\r"))
    if (!all(kc %in% key)) next             # a combination that does not exist
    sgn <- apply(corners, 1, function(z)
      prod(vapply(seq_along(vars), function(k)
        if (identical(z[[k]], sel[[k]][2])) 1 else -1, 1)))
    v <- numeric(nrow(d))
    v[match(kc, key)] <- sgn
    cols[[length(cols) + 1]] <- v
    nm <- c(nm, paste(vapply(seq_along(vars), function(k)
      sprintf("(%s - %s)", sel[[k]][2], sel[[k]][1]), ""), collapse = " x "))
  }
  if (!length(cols))
    stop("no interaction contrast is available: the design realizes no set of ",
         n_corner, " conditions with two levels of every one of `",
         paste(vars, collapse = "`, `"), "`. Within a partially nested design ",
         "some such sets do not exist; `at` a level of one variable, or name ",
         "fewer targets, to ask a comparison the design can answer.")
  M <- do.call(cbind, cols)
  colnames(M) <- nm
  M
}
