## nestimand: the weighting policy -----------------------------------------
## Every weighting scheme is a distribution p over the versions of a compound
## condition - the levels of the nested variable realized within a stratum -
## and each such p names a stochastic intervention. The estimand is therefore
## a causal contrast between two well-defined policies rather than an
## approximation to an atomic contrast.

policy_aliases <- c("equal", "proportional", "counterfactual", "hierarchical",
                    "nominated", "standardized", "within")

## Versions of the compound condition within each stratum of `target`.
## The stratum is a level of `target`; the versions are the realized
## combinations of the variables nested below it in the same family.
versions_of <- function(spec, target) {
  fam <- NULL
  for (f in spec$cat_families) if (target %in% f) fam <- f
  if (is.null(fam))
    stop("`", target, "` is not a categorical variable of any declared nesting ",
         "family (", paste(unlist(spec$cat_families), collapse = ", "), ").")
  pos <- match(target, fam)
  below <- fam[-seq_len(pos)]
  tab <- spec$cells
  key <- if (length(below)) do.call(paste, c(unname(tab[below]), sep = ".")) else
    rep("", nrow(tab))
  split(key, as.character(tab[[target]]))
}

## p as a per-stratum named vector, from an alias or a supplied distribution.
nest_policy <- function(spec, target, policy = "equal", at = NULL, data = spec$data) {
  vs <- versions_of(spec, target)
  if (is.character(policy) && length(policy) == 1L) {
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
      proportional = ,
      counterfactual = {
        key <- do.call(paste, c(unname(lapply(spec$cell_vars, function(x)
          as.character(data[[x]]))), sep = "."))
        cellcount <- table(factor(key, levels = spec$cell_levels))
        vkey <- stats::setNames(
          do.call(paste, c(unname(lapply(setdiff(spec$cell_vars, target), function(x)
            as.character(spec$cells[[x]]))), sep = ".")),
          as.character(spec$cells[[spec$cell_name]]))
        lapply(names(vs), function(s) {
          rows <- names(vkey)[as.character(spec$cells[[target]]) == s]
          cnt <- as.numeric(cellcount[rows])
          stats::setNames(cnt / sum(cnt), vs[[s]])
        }) |> stats::setNames(names(vs))
      },
      hierarchical = lapply(vs, function(v) stats::setNames(
        hierarchical_weights(v), v)),
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
hierarchical_weights <- function(v) {
  parts <- strsplit(v, ".", fixed = TRUE)
  d <- max(lengths(parts))
  w <- rep(1, length(v))
  for (k in seq_len(d)) {
    key <- vapply(parts, function(x) paste(utils::head(x, k - 1), collapse = "."), "")
    lev <- vapply(parts, function(x) paste(utils::head(x, k), collapse = "."), "")
    nsib <- ave_unique(lev, key)
    w <- w / nsib
  }
  w / sum(w)
}
ave_unique <- function(lev, key)
  as.numeric(stats::ave(seq_along(lev), key,
                        FUN = function(i) length(unique(lev[i]))))

version_key <- function(spec, target, at, versions) {
  fam <- NULL
  for (f in spec$cat_families) if (target %in% f) fam <- f
  below <- fam[-seq_len(match(target, fam))]
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
counterfactual_grid <- function(spec, data = spec$data, policy = NULL) {
  cells <- spec$cells
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
  fam <- NULL
  for (f in spec$cat_families) if (target %in% f) fam <- f
  below <- fam[-seq_len(match(target, fam))]
  vkey <- if (length(below))
    do.call(paste, c(unname(lapply(below, function(v) as.character(grid[[v]]))), sep = "."))
  else rep("", nrow(grid))
  strat <- as.character(grid[[target]])
  w <- numeric(nrow(grid))
  for (s in unique(strat)) {
    i <- strat == s
    pv <- policy$p[[s]]
    w[i] <- as.numeric(pv[match(vkey[i], names(pv))])
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
