## nestimand: chain-mode declarations for brms -----------------------------
## The chain parameterization keeps the original variables as predictors, so
## its coefficients read as effects rather than cell means, and `emmeans` sees
## the nesting. The price is that the design carries columns the data cannot
## inform: those for conditions that do not exist, and a further set that can
## be reconstructed from the others. brms can be told to hold each of them at
## zero. This file derives which ones, from the design matrix rather than by
## assumption, and distinguishes the two kinds - they mean different things and
## a methods section should say so.

## Which columns of a design carry no information, and why.
##   structural     - identically zero: the condition does not exist
##   identification - reconstructable from the others: a coding choice, like
##                    picking a reference level, with no effect on estimands
zero_columns <- function(formula, data, drop_intercept = TRUE) {
  X <- stats::model.matrix(formula, data)
  icept <- colnames(X) == "(Intercept)"
  Xd <- if (drop_intercept) X[, !icept, drop = FALSE] else X
  empty <- colSums(Xd != 0) == 0
  structural <- colnames(Xd)[empty]
  ## The intercept is not a candidate for holding at zero, but it does consume a
  ## dimension, so it must be present for the rank analysis. Leaving it out
  ## understates the number of constraints whenever the sentinel is not the
  ## reference level, and the model would then be left unidentified.
  Xr <- Xd[, !empty, drop = FALSE]
  Xq <- if (drop_intercept && any(icept)) cbind(X[, icept, drop = FALSE], Xr) else Xr
  q <- qr(Xq)
  kept <- setdiff(colnames(Xq)[q$pivot[seq_len(q$rank)]], "(Intercept)")
  identification <- setdiff(colnames(Xr), kept)
  list(structural = structural, identification = identification, kept = kept,
       columns = ncol(Xd), rank = qr(X)$rank)
}

## The random-effects side uses the same design, so the same examination
## settles both. Returned per grouping factor, since each bar has its own.
chain_random_zeros <- function(spec, data = NULL) {
  if (is.null(data)) data <- sentinel_first(spec)
  bars <- spec$random_original
  if (is.null(bars)) return(list())
  bl <- gsub("^\\(|\\)$", "", regmatches(bars, gregexpr("\\(([^()]*)\\)", bars))[[1]])
  out <- list()
  for (b in bl) {
    parts <- strsplit(b, "|", fixed = TRUE)[[1]]
    grp <- trimws(parts[2])
    lhs <- trimws(parts[1])
    lab <- attr(stats::terms(stats::as.formula(paste("~", lhs))), "term.labels")
    if (!any(vapply(strsplit(lab, ":"), function(v)
          any(v %in% boundary_vars(spec)), TRUE)))
      next            # does not cross the boundary: identified, nothing to declare
    ## The declarations attach to coefficients, so they exist only in the chain
    ## form: under a crossed random term no coefficient belongs to an
    ## impossible condition, and there is nothing to hold at zero. The columns
    ## are therefore taken from the chain counterpart of the declared term,
    ## which is what chain mode fits.
    covs <- lab[!vapply(strsplit(lab, ":"), function(v)
      any(v %in% spec$cell_vars), TRUE)]
    rhs <- paste(c("0", chain_terms(spec), covs), collapse = " + ")
    z <- zero_columns(stats::as.formula(paste("~", rhs)), data,
                      drop_intercept = FALSE)
    out[[grp]] <- z
  }
  out
}

chain_priors <- function(spec, regularize = "normal(0, 5)", data = NULL) {
  if (is.null(data)) data <- sentinel_first(spec)
  if (!identical(spec$fit, "brms"))
    stop("chain-mode declarations are a brms facility: they hold coefficients ",
         "at zero through the prior, which frequentist engines have no way to ",
         "express. The declared engine is `", spec$fit, "`; use the cell ",
         "parameterization there.")
  fz <- zero_columns(cell_formula(spec, "effects"), data)
  rz <- chain_random_zeros(spec, data)
  rows <- rbind(
    if (length(fz$structural))
      data.frame(part = "fixed", group = NA_character_, coef = fz$structural,
                 kind = "structural zero", row.names = NULL),
    if (length(fz$identification))
      data.frame(part = "fixed", group = NA_character_, coef = fz$identification,
                 kind = "identification constraint", row.names = NULL),
    do.call(rbind, lapply(names(rz), function(g) rbind(
      if (length(rz[[g]]$structural))
        data.frame(part = "random", group = g, coef = rz[[g]]$structural,
                   kind = "structural zero", row.names = NULL),
      if (length(rz[[g]]$identification))
        data.frame(part = "random", group = g, coef = rz[[g]]$identification,
                   kind = "identification constraint", row.names = NULL)))))
  code <- c(
    "## Coefficients the design cannot inform, derived from the design matrix.",
    "## `structural zero`: the condition does not exist.",
    "## `identification constraint`: reconstructable from the others - a coding",
    "## choice, like a reference level, which leaves every estimand unchanged.",
    sprintf('prior_block <- c('),
    if (!is.null(regularize))
      sprintf('  set_prior("%s", class = "b"),   # regularization, not identification',
              regularize),
    unlist(lapply(seq_len(nrow(rows)), function(i) {
      r <- rows[i, ]
      last <- i == nrow(rows)
      if (identical(r$part, "fixed"))
        sprintf('  set_prior("constant(0)", class = "b", coef = "%s")%s  # %s',
                r$coef, if (last) "" else ",", r$kind)
      else
        sprintf('  set_prior("constant(0)", class = "sd", group = "%s", coef = "%s")%s  # %s',
                r$group, r$coef, if (last) "" else ",", r$kind)
    })), ")")
  structure(list(spec = spec, table = rows, regularize = regularize,
                 fixed = fz, random = rz, code = code),
            class = "nestimand_chain_priors")
}

## The brms object itself, built from the same table the code shows.
chain_prior_object <- function(x) {
  if (!requireNamespace("brms", quietly = TRUE))
    stop("brms is required to build the prior block.")
  pr <- if (!is.null(x$regularize)) brms::set_prior(x$regularize, class = "b") else NULL
  for (i in seq_len(nrow(x$table))) {
    r <- x$table[i, ]
    p <- if (identical(r$part, "fixed"))
      brms::set_prior("constant(0)", class = "b", coef = r$coef)
    else brms::set_prior("constant(0)", class = "sd", group = r$group, coef = r$coef)
    pr <- if (is.null(pr)) p else pr + p
  }
  pr
}

print.nestimand_chain_priors <- function(x, ...) {
  cat("nestimand chain-mode declarations (brms)\n")
  cat("  fixed design: ", x$fixed$columns, " columns, rank ", x$fixed$rank, "\n", sep = "")
  tb <- x$table
  for (p in unique(tb$part)) {
    s <- tb[tb$part == p, ]
    cat("  ", p, ": ", nrow(s), " coefficient(s) held at zero\n", sep = "")
    for (i in seq_len(nrow(s)))
      cat(sprintf("    %-34s %s%s\n", s$coef[i], s$kind[i],
                  if (!is.na(s$group[i])) paste0("  [", s$group[i], "]") else ""))
  }
  if (!is.null(x$regularize))
    cat("  remaining coefficients: ", x$regularize, "\n", sep = "")
  cat("  show_code() prints the block as brms code.\n")
  invisible(x)
}

show_code.nestimand_chain_priors <- function(x, ...) {
  cat(paste(x$code, collapse = "\n"), "\n", sep = "")
  invisible(structure(paste(x$code, collapse = "\n"), class = "nestimand_code"))
}
