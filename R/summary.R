## nestimand: the fit, summarized in the original parameterization ----------
## A cell fit has one coefficient per realized condition, which is what makes it
## well posed but not what a report states. The identified effect coefficients
## are an exact linear function of them - m = A^-1 mu - so the summary a reader
## expects can be produced from the cell fit without refitting anything, and
## with the same standard errors a chain fit would give.
##
## The coefficient names R derives from the design matrix are the ones Section
## 3.2 of the companion document warns about: `chord_typemaj` names a comparison
## at one inversion, not an average over them, and which one depends on the
## pivot. Each row therefore carries the combination of conditions it actually
## equals, so no label has to be read in isolation.

nest_summary <- function(model, space = c("effects", "cells"),
                         conf_level = 0.95, random = FALSE, data = NULL,
                         spec = NULL) {
  space <- match.arg(space)
  spec <- resolve_spec(model, spec)
  if (is.null(data)) data <- spec$data
  b <- coef_vector(model)
  V <- as.matrix(stats::vcov(model))
  ## A fit whose coef() returns per-group values rather than fixed effects - or
  ## any other mismatch between the coefficients and their covariance - would
  ## otherwise be caught only by an out-of-memory failure inside the matrix
  ## product, with nothing to say what went wrong.
  if (is.null(names(b)) || is.null(colnames(V)))
    stop("cannot read named coefficients and a matching covariance matrix from ",
         "this model, so the translation cannot be formed. Report the classes ",
         "involved: class(model) is ", paste(class(model), collapse = "/"),
         ", coef is ", length(b), " long, vcov is ",
         paste(dim(V), collapse = " x "), ".")
  keep <- intersect(names(b), colnames(V))
  if (length(keep) < 2)
    stop("the coefficients and the covariance matrix of this model share ",
         length(keep), " names, so they cannot be matched up. class(model) is ",
         paste(class(model), collapse = "/"), "; coef is ", length(b),
         " long, vcov is ", paste(dim(V), collapse = " x "), ".")
  if (length(b) > 5000)
    stop("this model reports ", length(b), " coefficients, which is more than a ",
         "cell parameterization produces: coef() has most likely returned ",
         "per-group values rather than fixed effects. Report class(model): ",
         paste(class(model), collapse = "/"), ".")
  b <- b[keep]; V <- V[keep, keep, drop = FALSE]
  cells <- as.character(spec$cells[[spec$cell_name]])
  cn <- spec$cell_name
  if (!any(paste0(cn, cells) %in% names(b))) {
    ## A fit in the effects parameterization needs no translation: its
    ## coefficients are already the effects this function exists to produce.
    ## That happens when the declaration asks for less than the saturated
    ## structure, which the cell factor cannot express.
    md <- attr(model, "nestimand_mode")
    if (identical(md, "effects") || identical(md, "reduced"))
      return(effects_fit_summary(model, spec, conf_level, random, space = md))
    stop("nest_summary() reads the cell coefficients, and this model has none: ",
         "it was not fitted in the cell parameterization, and does not carry ",
         "the declaration that says it was fitted as effects. Refit with ",
         "nest_fit(spec), whose coefficients this function translates, or ",
         "summarize the model with its own summary() method.")
  }
  A <- effect_basis(spec)

  ## The coefficients fall into blocks, each holding one quantity per realized
  ## cell: the cell means, and one set of slopes for every covariate that
  ## interacts with the cell factor. Each block translates by the same map, so a
  ## covariate slope becomes a reference-cell slope plus differences from it,
  ## which is what a chain fit reports.
  blocks <- list(list(label = character(0),
                      cols = intersect(c("(Intercept)", paste0(cn, cells)), names(b)),
                      order = c("(Intercept)", paste0(cn, cells))))
  ## A numeric covariate contributes one column per cell, named `cell<k>:x`. A
  ## factor covariate contributes one such set per non-reference level, named
  ## `cell<k>:x<level>` - so the column name is the variable plus a suffix, and
  ## matching the variable name exactly found none of them. Each level is its
  ## own block, since each is a separate quantity per cell. Longer covariate
  ## names are matched first, and a column already claimed is not re-used, so a
  ## covariate whose name is a prefix of another does not swallow its columns.
  claimed <- character(0)
  for (cv in spec$covariates[order(-nchar(spec$covariates))]) {
    hit <- names(b)[grepl(paste0("(^|:)", cn), names(b)) &
                    grepl(paste0("(^|:)", cv, "[^:]*($|:)"), names(b))]
    hit <- setdiff(hit, claimed)
    if (!length(hit)) next
    col_of <- vapply(hit, function(z) {
      pieces <- strsplit(z, ":", fixed = TRUE)[[1]]
      pieces <- pieces[startsWith(pieces, cv)]
      if (length(pieces)) pieces[1] else NA_character_ }, "")
    for (lv in unique(stats::na.omit(col_of))) {
      h <- hit[!is.na(col_of) & col_of == lv]
      ord <- vapply(cells, function(k) {
        m <- h[grepl(paste0(cn, k, "(:|$)"), h, fixed = FALSE)]
        if (length(m) == 1) m else NA_character_ }, "")
      if (anyNA(ord)) next                  # not one slope per cell: leave alone
      claimed <- c(claimed, unname(ord))
      blocks[[length(blocks) + 1]] <- list(label = lv, cols = unname(ord),
                                           order = unname(ord))
    }
  }
  used <- unlist(lapply(blocks, `[[`, "cols"))
  extras <- setdiff(names(b), used)

  rows <- list()
  for (bl in blocks) {
    g <- cell_grid(spec, data)
    for (cvz in spec$covariates) if (is.numeric(data[[cvz]])) g[[cvz]] <- 0
    if (!length(bl$label)) {
      rhs <- paste(deparse(cell_formula(spec)[[3]]), collapse = " ")
      X <- stats::model.matrix(stats::as.formula(paste("~", rhs)), g)
      X <- X[, intersect(colnames(X), bl$cols), drop = FALSE]
      M <- matrix(0, length(cells), length(b), dimnames = list(cells, names(b)))
      M[cells, colnames(X)] <- X
      nm <- cells
    } else {
      M <- matrix(0, length(cells), length(b), dimnames = list(cells, names(b)))
      for (i in seq_along(cells)) M[i, bl$cols[i]] <- 1
      nm <- paste0(cells, " slope on ", bl$label)
    }
    if (space == "effects") {
      M <- solve(A) %*% M
      nm <- if (!length(bl$label)) colnames(A)
            else ifelse(colnames(A) == "(Intercept)", bl$label,
                        paste0(colnames(A), ":", bl$label))
    }
    rownames(M) <- nm
    rows[[length(rows) + 1]] <- list(M = M, label = bl$label)
  }
  if (length(extras)) {
    E <- matrix(0, length(extras), length(b), dimnames = list(extras, names(b)))
    for (nm in extras) E[nm, nm] <- 1
    rows[[length(rows) + 1]] <- list(M = E, label = NA_character_)
  }
  Tm <- do.call(rbind, lapply(rows, `[[`, "M"))

  ## A Bayesian fit is summarized from its draws, as the estimand functions are.
  ## Mapping the posterior covariance through the delta method would give a
  ## normal approximation to the posterior, and a p-value computed from it has
  ## no Bayesian reading: it answers a question about repeated sampling that the
  ## model was not fitted to ask. The draws give the posterior of each row
  ## directly - its mean, its standard deviation, a quantile interval, and the
  ## probability of direction.
  if (inherits(model, "brmsfit")) {
    D <- as.matrix(brms::as_draws_matrix(model))
    nm <- draw_names(colnames(Tm), colnames(D))
    out <- draws_summary(D[, nm, drop = FALSE] %*% t(Tm), rownames(Tm), conf_level)
  } else {
    est <- as.numeric(Tm %*% b)
    se  <- sqrt(quad_form_diag(Tm, V))    # diag(T V T') without forming it
    z   <- stats::qnorm(1 - (1 - conf_level) / 2)
    out <- data.frame(term = rownames(Tm), estimate = est, std.error = se,
                      statistic = est / se,
                      p.value = 2 * stats::pnorm(-abs(est / se)),
                      conf.low = est - z * se, conf.high = est + z * se,
                      row.names = NULL)
  }
  out$meaning <- unlist(lapply(rows, function(r)
    block_meaning(spec, r$label, cells, A, space, rownames(r$M))))
  if (has_thresholds(spec) && nrow(out) && isTRUE(out$std.error[1] == 0))
    out$meaning[1] <- paste0(out$meaning[1], " (fixed at 0: absorbed into the thresholds)")
  if (isTRUE(random))
    attr(out, "nestimand_random") <-
      ## the reduced form leaves the random side in the cell parameterization
      tryCatch(random_covariance(model, spec,
                                 if (identical(space, "reduced")) "cells" else space),
               error = function(e) structure(list(), note = conditionMessage(e)))
  attr(out, "nestimand_space") <- space
  ## which engine produced the coefficients being translated: the same table
  ## can come from any of them, and the reader should not have to remember
  attr(out, "nestimand_fit") <- spec$fit
  attr(out, "nestimand_call") <- tryCatch(
    paste(deparse(stats::formula(if (inherits(stats::formula(model), "bform"))
                                   stats::formula(stats::formula(model))
                                 else stats::formula(model))), collapse = " "),
    error = function(e) NULL)
  attr(out, "nestimand_map") <- Tm
  class(out) <- c("nestimand_summary", class(out))
  out
}

## A fit whose coefficients are already effects. Nothing is translated: the
## aliased columns the chain form carries are dropped, since they are held at
## zero and say nothing, and each row is labelled with the term it belongs to.
effects_fit_summary <- function(model, spec, conf_level = 0.95, random = FALSE,
                               space = "effects") {
  b <- coef_vector(model)
  keep <- names(b)[!is.na(b)]
  V <- vcov_beta(model, keep)
  keep <- intersect(keep, colnames(V))
  b <- b[keep]; V <- V[keep, keep, drop = FALSE]
  se <- sqrt(diag(V))
  if (inherits(model, "brmsfit")) {
    D <- as.matrix(brms::as_draws_matrix(model))
    S <- D[, draw_names(keep, colnames(D)), drop = FALSE]
    out <- draws_summary(S, keep, conf_level)
  } else {
    z <- stats::qnorm(1 - (1 - conf_level) / 2)
    out <- data.frame(term = keep, estimate = as.numeric(b), std.error = se,
                      statistic = b / se, p.value = 2 * stats::pnorm(-abs(b / se)),
                      conf.low = b - z * se, conf.high = b + z * se,
                      row.names = NULL)
  }
  ## A coefficient held at zero says which kind of constraint holds it, since a
  ## row of zeros with no direction to report is otherwise unreadable.
  held <- if (identical(space, "reduced")) character(0) else tryCatch({
    tb <- chain_priors(spec)$table
    stats::setNames(tb$kind[tb$part == "fixed"], tb$coef[tb$part == "fixed"])
  }, error = function(e) character(0))
  out$meaning <- ifelse(out$term %in% names(held),
                        paste("held at zero:", held[out$term]), "as fitted")
  ## the reduced form's columns are named syntactically so that every engine
  ## carries them unaltered; the report names the effect each one stands for
  if (identical(space, "reduced"))
    out$term <- tryCatch(reduced_labels(spec, out$term),
                         error = function(e) out$term)
  attr(out, "nestimand_space") <- space
  attr(out, "nestimand_fit") <- spec$fit
  ## a brmsformula deparses to its whole object; one more formula() call reduces
  ## it to the formula a reader wants to see
  attr(out, "nestimand_call") <- tryCatch({
    fm <- stats::formula(model)
    if (inherits(fm, "bform") || inherits(fm, "brmsformula"))
      fm <- tryCatch(stats::formula(fm), error = function(e) fm)
    paste(deparse(fm), collapse = " ")
  }, error = function(e) NULL)
  if (isTRUE(random))
    attr(out, "nestimand_random") <-
      ## the reduced form leaves the random side in the cell parameterization
      tryCatch(random_covariance(model, spec,
                                 if (identical(space, "reduced")) "cells" else space),
               error = function(e) structure(list(), note = conditionMessage(e)))
  class(out) <- c("nestimand_summary", class(out))
  out
}

## One row per column of `S`, the draws of that row's quantity: the posterior
## mean, its standard deviation, a quantile interval, and the probability of
## direction - the posterior mass on whichever side of zero holds more of it. A
## row the parameterization holds at zero has no direction to report.
draws_summary <- function(S, terms, conf_level = 0.95) {
  lo <- (1 - conf_level) / 2
  sdv <- apply(S, 2, stats::sd)
  out <- data.frame(
    term = terms, estimate = apply(S, 2, mean), std.error = sdv,
    conf.low = apply(S, 2, stats::quantile, probs = lo),
    conf.high = apply(S, 2, stats::quantile, probs = 1 - lo),
    pd = ifelse(sdv == 0, NA_real_,
                apply(S, 2, function(z) max(mean(z > 0), mean(z < 0)))),
    row.names = NULL)
  attr(out, "nestimand_draws") <- S
  out
}

## What each row equals, as a combination of realized conditions. For a slope
## block the same combination applies, of slopes rather than means.
block_meaning <- function(spec, label, cells, A, space, rownames_M) {
  if (length(label) && !is.null(label) && is.na(label)) {
    ## Coefficients that are not one quantity per cell. The translation leaves
    ## them as fitted, and their usual reading is unaffected - but only a
    ## threshold may be called one: labelling whatever is left over as a
    ## threshold turns an unrecognized coefficient into a confident wrong
    ## answer, which is how a covariate crossed with the cells once read.
    is_threshold <- grepl("^Intercept\\[[0-9]+\\]$", rownames_M) |
      grepl("^[^|]+\\|[^|]+$", rownames_M)
    return(ifelse(rownames_M %in% spec$covariates, "common slope",
                  ifelse(is_threshold, "threshold", "as fitted, not translated")))
  }
  suffix <- if (is.null(label) || !length(label)) ""
            else paste0(" (slope on ", label, ")")
  if (space == "cells")
    return(paste0(cells, suffix))
  W <- solve(A); dimnames(W) <- list(colnames(A), cells)
  vapply(seq_len(nrow(W)), function(i) {
    w <- W[i, ]; w <- w[abs(w) > 1e-8]
    fmt <- function(v, nm) {
      if (isTRUE(all.equal(unname(v), 1))) nm
      else if (isTRUE(all.equal(unname(v), -1))) paste0("-", nm)
      else sprintf("%.3g*%s", v, nm)
    }
    paste0(paste(mapply(fmt, w, names(w)), collapse = " + "), suffix)
  }, "")
}

print.nestimand_summary <- function(x, digits = 4, ...) {
  space <- attr(x, "nestimand_space")
  engine <- attr(x, "nestimand_fit")
  cat("nestimand summary: ", space, " parameterization",
      if (!is.null(engine)) paste0(", fitted with ", engine_call(engine)), "\n",
      sep = "")
  if (!is.null(attr(x, "nestimand_call")))
    cat("  ", attr(x, "nestimand_call"), "\n", sep = "")
  d <- as.data.frame(x)
  d$estimate <- round(d$estimate, digits)
  d$std.error <- round(d$std.error, digits)
  ## a posterior summary reports the probability of direction, a frequentist one
  ## a p-value: the column shown is the one the fit can support
  bayes <- "pd" %in% names(d)
  if (bayes) d$pd <- round(d$pd, max(3, digits - 1))
  else d$p.value <- format.pval(d$p.value, digits = max(2, digits - 1),
                                eps = 10^-digits)
  cols <- c("term", "estimate", "std.error", if (bayes) "pd" else "p.value")
  if (identical(space, "effects")) cols <- c(cols, "meaning")
  print_aligned(d[, cols])
  rc <- attr(x, "nestimand_random")
  if (!is.null(rc)) {
    cat("\nRandom effects\n")
    if (!length(rc)) cat("  ", attr(rc, "note"), "\n", sep = "") else print(rc)
  }
  invisible(x)
}


## The function the engine name stands for, as it appears in a call.
engine_call <- function(fit) {
  switch(fit, lm = "lm()", glm = "glm()", lmer = "lme4::lmer()",
         glmer = "lme4::glmer()", clm = "ordinal::clm()",
         clmm = "ordinal::clmm()", brm = "brms::brm()", paste0(fit, "()"))
}
