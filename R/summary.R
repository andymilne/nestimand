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

nest_summary <- function(model, spec = NULL, space = c("effects", "cells"),
                         conf_level = 0.95, data = NULL) {
  space <- match.arg(space)
  spec <- resolve_spec(model, spec)
  if (is.null(data)) data <- spec$data
  b <- coef_vector(model)
  V <- stats::vcov(model)
  keep <- intersect(names(b), colnames(V))
  b <- b[keep]; V <- V[keep, keep, drop = FALSE]
  cells <- as.character(spec$cells[[spec$cell_name]])
  cn <- spec$cell_name
  A <- effect_basis(spec)

  ## The coefficients fall into blocks, each holding one quantity per realized
  ## cell: the cell means, and one set of slopes for every covariate that
  ## interacts with the cell factor. Each block translates by the same map, so a
  ## covariate slope becomes a reference-cell slope plus differences from it,
  ## which is what a chain fit reports.
  blocks <- list(list(label = character(0),
                      cols = intersect(c("(Intercept)", paste0(cn, cells)), names(b)),
                      order = c("(Intercept)", paste0(cn, cells))))
  for (cv in spec$covariates) {
    hit <- names(b)[grepl(paste0("(^|:)", cn), names(b)) &
                    grepl(paste0("(^|:)", cv, "($|:)"), names(b))]
    if (!length(hit)) next
    ord <- vapply(cells, function(k) {
      m <- hit[grepl(paste0(cn, k, "(:|$)"), hit, fixed = FALSE)]
      if (length(m) == 1) m else NA_character_ }, "")
    if (anyNA(ord)) next                    # not one slope per cell: leave alone
    blocks[[length(blocks) + 1]] <- list(label = cv, cols = unname(ord),
                                         order = unname(ord))
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

  est <- as.numeric(Tm %*% b)
  se  <- sqrt(diag(Tm %*% V %*% t(Tm)))
  z   <- stats::qnorm(1 - (1 - conf_level) / 2)
  out <- data.frame(term = rownames(Tm), estimate = est, std.error = se,
                    statistic = est / se,
                    p.value = 2 * stats::pnorm(-abs(est / se)),
                    conf.low = est - z * se, conf.high = est + z * se,
                    row.names = NULL)
  out$meaning <- unlist(lapply(rows, function(r)
    block_meaning(spec, r$label, cells, A, space, rownames(r$M))))
  if (has_thresholds(spec) && nrow(out) && isTRUE(out$std.error[1] == 0))
    out$meaning[1] <- paste0(out$meaning[1], " (fixed at 0: absorbed into the thresholds)")
  attr(out, "nestimand_space") <- space
  attr(out, "nestimand_map") <- Tm
  class(out) <- c("nestimand_summary", class(out))
  out
}

## What each row equals, as a combination of realized conditions. For a slope
## block the same combination applies, of slopes rather than means.
block_meaning <- function(spec, label, cells, A, space, rownames_M) {
  if (length(label) && !is.null(label) && is.na(label))
    ## coefficients that are not one quantity per cell: a covariate slope common
    ## to every condition, or an ordinal threshold. The translation leaves them
    ## as fitted, and their usual reading is unaffected.
    return(ifelse(rownames_M %in% spec$covariates, "common slope", "threshold"))
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
  cat("nestimand summary: ", space, " parameterization\n", sep = "")
  d <- as.data.frame(x)
  d$estimate <- round(d$estimate, digits)
  d$std.error <- round(d$std.error, digits)
  d$p.value <- format.pval(d$p.value, digits = 3, eps = 1e-4)
  cols <- c("term", "estimate", "std.error", "p.value")
  if (identical(space, "effects")) cols <- c(cols, "meaning")
  print(d[, cols], row.names = FALSE, right = FALSE)
  invisible(x)
}
