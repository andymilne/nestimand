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

  ## the map from fitted coefficients to cell means, taken from the design
  ## matrix so that it holds under either coding, ordinal thresholds included
  g <- cell_grid(spec, data)
  for (cv in spec$covariates) if (is.numeric(data[[cv]])) g[[cv]] <- 0
  rhs <- paste(deparse(cell_formula(spec)[[3]]), collapse = " ")
  X <- stats::model.matrix(stats::as.formula(paste("~", rhs)), g)
  cells <- as.character(spec$cells[[spec$cell_name]])
  ## Columns that carry a cell mean: the cell dummies themselves, plus the
  ## intercept where the coding has one. Everything else - covariate slopes, and
  ## any covariate-by-cell interaction - is a slope rather than a mean, and is
  ## reported untranslated.
  mean_cols <- intersect(c("(Intercept)", paste0(spec$cell_name, cells)), names(b))
  X <- X[, intersect(colnames(X), mean_cols), drop = FALSE]
  extras <- setdiff(names(b), mean_cols)

  Tm <- matrix(0, nrow = length(cells) + length(extras), ncol = length(b),
               dimnames = list(c(cells, extras), names(b)))
  Tm[cells, colnames(X)] <- X
  for (nm in extras) Tm[nm, nm] <- 1

  if (space == "effects") {
    A <- effect_basis(spec)
    E <- matrix(0, nrow = nrow(Tm), ncol = nrow(Tm),
                dimnames = list(c(colnames(A), extras), rownames(Tm)))
    E[colnames(A), cells] <- solve(A)
    for (nm in extras) E[nm, nm] <- 1
    Tm <- E %*% Tm
  }

  est <- as.numeric(Tm %*% b)
  se  <- sqrt(diag(Tm %*% V %*% t(Tm)))
  z   <- stats::qnorm(1 - (1 - conf_level) / 2)
  out <- data.frame(term = rownames(Tm), estimate = est, std.error = se,
                    statistic = est / se,
                    p.value = 2 * stats::pnorm(-abs(est / se)),
                    conf.low = est - z * se, conf.high = est + z * se,
                    row.names = NULL)
  out$meaning <- coefficient_meaning(spec, Tm, cells, extras, space)
  ## Under a threshold family the reference condition has no free coefficient:
  ## its latent mean is absorbed into the thresholds and reads as exactly zero.
  if (has_thresholds(spec) && nrow(out) && isTRUE(out$std.error[1] == 0))
    out$meaning[1] <- paste0(out$meaning[1], " (fixed at 0: absorbed into the thresholds)")
  attr(out, "nestimand_space") <- space
  attr(out, "nestimand_map") <- Tm
  class(out) <- c("nestimand_summary", class(out))
  out
}

## What each row equals, written as a combination of realized conditions. The
## inherited names are unreliable; this is not.
coefficient_meaning <- function(spec, Tm, cells, extras, space) {
  A <- effect_basis(spec)
  W <- if (space == "effects") solve(A) else diag(length(cells))
  dimnames(W) <- list(if (space == "effects") colnames(A) else cells, cells)
  vapply(rownames(Tm), function(r) {
    if (r %in% extras) return("covariate, untranslated")
    w <- W[r, ]
    w <- w[abs(w) > 1e-8]
    if (!length(w)) return("")
    fmt <- function(v, nm) {
      if (isTRUE(all.equal(unname(v), 1))) nm
      else if (isTRUE(all.equal(unname(v), -1))) paste0("-", nm)
      else sprintf("%.3g*%s", v, nm)
    }
    paste(mapply(fmt, w, names(w)), collapse = " + ")
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
