## nestimand: the random-effects covariance, translated -----------------------
## A cell random effect carries one deviation per realized condition, so its
## covariance is already in the design's own terms - but labelled by the fitted
## factor rather than by the conditions, and printed by the engines in a form
## that truncates or wraps once the structure is large. This reads it out in
## full, in the original labels, and translates it into effect space where that
## is meaningful.

random_covariance <- function(model, spec = NULL,
                              space = c("cells", "effects")) {
  space <- match.arg(space)
  spec <- resolve_spec(model, spec)
  cn <- spec$cell_name
  cells <- as.character(spec$cells[[cn]])

  vc <- if (methods::is(model, "merMod") || inherits(model, "merMod"))
    lme4::VarCorr(model)
  else if (inherits(model, "brmsfit")) brms::VarCorr(model)
  else if (inherits(model, "clmm")) ordinal::VarCorr(model)
  else stop("no random-effects covariance can be read from a model of class ",
            paste(class(model), collapse = "/"), ". This is for mixed fits.")

  out <- list()
  for (g in names(vc)) {
    S <- if (inherits(model, "brmsfit")) {
      sd <- vc[[g]]$sd[, "Estimate"]
      cr <- if (!is.null(vc[[g]]$cor)) vc[[g]]$cor[, "Estimate", ] else diag(length(sd))
      cr <- matrix(cr, length(sd), length(sd),
                   dimnames = list(names(sd), names(sd)))
      diag(sd) %*% cr %*% diag(sd)
    } else as.matrix(vc[[g]])
    nm <- rownames(S)
    if (is.null(nm)) next
    ## label by the conditions rather than by the fitted factor
    plain <- sub(paste0("^", cn), "", nm)
    is_cells <- setequal(plain, cells)
    dimnames(S) <- list(plain, plain)
    if (identical(space, "effects")) {
      if (!is_cells) {
        out[[g]] <- structure(S, translated = FALSE,
          note = paste("this term is not a covariance over the realized cells,",
                       "so it has no effect-space counterpart; shown as fitted"))
        next
      }
      A <- effect_basis(spec)
      S <- S[rownames(A), rownames(A), drop = FALSE]
      Ai <- solve(A)
      S <- Ai %*% S %*% t(Ai)
      dimnames(S) <- list(colnames(A), colnames(A))
    } else if (is_cells) {
      S <- S[cells, cells, drop = FALSE]
    }
    out[[g]] <- structure(S, translated = is_cells)
  }
  structure(out, class = "nestimand_random", space = space)
}

## The heterogeneity of a named effect: c'Sigma c, with c the same weight vector
## that defines the corresponding fixed-effect contrast. This is the quantity a
## report needs, and it is not a coefficient of the fitted model.
random_heterogeneity <- function(model, target, policy = "equal", at = NULL,
                                 contrast = "pairwise", spec = NULL) {
  spec <- resolve_spec(model, spec)
  cv <- attr(latent_estimand(model, target, policy, at, contrast, spec = spec),
             "nestimand_cvecs")
  S <- random_covariance(model, spec, "cells")
  cells <- as.character(spec$cells[[spec$cell_name]])
  do.call(rbind, lapply(names(S), function(g) {
    M <- S[[g]]
    if (!isTRUE(attr(M, "translated"))) return(NULL)
    data.frame(group = g, term = names(cv),
               sd = vapply(cv, function(cc) {
                 ## the contrast vector is named by design-matrix columns; the
                 ## covariance by conditions
                 names(cc) <- sub(paste0("^", spec$cell_name), "", names(cc))
                 w <- cc[match(cells, names(cc))]
                 w[is.na(w)] <- 0
                 sqrt(as.numeric(t(w) %*% M[cells, cells] %*% w)) }, 1),
               row.names = NULL)
  }))
}

print.nestimand_random <- function(x, digits = 3, ...) {
  sp <- attr(x, "space")
  for (g in names(x)) {
    S <- x[[g]]
    cat("Group `", g, "`  (", sp, " parameterization, ", nrow(S),
        " dimensions)\n", sep = "")
    if (!is.null(attr(S, "note"))) cat("  ", attr(S, "note"), "\n", sep = "")
    sd <- sqrt(diag(S))
    cr <- stats::cov2cor(S)
    tab <- cbind(`Std.Dev.` = round(sd, digits), round(cr, digits))
    colnames(tab)[-1] <- abbreviate(colnames(S), minlength = 6)
    print(tab, quote = FALSE, right = TRUE)
    cat("\n")
  }
  invisible(x)
}
