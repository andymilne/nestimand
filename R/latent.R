## nestimand: estimands as linear functionals of the cell coefficients -----
## On a linear scale - the latent scale of an ordinal model, the link scale of
## a generalized linear model - a policy contrast is exactly c'b, with c a
## weighted difference of design-matrix rows. This is what the translation
## layer is for: the map from cells to named effects is a matrix, so the
## estimand needs no prediction machinery, is immune to the per-category
## expansion that ordinal fits provoke, and costs one matrix product.

## The contrast vector c for a policy, formed on the counterfactual grid so
## that the covariate distribution is common to every stratum. On a linear
## scale, averaging predictions over that grid and averaging its design-matrix
## rows are the same operation, so this is G-computation, not an approximation.
policy_contrast_matrix <- function(spec, target, policy, data = spec$data,
                                   model = NULL) {
  g <- counterfactual_grid(spec, data, policy)
  rhs <- paste(deparse(cell_formula(spec)[[3]]), collapse = " ")
  X <- stats::model.matrix(stats::as.formula(paste("~", rhs)), g)
  if (colnames(X)[1] == "(Intercept)") X <- X[, -1, drop = FALSE]
  if (!is.null(model)) {
    b <- coef_vector(model)
    missing_cols <- setdiff(colnames(X), names(b))
    if (length(missing_cols))
      stop("the fitted model has no coefficient for ",
           length(missing_cols), " column(s) of the cell design matrix (",
           paste(utils::head(missing_cols, 3), collapse = ", "),
           if (length(missing_cols) > 3) ", ...", "), so the linear map cannot ",
           "be formed. This happens when the model was not fitted from this ",
           "spec; refit with nest_fit().")
    X <- X[, colnames(X), drop = FALSE]
  }
  w <- g$.w / sum(g$.w[!duplicated(g$.row)]) # per-row weights, normalized below
  lev <- as.character(g[[target]])
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

coef_vector <- function(model) {
  if (inherits(model, "clm") || inherits(model, "clmm")) {
    b <- model$beta
    if (is.null(b)) b <- coef(model)
    return(b)
  }
  if (inherits(model, "merMod")) return(lme4::fixef(model))
  if (inherits(model, "brmsfit")) {
    b <- brms::fixef(model)[, "Estimate"]
    names(b) <- rownames(brms::fixef(model))
    return(b)
  }
  stats::coef(model)
}

vcov_beta <- function(model, nm) {
  V <- stats::vcov(model)
  if (inherits(model, "brmsfit")) {
    dn <- colnames(V)
    V <- V[dn %in% nm, dn %in% nm, drop = FALSE]
    return(V)
  }
  keep <- intersect(nm, colnames(V))
  V[keep, keep, drop = FALSE]
}

## Pairwise, reference, or sequential differences of the stratum rows.
contrast_pairs <- function(levs, contrast = "pairwise") {
  switch(contrast,
    pairwise = { ij <- utils::combn(seq_along(levs), 2)
                 lapply(seq_len(ncol(ij)), function(k) ij[, k]) },
    reference = lapply(seq_along(levs)[-1], function(j) c(1L, j)),
    sequential = lapply(seq_along(levs)[-1], function(j) c(j - 1L, j)),
    stop("contrast must be pairwise, reference, or sequential here."))
}

latent_estimand <- function(model, spec, target, policy = "equal", at = NULL,
                            contrast = "pairwise", data = spec$data,
                            conf_level = 0.95) {
  pol <- if (inherits(policy, "nestimand_policy")) policy
         else nest_policy(spec, target, policy, at, data)
  M <- policy_contrast_matrix(spec, target, pol, data, model)
  b <- coef_vector(model)[colnames(M)]
  V <- vcov_beta(model, colnames(M))
  levs <- rownames(M)
  prs <- contrast_pairs(levs, contrast)
  cvecs <- lapply(prs, function(p) M[p[1], ] - M[p[2], ])
  est <- vapply(cvecs, function(cc) sum(cc * b), 1)
  se  <- vapply(cvecs, function(cc)
    sqrt(as.numeric(t(cc) %*% V[colnames(M), colnames(M)] %*% cc)), 1)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  out <- data.frame(
    term = vapply(prs, function(p) paste(levs[p[1]], "-", levs[p[2]]), ""),
    estimate = est, std.error = se, statistic = est / se,
    p.value = 2 * stats::pnorm(-abs(est / se)),
    conf.low = est - z * se, conf.high = est + z * se,
    row.names = NULL)
  attr(out, "nestimand_scale") <- "linear predictor (latent / link)"
  attr(out, "nestimand_cvecs") <- stats::setNames(cvecs, out$term)
  out
}

## The Bayesian counterpart: the same c, applied draw by draw, so the result is
## a posterior rather than a delta-method interval.
latent_draws <- function(model, spec, target, policy = "equal", at = NULL,
                         contrast = "pairwise", data = spec$data) {
  if (!inherits(model, "brmsfit"))
    stop("draw-wise translation needs a posterior; this is a frequentist fit. ",
         "Use latent_estimand() for the delta-method interval.")
  pol <- if (inherits(policy, "nestimand_policy")) policy
         else nest_policy(spec, target, policy, at, data)
  M <- policy_contrast_matrix(spec, target, pol, data, model)
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
