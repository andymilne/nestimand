## nestimand: prior translation --------------------------------------------
## Priors are stated where the user thinks - on named effects, or on cell means
## - and fitted where estimation is well posed. For an elliptical prior the
## translation is exact and basis-free: if the effects carry m ~ N(m0, D), then
## the cell means carry mu ~ N(A m0, A D A'), and Stan accepts a correlated
## prior directly, so no reparameterization or Jacobian is involved.

## Resolve a prior stated before the declaration existed. Called by nest_fit().
resolve_prior <- function(priors, spec, .env = parent.frame()) {
  if (!inherits(priors, "nestimand_prior_pending")) return(priors)
  do.call(nest_prior, c(list(spec = spec), lapply(priors$call, eval, .env)))
}

nest_prior <- function(spec, mean, sd = NULL, cov = NULL,
                       on = c("effects", "cells"),
                       family = c("normal", "student_t"), df = 3,
                       covariate_mean = 0, covariate_sd = NULL) {
  ## A prior is stated against a declaration - it needs to know what the cells
  ## are - but the declaration need not exist yet. Called without one, this
  ## records the statement and hands back something nest_fit() resolves against
  ## the declaration it builds, so a prior can be written in the same call as
  ## the fit. What is stated is unchanged either way; only when it is dimensioned
  ## differs.
  if (missing(spec) || is.null(spec)) {
    args <- as.list(match.call())[-1]
    args$spec <- NULL
    return(structure(list(call = args), class = "nestimand_prior_pending"))
  }
  on <- match.arg(on)
  family <- match.arg(family)

  A <- effect_basis(spec)
  nm <- if (on == "effects") colnames(A) else rownames(A)
  mean <- expand_named(mean, nm, "mean")
  if (is.null(cov)) {
    if (is.null(sd))
      stop("state either `sd` (independent priors) or `cov` (a joint prior).")
    sd <- expand_named(sd, nm, "sd")
    if (any(sd <= 0)) stop("every `sd` must be positive.")
    cov <- diag(sd^2, nrow = length(sd))
    dimnames(cov) <- list(nm, nm)
  } else {
    if (!identical(dim(cov), c(length(nm), length(nm))))
      stop("`cov` must be ", length(nm), " x ", length(nm), ", one row and ",
           "column per ", on, ".")
    dimnames(cov) <- list(nm, nm)
  }
  ## translate into the fitting space, and back, so both are always available
  if (on == "effects") {
    eff_mean <- mean; eff_cov <- cov
    cell_mean <- as.numeric(A %*% mean)
    cell_cov  <- A %*% cov %*% t(A)
  } else {
    Ai <- solve(A)
    cell_mean <- mean; cell_cov <- cov
    eff_mean <- as.numeric(Ai %*% mean)
    eff_cov  <- Ai %*% cov %*% t(Ai)
  }
  names(cell_mean) <- rownames(A); names(eff_mean) <- colnames(A)
  dimnames(cell_cov) <- list(rownames(A), rownames(A))
  dimnames(eff_cov)  <- list(colnames(A), colnames(A))
  ## brms applies a `class = "b"` prior to EVERY population-level coefficient,
  ## covariates included, so the statement must span the whole vector or Stan
  ## receives mismatched dimensions and the chains fail without a usable
  ## message. The covariate block is appended here, independent of the cells,
  ## in the order the design matrix produces.
  ## Under a threshold family the fitted coefficients are contrasts of cell
  ## means against the reference condition, not the means themselves: the
  ## overall level has no coefficient, being absorbed into the thresholds. The
  ## translation is still linear - b = C mu - so the prior carries across in the
  ## same way, with one qualification recorded below.
  Cmap <- prior_cell_map(spec)
  covs <- fitted_covariate_names(spec)
  mean_mean <- as.numeric(Cmap %*% cell_mean)
  mean_cov  <- Cmap %*% cell_cov %*% t(Cmap)
  mean_nm   <- rownames(Cmap)
  if (length(covs)) {
    if (is.null(covariate_sd))
      stop("the model has population-level covariate coefficient(s) `",
           paste(covs, collapse = "`, `"), "`, and a brms prior of class \"b\" ",
           "covers every one of them, so a prior stated only on the cells ",
           "would be incomplete and Stan would fail on the dimension mismatch. ",
           "State `covariate_sd` (and `covariate_mean` if not zero); one value ",
           "applies to all, or name them individually.")
    cm <- expand_named(covariate_mean, covs, "covariate_mean")
    cs <- expand_named(covariate_sd, covs, "covariate_sd")
    full_mean <- c(stats::setNames(mean_mean, mean_nm), cm)
    full_cov <- rbind(cbind(mean_cov, matrix(0, nrow(mean_cov), length(covs))),
                      cbind(matrix(0, length(covs), ncol(mean_cov)),
                            diag(cs^2, nrow = length(cs))))
    dimnames(full_cov) <- list(names(full_mean), names(full_mean))
  } else {
    full_mean <- stats::setNames(mean_mean, mean_nm)
    full_cov <- mean_cov
    dimnames(full_cov) <- list(names(full_mean), names(full_mean))
  }
  structure(list(spec = spec, stated_on = on, family = family, df = df,
                 cell_mean = cell_mean, cell_cov = cell_cov,
                 eff_mean = eff_mean, eff_cov = eff_cov, A = A,
                 full_mean = full_mean, full_cov = full_cov,
                 covariates = covs, thresholds = has_thresholds(spec)),
            class = "nestimand_prior")
}

expand_named <- function(x, nm, what) {
  if (length(x) == 1L && is.null(names(x))) return(stats::setNames(rep(x, length(nm)), nm))
  if (is.null(names(x))) {
    if (length(x) != length(nm))
      stop("`", what, "` has ", length(x), " values but there are ", length(nm),
           " parameters; supply one value, a full vector, or a named vector.")
    return(stats::setNames(as.numeric(x), nm))
  }
  if (".default" %in% names(x)) {
    fill <- x[[".default"]]
    x <- x[names(x) != ".default"]
    add <- setdiff(nm, names(x))
    x <- c(x, stats::setNames(rep(fill, length(add)), add))
  }
  miss <- setdiff(nm, names(x))
  if (length(miss))
    stop("`", what, "` names no value for `", paste(miss, collapse = "`, `"),
         "`. A prior that is silent on a parameter is not a weaker prior, it is ",
         "an undefined one; state a value for every parameter, or supply a ",
         "single value to apply to all, or a `.default = ` element to cover ",
         "the parameters not named.")
  stats::setNames(as.numeric(x[nm]), nm)
}

## The audit table: whatever the prior was stated on, this is what it implies
## for the other space, and for any contrast of interest.
prior_audit <- function(prior, contrasts = NULL) {
  tab <- rbind(
    data.frame(space = "effects", parameter = names(prior$eff_mean),
               mean = as.numeric(prior$eff_mean),
               sd = sqrt(diag(prior$eff_cov)), row.names = NULL),
    data.frame(space = "cells", parameter = names(prior$cell_mean),
               mean = as.numeric(prior$cell_mean),
               sd = sqrt(diag(prior$cell_cov)), row.names = NULL))
  if (!is.null(contrasts)) {
    cs <- do.call(rbind, lapply(names(contrasts), function(k) {
      cc <- contrasts[[k]]
      data.frame(space = "contrast", parameter = k,
                 mean = sum(cc * prior$cell_mean),
                 sd = sqrt(as.numeric(t(cc) %*% prior$cell_cov %*% cc)),
                 row.names = NULL)
    }))
    tab <- rbind(tab, cs)
  }
  attr(tab, "family") <- prior$family
  tab
}

## The prior implied for a policy contrast, which is the quantity most reports
## actually state: "these cell priors imply maj - aug ~ normal(0, 2.1)".
prior_for_estimand <- function(prior, target, policy = "equal", at = NULL,
                               contrast = "pairwise") {
  spec <- prior$spec
  pol <- if (inherits(policy, "nestimand_policy")) policy
         else nest_policy(spec, target, policy, at, spec$data)
  g <- counterfactual_grid(spec, spec$data, pol)
  levs <- levels(factor(spec$data[[target]]))
  cells <- as.character(spec$cells[[spec$cell_name]])
  M <- t(vapply(levs, function(s) {
    i <- as.character(g[[target]]) == s
    w <- tapply(g$.w[i], factor(as.character(g[[spec$cell_name]])[i], levels = cells),
                sum)
    w[is.na(w)] <- 0
    as.numeric(w / sum(w))
  }, numeric(length(cells))))
  prs <- contrast_pairs(levs, contrast)
  cl <- stats::setNames(
    lapply(prs, function(p) M[p[1], ] - M[p[2], ]),
    vapply(prs, function(p) paste(levs[p[1]], "-", levs[p[2]]), ""))
  prior_audit(prior, contrasts = cl)[
    prior_audit(prior, contrasts = cl)$space == "contrast", ]
}

print.nestimand_prior <- function(x, ...) {
  cat("nestimand prior: ", x$family, ", stated on ", x$stated_on, "\n", sep = "")
  if (isTRUE(x$thresholds))
    cat("  ordinal family: the fitted coefficients are contrasts against `",
        rownames(x$A)[1], "`, and the overall\n  level is carried by the ",
        "thresholds, so the prior on that direction is not transmitted.\n",
        sep = "")
  tab <- prior_audit(x)
  for (sp in c("effects", "cells")) {
    cat("  ", sp, ":\n", sep = "")
    s <- tab[tab$space == sp, ]
    for (i in seq_len(nrow(s)))
      cat(sprintf("    %-34s mean %8.3f   sd %7.3f\n", s$parameter[i],
                  s$mean[i], s$sd[i]))
  }
  invisible(x)
}

## --- handing the prior to brms --------------------------------------------
## Stan accepts a correlated prior directly, so the whole translation reduces
## to two data blocks and one prior statement. Nothing is reparameterized, and
## no Jacobian correction arises.
prior_stanvars <- function(prior) {
  if (!requireNamespace("brms", quietly = TRUE))
    stop("brms is required to build stanvars.")
  brms::stanvar(as.numeric(prior$full_mean), name = "prior_mean") +
    brms::stanvar(as.matrix(prior$full_cov), name = "prior_cov")
}

## From cell means to the coefficients the engine will actually estimate. With a
## zero-intercept coding they are the same thing; under a threshold family they
## are contrasts against the reference condition, and the overall level is
## carried by the thresholds instead - so a prior on that direction is not
## transmitted, and must be stated on the thresholds separately.
prior_cell_map <- function(spec) {
  cells <- as.character(spec$cells[[spec$cell_name]])
  nm <- fitted_cell_names(spec)
  if (!has_thresholds(spec)) {
    M <- diag(length(cells)); dimnames(M) <- list(nm, cells); return(M)
  }
  ref <- cells[1]
  others <- setdiff(cells, ref)
  M <- matrix(0, length(others), length(cells),
              dimnames = list(paste0(spec$cell_name, others), cells))
  for (i in seq_along(others)) {
    M[i, ref] <- -1
    M[i, others[i]] <- 1
  }
  M
}

## The population-level coefficient names the fitted model will carry, taken
## from the design matrix rather than assumed, so the prior can be checked
## against them before sampling starts.
fitted_coef_names <- function(spec, mode = "cells") {
  X <- stats::model.matrix(cell_formula(spec, mode), spec$data)
  setdiff(colnames(X), "(Intercept)")
}
## the design matrix prefixes the factor name: `cellaug.none`, not `aug.none`
fitted_cell_names <- function(spec) {
  cells <- as.character(spec$cells[[spec$cell_name]])
  if (has_thresholds(spec)) cells <- cells[-1]     # the reference has no column
  paste0(spec$cell_name, cells)
}
fitted_covariate_names <- function(spec)
  setdiff(fitted_coef_names(spec), fitted_cell_names(spec))

## Checked before the fit, because Stan reports a dimension mismatch only as
## chains finishing unexpectedly, with nothing to say why.
check_prior_dimension <- function(prior, spec, mode = "cells") {
  want <- fitted_coef_names(spec, mode)
  have <- names(prior$full_mean)
  if (length(want) != length(have))
    stop("the prior covers ", length(have), " coefficient(s) but the model has ",
         length(want), " (", paste(utils::head(want, 12), collapse = ", "),
         "). A brms prior of class \"b\" must span every population-level ",
         "coefficient; Stan reports a mismatch only as chains failing.")
  if (!identical(want, have))
    warning("the prior's coefficient order (", paste(utils::head(have, 3), collapse = ", "),
            ", ...) differs from the model's (", paste(utils::head(want, 3), collapse = ", "),
            ", ...); the prior will be applied in the model's order.",
            call. = FALSE)
  invisible(TRUE)
}

prior_statement <- function(prior) {
  txt <- if (prior$family == "normal")
    "multi_normal(prior_mean, prior_cov)"
  else sprintf("multi_student_t(%s, prior_mean, prior_cov)", prior$df)
  if (!requireNamespace("brms", quietly = TRUE)) return(txt)
  brms::set_prior(txt, class = "b")
}

## Priors that pin the fitting basis: independent non-elliptical, or bounded
## support, are declarable only where the coordinates are the effects.
requires_effect_basis <- function(prior) isTRUE(attr(prior, "requires_effect_basis"))


## --- priors written in brms syntax -----------------------------------------
## A user states every prior as they would for brm(). The `class = "b"` block
## describes the mean structure and is translated into cell space, travelling to
## Stan through stanvars rather than as a prior argument; everything else - the
## thresholds, the random-effects standard deviations and correlations, sigma -
## is passed through untouched. Which space the b priors are written in is an
## argument rather than a guess: `effects` by default, since that is the space
## the questions are asked in.
split_brms_prior <- function(prior) {
  if (is.null(prior)) return(list(b = NULL, other = NULL))
  d <- as.data.frame(prior)
  list(b = d[d$class == "b", , drop = FALSE],
       other = d[d$class != "b", , drop = FALSE])
}

## normal(m, s) and student_t(df, m, s) are elliptical and translate exactly;
## anything else does not, and is refused rather than approximated.
parse_normal_prior <- function(txt) {
  txt <- trimws(txt)
  num <- function(x) suppressWarnings(as.numeric(trimws(x)))
  if (grepl("^normal\\(", txt)) {
    a <- num(strsplit(sub("^normal\\((.*)\\)$", "\\1", txt), ",")[[1]])
    if (length(a) == 2 && !anyNA(a))
      return(list(family = "normal", mean = a[1], sd = a[2], df = 3))
  }
  if (grepl("^student_t\\(", txt)) {
    a <- num(strsplit(sub("^student_t\\((.*)\\)$", "\\1", txt), ",")[[1]])
    if (length(a) == 3 && !anyNA(a))
      return(list(family = "student_t", mean = a[2], sd = a[3], df = a[1]))
  }
  NULL
}

## The whole thing: a brms prior in, a nestimand prior out, with the rest of the
## specification carried alongside it.
prior_from_brms <- function(spec, prior, space = c("effects", "cells"),
                            covariate_mean = 0, covariate_sd = NULL) {
  space <- match.arg(space)
  sp <- split_brms_prior(prior)
  if (!nrow(sp$b)) return(list(translated = NULL, other = prior))
  if (any(nzchar(sp$b$coef)))
    stop("a `class = \"b\"` prior naming individual coefficients cannot be ",
         "translated as it stands: the translation is a joint one, so the whole ",
         "mean structure has to be described. State one prior for the class - ",
         "set_prior(\"normal(0, 1)\", class = \"b\") - or build the prior with ",
         "nest_prior(), which takes a mean and a standard deviation per ",
         "parameter and names them in the space you choose.")
  if (nrow(sp$b) > 1)
    stop("several `class = \"b\"` priors were given; the translation needs one ",
         "description of the mean structure. Combine them, or use nest_prior().")
  pr <- parse_normal_prior(sp$b$prior[1])
  if (is.null(pr))
    stop("`", sp$b$prior[1], "` is not an elliptical prior, and only those ",
         "translate exactly: a normal or student_t prior on the effects induces ",
         "a normal or student_t prior on the cells, mu ~ N(A m, A D A'). For a ",
         "prior of another shape, fit in the chain parameterization, where the ",
         "coefficients are the effects themselves: nest_fit(spec, ",
         "mode = \"effects\", priors = chain_priors(spec)).")
  np <- nest_prior(spec, mean = pr$mean, sd = pr$sd, on = space,
                   family = pr$family, df = pr$df,
                   covariate_mean = covariate_mean,
                   covariate_sd = if (is.null(covariate_sd)) pr$sd else covariate_sd)
  list(translated = np,
       other = if (nrow(sp$other)) do.call(rbind, list(sp$other)) else NULL,
       stated = sp$b$prior[1], space = space)
}
