## nestimand: estimands, with the code that produced them ------------------
## The function computes; the code view is a by-product, not a reconstruction.
## Every estimand is assembled as text, then evaluated, so `show_code()` cannot
## drift from what was run: it is the same object. Non-core arguments in `...`
## are deparsed into that text and so reach the destination function - brms,
## marginaleffects, emmeans - unaltered, and appear in the saved code.

estimand <- function(model, target, policy = "equal", at = NULL,
                     contrast = c("pairwise", "reference", "sequential", "within"),
                     route = c("g_computation", "observed", "cells"),
                     weights = NULL,
                     scale = c("response", "latent"),
                     data = NULL, bounds = TRUE, self_check = TRUE,
                     dry_run = FALSE, ..., spec = NULL, .env = parent.frame()) {
  ## The declaration travels with a fit from nest_fit(); `spec` is needed only
  ## for a model fitted by calling the engine directly.
  spec_expr <- substitute(spec)
  recovered <- is.null(spec)
  spec <- resolve_spec(model, spec)
  check_model_spec(model, spec)
  contrast <- match.arg(contrast)
  route <- match.arg(route)
  wq <- substitute(weights)
  scale <- match.arg(scale)
  if (identical(scale, "latent") && identical(contrast, "within"))
    stop("contrast = \"within\" is computed through the prediction route and ",
         "has no linear-map form here; use scale = \"response\".")
  if (!inherits(spec, "nesting_spec"))
    stop("`spec` must be a nesting_spec object, as returned by nesting_spec().")
  ## A bare name is read as the variable itself - estimand(m, chord_type) - but
  ## a name that holds the target as a string is read for its value, so that
  ## the target can be supplied programmatically.
  tg <- substitute(target)
  if (is.name(tg)) {
    lit <- deparse(tg)
    target <- if (lit %in% spec$cell_vars) lit else
      tryCatch({ v <- eval(tg, .env)
                 if (is.character(v) && length(v) == 1L) v else lit },
               error = function(e) lit)
  }
  if (!target %in% spec$cell_vars)
    stop("`", target, "` is not one of the declared categorical nesting ",
         "variables (", paste(spec$cell_vars, collapse = ", "), "). Contrasts ",
         "of a covariate do not cross the structural boundary and need no ",
         "policy; compute them directly.")

  ## A marginal contrast of a nested variable is computed only over the strata
  ## in which that variable varies. Elsewhere it holds a single level, and a
  ## comparison against that level leaves the variable's own levels and compares
  ## strata instead - augmented chords against root-position triads, say, under
  ## an inversion label.
  deg <- if (!identical(contrast, "within")) degenerate_strata(spec, target)
  restricted <- !is.null(deg) && length(deg$drop)
  if (!restricted) deg <- NULL

  ## The model and the declaration may be supplied as expressions rather than
  ## as named objects - estimand(nest_fit(sp), ...) - and an expression cannot
  ## be assigned to in the emitted code.
  syntactic <- function(x) grepl("^[.A-Za-z][.A-Za-z0-9._]*$", x)
  model_name <- deparse(substitute(model))
  if (!syntactic(model_name)) model_name <- "model"
  ## When the declaration came from the fit, refer to it by the name it had
  ## when the model was fitted, so the emitted code reads as the user wrote it.
  spec_name  <- if (recovered)
    (if (is.null(attr(model, "nestimand_spec_name"))) "spec"
     else attr(model, "nestimand_spec_name")) else paste(deparse(spec_expr), collapse = "")
  if (!syntactic(spec_name)) spec_name <- "spec"
  data_name <- if (missing(data)) paste0(spec_name, "$data")
               else paste(deparse(substitute(data)), collapse = " ")
  if (is.null(data)) { data <- spec$data; data_name <- paste0(spec_name, "$data") }

  ## Unit weights standardize to a population other than the sample: survey
  ## weights, or post-stratification to a target distribution of the covariates.
  ## They are a separate question from the policy, which weights the versions of
  ## a condition, and they multiply it.
  weights_txt <- NULL
  if (!is.null(weights)) {
    if (identical(route, "cells"))
      stop("unit weights standardize over the rows of the data, and ",
           "route = \"cells\" evaluates one row per condition with the ",
           "covariates at their means, so there are no units to weight. Use ",
           "route = \"g_computation\" to standardize to a target population.")
    if (is.character(weights) && length(weights) == 1L) {
      if (!weights %in% names(data))
        stop("no column `", weights, "` in the data the estimand is computed ",
             "over. A weight column added after the model was fitted is not in ",
             "the declaration the fit carries; add it before declaring, pass ",
             "the data with `data =`, or supply the weights as a vector.")
      weights_txt <- sprintf('%s[["%s"]]', data_name, weights)
    } else {
      if (is.numeric(weights) && length(weights) != nrow(data))
        stop("`weights` has ", length(weights), " values but the data has ",
             nrow(data), " rows: one weight per row is needed.")
      weights_txt <- paste(deparse(wq), collapse = " ")
    }
  }

  ## non-core arguments, passed through verbatim to the destination function
  dots <- as.list(substitute(list(...)))[-1]
  ## A mixed fit predicts for a particular group unless told otherwise, and the
  ## estimands here are population-level. The exclusion is stated rather than
  ## assumed, and in the spelling the engine expects: brms and the frequentist
  ## engines disagree, and the wrong one is silently ignored.
  re_arg <- if (length(grouping_vars(spec)) &&
                spec$fit %in% c("lmer", "glmer", "clmm", "brms"))
    (if (identical(spec$fit, "brms")) "re_formula" else "re.form")
  re_note <- NULL
  if (!is.null(re_arg) && !re_arg %in% names(dots) && !identical(scale, "latent")) {
    dots[[re_arg]] <- quote(NA)
    re_note <- c(
      sprintf("## %s = NA: the estimand is population-level, describing a typical", re_arg),
      "## group rather than any sampled one. brms and the frequentist engines",
      "## spell this differently, and the wrong spelling is silently ignored.")
  }
  if (has_thresholds(spec) && identical(scale, "response") && !"type" %in% names(dots))
    stop("`", spec$fit, "` with an ordinal family has no single response scale, ",
         "so the scale of the contrast must be stated rather than assumed. On ",
         "the latent scale the contrast is one number per comparison, parallel ",
         "to a linear analysis; on the response scale it is one number per ",
         "outcome category, and the computation is correspondingly large. Pass ",
         "Either take the latent scale through the translation layer - ",
         "scale = \"latent\", which is exact, costs one matrix product, and ",
         "returns one number per contrast - or pass a scale through to the ",
         "engine: type = \"linear.predictor\" or ",
         "type = \"prob\" for clm and clmm, type = \"link\" or ",
         "type = \"response\" for brms. Note the support boundary: on ordinal ",
         "fits marginaleffects returns per-category output even on the latent ",
         "scale, and may refuse the pairwise comparison as too large.")
  dots_txt <- if (length(dots))
    paste0(", ", paste(mapply(function(nm, v) {
      v <- paste(deparse(v), collapse = " ")
      if (nzchar(nm)) paste0(nm, " = ", v) else v
    }, names(dots), dots), collapse = ", ")) else ""

  code <- estimand_code(spec, target, policy, at, contrast, dots_txt,
                        model_name, spec_name, data_name, bounds, scale,
                        deg, route, weights_txt, re_note)
  ## if the model was fitted by nest_fit(), its call travels with it, so the
  ## code view is the whole pipeline rather than its second half
  fit_code <- attr(model, "nestimand_code")
  if (!is.null(fit_code)) {
    fit_code <- sub("^m <- ", paste0(model_name, " <- "), fit_code)
    code <- c(fit_code, "", code)
  }
  ## `show_code(estimand(...))` computes first and prints after, because the
  ## inner call is evaluated first: show_code() is a printer, not a preview.
  ## dry_run returns the code without running anything, as nest_fit() does.
  if (isTRUE(dry_run))
    return(structure(paste(code, collapse = "\n"), class = "nestimand_code",
                     nestimand_code = code))
  env <- new.env(parent = .env)
  assign(model_name, model, envir = env)
  assign(spec_name,  spec,  envir = env)
  if (!exists(data_name, envir = env, inherits = TRUE))
    assign(data_name, data, envir = env)
  out <- eval(parse(text = paste(code, collapse = "\n")), envir = env)

  check <- if (isTRUE(self_check))
    reorder_check(model, spec, target, policy, at, contrast, dots_txt, data, scale,
                  route)
  else NULL

  structure(out,
            nestimand = list(build = nestimand_build, target = target,
                             policy = policy, at = at, contrast = contrast,
                             scale = scale, route = route,
                             bounds = attr(out, "nestimand_bounds"),
                             self_check = check,
                             code = code),
            class = c("nestimand_estimand", class(out)))
}

## --- code assembly ---------------------------------------------------------

estimand_code <- function(spec, target, policy, at, contrast, dots_txt,
                          model_name, spec_name, data_name, bounds,
                          scale = "response", deg = NULL,
                          route = "g_computation", weights_txt = NULL,
                          re_note = NULL) {
  cn <- spec$cell_name
  pol_txt <- if (is.character(policy) && length(policy) == 1L)
    sprintf('"%s"', policy)
  else paste0("c(", paste(sprintf('"%s" = %s', names(policy), policy),
                          collapse = ", "), ")")
  at_txt <- if (is.null(at)) "" else
    paste0(", at = c(", paste(sprintf('%s = "%s"', names(at), at),
                              collapse = ", "), ")")
  hdr <- c(sprintf("## nestimand %s -- estimand of `%s`", nestimand_build, target),
           sprintf("library(marginaleffects)"))
  if (contrast == "within") {
    fam <- NULL
    for (f in spec$cat_families) if (target %in% f) fam <- f
    pos <- match(target, fam)
    if (pos == 1)
      stop("contrast = \"within\" emits per-stratum contrasts of a *nested* ",
           "variable; `", target, "` is not nested within anything, so every ",
           "contrast of it crosses the structural boundary and requires a ",
           "policy. Use policy = instead.")
    parents <- fam[seq_len(pos - 1)]
    return(c(hdr,
      "## within-stratum contrasts: no boundary is crossed, so no policy applies",
      sprintf('levs <- levels(factor(%s$%s))', data_name, target),
      sprintf('parents <- c(%s)', paste(sprintf('"%s"', parents), collapse = ", ")),
      sprintf('parts <- split(%s, interaction(%s[, parents, drop = FALSE], drop = TRUE))',
              data_name, data_name),
      'est <- do.call(rbind, lapply(names(parts), function(s) {',
      '  d_s <- parts[[s]]',
      sprintf('  if (length(unique(d_s$%s)) < 2) return(NULL)  # degenerate stratum', target),
      sprintf('  cbind(stratum = s, mfx_canonical(as.data.frame('),
      sprintf('    avg_predictions(%s, newdata = d_s, by = "%s",', model_name, target),
      sprintf('                    hypothesis = %s%s)), levs))',
              mfx_hypothesis_txt("pairwise"), dots_txt),
      '}))', 'est'))
  }
  if (identical(scale, "latent")) {
    body <- c(hdr[1],
      "## latent-scale estimand: on a linear scale the policy contrast is c'b,",
      "## with c a weighted difference of design-matrix rows. Exact, one matrix",
      "## product, and free of the per-category expansion of ordinal fits.",
      sprintf('pol  <- nest_policy(%s, "%s", %s%s)', spec_name, target, pol_txt, at_txt),
      sprintf('est  <- latent_estimand(%s, "%s", pol, contrast = "%s", spec = %s%s)',
              model_name, target, contrast, spec_name, dots_txt))
    if (isTRUE(bounds))
      body <- c(body,
        "## partial-identification bounds over all admissible policies",
        sprintf('est  <- add_bounds(est, %s, %s, "%s", "%s", scale = "latent"%s%s)',
                model_name, spec_name, target, contrast,
                if (is.null(deg)) "" else ", cells = cells", dots_txt))
    return(c(body, "est"))
  }
  cells_txt <- if (is.null(deg)) sprintf("%s$cells", spec_name) else
    sprintf('subset(%s$cells, %s %%in%% c(%s))', spec_name, deg$vars[1],
            paste(sprintf('"%s"', deg$keep), collapse = ", "))
  restrict_note <- if (is.null(deg)) NULL else c(
    sprintf("## `%s` does not vary in every stratum. The contrast is pooled over", target),
    sprintf("## the strata in which it does - %s - since elsewhere a comparison",
            paste(deg$keep, collapse = ", ")),
    "## would leave its own levels and compare strata instead.",
    sprintf("cells <- %s", cells_txt))
  body <- c(hdr, restrict_note,
    sprintf("## the policy: a distribution over the versions of each compound condition"),
    sprintf('pol  <- nest_policy(%s, "%s", %s%s%s)', spec_name, target, pol_txt,
            at_txt, if (is.null(deg)) "" else ", cells = cells"),
    switch(route,
      g_computation = c(
        "## G-computation grid: every row of the data crossed with every realized",
        "## cell, so each condition is averaged over the same covariate distribution",
        sprintf('grid <- counterfactual_grid(%s, %s, pol%s)', spec_name, data_name,
                if (is.null(deg)) "" else ", cells = cells"),
        if (!is.null(weights_txt)) c(
          "## unit weights: the estimand is standardized to the population those",
          "## weights describe, rather than to the sample as observed",
          sprintf('grid$.w <- grid$.w * (%s)[grid$.row]', weights_txt))),
      observed = c(
        "## the observed rows as they stand, reweighted so that the versions enter",
        "## at the policy's proportions: no prediction is made for a row never run",
        sprintf('grid <- %s', data_name),
        sprintf('grid$.w <- observed_weights(%s, grid, pol)', spec_name),
        if (!is.null(deg))
          c("## rows in the excluded strata carry no weight, and are dropped so",
            "## that no empty group reaches the engine",
            'grid <- grid[grid$.w > 0, ]'),
        if (!is.null(weights_txt)) c(
          "## unit weights: standardized to the population they describe",
          sprintf('grid$.w <- grid$.w * (%s)', weights_txt))),
      cells = c(
        "## one row per realized cell, covariates at their means: the conditional",
        "## effect at an average covariate value",
        sprintf('grid <- cell_grid(%s)', spec_name),
        if (!is.null(deg))
          sprintf('grid <- grid[as.character(grid$%s) %%in%% as.character(cells$%s), ]',
                  cn, cn),
        sprintf('grid$.w <- policy_weights(%s, grid, pol)', spec_name))),
    "## estimand in the original variable space: `by =` names an original",
    sprintf("## factor, not the fitted `%s` predictor", cn), re_note,
    sprintf('est  <- avg_predictions(%s, newdata = grid, by = "%s", wts = grid$.w,',
            model_name, target),
    sprintf('                        hypothesis = %s%s)',
            mfx_hypothesis_txt(contrast), dots_txt),
    "## contrast direction is fixed by nestimand, not inherited from the engine:",
    "## declared factor-level order, earlier level minus later. Engine versions",
    "## differ on this, and an inherited convention would flip reported signs.",
    sprintf('est  <- mfx_canonical(est, levels(factor(%s$%s)))',
            data_name, target))
  if (isTRUE(bounds))
    body <- c(body,
      "## partial-identification bounds: every mixture estimand is a convex",
      "## combination of the single-version contrasts, so their range is the",
      "## region over all admissible policies (Manski 1990)",
      sprintf('est  <- add_bounds(est, %s, %s, "%s", "%s"%s%s)',
              model_name, spec_name, target, contrast,
              if (is.null(deg)) "" else ", cells = cells", dots_txt))
  c(body, "est")
}

## --- the bounds companion --------------------------------------------------

add_bounds <- function(est, model, spec, target, contrast = "pairwise",
                       scale = c("response", "latent"), cells = spec$cells, ...) {
  scale <- match.arg(scale)
  vs <- versions_of(spec, target, cells)
  vert <- expand.grid(lapply(vs, seq_along), KEEP.OUT.ATTRS = FALSE)
  if (nrow(vert) > 64) {
    attr(est, "nestimand_bounds") <- NULL
    message("partial-identification bounds not computed: the design admits ",
            nrow(vert), " single-version policies, beyond the enumeration ",
            "limit of 64. Request specific vertices with policy = \"nominated\".")
    return(est)
  }
  vals <- lapply(seq_len(nrow(vert)), function(i) {
    p <- structure(list(kind = "vertex", target = target,
      p = stats::setNames(lapply(names(vs), function(s)
        stats::setNames(as.numeric(seq_along(vs[[s]]) == vert[i, s]), vs[[s]])),
        names(vs)), at = NULL), class = "nestimand_policy")
    if (identical(scale, "latent"))
      return(latent_estimand(model, target, p, contrast = contrast, spec = spec)$estimate)
    g <- counterfactual_grid(spec, spec$data, p, cells = cells)
    e <- marginaleffects::avg_predictions(model, newdata = g, by = target,
           wts = g$.w, hypothesis = mfx_hypothesis(contrast), ...)
    mfx_canonical(as.data.frame(e), levels(factor(spec$data[[target]])))$estimate
  })
  M <- do.call(cbind, vals)
  ## `est` has already been put in declared level order; canonicalizing it a
  ## second time, without those levels, could relabel it against a different
  ## convention and pair each contrast with the wrong bounds.
  est_df <- as.data.frame(est)
  terms <- if ("term" %in% names(est_df)) est_df$term else
    mfx_canonical(est_df, levels(factor(spec$data[[target]])))$term
  attr(est, "nestimand_bounds") <- data.frame(
    term = terms,
    estimate = as.data.frame(est)$estimate,
    policy_low = apply(M, 1, min), policy_high = apply(M, 1, max))
  est
}

## --- the reorder self-check ------------------------------------------------
## Order instability is impossible under the cell parameterization, so this is
## a belt-and-braces check on the translation layer rather than on the fit.
reorder_check <- function(model, spec, target, policy, at, contrast, dots_txt, data,
                          scale = "response", route = "g_computation") {
  ## Order instability is a fixed-effects phenomenon: it arises from which
  ## columns a pivot drops, not from the random structure. Refitting a mixed
  ## model to test for it is therefore both expensive and unreliable - a large
  ## random structure can settle on a different optimum, and the estimands then
  ## differ by an amount that has nothing to do with level order. The check
  ## accordingly runs on the fixed-effects shadow model: the same cell formula
  ## fitted without the random terms.
  shadow <- switch(spec$fit,
    lmer = "lm", glmer = "glm", clmm = "clm",
    brms = if (has_thresholds(spec)) "clm" else if (is.null(spec$family)) "lm" else "glm",
    NULL)
  fit_shadow <- function(spx, dta) {
    f <- cell_formula(spx)
    switch(shadow,
      lm  = stats::lm(f, data = dta),
      glm = stats::glm(f, data = dta,
                       family = eval(parse(text = spx$family %||% "gaussian"))),
      clm = ordinal::clm(f, data = dta))
  }
  d2 <- data
  d2[[spec$cell_name]] <- NULL   # rebuilt from the permuted factors below
  for (v in spec$cell_vars) {
    lv <- levels(factor(d2[[v]]))
    d2[[v]] <- factor(d2[[v]], levels = c(lv[1], sample(lv[-1])))
  }
  out <- tryCatch({
    sp2 <- nesting_spec_quiet(spec, d2)
    if (is.null(shadow)) {
      m1 <- model; m2 <- stats::update(model, data = sp2$data)
    } else {
      m1 <- fit_shadow(spec, spec$data); m2 <- fit_shadow(sp2, sp2$data)
    }
    e1 <- unname(estimand_values(m1, spec, target, policy, at, contrast,
                                 spec$data, scale, route))
    e2 <- unname(estimand_values(m2, sp2, target, policy, at, contrast,
                                 sp2$data, scale, route))
    note <- if (is.null(shadow)) "estimand unchanged under level permutation"
            else paste0("estimand unchanged under level permutation (checked on ",
                        "the fixed-effects shadow model, since order instability ",
                        "is a fixed-effects phenomenon)")
    if (length(e1) != length(e2))
      list(status = "failed", note = "the two runs returned different numbers of contrasts")
    else if (!isTRUE(all.equal(sort(round(abs(e1), 8)), sort(round(abs(e2), 8)),
                               tolerance = 1e-8)))
      list(status = "failed",
           note = sprintf("max |change| = %s", format(max(abs(sort(abs(e1)) - sort(abs(e2)))))))
    else list(status = "passed", note = note)
  }, error = function(e) list(status = "error", note = conditionMessage(e)))
  if (identical(out$status, "failed"))
    warning("reorder self-check FAILED (", out$note, "): the estimate moved when ",
            "nested-factor levels were permuted, so this estimand depends on ",
            "unrealized-cell predictions and is not identified by the design. ",
            "Do not report it.", call. = FALSE)
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a

nesting_spec_quiet <- function(spec, data) {
  suppressMessages(nesting_spec(data, spec$formula_in,
    nests = spec_nests(spec), fit = spec$fit, family = spec$family,
    random = spec$random_original, cell_name = spec$cell_name))
}
spec_nests <- function(spec)
  unlist(lapply(spec$cat_families, function(f)
    if (length(f) > 1) sprintf("%s %%in%% %s", f[-1], f[-length(f)]) else character(0)))

estimand_values <- function(model, spec, target, policy, at, contrast, data,
                            scale = "response", route = "g_computation") {
  if (contrast == "within") {
    fam <- NULL; for (f in spec$cat_families) if (target %in% f) fam <- f
    parents <- fam[seq_len(match(target, fam) - 1)]
    parts <- split(data, interaction(data[, parents, drop = FALSE], drop = TRUE))
    return(unlist(lapply(parts, function(d_s) {
      if (length(unique(d_s[[target]])) < 2) return(NULL)
      mfx_canonical(as.data.frame(marginaleffects::avg_predictions(model,
        newdata = d_s, by = target, hypothesis = mfx_hypothesis("pairwise"))),
        levels(factor(data[[target]])))$estimate
    })))
  }
  deg <- if (!identical(contrast, "within")) degenerate_strata(spec, target)
  cells <- if (is.null(deg) || !length(deg$drop)) spec$cells else
    spec$cells[as.character(spec$cells[[deg$vars[1]]]) %in% deg$keep, , drop = FALSE]
  pol <- nest_policy(spec, target, policy, at, data, cells = cells)
  if (identical(scale, "latent"))
    return(latent_estimand(model, target, pol, contrast = contrast,
                           data = data, spec = spec)$estimate)
  g <- switch(route,
    g_computation = counterfactual_grid(spec, data, pol, cells = cells),
    observed = { d <- data; d$.w <- observed_weights(spec, d, pol)
                 d[d$.w > 0, , drop = FALSE] },
    cells = { d <- cell_grid(spec, data)
              d <- d[as.character(d[[spec$cell_name]]) %in%
                     as.character(cells[[spec$cell_name]]), , drop = FALSE]
              d$.w <- policy_weights(spec, d, pol); d })
  mfx_canonical(as.data.frame(marginaleffects::avg_predictions(model,
    newdata = g, by = target, wts = g$.w, hypothesis = mfx_hypothesis(contrast))),
    levels(factor(data[[target]])))$estimate
}

## --- the code view ---------------------------------------------------------

show_code <- function(x, ...) UseMethod("show_code")
show_code.default <- function(x, ...) {
  code <- attr(x, "nestimand_code")
  if (is.null(code))
    stop("no code is attached to this object. Code is attached by nest_fit() ",
         "and estimand(); a model fitted by calling the engine directly ",
         "carries none.")
  cat(paste(code, collapse = "\n"), "\n", sep = "")
  invisible(structure(paste(code, collapse = "\n"), class = "nestimand_code"))
}
show_code.nestimand_estimand <- function(x, ...) {
  code <- attr(x, "nestimand")$code
  cat(paste(code, collapse = "\n"), "\n", sep = "")
  invisible(structure(paste(code, collapse = "\n"), class = "nestimand_code"))
}

print.nestimand_estimand <- function(x, ...) {
  meta <- attr(x, "nestimand")
  y <- x
  attr(y, "nestimand") <- NULL
  class(y) <- setdiff(class(y), "nestimand_estimand")
  print(y, ...)
  pol <- if (is.character(meta$policy)) meta$policy else "supplied"
  cat("\nPolicy: ", pol, "   route: ", meta$route, "   contrast: ", meta$contrast,
      if (identical(meta$scale, "latent")) "   scale: latent" else "", sep = "")
  if (!is.null(meta$self_check))
    cat("   reorder check: ", meta$self_check$status, sep = "")
  cat("\n")
  if (!is.null(meta$bounds)) {
    cat("Bounds over all admissible policies:\n")
    b <- meta$bounds
    for (i in seq_len(nrow(b)))
      cat(sprintf("  %-14s %8.4f   [%8.4f, %8.4f]\n", b$term[i], b$estimate[i],
                  b$policy_low[i], b$policy_high[i]))
  }
  cat("Code: show_code() prints the ", length(meta$code),
      " lines that produced this.\n", sep = "")
  invisible(x)
}
