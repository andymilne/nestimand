## nestimand: estimands, with the code that produced them ------------------
## The function computes; the code view is a by-product, not a reconstruction.
## Every estimand is assembled as text, then evaluated, so `show_code()` cannot
## drift from what was run: it is the same object. Non-core arguments in `...`
## are deparsed into that text and so reach the destination function - brms,
## marginaleffects, emmeans - unaltered, and appear in the saved code.

estimand <- function(model, target, policy = "equal", at = NULL,
                     contrast = c("pairwise", "reference", "sequential", "within"),
                     scale = c("response", "latent"),
                     data = NULL, bounds = TRUE, self_check = TRUE,
                     ..., spec = NULL, .env = parent.frame()) {
  ## The declaration travels with a fit from nest_fit(); `spec` is needed only
  ## for a model fitted by calling the engine directly.
  spec_expr <- substitute(spec)
  recovered <- is.null(spec)
  spec <- resolve_spec(model, spec)
  check_model_spec(model, spec)
  contrast <- match.arg(contrast)
  scale <- match.arg(scale)
  if (identical(scale, "latent") && identical(contrast, "within"))
    stop("contrast = \"within\" is computed through the prediction route and ",
         "has no linear-map form here; use scale = \"response\".")
  if (!inherits(spec, "nesting_spec"))
    stop("`spec` must be a nesting_spec object, as returned by nesting_spec().")
  tg <- substitute(target)
  if (is.name(tg)) target <- deparse(tg)
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

  model_name <- deparse(substitute(model))
  ## When the declaration came from the fit, refer to it by the name it had
  ## when the model was fitted, so the emitted code reads as the user wrote it.
  spec_name  <- if (recovered)
    (if (is.null(attr(model, "nestimand_spec_name"))) "spec"
     else attr(model, "nestimand_spec_name")) else deparse(spec_expr)
  data_name <- if (missing(data)) paste0(spec_name, "$data")
               else paste(deparse(substitute(data)), collapse = " ")
  if (is.null(data)) { data <- spec$data; data_name <- paste0(spec_name, "$data") }

  ## non-core arguments, passed through verbatim to the destination function
  dots <- as.list(substitute(list(...)))[-1]
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
                        deg)
  ## if the model was fitted by nest_fit(), its call travels with it, so the
  ## code view is the whole pipeline rather than its second half
  fit_code <- attr(model, "nestimand_code")
  if (!is.null(fit_code)) {
    fit_code <- sub("^m <- ", paste0(model_name, " <- "), fit_code)
    code <- c(fit_code, "", code)
  }
  env <- new.env(parent = .env)
  assign(model_name, model, envir = env)
  assign(spec_name,  spec,  envir = env)
  if (!exists(data_name, envir = env, inherits = TRUE))
    assign(data_name, data, envir = env)
  out <- eval(parse(text = paste(code, collapse = "\n")), envir = env)

  check <- if (isTRUE(self_check))
    reorder_check(model, spec, target, policy, at, contrast, dots_txt, data, scale)
  else NULL

  structure(out,
            nestimand = list(build = nestimand_build, target = target,
                             policy = policy, at = at, contrast = contrast,
                             scale = scale,
                             bounds = attr(out, "nestimand_bounds"),
                             self_check = check,
                             code = code),
            class = c("nestimand_estimand", class(out)))
}

## --- code assembly ---------------------------------------------------------

estimand_code <- function(spec, target, policy, at, contrast, dots_txt,
                          model_name, spec_name, data_name, bounds,
                          scale = "response", deg = NULL) {
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
        sprintf('est  <- add_bounds(est, %s, %s, "%s", "%s", scale = "latent"%s)',
                model_name, spec_name, target, contrast, dots_txt))
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
    "## G-computation grid: every row of the data crossed with every realized",
    "## cell, with the policy attached as row weights",
    sprintf('grid <- counterfactual_grid(%s, %s, pol%s)', spec_name, data_name,
            if (is.null(deg)) "" else ", cells = cells"),
    "## estimand in the original variable space: `by =` names an original",
    sprintf("## factor, not the fitted `%s` predictor", cn),
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
      sprintf('est  <- add_bounds(est, %s, %s, "%s", "%s"%s)',
              model_name, spec_name, target, contrast, dots_txt))
  c(body, "est")
}

## --- the bounds companion --------------------------------------------------

add_bounds <- function(est, model, spec, target, contrast = "pairwise",
                       scale = c("response", "latent"), ...) {
  scale <- match.arg(scale)
  vs <- versions_of(spec, target)
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
    g <- counterfactual_grid(spec, spec$data, p)
    e <- marginaleffects::avg_predictions(model, newdata = g, by = target,
           wts = g$.w, hypothesis = mfx_hypothesis(contrast), ...)
    mfx_canonical(as.data.frame(e), levels(factor(spec$data[[target]])))$estimate
  })
  M <- do.call(cbind, vals)
  attr(est, "nestimand_bounds") <- data.frame(
    term = mfx_canonical(as.data.frame(est))$term,
    estimate = as.data.frame(est)$estimate,
    policy_low = apply(M, 1, min), policy_high = apply(M, 1, max))
  est
}

## --- the reorder self-check ------------------------------------------------
## Order instability is impossible under the cell parameterization, so this is
## a belt-and-braces check on the translation layer rather than on the fit.
reorder_check <- function(model, spec, target, policy, at, contrast, dots_txt, data,
                          scale = "response") {
  if (inherits(model, "brmsfit"))
    return(list(status = "skipped", note = paste(
      "reorder check skipped: refitting a brms model doubles sampling time.",
      "Re-run with self_check = TRUE on a cheaper fit of the same structure,",
      "or accept the check as unnecessary under the cell parameterization.")))
  d2 <- data
  d2[[spec$cell_name]] <- NULL   # rebuilt from the permuted factors below
  for (v in spec$cell_vars) {
    lv <- levels(factor(d2[[v]]))
    d2[[v]] <- factor(d2[[v]], levels = c(lv[1], sample(lv[-1])))
  }
  out <- tryCatch({
    sp2 <- nesting_spec_quiet(spec, d2)
    m2  <- stats::update(model, data = sp2$data)
    e1  <- unname(estimand_values(model, spec, target, policy, at, contrast,
                                  spec$data, scale))
    e2  <- unname(estimand_values(m2, sp2, target, policy, at, contrast,
                                  sp2$data, scale))
    if (length(e1) != length(e2))
      list(status = "failed", note = "the two runs returned different numbers of contrasts")
    else if (!isTRUE(all.equal(sort(round(abs(e1), 8)), sort(round(abs(e2), 8)),
                               tolerance = 1e-8)))
      list(status = "failed",
           note = sprintf("max |change| = %s", format(max(abs(sort(abs(e1)) - sort(abs(e2)))))))
    else list(status = "passed", note = "estimand unchanged under level permutation")
  }, error = function(e) list(status = "error", note = conditionMessage(e)))
  if (identical(out$status, "failed"))
    warning("reorder self-check FAILED (", out$note, "): the estimate moved when ",
            "nested-factor levels were permuted, so this estimand depends on ",
            "unrealized-cell predictions and is not identified by the design. ",
            "Do not report it.", call. = FALSE)
  out
}

nesting_spec_quiet <- function(spec, data) {
  suppressMessages(nesting_spec(data, spec$formula_in,
    nests = spec_nests(spec), fit = spec$fit, family = spec$family,
    random = spec$random_original, cell_name = spec$cell_name))
}
spec_nests <- function(spec)
  unlist(lapply(spec$cat_families, function(f)
    if (length(f) > 1) sprintf("%s %%in%% %s", f[-1], f[-length(f)]) else character(0)))

estimand_values <- function(model, spec, target, policy, at, contrast, data,
                            scale = "response") {
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
  g <- counterfactual_grid(spec, data, pol, cells = cells)
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
  cat("\nPolicy: ", pol, "   contrast: ", meta$contrast,
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
