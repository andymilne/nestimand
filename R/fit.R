## nestimand: the fitting side ---------------------------------------------
## Same discipline as the estimand side: the call is assembled as text and then
## evaluated, so the code attached to the returned fit is the call that made it.
## Non-core arguments in `...` are deparsed into that text and so reach the
## engine - brms, lme4, ordinal - unaltered, and appear in the saved code.

## --- random-effects translation -------------------------------------------
## This is the decisive argument for the cell parameterization, not a
## consequence of it. Under the chain form, a random slope over the nesting
## structure is structurally rank-deficient: the columns for impossible
## combinations are identically zero for every group, so the covariance matrix
## has more dimensions than the data can ever identify. lme4 does not refuse
## it - on the demonstration data it fits a 16-dimensional covariance of rank
## 10, with three all-zero columns and a degenerate Hessian, and reports only a
## convergence warning. The chain form is therefore confined in practice to
## nested grouping factors, `(1|p) + (1|p:chord) + (1|p:cell)`, which impose
## compound symmetry: one variance per level, equal correlation within a
## stratum. The cell form places the whole structure on a single factor, where
## every random-effects structure R offers is available and identified -
## unstructured, diagonal, or the grouping chain as a submodel.

## Which variables actually span the structural boundary: the nested ones, not
## every variable in a nesting family. A slope on the root of a family - chord
## type here - is encountered at every level by every unit, so it is identified
## as written and must be left exactly as declared. Translating it would
## silently enlarge the model the user asked for.
boundary_vars <- function(spec)
  unlist(lapply(spec$cat_families, function(f) f[-1]), use.names = FALSE)

## The parts of the declared random terms, read from the parsed expression
## rather than by matching parentheses in the text. A bar whose left side is
## itself bracketed - `(chord_type * (inversion + top) | participant)`, which is
## how anyone writes a slope over a set of variables - defeated the regex: it
## matched the inner group, and the term arrived with no grouping factor at all.
## A covariance-structure wrapper such as `diag(...)` is peeled off and kept, to
## be restored after translation.
bar_terms_of <- function(bars) {
  e <- str2lang(paste("~", bars))[[2]]
  flat <- function(z) if (is.call(z) && identical(as.character(z[[1]]), "+"))
    c(flat(z[[2]]), flat(z[[3]])) else list(z)
  lapply(flat(e), function(z) {
    wrapper <- ""
    while (is.call(z) && !as.character(z[[1]]) %in% c("|", "||")) {
      if (!identical(as.character(z[[1]]), "(")) wrapper <- as.character(z[[1]])
      z <- z[[2]]
    }
    if (!is.call(z) || !as.character(z[[1]]) %in% c("|", "||"))
      stop("`", paste(deparse(z), collapse = " "), "` is not a random-effects ",
           "term: a bar and a grouping factor were expected.")
    list(wrapper = wrapper, op = as.character(z[[1]]),
         lhs = paste(deparse(z[[2]]), collapse = " "),
         grp = paste(deparse(z[[3]]), collapse = " "))
  })
}

## The random side follows the declaration, exactly as the fixed side does.
## A bar that names no nesting variable is passed through as written. A bar that
## does reach across the structure is written on the same columns the mean
## structure is fitted on: over the original factors it would carry columns for
## conditions that do not exist, which no amount of data can inform.
random_terms <- function(spec) {
  bars <- spec$random_original
  if (is.null(bars)) return(NULL)
  ## A bar may be wrapped in a covariance-structure call - `diag(...)` in lme4
  ## 2.0-0 and later, for a diagonal covariance. The wrapper states what the
  ## covariance looks like and must survive translation: dropping it would turn
  ## a diagonal request into an unstructured one without saying so.
  bl <- bar_terms_of(bars)
  out <- unlist(lapply(seq_along(bl), function(k) {
    wrapper <- bl[[k]]$wrapper
    rewrap <- function(x) if (nzchar(wrapper))
      vapply(x, function(z) paste0(wrapper, z), "") else x
    grp <- bl[[k]]$grp
    op  <- bl[[k]]$op
    lhs_txt <- bl[[k]]$lhs
    lhs <- attr(stats::terms(stats::as.formula(paste("~", lhs_txt))), "term.labels")
    ## a term needs translation only if it reaches across the boundary
    structural <- vapply(strsplit(lhs, ":"), function(vs)
      any(vs %in% boundary_vars(spec)), TRUE)
    ## covariates travel with the term; anything built from the nesting
    ## variables is subsumed by the translated structure and must not be
    ## repeated alongside it
    covs <- lhs[!vapply(strsplit(lhs, ":"), function(vs)
      any(vs %in% spec$cell_vars), TRUE)]
    ## A covariate declared crossed with the structure keeps that crossing: it
    ## gets a slope per column of the translated design and one for the
    ## reference condition, exactly as on the fixed side.
    crossed <- vapply(covs, function(k) {
      kv <- strsplit(k, ":")[[1]]
      any(vapply(strsplit(lhs, ":"), function(vs)
        all(kv %in% vs) && any(vs %in% spec$cell_vars), TRUE))
    }, TRUE)
    if (!any(structural))
      return(rewrap(sprintf("(%s %s %s)", lhs_txt, op, grp)))
    ## The declared structure varying by group: the same columns the fixed side
    ## is fitted on, so the random effects are deviations in exactly the effects
    ## the bar names, and nothing larger.
    want <- declared_terms(spec, lhs[vapply(strsplit(lhs, ":"), function(vs)
      any(vs %in% spec$cell_vars), TRUE)])
    extra <- want[!term_key(want) %in% term_key(declared_terms(spec))]
    if (length(extra))
      stop("the random term `(", lhs_txt, " ", op, " ", grp, ")` asks for `",
           paste(extra, collapse = "`, `"), "`, which the model formula does ",
           "not contain. A random effect for a term the mean structure leaves ",
           "out says it averages to zero across groups but varies between them ",
           "- coherent, but not something the formula states. Add the term to ",
           "the model formula, or drop it from the bar.")
    cols <- reduced_columns(spec, want)
    cov_red <- if (!length(covs)) "" else paste(" +", paste(unlist(
      lapply(seq_along(covs), function(i) if (!crossed[[i]]) covs[[i]] else
        c(covs[[i]], paste0(cols, ":", covs[[i]])))), collapse = " + "))
    rewrap(sprintf("(1 + %s%s %s %s)", paste(cols, collapse = " + "),
                   cov_red, op, grp))
  }))
  paste(unique(out), collapse = " + ")
}

## --- the fit ---------------------------------------------------------------
nest_fit <- function(fit, formula, data, family = NULL, nests = NULL,
                     random = NULL, dry_run = FALSE, ..., .env = parent.frame()) {
  ## The declaration and the fit are one step. The constructor is called here
  ## with the user's own expressions rather than their values, so that
  ## nesting_spec() sees `nests` unevaluated - which is what lets
  ## `inversion %in% chord_type` be written unquoted.
  cl0 <- match.call()
  if (missing(fit) || missing(formula) || missing(data))
    stop("nest_fit() takes an engine, a formula and the data, in that order: ",
         "nest_fit(\"lm\", rating ~ chord_type * inversion, dat). `nests` may ",
         "be left out, in which case the structure is read from the data - a ",
         "variable's structurally undefined rows being marked NA.")
  args <- list(quote(nesting_spec), data = cl0$data, formula = cl0$formula,
               fit = cl0$fit)
  for (nm in c("nests", "family", "random"))
    if (!is.null(cl0[[nm]])) args[[nm]] <- cl0[[nm]]
  spec_call <- as.call(args)
  spec <- eval(spec_call, .env)
  spec_name <- "spec"
  data_name <- "spec$data"
  data <- spec$data
  fit <- spec$fit

  ## The parameterization is decided by the formula, never asked for. A formula
  ## that crosses the structure fully names one mean per realized condition, and
  ## the cell factor is the compact way to write that; anything less is written
  ## as the design the formula gives over the realized conditions.
  fm <- fitting_mode(spec)
  mode <- as.character(fm)
  mode_note <- sprintf("## parameterization: %s (%s)", mode, attr(fm, "reason"))
  f <- paste(deparse(cell_formula(spec, mode)), collapse = " ")
  ## The emitted code has to be runnable, so it names the fitted columns rather
  ## than the structure they encode. A line saying what `cell` is costs nothing
  ## and is the difference between readable and cryptic.
  cell_note <- if (identical(mode, "cells"))
    c(sprintf("## `%s` is a factor whose levels are the %d realized conditions, so its",
              spec$cell_name, nrow(spec$cells)),
      "##   coefficients are their means. The formula crosses the structure fully,",
      "##   so that is the model it names, written compactly.")
  re <- if (fit %in% c("lmer", "glmer", "clmm", "brm")) random_terms(spec)
  ## The reduced form states the declared structure as one column per identified
  ## effect, computed once per realized cell and looked up by cell. It is the
  ## same kind of object the cell factor is - a design the data can inform in
  ## full - so nothing has to be held at zero, and the columns are ordinary
  ## numeric variables that every engine treats alike.
  aug_code <- NULL
  if (identical(mode, "reduced") ||
      (!is.null(re) && grepl("dm_", re, fixed = TRUE))) {  # == uses_reduced()
    raw_data_name <- data_name
    raw_data <- data
    data <- with_reduced(spec, data)
    data_name <- ".nestimand_data"
    aug_code <- c(
      "## one column per identified effect of the declared structure, looked up",
      "## by realized cell. reduced_design(spec) shows the columns and the",
      "## effect each one stands for.",
      sprintf("%s <- with_reduced(%s, %s)", data_name, spec_name, raw_data_name))
  }
  if (!is.null(re)) f <- paste(f, "+", re)
  fam <- if (!is.null(spec$family)) sprintf(", family = %s", spec$family) else ""
  dots <- as.list(substitute(list(...)))[-1]
  ## brms's own `prior` is the route to everything the translation does not
  ## touch - the random-effects sd and correlations, the thresholds, sigma. It
  ## is passed through untouched. A `class = "b"` entry is a different matter:
  ## those coefficients are the realized cells, not the original effects, so a
  ## prior written for the effects would mean something else there, and is
  ## translated rather than handed over as written.
  priors <- NULL
  prior_name <- ".prior"
  user_prior <- NULL
  user_prior_txt <- NULL
  if ("prior" %in% names(dots)) {
    user_prior <- tryCatch(eval(dots$prior, .env), error = function(e) NULL)
    user_prior_txt <- paste(deparse(dots$prior), collapse = " ")
  }
  if (inherits(user_prior, "brmsprior") && any(user_prior$class == "b") &&
      identical(mode, "cells")) {
    split <- prior_from_brms(spec, user_prior)
    priors <- split$translated
    message("the `class = \"b\"` prior was stated for the original variables ",
            "and has been translated to the coefficients the model is fitted ",
            "on: mu ~ N(A m, A D A'), reaching Stan through stanvars. ",
            "prior_for_estimand() shows what it implies for a contrast.")
    user_prior_txt <- if (!is.null(split$other) && nrow(split$other))
      sprintf("brms::prior_string(%s, class = %s%s)",
              paste0("c(", paste(sprintf('"%s"', split$other$prior), collapse = ", "), ")"),
              paste0("c(", paste(sprintf('"%s"', split$other$class), collapse = ", "), ")"),
              "") else NULL
    dots$prior <- NULL
    user_prior <- split$other
  }
  if (!is.null(priors) && !is.null(user_prior) && "prior" %in% names(dots))
    dots$prior <- NULL
  dots_txt <- if (length(dots))
    paste0(", ", paste(mapply(function(nm, v) {
      v <- paste(deparse(v), collapse = " ")
      if (nzchar(nm)) paste0(nm, " = ", v) else v
    }, names(dots), dots), collapse = ", ")) else ""
  ## brms writes one prior draw per class when asked for prior samples, and
  ## declares it `real` - which a multivariate prior on the coefficients cannot
  ## satisfy, so the Stan program fails to compile. The alternative is a
  ## separate prior-only run.
  if (!is.null(priors) && "sample_prior" %in% names(dots)) {
    sv <- tryCatch(eval(dots$sample_prior, .env), error = function(e) NULL)
    if (isTRUE(sv) || identical(sv, "yes"))
      stop("`sample_prior = \"yes\"` cannot be combined with a prior on the ",
           "mean structure. The translation puts a multivariate normal on the ",
           "coefficients, and brms draws one prior sample per class into a ",
           "scalar, so Stan reports an ill-typed assignment - ",
           "`real prior_b = multi_normal_rng(...)`. Run the prior separately ",
           "with sample_prior = \"only\", which works and is what the prior ",
           "check in tests/test_brms.R does.")
  }
  prior_txt <- ""
  prior_note <- NULL
  if (!is.null(priors)) {
    check_prior_dimension(priors, spec, mode)
    ## a prior for the classes the translation does not cover is combined with
    ## it rather than emitted a second time, which brms would refuse
    other <- user_prior_txt
    prior_txt <- if (is.null(other))
      sprintf(", prior = prior_statement(%s), stanvars = prior_stanvars(%s)",
              prior_name, prior_name)
    else sprintf(", prior = c(prior_statement(%s), %s), stanvars = prior_stanvars(%s)",
                 prior_name, other, prior_name)
    prior_note <- c(
      sprintf("## prior stated on %s and translated to the cell means:", priors$stated_on),
      "## mu ~ N(A m, A D A'), exact and basis-free, since Stan accepts a",
      "## correlated prior directly. prior_audit() shows what it implies in",
      "## both spaces, and prior_for_estimand() what it implies for a contrast.")
  }
  fn <- switch(fit, lm = "lm", glm = "glm", lmer = "lme4::lmer", glmer = "lme4::glmer",
               clm = "ordinal::clm", clmm = "ordinal::clmm", brm = "brms::brm")
  lib <- switch(fit, lmer = "library(lme4)", glmer = "library(lme4)",
                clm = "library(ordinal)", clmm = "library(ordinal)",
                brm = "library(brms)")
  re_note <- if (!is.null(re) && !identical(re, spec$random_original))
    c(sprintf("## random structure %s", spec$random_original),
      sprintf("## translated to   %s", re),
      "## - a random slope written over the original factors carries columns the",
      "##   data cannot inform: those for conditions that do not exist, and a",
      "##   further set reconstructable from the others. The columns above are",
      "##   the ones the declared structure identifies, so every parameter of",
      "##   the covariance is estimable and nothing is held at zero.")
  ## the declaration is part of the code that reproduces the fit
  spec_code <- c("## the declaration, built from the data and the formula",
                 paste("spec <-", paste(deparse(spec_call), collapse = " ")))
  code <- c(sprintf("## nestimand %s -- fit", nestimand_build), lib, spec_code,
            mode_note,
            cell_note, aug_code, re_note, prior_note,
            sprintf("m <- %s(%s%s, data = %s%s%s)", fn, f, fam, data_name,
                    prior_txt, dots_txt))
  if (isTRUE(dry_run))
    return(structure(paste(code, collapse = "\n"), class = "nestimand_code",
                     nestimand_code = code))
  env <- new.env(parent = .env)
  assign(spec_name, spec, envir = env)
  if (!is.null(priors)) assign(prior_name, priors, envir = env)
  m <- eval(parse(text = paste(c(code, "m"), collapse = "\n")), envir = env)
  attr(m, "nestimand_code") <- code
  attr(m, "nestimand_mode") <- mode
  ## the fit carries its declaration, so the summary and estimand functions can
  ## be called on the model alone
  attr(m, "nestimand_spec") <- spec
  attr(m, "nestimand_spec_name") <- spec_name
  ## An extra class, so that update() can carry the declaration across a refit.
  ## It goes *after* the engine's own classes, not before: packages that
  ## dispatch or validate on class(model)[1] - marginaleffects among them -
  ## would otherwise see a class they do not know and warn about arguments that
  ## are perfectly valid for the underlying fit.
  if (!isS4(m)) class(m) <- unique(c(class(m), "nestimand_fit"))
  m
}

## update() rebuilds the model from its call. That call names the declaration
## the package built - `spec` - which exists only in the environment nest_fit()
## evaluated it in, so the default method, which evaluates in the caller's
## frame, cannot find it. The declaration is put back within reach here and
## reattached to the result, so a fit that has been simplified - a covariate
## dropped, a random term reduced - remains usable without anything having to be
## supplied again. Whether the updated model still corresponds to the
## declaration is checked where it matters, by check_model_spec().
update.nestimand_fit <- function(object, formula., ..., evaluate = TRUE) {
  cl <- stats::getCall(object)
  if (is.null(cl))
    stop("this fit carries no call, so it cannot be updated.")
  if (!missing(formula.))
    cl$formula <- stats::update.formula(stats::formula(object), formula.)
  extras <- match.call(expand.dots = FALSE)$...
  if (length(extras))
    for (nm in names(extras))
      if (is.null(extras[[nm]])) cl[[nm]] <- NULL else cl[[nm]] <- extras[[nm]]
  if (!evaluate) return(cl)
  e <- new.env(parent = parent.frame())
  nm <- attr(object, "nestimand_spec_name"); if (is.null(nm)) nm <- "spec"
  assign(nm, attr(object, "nestimand_spec"), envir = e)
  out <- eval(cl, e)
  for (a in c("nestimand_spec", "nestimand_spec_name", "nestimand_mode"))
    attr(out, a) <- attr(object, a)
  attr(out, "nestimand_code") <-
    c(sprintf("## nestimand %s -- fit, updated from the original call",
              nestimand_build),
      paste("m <-", paste(deparse(cl), collapse = " ")))
  if (!isS4(out)) class(out) <- unique(c(class(out), "nestimand_fit"))
  out
}

## Does the fit actually correspond to the declaration? A model that does not
## contain the declared structure cannot respond to it, so averaging over
## versions returns whatever weighting the fit already implies rather than the
## policy asked for - and on balanced data the two can coincide, which makes the
## error invisible exactly where it is easiest to make.
check_model_spec <- function(model, spec) {
  ## a brmsformula yields nothing to all.vars(); one more formula() call
  ## reduces it to the ordinary formula underneath
  fm <- tryCatch(stats::formula(model), error = function(e) NULL)
  if (inherits(fm, "bform") || inherits(fm, "brmsformula"))
    fm <- tryCatch(stats::formula(fm), error = function(e) fm)
  v <- tryCatch(all.vars(fm), error = function(e) NULL)
  if (!length(v)) v <- NULL
  if (is.null(v)) return(invisible(TRUE))          # engine without a formula
  if (spec$cell_name %in% v) return(invisible(TRUE))
  ## The reduced form names the design's own columns rather than the original
  ## factors, and each is a fixed function of the cell, so the fit does respond
  ## to the declared structure - by construction, and in exactly the terms the
  ## declaration asked for.
  if (identical(attr(model, "nestimand_mode"), "reduced") &&
      all(setdiff(colnames(reduced_design(spec)), "(Intercept)") %in% v))
    return(invisible(TRUE))
  miss <- setdiff(spec$cell_vars, v)
  if (length(miss))
    stop("the model does not contain `", paste(miss, collapse = "`, `"),
         "`, which the declaration nests. Its predictions cannot respond to ",
         "that variable, so averaging over versions would return the weighting ",
         "the fit already implies - not the policy requested, whatever the ",
         "output is labelled. On balanced data the two can agree, which is why ",
         "this is checked rather than left to be noticed. Fit with nest_fit(), ",
         "or include the declared structure in the model.")
  invisible(TRUE)
}

## Recover the declaration from a fit, when it was not passed explicitly.
resolve_spec <- function(model, spec = NULL) {
  if (inherits(spec, "nesting_spec")) return(spec)
  s <- attr(model, "nestimand_spec")
  if (is.null(s))
    stop("this model was not fitted by nest_fit(), so it carries no ",
         "declaration of the nesting structure and nothing here can be read ",
         "against one. Fit it with nest_fit().")
  s
}

print.nestimand_code <- function(x, ...) {
  cat("## code only; nothing below has been run\n")
  cat(unclass(x), "\n", sep = "")
  invisible(x)
}
