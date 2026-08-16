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

## Reports the rank deficiency of a declared structural random slope, so the
## constraint can be shown rather than asserted.
random_slope_rank <- function(spec, data = spec$data) {
  f <- stats::as.formula(paste("~ 0 +", paste(chain_terms(spec), collapse = " + ")))
  X <- stats::model.matrix(f, data)
  list(columns = ncol(X), rank = qr(X)$rank,
       zero_columns = sum(colSums(X != 0) == 0),
       identified = qr(X)$rank == ncol(X))
}
## Which variables actually span the structural boundary: the nested ones, not
## every variable in a nesting family. A slope on the root of a family - chord
## type here - is encountered at every level by every unit, so it is identified
## as written and must be left exactly as declared. Translating it would
## silently enlarge the model the user asked for.
boundary_vars <- function(spec)
  unlist(lapply(spec$cat_families, function(f) f[-1]), use.names = FALSE)

random_terms <- function(spec, structure = c("cells", "chain", "chain_slope",
                                             "as_declared")) {
  structure <- match.arg(structure)
  bars <- spec$random_original
  if (is.null(bars)) return(NULL)
  if (structure == "as_declared") {
    if (any(vapply(boundary_vars(spec), function(v)
          grepl(v, bars, fixed = TRUE), TRUE))) {
      rk <- random_slope_rank(spec)
      if (!rk$identified)
        warning("the declared random structure places slopes on the nesting ",
                "variables, whose design is rank-deficient by construction: ",
                rk$columns, " columns of rank ", rk$rank, ", with ",
                rk$zero_columns, " identically zero. The engine will fit it ",
                "without refusing, but ", rk$columns - rk$rank, " dimensions of ",
                "the covariance are not identified by any amount of data. ",
                "random_structure = \"cells\" expresses the same intent on the ",
                "realized cells, where every dimension is identified.",
                call. = FALSE)
    }
    return(bars)
  }
  ## A bar may be wrapped in a covariance-structure call - `diag(...)` in lme4
  ## 2.0-0 and later, for a diagonal covariance. The wrapper states what the
  ## covariance looks like and must survive translation: dropping it would turn
  ## a diagonal request into an unstructured one without saying so.
  pieces <- regmatches(bars, gregexpr("[A-Za-z.][A-Za-z0-9._]*\\([^()]*\\)|\\([^()]*\\)",
                                      bars))[[1]]
  wrap <- sub("\\(.*$", "", pieces)
  bl <- sub("^[A-Za-z.][A-Za-z0-9._]*\\(|^\\(", "", sub("\\)$", "", pieces))
  cn <- spec$cell_name
  out <- unlist(lapply(seq_along(bl), function(k) {
    b <- bl[[k]]; wrapper <- wrap[[k]]
    rewrap <- function(x) if (nzchar(wrapper))
      vapply(x, function(z) paste0(wrapper, z), "") else x
    parts <- strsplit(b, "|", fixed = TRUE)[[1]]
    grp <- trimws(parts[2])
    lhs <- attr(stats::terms(stats::as.formula(paste("~", parts[1]))), "term.labels")
    lhs_txt <- trimws(parts[1])
    ## a term needs translation only if it reaches across the boundary
    structural <- vapply(strsplit(lhs, ":"), function(vs)
      any(vs %in% boundary_vars(spec)), TRUE)
    ## covariates travel with the term; anything built from the nesting
    ## variables is subsumed by the translated structure and must not be
    ## repeated alongside it
    covs <- lhs[!vapply(strsplit(lhs, ":"), function(vs)
      any(vs %in% spec$cell_vars), TRUE)]
    if (!any(structural))
      return(rewrap(sprintf("(%s | %s)", lhs_txt, grp)))
    ## Does the declared term describe the whole categorical structure, or only
    ## part of it? Only the first has a direct expression over realized cells.
    struct_vars <- unique(unlist(strsplit(lhs[structural], ":")))
    struct_vars <- struct_vars[struct_vars %in% spec$cell_vars]
    fam <- NULL
    for (ff in spec$cat_families) if (any(struct_vars %in% ff)) fam <- ff
    partial <- !all(fam %in% struct_vars)
    if (partial && structure %in% c("cells", "chain_slope"))
      stop("the random term `(", lhs_txt, " | ", grp, ")` varies with `",
           paste(intersect(struct_vars, boundary_vars(spec)), collapse = "`, `"),
           "` but not with `", paste(setdiff(fam, struct_vars), collapse = "`, `"),
           "`. Over the realized conditions that is a covariance of reduced ",
           "rank, which no formula states directly. Three ways forward: declare ",
           "the whole structure (`", paste(fam, collapse = " * "),
           "`) for an unstructured covariance over the realized conditions; use ",
           "random_structure = \"chain\" for the grouping-factor submodel, which ",
           "is the parsimonious counterpart and always identified; or ",
           "random_structure = \"as_declared\" to fit it as written, noting that ",
           "the sentinel level then enters as though it were a level of the ",
           "nested variable.")
    cov_txt <- if (length(covs)) paste(" +", paste(covs, collapse = " + ")) else ""
    if (structure == "chain_slope")
      ## the chain counterpart of an unstructured cell covariance: the same
      ## columns as the fixed part, so one set of declarations covers both
      rewrap(sprintf("(0 + %s%s | %s)", paste(chain_terms(spec), collapse = " + "),
                     cov_txt, grp))
    else if (structure == "cells")
      rewrap(sprintf("(0 + %s%s | %s)", cn, cov_txt, grp))
    else {
      ## the grouping-chain submodel: parsimonious, but it constrains conditions
      ## sharing a stratum to equal correlation, where the cell form does not.
      ## It is always identified, and is the remedy offered for a partial
      ## structure, so it accepts one: the chain is truncated at the deepest
      ## declared level.
      depth <- if (partial) max(match(intersect(fam, struct_vars), fam)) else length(fam)
      c(rewrap(sprintf("(1%s | %s)", cov_txt, grp)),
        vapply(seq_len(depth), function(k)
          sprintf("(1 | %s)", paste(c(grp, fam[seq_len(k)]), collapse = ":")), ""))
    }
  }))
  paste(unique(out), collapse = " + ")
}

## --- the fit ---------------------------------------------------------------
nest_fit <- function(spec, mode = NULL, random_structure = c("cells", "chain", "as_declared"),
                     data = NULL, engine = "marginaleffects",
                     dry_run = FALSE, ..., priors = NULL, .env = parent.frame()) {
  random_structure <- match.arg(random_structure)
  spec_name  <- deparse(substitute(spec))
  prior_name <- deparse(substitute(priors))
  data_name <- if (missing(data) || is.null(substitute(data)))
    paste0(spec_name, "$data") else paste(deparse(substitute(data)), collapse = " ")
  if (is.null(data)) { data <- spec$data; data_name <- paste0(spec_name, "$data") }
  fit <- spec$fit
  if (is.null(mode)) {
    fm <- fitting_mode(spec, engine = engine, priors = priors)
    mode <- as.character(fm)
    mode_note <- sprintf("## parameterization: %s (%s)", mode, attr(fm, "reason"))
  } else {
    mode <- match.arg(mode, c("cells", "effects"))
    mode_note <- sprintf("## parameterization: %s (requested)", mode)
    if (identical(mode, "effects") && identical(fit, "brms") && is.null(priors))
      message("chain parameterization on brms without declarations: the design ",
              "carries columns the data cannot inform, and the posterior will ",
              "be improper along them. chain_priors(spec) derives the ",
              "constant(0) block; pass it as priors =.")
  }
  f <- paste(deparse(cell_formula(spec, mode)), collapse = " ")
  if (identical(mode, "effects") && identical(random_structure, "cells"))
    random_structure <- "chain_slope"   # match the random side to the fixed side
  re <- if (fit %in% c("lmer", "glmer", "clmm")) random_terms(spec, random_structure)
        else if (identical(fit, "brms")) random_terms(spec, random_structure)
  if (!is.null(re)) f <- paste(f, "+", re)
  fam <- if (!is.null(spec$family)) sprintf(", family = %s", spec$family) else ""
  dots <- as.list(substitute(list(...)))[-1]
  dots_txt <- if (length(dots))
    paste0(", ", paste(mapply(function(nm, v) {
      v <- paste(deparse(v), collapse = " ")
      if (nzchar(nm)) paste0(nm, " = ", v) else v
    }, names(dots), dots), collapse = ", ")) else ""
  prior_txt <- ""
  prior_note <- NULL
  ## brms's own `prior` argument must reach the engine untouched. With `priors`
  ## declared before `...` R would partial-match it here and refuse the call, so
  ## `priors` sits after the dots and has to be named in full.
  if (inherits(priors, "brmsprior"))
    stop("`priors` takes a nestimand prior object - nest_prior() for a translated ",
         "prior, chain_priors() for the chain-mode declarations. A brms prior ",
         "goes to the engine under its own name: pass `prior = ` and it will be ",
         "handed through unaltered.")
  if (inherits(priors, "nestimand_chain_priors")) {
    if (!identical(mode, "effects"))
      stop("chain-mode declarations belong to the chain parameterization, but ",
           "the fitting mode is `", mode, "`. Pass mode = \"effects\", or use ",
           "nest_prior() for a translated prior on the cell parameterization.")
    prior_txt <- sprintf(", prior = chain_prior_object(%s)", prior_name)
    prior_note <- c(priors$code[1:4],
      "## Held at zero below; every other coefficient is estimated as usual.")
    priors_obj <- priors
    priors <- NULL
  } else priors_obj <- NULL
  if (!is.null(priors)) {
    if (!inherits(priors, "nestimand_prior"))
      stop("`priors` must be a nestimand_prior object, as returned by nest_prior().")
    if (!identical(fit, "brms"))
      stop("translated priors apply to the brms engine; the declared engine is `",
           fit, "`, which has no prior to state.")
    check_prior_dimension(priors, spec, mode)
    prior_txt <- sprintf(", prior = prior_statement(%s), stanvars = prior_stanvars(%s)",
                         prior_name, prior_name)
    prior_note <- c(
      sprintf("## prior stated on %s and translated to the cell means:", priors$stated_on),
      "## mu ~ N(A m, A D A'), exact and basis-free, since Stan accepts a",
      "## correlated prior directly. prior_audit() shows what it implies in",
      "## both spaces, and prior_for_estimand() what it implies for a contrast.")
  }
  fn <- switch(fit, lm = "lm", glm = "glm", lmer = "lme4::lmer", glmer = "lme4::glmer",
               clm = "ordinal::clm", clmm = "ordinal::clmm", brms = "brms::brm")
  lib <- switch(fit, lmer = "library(lme4)", glmer = "library(lme4)",
                clm = "library(ordinal)", clmm = "library(ordinal)",
                brms = "library(brms)")
  re_note <- if (!is.null(re) && !identical(re, spec$random_original))
    c(sprintf("## random structure %s", spec$random_original),
      sprintf("## translated to   %s", re),
      "## - a random slope written over the original factors carries columns the",
      "##   data cannot inform: those for conditions that do not exist, and a",
      "##   further set reconstructable from the others. The cell form has one",
      "##   column per realized condition, so every parameter is estimable and",
      "##   nothing has to be held at zero. See chain_priors() for the",
      "##   alternative, which keeps the original factors and declares the",
      "##   uninformative columns instead.")
  relevel_code <- if (identical(mode, "effects"))
    sentinel_relevel_code(spec, data_name)
  code <- c(sprintf("## nestimand %s -- fit", nestimand_build), lib, mode_note,
            relevel_code, re_note, prior_note,
            sprintf("m <- %s(%s%s, data = %s%s%s)", fn, f, fam, data_name,
                    prior_txt, dots_txt))
  if (isTRUE(dry_run))
    return(structure(paste(code, collapse = "\n"), class = "nestimand_code",
                     nestimand_code = code))
  env <- new.env(parent = .env)
  assign(spec_name, spec, envir = env)
  if (!is.null(priors)) assign(prior_name, priors, envir = env)
  if (!is.null(priors_obj)) assign(prior_name, priors_obj, envir = env)
  if (!exists(data_name, envir = env, inherits = TRUE))
    assign(data_name, data, envir = env)
  m <- eval(parse(text = paste(c(code, "m"), collapse = "\n")), envir = env)
  attr(m, "nestimand_code") <- code
  attr(m, "nestimand_mode") <- mode
  ## the fit carries its declaration, so the summary and estimand functions can
  ## be called on the model alone
  attr(m, "nestimand_spec") <- spec
  attr(m, "nestimand_spec_name") <- spec_name
  ## An extra class, so that update() can carry the declaration across a refit.
  ## S4 fits - brmsfit among them - are left alone: their class is part of the
  ## object's formal definition, and the attributes travel on their own.
  if (!isS4(m)) class(m) <- unique(c("nestimand_fit", class(m)))
  m
}

## update() rebuilds the model from its call, which drops everything attached to
## it. The declaration is reattached here, so a fit that has been simplified -
## a covariate dropped, a random term reduced - remains usable without the spec
## having to be supplied again. Whether the updated model still corresponds to
## the declaration is checked where it matters, by check_model_spec().
update.nestimand_fit <- function(object, ...) {
  out <- NextMethod()
  for (a in c("nestimand_spec", "nestimand_spec_name", "nestimand_mode"))
    attr(out, a) <- attr(object, a)
  cl <- tryCatch(paste(deparse(stats::getCall(out)), collapse = " "),
                 error = function(e) NULL)
  attr(out, "nestimand_code") <- if (is.null(cl)) attr(object, "nestimand_code")
    else c(sprintf("## nestimand %s -- fit, updated from the original call",
                   nestimand_build),
           paste("m <-", cl))
  if (!isS4(out)) class(out) <- unique(c("nestimand_fit", class(out)))
  out
}

## Does the fit actually correspond to the declaration? A model that does not
## contain the declared structure cannot respond to it, so averaging over
## versions returns whatever weighting the fit already implies rather than the
## policy asked for - and on balanced data the two can coincide, which makes the
## error invisible exactly where it is easiest to make.
check_model_spec <- function(model, spec) {
  v <- tryCatch(all.vars(stats::formula(model)), error = function(e) NULL)
  if (is.null(v)) return(invisible(TRUE))          # engine without a formula
  if (spec$cell_name %in% v) return(invisible(TRUE))
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
    stop("no `spec` supplied, and this model does not carry one. A model fitted ",
         "by nest_fit() carries its declaration; one fitted by calling the ",
         "engine directly does not, so pass the nesting_spec explicitly.")
  s
}

print.nestimand_code <- function(x, ...) {
  cat("## code only; nothing below has been run\n")
  cat(unclass(x), "\n", sep = "")
  invisible(x)
}
