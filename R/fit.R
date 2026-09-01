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

random_terms <- function(spec, structure = c("cells", "reduced", "chain",
                                             "chain_slope", "as_declared")) {
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
  bl <- bar_terms_of(bars)
  cn <- spec$cell_name
  out <- unlist(lapply(seq_along(bl), function(k) {
    wrapper <- bl[[k]]$wrapper
    rewrap <- function(x) if (nzchar(wrapper))
      vapply(x, function(z) paste0(wrapper, z), "") else x
    grp <- bl[[k]]$grp
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
    ## A covariate declared crossed with the structure keeps that crossing. The
    ## structural terms it appeared in are dropped as subsumed by the
    ## translation, but `cell` alone carries no slope that varies by condition,
    ## so the crossing has to be restated against the translated factor - the
    ## same distinction the fixed side draws with `cov_by_cell`.
    crossed <- vapply(covs, function(k) {
      kv <- strsplit(k, ":")[[1]]
      any(vapply(strsplit(lhs, ":"), function(vs)
        all(kv %in% vs) && any(vs %in% spec$cell_vars), TRUE))
    }, TRUE)
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
    ## `with` is the translated structure the crossing is restated against -
    ## the cell factor, or the chain terms, which distribute over the covariate
    ## one term at a time. NULL crosses nothing, for the compound-symmetric
    ## submodel, where a crossed slope would not be the parsimonious
    ## counterpart it is offered as.
    cov_term <- function(with) {
      if (!length(covs)) return("")
      out <- unlist(lapply(seq_along(covs), function(i)
        if (is.null(with) || !crossed[[i]]) covs[[i]]
        else paste0(with, ":", covs[[i]])))
      paste(" +", paste(out, collapse = " + "))
    }
    cov_txt <- cov_term(cn)
    if (structure == "reduced") {
      ## The declared structure varying by group: the same columns the fixed
      ## side is fitted on, so the random effects are deviations in exactly the
      ## effects the mean structure states, and the two have the same dimension.
      ## The cell form would enlarge it to the saturated covariance - which is
      ## more than the declaration asks for, and the same fault the reduced
      ## fixed design exists to avoid.
      ## every declared term that names a design variable, not only those that
      ## reach across a structural boundary: `chord_type` and `top` are part of
      ## the structure the group varies in just as `chord_type:inversion` is,
      ## and the cell factor subsumed them only because it spans everything
      want <- declared_terms(spec, lhs[vapply(strsplit(lhs, ":"), function(vs)
        any(vs %in% spec$cell_vars), TRUE)])
      extra <- setdiff(want, declared_terms(spec))
      if (length(extra))
        stop("the random term `(", lhs_txt, " | ", grp, ")` asks for `",
             paste(extra, collapse = "`, `"), "`, which the mean structure does ",
             "not contain. A random effect for something held out of the mean ",
             "says it averages to zero across groups but varies between them - ",
             "coherent, but not what the reduced random structure expresses, ",
             "which is the declared mean structure varying by group. Add the ",
             "term to the model formula, or use random_structure = \"cells\" for ",
             "the unstructured covariance over realized conditions.")
      cols <- reduced_columns(spec, want)
      ## a covariate crossed with the structure gets a slope per column and one
      ## for the reference condition, exactly as on the fixed side
      cov_red <- if (!length(covs)) "" else paste(" +", paste(unlist(
        lapply(seq_along(covs), function(i) if (!crossed[[i]]) covs[[i]] else
          c(covs[[i]], paste0(cols, ":", covs[[i]])))), collapse = " + "))
      rewrap(sprintf("(1 + %s%s | %s)", paste(cols, collapse = " + "),
                     cov_red, grp))
    }
    else if (structure == "chain_slope")
      ## the chain counterpart of an unstructured cell covariance: the same
      ## columns as the fixed part, so one set of declarations covers both
      rewrap(sprintf("(0 + %s%s | %s)", paste(chain_terms(spec), collapse = " + "),
                     cov_term(chain_terms(spec)), grp))
    else if (structure == "cells")
      rewrap(sprintf("(0 + %s%s | %s)", cn, cov_txt, grp))
    else {
      ## the grouping-chain submodel: parsimonious, but it constrains conditions
      ## sharing a stratum to equal correlation, where the cell form does not.
      ## It is always identified, and is the remedy offered for a partial
      ## structure, so it accepts one: only the declared variables and their
      ## ancestors contribute a rung.
      ##
      ## A rung is the ancestor path of a variable, not a prefix of the family
      ## vector. The two coincide while the family is a chain - p:a, p:a:b,
      ## p:a:b:c - but where a parent holds two children the prefix form would
      ## make one of them the coarser division and the other the finer one,
      ## which the design does not say and which the order of the declaration
      ## would then decide. Siblings therefore enter crossed, one variance
      ## each, with the whole structure as the finest rung.
      vars <- if (partial) intersect(fam, struct_vars) else fam
      vars <- unique(unlist(lapply(vars, function(v)
        c(nest_ancestors(spec, v), v))))
      ## Order the rungs by depth and then by name rather than by the order the
      ## nests were declared: siblings commute inside a grouping factor, so the
      ## same design should emit the same formula however it was written.
      depth_of <- function(v) length(nest_ancestors(spec, v))
      vars <- vars[order(vapply(vars, depth_of, 1L), vars)]
      rungs <- vapply(vars, function(v)
        paste(c(grp, nest_ancestors(spec, v), v), collapse = ":"), "")
      if (length(vars) > 1)
        rungs <- c(rungs, paste(c(grp, vars), collapse = ":"))
      rungs <- unique(rungs)
      rungs <- rungs[order(lengths(strsplit(rungs, ":")), rungs)]
      c(rewrap(sprintf("(1%s | %s)", cov_term(NULL), grp)),
        vapply(rungs, function(z) sprintf("(1 | %s)", z), ""))
    }
  }))
  paste(unique(out), collapse = " + ")
}

## --- the fit ---------------------------------------------------------------
nest_fit <- function(spec, mode = NULL,
                     random_structure = c("cells", "reduced", "chain", "as_declared"),
                     data = NULL, engine = "marginaleffects",
                     dry_run = FALSE, ..., priors = NULL,
                     prior_space = c("effects", "cells"), .env = parent.frame()) {
  ## `priors` and `prior_space` sit after the dots: both would otherwise be
  ## partial-matched by brms's own `prior`, which must reach the engine.
  prior_space <- match.arg(prior_space)
  rs_default <- missing(random_structure)
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
    mode <- match.arg(mode, c("cells", "reduced", "effects"))
    mode_note <- sprintf("## parameterization: %s (requested)", mode)
  }
  f <- paste(deparse(cell_formula(spec, mode)), collapse = " ")
  ## The emitted code has to be runnable, so it names the fitted columns rather
  ## than the structure they encode. A line saying what `cell` is costs nothing
  ## and is the difference between readable and cryptic.
  cell_note <- if (identical(mode, "cells"))
    c(sprintf("## `%s` is a factor whose levels are the %d realized conditions, so its",
              spec$cell_name, nrow(spec$cells)),
      "##   coefficients are their means. The formula crosses the structure fully,",
      "##   so that is the model it names, written compactly.")
  if (identical(mode, "effects") && identical(random_structure, "cells"))
    random_structure <- "chain_slope"   # match the random side to the fixed side
  ## Same argument, the other half of the model: a random structure that asks
  ## for less than the saturated one is restricted whatever the mean structure
  ## does, and the cell factor can express it no better there than in the mean.
  ## An explicit `random_structure` is left alone - the saturated covariance is
  ## a reasonable thing to want - so only the default follows the declaration.
  if (rs_default && (identical(mode, "reduced") || random_restricted(spec)))
    random_structure <- "reduced"
  ## The reduced form states the declared structure as one column per identified
  ## effect, computed once per realized cell and looked up by cell. It is the
  ## same kind of object the cell factor is - a design the data can inform in
  ## full - so nothing has to be held at zero, and the columns are ordinary
  ## numeric variables that every engine treats alike.
  aug_code <- NULL
  if (identical(mode, "reduced") || identical(random_structure, "reduced")) {
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
  re <- if (fit %in% c("lmer", "glmer", "clmm")) random_terms(spec, random_structure)
        else if (identical(fit, "brm")) random_terms(spec, random_structure)
  if (!is.null(re)) f <- paste(f, "+", re)
  fam <- if (!is.null(spec$family)) sprintf(", family = %s", spec$family) else ""
  dots <- as.list(substitute(list(...)))[-1]
  ## brms's own `prior` is the route to everything the translation does not
  ## touch - the random-effects sd and correlations, the thresholds, sigma. It
  ## is passed through untouched. A `class = "b"` entry is a different matter:
  ## those coefficients are the cells, not the original effects, so a prior
  ## written for the effects would mean something else there.
  user_prior <- NULL
  user_prior_txt <- NULL
  if ("prior" %in% names(dots)) {
    user_prior <- tryCatch(eval(dots$prior, .env), error = function(e) NULL)
    user_prior_txt <- paste(deparse(dots$prior), collapse = " ")
  }
  ## A `class = "b"` prior describes the mean structure. Written in the original
  ## variables - the default - it is translated into cell space and travels to
  ## Stan through stanvars, so the user states priors once, in brms syntax, and
  ## the translation is the package's business. `prior_space = "cells"` says the
  ## prior was written for the fitted coefficients and leaves it alone.
  if (inherits(user_prior, "brmsprior") && any(user_prior$class == "b") &&
      identical(mode, "cells") && identical(prior_space, "effects")) {
    if (!is.null(priors))
      stop("`prior` carries an entry for class \"b\" and `priors` supplies one ",
           "too, so the mean structure would be given two priors. Use one or ",
           "the other: `prior` in brms syntax, or a nestimand prior built with ",
           "nest_prior().")
    split <- prior_from_brms(spec, user_prior, space = prior_space)
    priors <- split$translated
    prior_name <- ".prior"
    message("the `class = \"b\"` prior was stated for the original variables ",
            "and has been translated to the ", spec$cell_name, " coefficients: ",
            "mu ~ N(A m, A D A'), reaching Stan through stanvars. ",
            "prior_for_estimand() shows what it implies for a contrast; ",
            "prior_space = \"cells\" leaves it as written instead.")
    user_prior_txt <- if (!is.null(split$other) && nrow(split$other))
      sprintf("brms::prior_string(%s, class = %s%s)",
              paste0("c(", paste(sprintf('"%s"', split$other$prior), collapse = ", "), ")"),
              paste0("c(", paste(sprintf('"%s"', split$other$class), collapse = ", "), ")"),
              "") else NULL
    dots$prior <- NULL
    user_prior <- split$other
  }
  ## The effects parameterization carries columns the data cannot inform, and
  ## brms cannot drop them the way `lm` does: it samples them, and the posterior
  ## is improper along them. They are known from the design, not chosen, so they
  ## are held at zero rather than left for the user to discover from a
  ## divergence. Which coefficients, and of which kind, is said aloud - passing
  ## `priors` explicitly overrides this.
  if (identical(fit, "brm") && identical(mode, "effects") && is.null(priors)) {
    user_b <- inherits(user_prior, "brmsprior") && any(user_prior$class == "b")
    cp <- chain_priors(spec, regularize = if (user_b) NULL else "normal(0, 5)")
    if (nrow(cp$table)) {
      ## Say where the coefficients are as well as how many, since the count is
      ## of design columns and the declaration was reported in cell dimensions:
      ## the two differ, and comparing them without that is confusing. A
      ## structural term is crossed with each covariate it interacts with, so
      ## one uninformative condition can carry several columns, and the random
      ## side repeats the whole set for every grouping factor.
      tb <- cp$table
      n_str <- sum(tb$kind == "structural zero")
      n_id  <- sum(tb$kind == "identification constraint")
      where <- c(
        if (any(tb$part == "fixed"))
          sprintf("%d in the mean structure", sum(tb$part == "fixed")),
        unlist(lapply(unique(stats::na.omit(tb$group)), function(g)
          sprintf("%d in the random structure for `%s`",
                  sum(tb$part == "random" & tb$group %in% g), g))))
      message("the effects parameterization states the structure as ",
              "coefficients rather than as cells, and carries ", nrow(tb),
              " of them the data cannot inform - ",
              paste(where, collapse = ", "), ". Of these, ", n_str,
              " are structural zeros - conditions the design does not realize - ",
              "and ", n_id, " are identification constraints, a coding choice ",
              "like a reference level, which leaves every estimand unchanged. ",
              "All are held at zero by constant(0) priors, so that the ",
              "posterior is proper. This counts coefficients, not cells: a ",
              "structural term is crossed with each covariate it interacts ",
              "with, so one condition that does not exist can carry several ",
              "columns. chain_priors(", spec_name, ") lists them and ",
              "show_code() prints the block; the rest take ",
              if (user_b) "the prior you supplied" else "normal(0, 5)",
              ", and passing `priors =` yourself replaces all of this.")
      priors <- cp
      ## The emitted code re-derives the block rather than naming an object the
      ## fitting environment held, so it has to carry the arguments it was
      ## derived with: with the default regularizer it would state a second
      ## `class = "b"` prior beside the user's, which brms refuses as a
      ## duplicate - and the code that was run would not be the code shown.
      prior_name <- sprintf("chain_priors(%s%s)", spec_name,
                            if (user_b) ", regularize = NULL" else "")
    }
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
  ## satisfy, so the Stan program fails to compile. The alternatives are a
  ## separate prior-only run, or leaving the prior in cell space.
  if (!is.null(priors) && "sample_prior" %in% names(dots)) {
    sv <- tryCatch(eval(dots$sample_prior, .env), error = function(e) NULL)
    if (isTRUE(sv) || identical(sv, "yes"))
      stop("`sample_prior = \"yes\"` cannot be combined with a translated prior. ",
           "The translation puts a multivariate normal on the coefficients, and ",
           "brms draws one prior sample per class into a scalar, so Stan reports ",
           "an ill-typed assignment - `real prior_b = multi_normal_rng(...)`. ",
           "Two ways round it: run the prior separately with ",
           "sample_prior = \"only\", which works and is what the prior check in ",
           "tests/test_brms.R does; or state the prior in cell space with ",
           "prior_space = \"cells\", where it is univariate per coefficient and ",
           "brms can sample it.")
  }
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
    other <- user_prior_txt
    prior_txt <- if (is.null(other))
      sprintf(", prior = chain_prior_object(%s)", prior_name)
    else sprintf(", prior = c(chain_prior_object(%s), %s)", prior_name, other)
    prior_note <- c(priors$code[1:4],
      "## Held at zero below; every other coefficient is estimated as usual.")
    priors_obj <- priors
    priors <- NULL
  } else priors_obj <- NULL
  if (!is.null(priors)) {
    if (!inherits(priors, "nestimand_prior"))
      stop("`priors` must be a nestimand_prior object, as returned by nest_prior().")
    if (!identical(fit, "brm"))
      stop("translated priors apply to the brms engine; the declared engine is `",
           fit, "`, which has no prior to state.")
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
      if (identical(random_structure, "reduced"))
        c(sprintf(paste("## - the structure declared after the bar, varying by group:",
                        "%d columns, one\n##   for each effect that structure",
                        "names over the realized conditions."),
                  length(reduced_columns(spec, declared_terms(spec,
                    unlist(lapply(bar_terms_of(spec$random_original), function(b)
                      tryCatch(attr(stats::terms(stats::as.formula(paste("~", b$lhs))),
                                    "term.labels"), error = function(e) character(0))))))) + 1L),
          "##   Written over the original factors it would carry columns the data",
          "##   cannot inform; written as the cell factor it would let everything",
          "##   vary, including terms the bar does not name.",
          "##   random_structure = \"cells\" gives that larger covariance if wanted.")
      else
        c("## - a random slope written over the original factors carries columns the",
          "##   data cannot inform: those for conditions that do not exist, and a",
          "##   further set reconstructable from the others. The cell form has one",
          "##   column per realized condition, so every parameter is estimable and",
          "##   nothing has to be held at zero. See chain_priors() for the",
          "##   alternative, which keeps the original factors and declares the",
          "##   uninformative columns instead."))
  relevel_code <- if (identical(mode, "effects"))
    sentinel_relevel_code(spec, data_name)
  code <- c(sprintf("## nestimand %s -- fit", nestimand_build), lib, mode_note,
            cell_note, relevel_code, aug_code, re_note, prior_note,
            sprintf("m <- %s(%s%s, data = %s%s%s)", fn, f, fam, data_name,
                    prior_txt, dots_txt))
  if (isTRUE(dry_run))
    return(structure(paste(code, collapse = "\n"), class = "nestimand_code",
                     nestimand_code = code))
  env <- new.env(parent = .env)
  assign(spec_name, spec, envir = env)
  ## `prior_name` may be a call - `chain_priors(sp)` - so that the emitted code
  ## stands on its own; there is then nothing to assign, the call being
  ## evaluated where it stands.
  if (grepl("^[.A-Za-z][.A-Za-z0-9._]*$", prior_name)) {
    if (!is.null(priors)) assign(prior_name, priors, envir = env)
    if (!is.null(priors_obj)) assign(prior_name, priors_obj, envir = env)
  }
  if (!exists(data_name, envir = env, inherits = TRUE))
    assign(data_name, data, envir = env)
  ## the reduced form builds its data in the emitted code, so the frame it
  ## builds it from has to be reachable there too
  if (identical(mode, "reduced") && grepl("^[.A-Za-z][.A-Za-z0-9._]*$", raw_data_name) &&
      !exists(raw_data_name, envir = env, inherits = TRUE))
    assign(raw_data_name, raw_data, envir = env)
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
