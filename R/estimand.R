## nestimand: estimands, with the code that produced them ------------------
## The function computes; the code view is a by-product, not a reconstruction.
## Every estimand is assembled as text, then evaluated, so `show_code()` cannot
## drift from what was run: it is the same object. Non-core arguments in `...`
## are deparsed into that text and so reach the destination function - brms,
## marginaleffects, emmeans - unaltered, and appear in the saved code.

estimand <- function(model, target, policy = "equal", at = NULL,
                     contrast = c("pairwise", "reference", "sequential", "within",
                                  "interaction"),
                     route = c("g_computation", "cells"),
                     weights = NULL, type = NULL, subsample = NULL,
                     data = NULL, bounds = TRUE, self_check = TRUE,
                     dry_run = FALSE, by = NULL, ..., spec = NULL,
                     .env = parent.frame(), .restrict = NULL) {
  ## The declaration travels with a fit from nest_fit(); `spec` is needed only
  ## for a model fitted by calling the engine directly.
  spec_expr <- substitute(spec)
  recovered <- is.null(spec)
  spec <- resolve_spec(model, spec)
  check_model_spec(model, spec)
  contrast <- if (identical(contrast, "interaction")) "interaction"
              else match.arg(contrast)
  route <- match.arg(route)
  wq <- substitute(weights)
  ## One argument names the quantity, in the engines' own vocabulary; which
  ## machinery computes it is the package's business. On the linear predictor
  ## the contrast is c'b and is taken from the coefficients: exact, and it works
  ## where an engine's ordinal support does not. Every other quantity goes to
  ## the prediction function. An ordinal family defaults to the linear
  ## predictor, giving one number per contrast rather than one per outcome
  ## category; elsewhere the two coincide.
  ## `eta` is nestimand's own name for the linear predictor - the mean
  ## structure, thresholds and link excluded - computed from the coefficients as
  ## c'b. It is deliberately not a name marginaleffects uses for anything:
  ## "link", "latent", "linpred", "linear.predictor", "lp" and "xb" all appear
  ## in its dictionary for one class or another, so any of them would be a
  ## claim on its vocabulary. Every other value is passed through untouched.
  linear_types <- "eta"
  ## A type the user typed is theirs: if the engine does not accept it for this
  ## model, it is passed through and the engine says so, rather than the package
  ## quietly recognising a name marginaleffects would reject. A type the package
  ## chooses for itself needs no such check - and must not have one, since some
  ## classes the engine cannot predict from at all.
  type_passthrough <- FALSE
  type_supplied <- !is.null(type)
  if (is.null(type)) type <- if (has_thresholds(spec)) "eta" else "response"
  scale <- if (type %in% linear_types) "latent" else "response"

  if (identical(scale, "latent") && identical(contrast, "within"))
    stop("contrast = \"within\" is computed through the prediction route and ",
         "has no linear-map form here; use scale = \"response\".")
  if (!inherits(spec, "nesting_spec"))
    stop("`spec` must be a nesting_spec object, as returned by nesting_spec().")
  ## Several targets may be named at once, bare or quoted: each is computed in
  ## turn and returned together, so one call covers a table of estimands.
  tg <- substitute(target)
  multi <- NULL
  ## The target is expanded as a formula expands its right-hand side, by the
  ## formula machinery itself rather than by walking the parse tree here: `a * b`
  ## is a, b, a:b; `a * (b + c)` is a, b, c, a:b, a:c; `a:b` is the interaction
  ## alone. Every bracketing and combination of operators then comes out right
  ## because R settles it, which a hand-rolled walk kept failing to do.
  if (is.call(tg) && as.character(tg[[1]]) %in% c("*", ":", "+", "(")) {
    labs <- attr(stats::terms(stats::as.formula(
      paste("~", paste(deparse(tg), collapse = " ")))), "term.labels")
    vs <- labs[!grepl(":", labs, fixed = TRUE)]
    inter <- lapply(labs[grepl(":", labs, fixed = TRUE)],
                    function(z) strsplit(z, ":", fixed = TRUE)[[1]])
    if (!length(vs) && length(inter) == 1L) {
      target <- inter[[1]]
      contrast <- "interaction"
      tg <- NULL
    } else {
      cl <- match.call()
      has_hyp <- "hypothesis" %in% names(as.list(substitute(list(...)))[-1])
      ## Each part is a separate call and would repeat the same notes, so they
      ## are collected and said once.
      seen <- character(0)
      once <- function(expr) withCallingHandlers(expr, message = function(m) {
        txt <- conditionMessage(m)
        if (txt %in% seen) invokeRestart("muffleMessage")
        seen[[length(seen) + 1L]] <<- txt
      })
      assign("star", TRUE, envir = nestimand_env)
      on.exit(assign("star", FALSE, envir = nestimand_env), add = TRUE)
      out <- lapply(vs, function(k) { c2 <- cl; c2$target <- k; once(eval(c2, .env)) })
      ## The interactions are defined by their own comparison matrices, so a
      ## `hypothesis` meant for the marginal contrasts does not apply to them -
      ## but it need not be dropped either: on grouped output an interaction
      ## is formed within each group, which is the same restriction the
      ## hypothesis was asking for.
      parts <- c(
        if (has_hyp)
          "your `hypothesis` sets the comparisons for the marginal contrasts, and their labels and direction come from marginaleffects",
        if (!identical(policy, "equal") && length(inter))
          "the policy weights the marginal contrasts and not the interactions, which use only cells that exist",
        if (length(inter))
          "each interaction has its own comparison matrix, formed within each group where the output is grouped")
      message("`", paste(deparse(tg), collapse = " "), "` expands as a formula ",
              "does: ", length(labs), " results - ",
              paste(labs, collapse = ", "),
              if (length(parts)) paste0(": ", paste(parts, collapse = "; ")), ".")
      out <- c(out, lapply(inter, function(k)
        once({ c2 <- cl; c2$target <- k
               c2$contrast <- "interaction"
               c2$hypothesis <- NULL; eval(c2, .env) })))
      names(out) <- c(vs, vapply(inter, paste, "", collapse = ":"))
      return(collect_estimands(out))
    }
  }
  if (!is.null(tg) && is.call(tg) && identical(tg[[1]], as.name("c")))
    multi <- vapply(as.list(tg)[-1], function(z)
      if (is.character(z)) z else deparse(z), "")
  else if (tryCatch(is.character(target) && length(target) > 1L &&
                    !identical(contrast, "interaction"), error = function(e) FALSE))
    multi <- target
  if (!is.null(multi)) {
    cl <- match.call()
    out <- lapply(multi, function(k) { cl$target <- k; eval(cl, .env) })
    names(out) <- multi
    return(collect_estimands(out))
  }
  if (is.name(tg)) {
    lit <- deparse(tg)
    target <- if (lit %in% spec$cell_vars) lit else
      tryCatch({ v <- eval(tg, .env)
                 if (is.character(v) && length(v) == 1L) v else lit },
               error = function(e) lit)
  }
  if (!all(target %in% spec$cell_vars)) {
    miss <- setdiff(target, spec$cell_vars)
    in_data <- miss[miss %in% names(spec$data)]
    categorical <- in_data[vapply(in_data, function(v)
      is.factor(spec$data[[v]]) || is.character(spec$data[[v]]), TRUE)]
    ## A categorical variable the formula crosses with the structure is folded
    ## into the cell factor without being declared, so one that is still a
    ## covariate is either additive - its effect declared common to every
    ## condition - or ordered, or numeric. Which it is decides the remedy.
    additive <- categorical[vapply(categorical, function(v)
      !any(vapply(strsplit(spec$term_labels, ":"), function(vs)
        v %in% vs && any(spec$cell_vars %in% vs), TRUE)), TRUE)]
    ordered_f <- setdiff(categorical[vapply(categorical, function(v)
      is.ordered(spec$data[[v]]), TRUE)], additive)
    stop("`", paste(miss, collapse = "`, `"),
         "` is not among the declared categorical design ",
         "variables (", paste(spec$cell_vars, collapse = ", "), "), so it has ",
         "no conditions to average over and no policy to apply.",
         if (length(additive))
           paste0(" `", paste(additive, collapse = "`, `"), "` enters ",
                  "additively, which declares its effect common to every ",
                  "condition; cross it with the structure in the formula - ",
                  "`", spec$cell_vars[1], " * ", additive[1],
                  "` - and it becomes part of the categorical design.")
         else if (length(ordered_f))
           paste0(" `", paste(ordered_f, collapse = "`, `"), "` is an ordered ",
                  "factor, which is kept as a covariate because its contrasts ",
                  "say it is meant as a quantity. Name it in `nests` - ",
                  "nests = c(\"", paste(c(spec_nests(spec), ordered_f[1]),
                                        collapse = "\", \""),
                  "\") - to fold it into the cell factor instead.")
         else " Contrasts of a covariate do not cross the structural boundary and need no policy; compute them directly.")
  }

  ## A marginal contrast of a nested variable is computed only over the strata
  ## in which that variable varies. Elsewhere it holds a single level, and a
  ## comparison against that level leaves the variable's own levels and compares
  ## strata instead - augmented chords against root-position triads, say, under
  ## an inversion label.
  ## an interaction exists only where every one of its variables varies, so the
  ## restrictions of the targets are intersected rather than the deepest of them
  ## taken as standing for the rest
  deg <- if (!identical(contrast, "within"))
    degenerate_strata_multi(spec, target)
  restricted <- !is.null(deg) && length(deg$drop)
  if (!restricted) deg <- NULL
  ## A recursive call from the `by` branch below carries the cells of its group.
  if (!is.null(.restrict)) deg <- .restrict

  ## `by` names design variables whose levels the estimand is computed within,
  ## as it names the grouping of an average in marginaleffects: the target's
  ## contrasts are formed inside each group, with the policy weighted over that
  ## group's conditions alone, and the groups are returned stacked and labelled.
  ## `contrast = "within"` is the same thing with the strata taken from a nested
  ## target's own ancestors, which a variable crossed with the structure has
  ## none of - hence naming them.
  bq <- substitute(by)
  byv <- if (is.null(bq)) character(0) else {
    bv <- tryCatch(eval(bq, .env), error = function(e) NULL)
    if (is.character(bv)) bv else all.vars(bq)
  }
  if (length(byv)) {
    bad <- setdiff(byv, spec$cell_vars)
    if (length(bad))
      stop("`by` groups the estimand by declared design variables; `",
           paste(bad, collapse = "`, `"), "` is not one (the design has ",
           paste(spec$cell_vars, collapse = ", "), "). A covariate has no ",
           "conditions to group by; cross it with the structure, or declare it.")
    if (length(intersect(byv, target)))
      stop("`by` names `", paste(intersect(byv, target), collapse = "`, `"),
           "`, which is also the target: a variable cannot be contrasted within ",
           "its own levels. For the two together, name them as one target - ",
           "`", paste(unique(c(target, byv)), collapse = " * "), "`.")
    base <- if (is.null(deg)) spec$cells else
      spec$cells[deg_key(spec$cells, deg$vars) %in% deg$keep, , drop = FALSE]
    gkey <- deg_key(base, byv)
    tgv <- as.character(base[[target[length(target)]]])
    groups <- unique(gkey)
    varies <- vapply(groups, function(g)
      length(unique(tgv[gkey == g])) > 1, TRUE)
    if (!any(varies))
      stop("`", target[length(target)], "` takes one level in every group of `",
           paste(byv, collapse = "`, `"), "`, so there is no contrast to form ",
           "in any of them.")
    if (any(!varies))
      message("`", target[length(target)], "` does not vary in ",
              sum(!varies), " group(s) - ", paste(groups[!varies], collapse = ", "),
              " - which are left out: a comparison there would leave the ",
              "variable's own levels.")
    cl <- match.call()
    parts <- lapply(groups[varies], function(g) {
      c2 <- cl
      c2$target <- target
      c2$by <- NULL
      vars_g <- union(if (is.null(deg)) character(0) else deg$vars, byv)
      c2$.restrict <- list(
        vars = vars_g,
        keep = unique(deg_key(base[gkey == g, , drop = FALSE], vars_g)),
        drop = "", by = byv, group = g)
      e <- eval(c2, .env)
      list(g = g, e = e, lab = base[gkey == g, byv, drop = FALSE][1, , drop = FALSE])
    })
    ## a dry run returns the scripts, one block per group, each of which stands
    ## on its own: the groups differ only in the cells they restrict to
    code <- unlist(lapply(parts, function(z)
      c(sprintf("## group: %s", z$g),
        if (inherits(z$e, "nestimand_code")) attr(z$e, "nestimand_code")
        else attr(z$e, "nestimand")$code)))
    if (isTRUE(dry_run))
      return(structure(paste(code, collapse = "\n"), class = "nestimand_code",
                       nestimand_code = code))
    res <- do.call(rbind, lapply(parts, function(z) {
      d <- as.data.frame(z$e)
      cbind(z$lab[rep(1L, nrow(d)), , drop = FALSE], d, row.names = NULL)
    }))
    m1 <- attr(parts[[1]]$e, "nestimand")
    bnd <- lapply(parts, function(z) {
      b <- attr(z$e, "nestimand")$bounds
      if (is.null(b)) NULL else
        cbind(z$lab[rep(1L, nrow(b)), , drop = FALSE], b, row.names = NULL)
    })
    bnd <- if (any(vapply(bnd, is.null, TRUE))) NULL else do.call(rbind, bnd)
    st <- vapply(parts, function(z)
      attr(z$e, "nestimand")$self_check$status %||% NA_character_, "")
    return(structure(res,
      nestimand = list(build = nestimand_build, target = target, by = byv,
                       policy = policy, at = at, contrast = contrast,
                       scale = m1$scale, type = m1$type, route = route,
                       groups = groups[varies], bounds = bnd,
                       self_check = if (all(is.na(st))) NULL else
                         list(status = if (all(st %in% "passed")) "passed"
                              else paste(unique(st[!is.na(st)]), collapse = "/")),
                       code = code),
      nestimand_bounds = bnd,
      class = c("nestimand_estimand", class(res))))
  }

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
  ## Under an identity link the response scale and the linear predictor are the
  ## same quantity, so the contrast can be taken from the coefficients: the
  ## averaging may be done on the design matrix rather than on one prediction
  ## per row per draw. Same answer, far less work. Gated on the link rather
  ## than the family - gaussian(link = "log") is not identity - and set aside
  ## when a hypothesis of the user's own needs the engine's own output.
  identity_link <- identical(model_link(model, spec), "identity")
  ## Any argument meant for the prediction function is a request to go through
  ## it, so the shortcut stands aside rather than refusing them.
  dot_names <- names(as.list(substitute(list(...)))[-1])
  ## The shortcut applies to the quantity the package itself asked for. A type
  ## the user named is theirs: it goes to the engine, which accepts it or does
  ## not.
  linear_shortcut <- identical(type, "response") && identity_link &&
    all(dot_names %in% c("conf_level", "ndraws")) &&
    !identical(contrast, "within") && linear_map_available(model, spec, data)
  if (linear_shortcut) scale <- "latent"
  ## G-computation averages over the units in the data. Where that is too many
  ## to evaluate - a nonlinear scale multiplies the rows by the draws and the
  ## outcome categories - a random sample of them estimates the same average,
  ## with Monte Carlo error in place of exactness. The rows are drawn once and
  ## written into the code, so the result can be reproduced.
  subsample_txt <- NULL
  if (!is.null(subsample) && identical(route, "g_computation")) {
    if (subsample < nrow(data)) {
      ## the rows are drawn from a recorded seed rather than listed: a literal
      ## index vector reproduces just as well but becomes unusable as code
      ## once the sample is large
      sub_seed <- sample.int(.Machine$integer.max, 1)
      set.seed(sub_seed)
      idx <- sort(sample.int(nrow(data), subsample))
      data <- data[idx, , drop = FALSE]
      subsample_txt <- sub_seed
      message("averaging over a random ", subsample, " of ",
              format(nrow(spec$data), big.mark = ","), " units rather than all ",
              "of them: the same population average, with Monte Carlo error. ",
              "The rows drawn are written into the code.")
    }
  }


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
  ## `ndraws` is the engine's own name for thinning a posterior, and it works
  ## on either route: passed through to the prediction function, or applied to
  ## the draws matrix here. The other dimension of the same problem as
  ## `subsample`, and usually the cheaper one, since a posterior has far more
  ## draws than the quantities they inform.
  draws_txt <- ""
  if (!is.null(dots$ndraws) && identical(scale, "latent")) {
    draws_txt <- sprintf(", ndraws = %s", deparse(dots$ndraws))
    dots$ndraws <- NULL
  }

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

  ## The other dimension of the same problem: a posterior has as many draws as
  ## it has, and the response scale evaluates the model once per draw. Thinning
  ## them estimates the same posterior summaries with Monte Carlo error, and is
  ## usually the cheaper of the two levers, since draws are far more numerous
  ## than the quantities they inform.
  if (!is.null(dots$ndraws) && inherits(model, "brmsfit")) {
    nd <- tryCatch(brms::ndraws(model), error = function(e) NA_integer_)
    k <- tryCatch(eval(dots$ndraws, .env), error = function(e) NA_integer_)
    if (!is.na(nd) && !is.na(k) && k < nd)
      message("using ", k, " of ", format(nd, big.mark = ","),
              " posterior draws. The posterior is unchanged; only how precisely ",
              "it is summarized - the Monte Carlo error on a posterior mean is ",
              "about 1/sqrt(", k, ") = ", signif(1 / sqrt(k), 2),
              " posterior standard deviations.")
  }
  ## A mixed fit predicts for a particular group unless told otherwise, and the
  ## estimands here are population-level. The exclusion is stated rather than
  ## assumed, and in the spelling the engine expects: brms and the frequentist
  ## engines disagree, and the wrong one is silently ignored.
  ## Whether the group-level effects enter the predictions is the engine's
  ## default to make, not the package's: `re.form` and `re_formula` are passed
  ## through untouched, as every other engine argument is. What the package
  ## does is say which default is in force, since it decides whether the
  ## estimand describes a typical group or the sampled ones - and the answer is
  ## not obvious from the call.
  re_arg <- if (length(grouping_vars(spec)) &&
                spec$fit %in% c("lmer", "glmer", "clmm", "brm"))
    (if (identical(spec$fit, "brm")) "re_formula" else "re.form")
  re_note <- NULL
  if (!is.null(re_arg) && !re_arg %in% names(dots) && !identical(scale, "latent")) {
    re_note <- c(
      sprintf("## no %s given, so the engine's default stands: the group", re_arg),
      "## deviations enter the predictions, and the estimand averages over the",
      sprintf("## groups in the data. %s = NA holds them at zero instead - the", re_arg),
      "## average group - as a coefficient is read with the other terms at zero.")
    if (!in_star())
      message("no `", re_arg, "` given, so ",
              if (identical(spec$fit, "brm")) "brms" else "the engine",
              "'s default stands: the group deviations enter the predictions, ",
              "so the estimand averages over the groups in the data. `",
              re_arg, " = NA` holds them at zero instead - the average group - ",
              "as a coefficient is read with the other terms at zero. On a ",
              "linear scale the two give the same estimate, the deviations ",
              "averaging to zero; under a link they do not.")
  }

  ## An ordinal family has no single response scale, so `scale` chooses it and
  ## the engine's own spelling is supplied here: the user should not have to
  ## know that clm calls it "prob" and brms calls it "response".
  if (identical(scale, "response")) {
    ## the engines spell the expected-outcome scale differently, and the
    ## accepted set has changed between releases, so "response" is resolved
    ## against this fit rather than assumed
    dots$type <- if (identical(type, "response") && has_thresholds(spec))
      ordinal_response_type(model, spec, data) else type
    ## `prediction` draws from the posterior predictive. On an ordinal family
    ## those draws are category codes, so averaging them treats the rating
    ## scale as an interval one - a modelling assumption the ordinal family was
    ## chosen to avoid.
    if (has_thresholds(spec) && identical(type, "prediction"))
      message("ordinal fit: type = \"prediction\" draws from the posterior ",
              "predictive, and on this family those draws are category codes. ",
              "Averaging them treats the rating scale as an interval one, ",
              "which is the assumption the ordinal family was chosen to avoid. ",
              "type = \"response\" gives the probability of each category, and ",
              "type = \"eta\" the linear predictor.")
    ## The response scale evaluates the model at every grid row, for every
    ## posterior draw, for every outcome category. On the G-computation grid
    ## that is the whole data set times the cells, and the work can run to
    ## hours; the cells route asks the same question of one row per condition.
    nrows <- if (identical(route, "g_computation"))
      nrow(data) * nrow(if (is.null(deg)) spec$cells else
        spec$cells[deg_key(spec$cells, deg$vars) %in% deg$keep, ,
                   drop = FALSE]) else nrow(spec$cells)
    ## the draws actually used, not the draws the model holds
    ndraw <- if (!is.null(dots$ndraws)) dots$ndraws
             else if (inherits(model, "brmsfit"))
               tryCatch(brms::ndraws(model), error = function(e) NA_integer_)
             else NA_integer_
    ncat0 <- tryCatch(nlevels(data[[spec$outcome]]), error = function(e) NA_integer_)
    work <- nrows * (if (is.na(ndraw)) 1 else ndraw) * (if (is.na(ncat0)) 1 else ncat0)
    if (is.finite(work) && work > 5e6)
      message("this will take a while: the response scale evaluates the model ",
              "at ", format(nrows, big.mark = ","), " grid rows",
              if (!is.na(ndraw)) paste0(" for each of ", format(ndraw, big.mark = ","),
                                        " posterior draws") else "",
              if (!is.na(ncat0)) paste0(" and each of ", ncat0, " outcome categories")
              else "", ". The averaging cannot be done first here, as it can on ",
              "the link scale, because the link stands between the model and ",
              "the mean.",
              if (is.null(subsample))
                paste0(" `subsample = 500` averages over a random 500 units",
                       " instead, which estimates the same quantity;")
              else "",
              if (is.null(dots$ndraws) && !is.na(ndraw))
                paste0(" `ndraws = 500` uses fewer posterior draws, likewise;")
              else "",
              " route = \"cells\" answers a different question - the effect at ",
              "average covariates - which on a nonlinear scale is not the ",
              "population average.")

    ## Where the engine's own quantity is the linear predictor, `eta` gives the
    ## same numbers from the coefficients, without a prediction per row - and
    ## per draw, on a posterior. Worth saying once; the choice stays the user's.
    ## Only where the engine will actually produce it: if the type is about to
    ## be refused, an equivalence is beside the point and the refusal is the
    ## information.
    same_as_eta <- identical(type, "link") &&
      (identity_link || (has_thresholds(spec) && inherits(model, "brmsfit"))) &&
      engine_accepts(model, spec, type, data)
    if (isTRUE(same_as_eta))
      message("type = \"eta\" gives these same numbers from the coefficients, ",
              "without evaluating the model at every grid row",
              if (inherits(model, "brmsfit")) " for every draw" else "",
              ". The linear predictor here is the mean structure, which is what ",
              "eta computes directly.")

    ## One note, and only where it can still change what the user does. With a
    ## hypothesis of their own the comparison set is already theirs.
    ## Only for a plain pairwise comparison: an interaction already forms its
    ## contrasts within each group, so there is nothing to recommend.
    if (has_thresholds(spec) && identical(type, "response") &&
        !"hypothesis" %in% names(dots) && identical(contrast, "pairwise")) {
      ncat <- tryCatch(nlevels(data[[spec$outcome]]), error = function(e) NA_integer_)
      nlev <- length(unique(as.character(data[[target[1]]])))
      message("ordinal fit on the response scale: every condition becomes ",
              if (is.na(ncat)) "one number per outcome category"
              else paste(ncat, "numbers, one per outcome category"),
              ", so comparing them all gives ",
              if (is.na(ncat)) "many more contrasts than expected"
              else paste0(choose(nlev * ncat, 2), " contrasts rather than ",
                          choose(nlev, 2)),
              ", most of them across categories rather than within one. For ",
              "comparisons within each category, hypothesis = ~ pairwise | ",
              "group; for one number per contrast, type = \"link\".")
    }
  }

  dots_txt <- if (length(dots))
    paste0(", ", paste(mapply(function(nm, v) {
      v <- paste(deparse(v), collapse = " ")
      if (nzchar(nm)) paste0(nm, " = ", v) else v
    }, names(dots), dots), collapse = ", ")) else ""

  if (identical(contrast, "interaction")) {
    if (length(target) < 2)
      stop("an interaction contrast needs at least two targets, e.g. ",
           "estimand(model, chord_type:inversion).")
    ## A policy weights the versions of a compound condition. An interaction
    ## uses only cells that exist, so it crosses no boundary and there is
    ## nothing for a policy to weight. Rather than refuse the call - which would
    ## also refuse the marginal contrasts of `a * b`, where the policy does
    ## apply - it is dropped here, and said so.
    if (!identical(policy, "equal")) {
      if (!in_star())
      message("`policy = \"", if (is.character(policy)) policy else "supplied",
              "\"` does not apply to an interaction contrast: every cell it uses ",
              "exists, so no boundary is crossed and there are no versions to ",
              "weight. The interaction was computed without it; any marginal ",
              "contrasts requested alongside it keep the policy.")
      policy <- "equal"
      policy_dropped <- TRUE
    }
    bounds <- FALSE
  }
  ## A hypothesis of the user's own replaces the one the contrast implies:
  ## they are two answers to the same question, and the engine would refuse
  ## both. The direction convention and the contrast argument go with it.
  ## The linear predictor gives one value per condition, so the comparisons are
  ## the ones `contrast` names and there is nothing for a hypothesis to group.
  if (identical(scale, "latent")) {
    if ("hypothesis" %in% names(dots))
      stop("the linear predictor gives one value per condition, so there are ",
           "no groups to compare within and a `hypothesis` cannot be applied ",
           "here: the comparisons come from `contrast` - pairwise, reference, ",
           "sequential, or interaction. A grouped hypothesis such as ",
           "~ pairwise | group belongs with a quantity that is grouped, which ",
           "on an ordinal fit means type = \"response\".")
    ## eta is the contrast with the group deviations at zero. Asking for the
    ## sampled-group quantity is asking for a different spread, which this
    ## route cannot produce: it never touches the group-level parameters.
    re_named <- intersect(names(dots), c("re.form", "re_formula"))
    ## On a brms fit the group-level parameters are in the draws, so eta can
    ## average over the sampled groups after all: the estimate is unchanged and
    ## the posterior widens. Elsewhere there is nothing to average over.
    eta_re <- ""
    if (length(re_named) && inherits(model, "brmsfit")) {
      val <- tryCatch(eval(dots[[re_named[1]]], .env), error = function(e) NA)
      if (is.null(val)) {
        eta_re <- ", re_formula = NULL"
        dots[[re_named[1]]] <- NULL
        re_named <- character(0)
        message("averaging over the groups in the data on the linear ",
                "predictor: the group-level draws are added before the ",
                "contrast is taken, so the estimate is unchanged and the ",
                "posterior is wider than for the average group.")
      }
    }
    if (length(re_named))
      stop("`", re_named[1], "` and type = \"eta\" ask for different things: ",
           "eta holds the group deviations at zero - the average group - and ",
           "is computed from the fixed effects alone. ",
           if (inherits(model, "brmsfit"))
             paste0("On a brms fit that distinction is real: the draws carry ",
                    "the group-level parameters, so ", re_named[1], " = NULL ",
                    "gives a wider posterior. Ask for a predicted quantity - ",
                    "type = \"response\" - to get it.")
           else
             paste0("On this engine marginaleffects builds the interval from ",
                    "the fixed-effect covariance whichever ", re_named[1],
                    " is given, so the sampled-group spread is not available ",
                    "on any route: dropping the argument loses nothing."))
    other <- setdiff(names(dots), "conf_level")
    if (length(other))
      stop("`", paste(other, collapse = "`, `"), "` ",
           if (length(other) > 1) "are arguments" else "is an argument",
           " of marginaleffects::avg_predictions, which the linear predictor ",
           "does not go through: the contrast is c'b, taken from the ",
           "coefficients. `conf_level` and `ndraws` work on either route; for ",
           "the rest, ask for a quantity that is computed by prediction, such ",
           "as type = \"response\".")
  }
  user_hyp <- "hypothesis" %in% names(dots)
  if (user_hyp) {
    if (identical(contrast, "interaction"))
      stop("`hypothesis` and contrast = \"interaction\" both define which ",
           "comparisons are formed. The interaction is built as a matrix of ",
           "differences of differences; supply your own hypothesis instead, or ",
           "drop it and let the interaction stand.")
    if (!in_star())
    message("using your `hypothesis` instead of contrast = \"", contrast,
            "\". The comparisons, their names, and which way round each ",
            "subtraction goes are all as marginaleffects returns them.")
  }
  if (!exists("policy_dropped", inherits = FALSE)) policy_dropped <- FALSE
  sub_code <- NULL
  if (!is.null(subsample_txt)) {
    sub_code <- c(
      sprintf("## averaging over a random %d of %s units, drawn from a recorded",
              subsample, format(nrow(spec$data), big.mark = ",")),
      "## seed so that the same rows come back",
      sprintf("set.seed(%d)", subsample_txt),
      sprintf(".rows <- sort(sample.int(nrow(%s), %d))", data_name, subsample),
      sprintf(".data <- %s[.rows, ]", data_name))
    data_name <- ".data"
  }
  code <- estimand_code(spec, target, policy, at, contrast, dots_txt,
                        model_name, spec_name, data_name, bounds, scale,
                        deg, route, weights_txt, re_note, user_hyp, draws_txt,
                        linear_shortcut, if (exists("eta_re", inherits = FALSE)) eta_re else "")
  if (!is.null(sub_code)) {
    i <- grep("^library\\(", code)
    at_line <- if (length(i)) max(i) else 1L
    code <- append(code, sub_code, after = at_line)
  }
  ## if the model was fitted by nest_fit(), its call travels with it, so the
  ## code view is the whole pipeline rather than its second half
  ## The fit belongs in the code view, so that what is shown is a whole
  ## analysis - but it must not be re-run here: the model is already fitted and
  ## has been passed in. Only the estimand lines are evaluated.
  run_code <- code
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
  if (grepl("^[.A-Za-z][.A-Za-z0-9._]*$", data_name) &&
      !exists(data_name, envir = env, inherits = TRUE))
    assign(data_name, data, envir = env)
  out <- tryCatch(eval(parse(text = paste(run_code, collapse = "\n")), envir = env),
    error = function(err) {
      m <- conditionMessage(err)
      ## the engine's refusal here is opaque, and the cause is the per-category
      ## output rather than anything about the contrast that was asked for
      if (has_thresholds(spec) && identical(scale, "response") &&
          grepl("pairwise|reference|sequential", m))
        stop("this version of marginaleffects will not form a `", contrast,
             "` comparison of per-category probabilities: on an ordinal fit ",
             "type = \"", dots$type, "\" returns one value per outcome ",
             "category, and the comparison is defined over a single value per ",
             "condition. Use type = \"link\" for one number per contrast, or ",
             "add `by` and a hypothesis of your own to compare particular ",
             "categories.\n  The engine said: ", m, call. = FALSE)
      stop(err)
    })

  ## A Bayesian fit is summarized from its posterior on this route too: the
  ## engine returns draws, so the probability of direction is available, and a
  ## Wald statistic computed from a posterior covariance is not meaningful.
  out <- add_posterior_summary(out, model)

  check <- if (isTRUE(self_check))
    reorder_check(model, spec, target, policy, at, contrast, dots_txt, data, scale,
                  route)
  else NULL

  structure(out,
            nestimand = list(build = nestimand_build, target = target,
                             policy = policy, at = at, contrast = contrast,
                             scale = scale, type = type, route = route,
                             policy_dropped = policy_dropped,
                             bounds = attr(out, "nestimand_bounds"),
                             self_check = check,
                             code = code),
            class = c("nestimand_estimand", class(out)))
}

## Which name the installed engine gives the response scale for this fit. The
## accepted set has differed between marginaleffects versions and between model
## classes, so it is probed on two rows rather than assumed.
ordinal_response_type <- function(model, spec, data) {
  cand <- if (identical(spec$fit, "brm")) c("response", "prob")
          else c("prob", "response")
  g <- utils::head(data, 2)
  for (ty in cand) {
    ok <- tryCatch({
      marginaleffects::predictions(model, newdata = g, type = ty); TRUE
    }, error = function(e) FALSE)
    if (ok) return(ty)
  }
  ## Some classes the prediction machinery does not handle at all - clmm among
  ## them - and no choice of type will help. The linear predictor is unaffected,
  ## being computed from the coefficients.
  unsupported <- tryCatch({
    marginaleffects::predictions(model, newdata = g); FALSE
  }, error = function(e) grepl("not supported", conditionMessage(e)))
  if (isTRUE(unsupported))
    stop("marginaleffects does not support models of class ",
         paste(class(model), collapse = "/"), ", so no quantity computed by ",
         "prediction is available for this fit - the response scale included. ",
         "type = \"eta\" works: the contrast is taken from the coefficients ",
         "and their covariance, which this model does provide. It is also the ",
         "default for an ordinal family.")
  stop("no response scale of this model was accepted by marginaleffects: ",
       paste(cand, collapse = " or "), " were tried. Use type = \"eta\", ",
       "which does not go through the prediction machinery, or pass a `type` ",
       "this version accepts.")
}

## Can the contrast be formed from the coefficients at all? Only if the fit
## carries the cell design's columns - a model fitted from this declaration.
## One written by hand in crossed or chain form has different coefficients, and
## must go through the prediction function.
linear_map_available <- function(model, spec, data) {
  tryCatch({
    g <- cell_grid(spec, data)
    rhs <- paste(deparse(cell_formula(spec)[[3]]), collapse = " ")
    X <- stats::model.matrix(stats::as.formula(paste("~", rhs)), g)
    need <- setdiff(colnames(X), "(Intercept)")
    all(need %in% names(coef_vector(model)))
  }, error = function(e) FALSE)
}

## Does the engine accept this `type` for this model? Probed on two rows,
## because the accepted set differs between model classes and between releases.
engine_accepts <- function(model, spec, type, data = spec$data) {
  tryCatch({
    marginaleffects::predictions(model, newdata = utils::head(data, 2),
                                 type = type)
    TRUE
  }, error = function(e) FALSE)
}

## The link function of a fit, where it has one. A model with no family - lm,
## lmer - is fitted on the scale of the outcome, which is the identity link.
model_link <- function(model, spec = NULL) {
  ## An ordinal fit has a link of its own and no family() method, so asking
  ## family() would fail and the default would be wrong in the one place it
  ## matters most.
  if (!is.null(spec) && has_thresholds(spec))
    return(tryCatch(model$link, error = function(e) "probit") %||% "probit")
  if (inherits(model, "clm") || inherits(model, "clmm"))
    return(tryCatch(model$link, error = function(e) "probit") %||% "probit")
  f <- tryCatch(stats::family(model), error = function(e) NULL)
  if (is.null(f)) return("identity")
  lk <- tryCatch(f$link, error = function(e) NULL)
  if (is.null(lk) || !nzchar(lk)) "identity" else lk
}

## --- code assembly ---------------------------------------------------------

estimand_code <- function(spec, target, policy, at, contrast, dots_txt,
                          model_name, spec_name, data_name, bounds,
                          scale = "response", deg = NULL,
                          route = "g_computation", weights_txt = NULL,
                          re_note = NULL, user_hyp = FALSE, draws_txt = "",
                          shortcut = FALSE, eta_re = "") {
  cn <- spec$cell_name
  pol_txt <- if (is.character(policy) && length(policy) == 1L)
    sprintf('"%s"', policy)
  else paste0("c(", paste(sprintf('"%s" = %s', names(policy), policy),
                          collapse = ", "), ")")
  at_txt <- if (is.null(at)) "" else
    paste0(", at = c(", paste(sprintf('%s = "%s"', names(at), at),
                              collapse = ", "), ")")
  hdr <- c(sprintf("## nestimand %s -- estimand of `%s`", nestimand_build,
                   paste(target, collapse = ":")),
           sprintf("library(marginaleffects)"))
  cells_txt <- if (is.null(deg)) sprintf("%s$cells", spec_name) else
    sprintf('subset(%s$cells, paste(%s, sep = ".") %%in%% c(%s))', spec_name,
            paste(deg$vars, collapse = ", "),
            paste(sprintf('"%s"', deg$keep), collapse = ", "))
  restrict_note <- if (is.null(deg)) NULL else if (length(deg$by)) c(
    sprintf("## within `%s` = %s: the contrast is formed over that group's",
            paste(deg$by, collapse = "`, `"), deg$group),
    "## conditions alone, and the policy is weighted over them.",
    sprintf("cells <- %s", cells_txt)) else c(
    sprintf("## `%s` does not vary in every stratum. The contrast is pooled over", target),
    sprintf("## the strata in which it does - %s - since elsewhere a comparison",
            paste(deg$keep, collapse = ", ")),
    "## would leave its own levels and compare strata instead.",
    sprintf("cells <- %s", cells_txt))

  if (contrast == "within") {
    fam <- NULL
    for (f in spec$cat_families) if (target %in% f) fam <- f
    pos <- match(target, fam)
    if (pos == 1)
      stop("contrast = \"within\" emits per-stratum contrasts of a *nested* ",
           "variable; `", target, "` is not nested within anything, so every ",
           "contrast of it crosses the structural boundary and requires a ",
           "policy. Use policy = instead.")
    ## the strata of a nested variable are its ancestors, not everything
    ## declared before it: a sibling is not one of its parents
    parents <- nest_ancestors(spec, target)
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
  if (identical(scale, "latent") && identical(contrast, "interaction")) {
    tv <- paste(sprintf('"%s"', target), collapse = ", ")
    return(c(hdr[1], restrict_note,
      "## interaction on the latent scale: a difference of differences among",
      "## cells, computed as c'b. Exact, and free of the per-category expansion",
      "## that ordinal fits provoke.",
      sprintf('est  <- latent_estimand(%s, c(%s), contrast = "interaction", spec = %s, data = %s%s%s)',
              model_name, tv, spec_name, data_name, dots_txt, draws_txt),
      "est"))
  }
  if (identical(scale, "latent")) {
    body <- c(hdr[1], restrict_note,
      if (!is.null(weights_txt))
        c("## unit weights: the estimand is standardized to the population those",
          "## weights describe. They multiply the row weights, so the contrast is",
          "## a weighted average of design rows as before."),
      if (isTRUE(shortcut))
        c("## the link is the identity, so the response scale and the linear",
          "## predictor are the same quantity: the contrast is taken from the",
          "## coefficients, which averages the design matrix instead of one",
          "## prediction per row. avg_predictions() gives the same answer."),
      "## latent-scale estimand: on a linear scale the policy contrast is c'b,",
      "## with c a weighted difference of design-matrix rows. Exact, one matrix",
      "## product, and free of the per-category expansion of ordinal fits.",
      sprintf('pol  <- nest_policy(%s, "%s", %s%s%s)', spec_name, target, pol_txt,
              at_txt, if (is.null(deg)) "" else ", cells = cells"),
      sprintf('est  <- latent_estimand(%s, "%s", pol, contrast = "%s", spec = %s, data = %s%s%s%s%s%s%s)',
              model_name, target, contrast, spec_name, data_name,
              if (identical(route, "g_computation")) "" else
                sprintf(', route = "%s"', route),
              if (is.null(deg)) "" else ", cells = cells", eta_re,
              if (is.null(weights_txt)) "" else sprintf(", weights = %s", weights_txt),
              dots_txt, draws_txt))
    if (isTRUE(bounds))
      body <- c(body,
        "## partial-identification bounds over all admissible policies",
        sprintf('est  <- add_bounds(est, %s, %s, "%s", "%s", scale = "latent"%s%s)',
                model_name, spec_name, target, contrast,
                if (is.null(deg)) "" else ", cells = cells", dots_txt))
    return(c(body, "est"))
  }
  if (identical(contrast, "interaction")) {
    tv <- paste(sprintf('"%s"', target), collapse = ", ")
    return(c(hdr, restrict_note,
      "## an interaction contrast: a difference of differences. Every cell it",
      "## uses exists, so no boundary is crossed and no policy applies.",
      sprintf('pol  <- nest_policy(%s, "%s", "equal"%s)', spec_name,
              target[length(target)],
              if (is.null(deg)) "" else ", cells = cells"),
      if (identical(route, "cells"))
        c("## one row per realized cell, covariates at their means",
          sprintf('grid <- cell_grid(%s)', spec_name),
          if (!is.null(deg))
            sprintf('grid <- grid[as.character(grid$%s) %%in%% as.character(cells$%s), ]',
                    cn, cn),
          sprintf('grid$.w <- policy_weights(%s, grid, pol)', spec_name))
      else
        sprintf('grid <- counterfactual_grid(%s, %s, pol%s)', spec_name, data_name,
                if (is.null(deg)) "" else ", cells = cells"),
      "## the cell means first, to learn the order the engine returns them in",
      sprintf('g0   <- avg_predictions(%s, newdata = grid, by = c(%s), wts = grid$.w%s)',
              model_name, tv, dots_txt),
      sprintf('H    <- interaction_matrix(g0, c(%s))', tv),
      sprintf('est  <- avg_predictions(%s, newdata = grid, by = c(%s), wts = grid$.w,',
              model_name, tv),
      sprintf('                        hypothesis = H%s)', dots_txt),
      "est"))
  }
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
    if (user_hyp)
      c(sprintf('est  <- avg_predictions(%s, newdata = grid, by = "%s", wts = grid$.w%s)',
                model_name, target, dots_txt))
    else
      c(sprintf('est  <- avg_predictions(%s, newdata = grid, by = "%s", wts = grid$.w,',
                model_name, target),
        sprintf('                        hypothesis = %s%s)',
                mfx_hypothesis_txt(contrast), dots_txt)),
    "## contrast direction is fixed by nestimand, not inherited from the engine:",
    "## later declared level minus earlier, so that it reads as a departure from",
    "## the reference condition. Engine versions differ on this, and an inherited",
    "## convention would flip reported signs.",
    sprintf('est  <- mfx_canonical(est, levels(factor(%s$%s)))',
            data_name, target))
  if (isTRUE(bounds))
    body <- c(body,
      "## partial-identification bounds: every mixture estimand is a convex",
      "## combination of the single-version contrasts, so their range is the",
      "## region over all admissible policies (Manski 1990)",
      sprintf('est  <- add_bounds(est, %s, %s, "%s", "%s"%s%s%s)',
              model_name, spec_name, target, contrast,
              if (is.null(deg)) "" else ", cells = cells",
              if (identical(route, "g_computation")) "" else
                sprintf(', route = "%s"', route), dots_txt))
  c(body, "est")
}

## --- the bounds companion --------------------------------------------------

add_bounds <- function(est, model, spec, target, contrast = "pairwise",
                       scale = c("response", "latent"), cells = spec$cells,
                       route = c("g_computation", "cells"), data = spec$data, ...) {
  ## The latent route computes from the coefficients, so arguments meant for
  ## the prediction function have nowhere to go. Said plainly, rather than
  ## surfacing as "unused argument" from a function the user did not call.
  engine_only <- c("type", "vcov", "p_adjust", "transform", "df", "byfun",
                   "equivalence", "numderiv", "re.form", "re_formula",
                   "variables", "newdata")
  scale <- match.arg(scale)
  route <- match.arg(route)
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
      return(latent_estimand(model, target, p, contrast = contrast, spec = spec,
                             data = data, cells = cells)$estimate)
    ## the bounds are the same estimand under other policies, so they are
    ## computed the same way: on whichever grid the estimand itself used
    g <- if (identical(route, "cells")) {
      d <- cell_grid(spec, data)
      d <- d[as.character(d[[spec$cell_name]]) %in%
             as.character(cells[[spec$cell_name]]), , drop = FALSE]
      d$.w <- policy_weights(spec, d, p); d
    } else counterfactual_grid(spec, data, p, cells = cells)
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
    brm = if (has_thresholds(spec)) "clm" else if (is.null(spec$family)) "lm" else "glm",
    NULL)
  converged <- function(fit) {
    if (inherits(fit, "clm") || inherits(fit, "clmm"))
      return(isTRUE(fit$convergence$code == 0))
    if (inherits(fit, "glm")) return(isTRUE(fit$converged))
    TRUE
  }
  fit_shadow <- function(spx, dta) {
    f <- cell_formula(spx)
    switch(shadow,
      lm  = stats::lm(f, data = dta),
      glm = stats::glm(f, data = dta,
                       family = eval(parse(text = spx$family %||% "gaussian"))),
      clm = ordinal::clm(f, data = dta))
  }
  ## The check permutes factor levels and asks whether the estimand moves, so
  ## it works from the declaration's own data: a subsample taken for the
  ## estimand may not contain every cell, and rebuilding a spec from it would
  ## compare a smaller design against the full one.
  d2 <- spec$data
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
    ## A shadow that did not converge cannot settle anything: two runs may
    ## differ because the optimizer stopped in different places, which says
    ## nothing about the parameterization. Saying so is better than telling
    ## someone not to report a result that is probably fine.
    if (!converged(m1) || !converged(m2))
      return(list(status = "inconclusive",
                  note = paste("the fixed-effects shadow model did not",
                               "converge, so a difference between the two",
                               "orderings would say more about the optimizer",
                               "than about the parameterization")))
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
           note = sprintf(paste("max |change| = %s over %d contrasts, checked on",
                                "the %s scale over the %s grid%s"),
                          format(max(abs(sort(abs(e1)) - sort(abs(e2))))),
                          length(e1), chk_scale, chk_route,
                          if (!is.null(shadow))
                            paste0(" of a ", shadow, " shadow model") else ""))
    else list(status = "passed", note = note)
  }, error = function(e) list(status = "error", note = conditionMessage(e)))
  if (identical(out$status, "inconclusive"))
    message("reorder self-check inconclusive: ", out$note, ".")
  if (identical(out$status, "failed"))
    warning("reorder self-check FAILED (", out$note, "): the estimate moved when ",
            "nested-factor levels were permuted, so this estimand depends on ",
            "unrealized-cell predictions and is not identified by the design. ",
            "Do not report it.", call. = FALSE)
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a

## Set while the parts of `a * b` are computed: each would otherwise repeat
## what the whole says once, from its own point of view.
nestimand_env <- new.env(parent = emptyenv())
in_star <- function() isTRUE(get0("star", envir = nestimand_env, ifnotfound = FALSE))

nesting_spec_quiet <- function(spec, data) {
  suppressMessages(nesting_spec(data, spec$formula_in,
    nests = spec_nests(spec), fit = spec$fit, family = spec$family,
    random = spec$random_original, cell_name = spec$cell_name))
}
## The declarations, rebuilt from the parent map rather than from adjacency in
## the family vector, which is not ancestry once a parent holds two children.
spec_nests <- function(spec)
  c(if (length(spec$parent))
      unname(sprintf("%s %%in%% %s", names(spec$parent), spec$parent)),
    ## a variable nested in nothing is declared by name, and has to be restated
    ## or the rebuilt spec would lose it from the cell factor
    unname(spec$crossed))

estimand_values <- function(model, spec, target, policy, at, contrast, data,
                            scale = "response", route = "g_computation") {
  if (identical(contrast, "interaction") && identical(scale, "latent"))
    return(latent_estimand(model, target, contrast = "interaction",
                           data = data, spec = spec)$estimate)
  if (identical(contrast, "interaction")) {
    deg <- degenerate_strata_multi(spec, target)
    cells <- if (is.null(deg) || !length(deg$drop)) spec$cells else
      spec$cells[deg_key(spec$cells, deg$vars) %in% deg$keep, , drop = FALSE]
    pol <- nest_policy(spec, target[length(target)], "equal", NULL, spec$data,
                       cells = cells)
    g <- if (identical(route, "cells")) {
      d <- cell_grid(spec, data)
      d <- d[as.character(d[[spec$cell_name]]) %in%
             as.character(cells[[spec$cell_name]]), , drop = FALSE]
      d$.w <- policy_weights(spec, d, pol); d
    } else counterfactual_grid(spec, data, pol, cells = cells)
    g0 <- marginaleffects::avg_predictions(model, newdata = g, by = target,
                                           wts = g$.w)
    H <- interaction_matrix(g0, target)
    return(as.data.frame(marginaleffects::avg_predictions(model, newdata = g,
             by = target, wts = g$.w, hypothesis = H))$estimate)
  }
  if (contrast == "within") {
    fam <- NULL; for (f in spec$cat_families) if (target %in% f) fam <- f
    parents <- nest_ancestors(spec, target)
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
    spec$cells[deg_key(spec$cells, deg$vars) %in% deg$keep, , drop = FALSE]
  ## The policy describes the design - how often each version is realized -
  ## and is a property of the whole data, not of whichever rows were sampled to
  ## average over. Computing it from a subsample would make the weights noisy,
  ## or undefined where a version is missing.
  pol <- nest_policy(spec, target, policy, at, spec$data, cells = cells)
  if (identical(scale, "latent"))
    return(latent_estimand(model, target, pol, contrast = contrast,
                           data = data, spec = spec, cells = cells)$estimate)
  g <- switch(route,
    g_computation = counterfactual_grid(spec, data, pol, cells = cells),
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

## Numbers right-aligned under right-aligned headers, as the neighbouring
## packages print them; the label columns left-aligned, header included, by
## padding both to a common width.
print_aligned <- function(d, ...) {
  labels <- c("term", "stratum", "meaning", "group", "parameter", "space")
  for (k in names(d)) {
    v <- if (is.numeric(d[[k]])) format(d[[k]], trim = TRUE, scientific = FALSE)
         else as.character(d[[k]])
    v[is.na(d[[k]])] <- "NA"
    w <- max(nchar(c(v, k)), na.rm = TRUE)
    ## label columns left-justified, header included; everything else right,
    ## so that a header sits over its own digits whatever the column's type
    just <- if (k %in% labels) -w else w
    d[[k]] <- formatC(v, width = just)
    names(d)[names(d) == k] <- formatC(k, width = just)
  }
  print(d, row.names = FALSE, right = TRUE, ...)
}

print.nestimand_estimand <- function(x, digits = 4, ...) {
  meta <- attr(x, "nestimand")
  d <- as.data.frame(x)
  ## `group` carries the outcome category on a grouped fit: without it the
  ## rows cannot be told apart
  keep <- intersect(c(meta$by, "stratum", "group", "term", "estimate", "std.error",
                      "conf.low", "conf.high", "statistic", "p.value", "pd"),
                    names(d))
  d <- d[, keep, drop = FALSE]
  for (k in intersect(c("estimate", "std.error", "conf.low", "conf.high",
                        "statistic", "pd"), names(d)))
    d[[k]] <- round(d[[k]], digits)
  if ("p.value" %in% names(d))
    d$p.value <- format.pval(d$p.value, digits = max(2, digits - 1), eps = 10^-digits)
  print_aligned(d, ...)
  pol <- if (identical(meta$contrast, "interaction")) "not applicable"
         else if (is.character(meta$policy)) meta$policy else "supplied"
  cat("\nPolicy: ", pol, "   route: ", meta$route, "   contrast: ", meta$contrast,
      "   type: ", meta$type, sep = "")
  if (!is.null(meta$self_check))
    cat("   reorder check: ", meta$self_check$status, sep = "")
  cat("\n")
  if (!is.null(meta$bounds)) {
    cat("Bounds over all admissible policies:\n")
    b <- meta$bounds
    bt <- data.frame(term = b$term, estimate = round(b$estimate, digits),
                     low = round(b$policy_low, digits),
                     high = round(b$policy_high, digits))
    ## a grouped estimand has one set of bounds per group, and they are not
    ## readable without the group beside them
    if (length(intersect(meta$by, names(b))))
      bt <- cbind(b[, intersect(meta$by, names(b)), drop = FALSE], bt)
    print_aligned(bt)
  }
  cat("Code: show_code() prints the ", length(meta$code),
      " lines that produced this.\n", sep = "")
  invisible(x)
}


print.nestimand_estimands <- function(x, digits = 4, ...) {
  if (!length(x)) { cat("no estimands\n"); return(invisible(x)) }
  nm <- names(x)
  if (is.null(nm)) nm <- paste0("[[", seq_along(x), "]]")
  for (i in seq_along(x)) {
    cat("== ", nm[i], " ==\n", sep = "")
    print(x[[i]], digits = digits, ...)
    cat("\n")
  }
  invisible(x)
}

show_code.nestimand_estimands <- function(x, ...) {
  nm <- names(x)
  if (is.null(nm)) nm <- paste0("[[", seq_along(x), "]]")
  code <- unlist(lapply(seq_along(x), function(i)
    c(sprintf("## --- %s ---", nm[i]), attr(x[[i]], "nestimand")$code, "")))
  cat(paste(code, collapse = "\n"), "\n", sep = "")
  invisible(structure(paste(code, collapse = "\n"), class = "nestimand_code"))
}


## Several targets in one call. With `dry_run` each part is code rather than a
## result, and the collection is then one script rather than a list of tables:
## a list of code objects would have no useful print or show_code.
collect_estimands <- function(out) {
  if (all(vapply(out, inherits, TRUE, "nestimand_code"))) {
    nm <- names(out)
    code <- unlist(lapply(seq_along(out), function(i)
      c(sprintf("## --- %s ---", nm[i]),
        attr(out[[i]], "nestimand_code"), "")))
    return(structure(paste(code, collapse = "\n"), class = "nestimand_code",
                     nestimand_code = code))
  }
  structure(out, class = "nestimand_estimands")
}


## Replace a frequentist summary with a posterior one where the draws exist.
add_posterior_summary <- function(out, model) {
  if (!inherits(model, "brmsfit")) return(out)
  d <- tryCatch(as.data.frame(out), error = function(e) NULL)
  if (is.null(d) || "pd" %in% names(d)) return(out)
  drw <- tryCatch(marginaleffects::posterior_draws(out), error = function(e) NULL)
  if (is.null(drw) || !all(c("drawid", "draw") %in% names(drw))) return(out)
  key <- if ("term" %in% names(drw)) "term" else NULL
  grp <- intersect(c("group", "by"), names(drw))
  idx <- do.call(paste, c(unname(drw[c(key, grp)]), sep = "\r"))
  pd <- tapply(drw$draw, idx, function(z) max(mean(z > 0), mean(z < 0)))
  own <- do.call(paste, c(unname(d[intersect(c(key, grp), names(d))]), sep = "\r"))
  d$pd <- as.numeric(pd[match(own, names(pd))])
  d$statistic <- NULL
  d$p.value <- NULL
  for (a in names(attributes(out)))
    if (!a %in% c("names", "row.names", "class")) attr(d, a) <- attr(out, a)
  class(d) <- unique(c(setdiff(class(out), "data.frame"), class(d)))
  d
}

## --- what the engine will accept -------------------------------------------
## The type vocabulary belongs to marginaleffects and differs by model class -
## and between its releases - so it is read from the installed version rather
## than written down here. `eta` is nestimand's own and is always available.
engine_types <- function(model = NULL) {
  d <- tryCatch(marginaleffects:::type_dictionary, error = function(e) NULL)
  if (is.null(d))
    stop("the type dictionary of the installed marginaleffects cannot be read; ",
         "see ?marginaleffects::avg_predictions for what it accepts.")
  cls <- if (is.null(model)) c("lm", "glm", "lmerMod", "glmerMod", "clm", "clmm",
                               "brmsfit")
         else class(model)
  out <- do.call(rbind, lapply(cls, function(k) {
    ty <- d$type[d$class == k]
    if (!length(ty)) return(NULL)
    data.frame(class = k, types = paste(ty, collapse = ", "), row.names = NULL)
  }))
  if (is.null(out) || !nrow(out))
    out <- data.frame(class = paste(cls, collapse = "/"),
                      types = "none: this class has no entry, so no quantity computed by prediction is available",
                      row.names = NULL)
  out <- rbind(data.frame(class = "any", types = "eta (nestimand's own: the linear predictor, from the coefficients)"),
               out)
  structure(out, class = c("nestimand_types", "data.frame"))
}

print.nestimand_types <- function(x, ...) {
  cat("types accepted by marginaleffects ",
      as.character(utils::packageVersion("marginaleffects")), "\n", sep = "")
  print_aligned(as.data.frame(x))
  invisible(x)
}
