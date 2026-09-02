## nestimand: declaration surface -----------------------------------------
## The declaration surface is carried over from the chain-based prototype
## unchanged in its user-facing form (bare or quoted arguments, validations,
## teaching errors, saturation and outcome-type guards). What changed beneath
## it is the parameterization: the declared structure is now compiled into a
## realized-cell factor rather than into a chain of interaction terms.

nestimand_build <- "2026-08-17.1"

## The variables named on the left of `%in%`. Read from the parsed expression
## rather than by splitting the text, so that `a + b`, `c(a, b)`, `(a + b)` and
## any nesting of those mean the same thing - the parentheses a reader is
## likely to add for the precedence of `%in%` among them. Anything that is not
## a list of names is refused: an interaction is not nested in anything.
nested_names <- function(txt, parent) {
  refuse <- function(e)
    stop("the left of `%in%` names the variables nested in `", parent,
         "`, separated by `+`, gathered with `c()`, or bracketed; `",
         paste(deparse(e), collapse = " "), "` is not a variable name. An ",
         "interaction is not nested in anything - declare each variable of ",
         "it separately.")
  walk <- function(e) {
    if (is.name(e)) return(as.character(e))
    ## a quoted declaration parses its names as string constants
    if (is.character(e) && length(e) == 1L &&
        grepl("^[A-Za-z.][A-Za-z0-9._]*$", e)) return(e)
    if (is.call(e) && as.character(e[[1]]) %in% c("(", "c", "+"))
      return(unlist(lapply(as.list(e)[-1], walk), use.names = FALSE))
    refuse(e)
  }
  e <- tryCatch(str2lang(trimws(txt)),
                error = function(err) stop("`", trimws(txt), "` on the left of ",
                  "`%in%` is not an R expression."))
  unique(walk(e))
}

## --- inferring the structure from the data ---------------------------------
## A partially nested design writes itself into the data: where a variable is
## structurally undefined it has no value, and the rows where that happens are
## exactly the rows some other variable's levels pick out. So the declaration
## can be read off rather than asked for.
##
## The rule: `v` is nested within the smallest set of design variables such that,
## within every combination of their levels, `v` is either always absent or never
## absent. That is what "structurally undefined here" means, and it is decidable.
## Where no such set exists the absences are not structural - they are missing
## data - and saying so is the point: telling structure from missingness is the
## one judgement here that cannot be got wrong quietly.
##
## Two cases the first version of this got wrong, both worth keeping in view:
## a parent may itself be absent (that is what a chain is), so its own absence
## is one of its strata rather than a disqualification; and two variables absent
## on the same rows each explain the other, which is what siblings look like, so
## a candidate that is never absent is preferred as the coarser division.
infer_nests <- function(data, vars, max_parents = 2L) {
  absent <- lapply(vars, function(v) is.na(data[[v]]))
  names(absent) <- vars
  nests <- character(0); notes <- character(0); ambiguous <- character(0)
  for (v in vars) {
    if (!any(absent[[v]])) next                     # never absent: crossed
    if (all(absent[[v]])) {
      notes <- c(notes, sprintf("`%s` is missing everywhere", v)); next
    }
    cand <- setdiff(vars, v); found <- NULL; ties <- list()
    for (k in seq_len(min(max_parents, length(cand)))) {
      sets <- utils::combn(cand, k, simplify = FALSE)
      ok <- Filter(function(P) {
        key <- do.call(paste, c(lapply(P, function(p)
          ifelse(absent[[p]], "\r<absent>", as.character(data[[p]]))), sep = "\r"))
        all(tapply(absent[[v]], key, function(z) length(unique(z)) == 1L))
      }, sets)
      if (length(ok)) {
        solid <- Filter(function(P)
          !any(vapply(P, function(p) any(absent[[p]]), TRUE)), ok)
        if (length(solid)) ok <- solid
        found <- ok[[1]]; ties <- ok; break
      }
    }
    if (is.null(found)) {
      notes <- c(notes, sprintf(paste(
        "`%s` is missing on rows that no combination of the other design",
        "variables picks out, so the gaps are not structural: they are missing",
        "data, which R would delete casewise and silently"), v))
    } else {
      if (length(ties) > 1)
        ambiguous <- c(ambiguous, sprintf("`%s` (%s)", v,
          paste(vapply(ties, paste, "", collapse = ":"), collapse = " or ")))
      nests <- c(nests, sprintf("%s %%in%% %s", v, paste(found, collapse = ":")))
    }
  }
  list(nests = nests, notes = notes, ambiguous = ambiguous)
}

nesting_spec <- function(data, formula, nests = NULL,
                         fit = c("lm", "glm", "lmer", "glmer", "clm", "clmm", "brm"),
                         family = NULL, random = NULL) {
  ## The engines are named for their fitting functions, so brms's is `brm`.
  ## Its package name is an easy slip and is accepted as the same thing.
  if (identical(fit, "brms")) fit <- "brm"
  fit <- match.arg(fit)
  fq <- substitute(family)
  family <- if (is.null(fq)) NULL
            else if (tryCatch(is.character(family), error = function(e) FALSE)) family
            else paste(deparse(fq), collapse = "")  # bare call, e.g. cumulative("probit")
  has_bars <- grepl("|", paste(deparse(formula), collapse = ""), fixed = TRUE)
  if (fit %in% c("lmer", "glmer", "clmm") && is.null(random) && !has_bars)
    stop("fit = '", fit, "' needs a `random` term, e.g. random = \"(1 | participant)\"")
  if (fit == "glmer" && is.null(family))
    stop("fit = 'glmer' needs a `family`, e.g. family = \"binomial\"")
  if (inherits(random, "formula")) random <- paste(deparse(random[[2]]), collapse = " ")
  data_name <- { dq <- substitute(data); if (is.name(dq)) deparse(dq) else "dat" }
  ## The declaration may arrive as text, as a character vector, or unquoted -
  ## `inversion %in% chord_type`, `c(inversion, X1) %in% chord_type`, the
  ## `inversion + X1 %in% chord_type` sugar, whose parse tree does not match its
  ## reading, or a bare name for a variable nested in nothing. Everything
  ## unquoted is carried as its own text and read below; a value that is already
  ## a vector of declarations is taken as it stands.
  ## Omitted, the declaration is read off the data: the variables the formula
  ## names, and which of them is undefined where. What was inferred is said, in
  ## the syntax the user would have written, so it can be checked and overridden.
  ## Inference sees the sample, not the design: gaps that are regular but not
  ## structural - a block of trials lost for one condition - look exactly like a
  ## nesting, which is why this is announced rather than done quietly.
  inferred <- NULL
  ## tested on the expression, not the value: `nests` is unquoted more often
  ## than not, and `dose %in% arm` names columns rather than objects, so forcing
  ## the promise merely to ask whether it is NULL would evaluate it and fail
  ne0 <- substitute(nests)
  if (missing(nests) || is.null(ne0)) {
    cvars <- unique(unlist(strsplit(
      attr(stats::terms(formula), "term.labels")[
        !grepl("|", attr(stats::terms(formula), "term.labels"), fixed = TRUE)],
      ":", fixed = TRUE)))
    cvars <- cvars[cvars %in% names(data)]
    cvars <- cvars[vapply(cvars, function(v)
      is.factor(data[[v]]) || is.character(data[[v]]) || is.logical(data[[v]]), TRUE)]
    inf <- infer_nests(data, cvars)
    if (length(inf$notes)) stop(paste(inf$notes, collapse = " "),
      ". Supply `nests` to say what the structure is, and those rows will be ",
      "checked against it; anything left over is missing data, and should be ",
      "resolved rather than coded as a condition.")
    if (length(inf$ambiguous))
      stop("the structure cannot be read off the data: ",
           paste(inf$ambiguous, collapse = "; "), " each explain the same gaps ",
           "equally well, so which is the parent is not something the data ",
           "says. Supply `nests`.")
    nests <- inf$nests
    inferred <- nests
    message(if (!length(nests))
      "no structurally undefined values found, so nothing is nested: the design is fully crossed."
      else paste0("nesting structure read from the data: nests = ",
                  if (length(nests) == 1) paste0('"', nests, '"') else
                    paste0("c(", paste(sprintf('"%s"', nests), collapse = ", "), ")"),
                  ". Every variable is undefined exactly where its parent's levels say ",
                  "it should be, which is what makes this readable off the data rather ",
                  "than guessed. Supply `nests` yourself to declare something else."))
    ## the undefined rows are now known to be structural, so they can be given
    ## the sentinel the rest of the package is built on
    for (z in nests) {
      v <- trimws(strsplit(z, "%in%")[[1]][1])
      ## the inference has just established that these absences are structural,
      ## so the check can only pass; it runs all the same, since one path in and
      ## one path out is worth more than the call it saves
      if (anyNA(data[[v]]))
        data <- set_sentinel(data, v, where = is.na(data[[v]]))
    }
  }
  ne <- substitute(nests)
  ne_txt <- paste(deparse(ne), collapse = " ")
  ne_val <- tryCatch(eval(ne, parent.frame()), error = function(e) NULL)
  if (!is.null(inferred)) { ne <- NULL; ne_val <- inferred }
  nests <- if (is.character(ne_val)) ne_val
  else if (is.null(ne)) ne_val
  else if (is.name(ne)) ne_txt
  else if (is.call(ne) && identical(ne[[1]], as.name("c")))
    vapply(as.list(ne)[-1], function(x)
      if (is.character(x)) x else paste(deparse(x), collapse = " "), "")
  else if (is.call(ne) && grepl("%in%", ne_txt, fixed = TRUE))
    ne_txt
  else eval(ne, parent.frame())

  ## One declaration may name several variables nested in the same parent,
  ## written `c(inversion, X1) %in% chord_type` or `inversion + X1 %in%
  ## chord_type`. The declaration is carried as text and split here rather than
  ## evaluated, which is what makes the second form work: `%in%` binds tighter
  ## than `+`, so R parses it as `inversion + (X1 %in% chord_type)`, and only
  ## the text says what was meant.
  parse1 <- function(s) {
    p <- trimws(strsplit(s, "%in%")[[1]])
    if (length(p) > 2)
      stop("`", trimws(s), "` chains `%in%` more than once. The operator is ",
           "left-associative, so this reads as `(", p[1], " %in% ", p[2],
           ") %in% ", p[3], "`, which names no design. Declare one level per ",
           "entry: c(\"", p[2], " %in% ", p[3], "\", \"", p[1], " %in% ",
           p[2], "\").")
    if (length(p) != 2)
      stop("a nesting declaration reads `child %in% parent`; `", trimws(s),
           "` has no `%in%`.")
    kids <- nested_names(p[1], p[2])
    stats::setNames(rep(p[2], length(kids)), kids)    # child -> parent
  }
  ## An entry with no `%in%` names a categorical variable of the design that is
  ## nested in nothing - crossed with the rest of the structure. It belongs in
  ## the cell factor all the same: the cells are the realized categorical
  ## design, and a variable left out of them is a covariate, which is a
  ## different thing to be and cannot be the target of an estimand.
  is_nest <- grepl("%in%", nests, fixed = TRUE)
  crossed <- unique(unlist(lapply(nests[!is_nest], function(z)
    nested_names(z, "the design"))))
  parent <- if (any(is_nest))
    do.call(c, lapply(nests[is_nest], parse1)) else character(0)
  if (anyDuplicated(names(parent)))
    stop("`", paste(unique(names(parent)[duplicated(names(parent))]),
                    collapse = "`, `"), "` is declared inside more than one ",
         "parent. A variable holds one position in the structure: nest it in ",
         "the deeper parent, which carries the other with it.")

  ## The fixed term labels are needed before the families are formed, since a
  ## categorical variable the formula crosses with the structure joins them.
  labs <- attr(stats::terms(formula), "term.labels")
  labs <- labs[!grepl("|", labs, fixed = TRUE)]

  ## --- validation --------------------------------------------------------
  nest_vars <- unique(c(names(parent), parent, crossed))
  for (v in nest_vars)
    if (!v %in% names(data)) stop("variable not in data: ", v)
  ## A declared structure says where each nested variable is undefined: in the
  ## strata of its parents in which it takes no value at all. That is checkable
  ## against the data, so the sentinel is applied here rather than asked for -
  ## and an NA the declaration does not account for is still refused, which is
  ## the whole of what the user used to do by hand.
  for (v in nest_vars) {
    if (!anyNA(data[[v]])) next
    anc <- v; repeat { p <- unname(parent[anc[1]]); if (is.na(p) || is.null(p)) break
                       anc <- c(p, anc) }
    up <- setdiff(anc, v)
    if (!length(up))
      stop("`", v, "` contains NA but is nested in nothing, so no structure ",
           "says where it should be undefined. Those NAs are missing data - R ",
           "would delete them casewise and silently. Declare `", v, "` inside ",
           "the variable that determines where it applies, or resolve them.")
    key <- do.call(paste, c(lapply(up, function(z)
      ifelse(is.na(data[[z]]), "\r<absent>", as.character(data[[z]]))), sep = "\r"))
    absent <- is.na(data[[v]])
    ok <- tapply(absent, key, function(z) length(unique(z)) == 1L)
    if (!all(ok))
      stop("`", v, "` is NA in some rows of a stratum of `",
           paste(up, collapse = ":"), "` and not others (", sum(!ok),
           " stratum/strata). A variable is undefined for a whole stratum or ",
           "for none of it, so these are missing data rather than structure. ",
           "Resolve them, or declare the structure that accounts for them.")
    data <- set_sentinel(data, v, where = absent)
  }
  for (ch in names(parent))
    if (is.numeric(data[[parent[[ch]]]]))
      stop("a continuous variable cannot nest another variable: ", parent[[ch]])
  ## The realized conditions are carried internally as a factor. It is an
  ## implementation detail, so a column of that name already in the data is a
  ## clash to be stepped around rather than a question to put to the user.
  cell_name <- "cell"
  while (cell_name %in% names(data)) cell_name <- paste0(".", cell_name)

  ## --- families (roots and their chains) --------------------------------
  ## A categorical variable the formula crosses with the declared structure is
  ## part of the categorical design, and does not have to be declared to be
  ## treated as one: `chord_type * inversion * top_note` says that the
  ## conditions are the combinations of all three. It is folded into the cell
  ## factor, which changes no fitted model - `cell + cell:top` and a cell factor
  ## over the enlarged design span the same space - but makes it a legal target
  ## of an estimand, and drops the columns of any combination the design does
  ## not realize, which the interaction form would carry as aliased
  ## coefficients. A variable entering additively is left alone: its effect is
  ## declared common to every condition, and folding it in would silently
  ## saturate it. So is an ordered factor, whose contrasts say it is meant as a
  ## quantity; name it in `nests` to fold that in deliberately.
  declared_vars <- unique(c(names(parent), parent, crossed))
  is_design_cat <- function(v)
    v %in% names(data) && !is.numeric(data[[v]]) && !is.ordered(data[[v]])
  auto_crossed <- setdiff(unique(unlist(lapply(strsplit(labs, ":"), function(vs)
    if (any(vs %in% declared_vars)) vs))), declared_vars)
  auto_crossed <- auto_crossed[vapply(auto_crossed, is_design_cat, TRUE)]
  if (length(auto_crossed))
    message("`", paste(auto_crossed, collapse = "`, `"), "` ",
            if (length(auto_crossed) > 1) "are factors" else "is a factor",
            " crossed with the declared structure, so ",
            if (length(auto_crossed) > 1) "they are " else "it is ",
            "part of the categorical design and ",
            if (length(auto_crossed) > 1) "join " else "joins ",
            "the cell factor - the fitted model is unchanged, and ",
            if (length(auto_crossed) > 1) "they can " else "it can ",
            "be the target of an estimand. A variable entering additively, a ",
            "numeric one, and an ordered factor stay covariates.")
  crossed <- unique(c(crossed, auto_crossed))

  ## A crossed variable is a family of its own, one variable deep.
  for (v in crossed) {
    if (is.numeric(data[[v]]))
      stop("`", v, "` is numeric, and a numeric variable crossed with the ",
           "structure is a covariate, which is what it already is if it is ",
           "left out of the declaration: write it in the formula and it enters ",
           "as a slope within the realized cells. Declare it here only if its ",
           "values are labels, in which case make it a factor first.")
    if (v %in% names(parent))
      stop("`", v, "` is declared both on its own and inside `", parent[[v]],
           "`. A variable holds one position in the structure: nested, or ",
           "crossed with everything else.")
  }
  roots <- unique(c(setdiff(parent, names(parent)), crossed))
  ## A family is a tree, not only a chain: one parent may hold several nested
  ## variables - inversion and X1 both inside chord_type - and each of those
  ## may hold further ones. The family is returned depth-first in declaration
  ## order, so an ancestor always precedes its descendants and the vector can
  ## still be read as a sequence where that is all a caller needs. Real
  ## ancestry is `nest_ancestors()`; a position in the vector is not it.
  family_of <- function(root) {
    walk <- function(v, seen) {
      if (v %in% seen)
        stop("the declared nesting is circular at `", v, "`.")
      c(v, unlist(lapply(names(parent)[parent == v], walk, seen = c(seen, v)),
                  use.names = FALSE))
    }
    walk(root, character(0))
  }
  families <- lapply(unique(roots), family_of)

  ## --- outcome and fit guards (unchanged) --------------------------------
  outcome <- deparse(formula[[2]])
  if (fit %in% c("clm", "clmm") && !(outcome %in% names(data) && is.ordered(data[[outcome]])))
    stop("fit = '", fit, "' needs an ordered-factor outcome; `", outcome, "` is not one.")
  if (fit == "brm" && outcome %in% names(data) && is.ordered(data[[outcome]]) &&
      is.null(family))
    stop("outcome `", outcome, "` is an ordered factor, and brms defaults to ",
         "family = gaussian, which requires a numeric response. For ordinal ",
         "ratings use family = \"cumulative()\" (the brms counterpart of clm); ",
         "the cell parameterization applies unchanged.")
  if (fit == "glm" && outcome %in% names(data) && is.ordered(data[[outcome]]))
    stop("outcome `", outcome, "` is an ordered factor, and glm has no ordinal ",
         "family. Use fit = \"clm\", or fit = \"brms\" with ",
         "family = \"cumulative()\".")
  if (fit %in% c("lm", "lmer") && outcome %in% names(data) && is.factor(data[[outcome]]))
    stop("outcome `", outcome, "` is a", if (is.ordered(data[[outcome]])) "n ordered",
         " factor. `", fit, "()` on a factor response computes coefficients from ",
         "the internal integer codes but cannot form residuals, so standard ",
         "errors and tests are lost. Either use fit = \"clm\" (same cell ",
         "parameterization), or convert explicitly - `as.integer(", outcome, ")` - ",
         "if a linear analysis of the codes is intended.")

  ## --- fixed and random parts of the declared formula --------------------
  all_labs <- attr(stats::terms(formula), "term.labels")
  bar_labs <- c(all_labs[grepl("|", all_labs, fixed = TRUE)],
                if (!is.null(random))
                  grep("|", attr(stats::terms(stats::as.formula(
                    paste("~", random))), "term.labels"), value = TRUE, fixed = TRUE))
  check_double_bar(bar_labs, data)
  ## `diag(...)`, and any other covariance-structure wrapper, is already a call:
  ## wrapping it again in parentheses would change the term.
  random_original <- if (length(bar_labs))
    paste(vapply(bar_labs, function(b)
      if (grepl("^[A-Za-z.][A-Za-z0-9._]*\\(", b)) b else paste0("(", b, ")"), ""),
      collapse = " + ")

  ## --- the cell structure ------------------------------------------------
  ## One factor for the full realized categorical structure, at the deepest
  ## level of every declared family. Continuous leaves are not part of the
  ## cell factor: they enter as covariates crossed with it.
  is_cont <- function(v) is.numeric(data[[v]])
  cat_families  <- lapply(families, function(f) f[!vapply(f, is_cont, TRUE)])
  cont_nested   <- unlist(lapply(families, function(f) f[vapply(f, is_cont, TRUE)]))
  cell_vars     <- unlist(cat_families, use.names = FALSE)
  if (!length(cell_vars))
    stop("every declared nesting variable is continuous; there is no ",
         "categorical structure to form cells from.")
  key <- do.call(paste, c(unname(data[cell_vars]), sep = "."))
  ord <- do.call(order, unname(data[cell_vars]))
  cell_levels <- unique(key[ord])
  cells <- unique(data[ord, cell_vars, drop = FALSE])
  rownames(cells) <- NULL
  cells[[cell_name]] <- factor(cell_levels, levels = cell_levels)
  data[[cell_name]] <- factor(key, levels = cell_levels)

  ## per-family realized tables, for reporting
  cells_by_family <- lapply(cat_families, function(f) unique(data[, f, drop = FALSE]))

  ## `all.vars` sees through a transformation - I(x^2), log(x), poly(x, 2) -
  ## and reports the variable underneath, so a formula written with one would
  ## be silently reduced to the untransformed column. Refused rather than
  ## reduced: the fitted model would not be the one that was declared.
  fixed_terms <- attr(stats::terms(stats::as.formula(
    paste("~", paste(deparse(formula[[3]]), collapse = " ")))), "term.labels")
  fixed_terms <- fixed_terms[!grepl("\\|", fixed_terms)]
  transformed <- fixed_terms[grepl("[A-Za-z._][A-Za-z0-9._]*\\(", fixed_terms)]
  if (length(transformed))
    stop("the formula transforms a variable in place: ",
         paste(unique(sub(".*:", "", transformed)), collapse = ", "),
         ". The declaration is translated term by term, and a transformation ",
         "written inside it would be lost. Compute the transformed variable as ",
         "a column of the data - dat$x2 <- dat$x^2 - and name that column in ",
         "the formula.")
  ## The covariates are the variables of the *fixed* terms that are not declared
  ## nesting variables, read from the fixed term labels, which already exclude
  ## the bars. Reading them from `all.vars(formula[[3]])` and then subtracting
  ## everything named in a bar took out the slope variables along with the
  ## grouping factors, so a variable given a random slope lost its fixed effect
  ## - a random slope with no fixed counterpart, and nothing said.
  covariates <- setdiff(unique(unlist(strsplit(labs, ":"))), unlist(families))
  covariates <- unique(c(covariates, cont_nested))

  ## does the declared formula ask for covariate x structure interaction?
  cov_by_cell <- vapply(covariates, function(cv)
    any(vapply(strsplit(labs, ":"), function(vs)
      cv %in% vs && any(unlist(families) %in% vs), TRUE)), TRUE)
  names(cov_by_cell) <- covariates

  ## Does the declared formula ask for less than the saturated structure over
  ## the realized cells? Cells mode fits `~ 0 + cell` whatever was written, so
  ## `+` between two nesting variables restricts nothing. Saying so is better
  ## than leaving it to be inferred from a residual degrees of freedom.
  ## A numeric nested variable is not part of the cell factor: it enters as a
  ## slope within the realized cells, which is what a continuous nesting means.
  ## The difference from a factor is large enough to be worth stating, since a
  ## variable whose values are labels is easily left numeric by accident.
  if (length(cont_nested))
    message("`", paste(cont_nested, collapse = "`, `"), "` ",
            if (length(cont_nested) > 1) "are numeric" else "is numeric",
            ", so ", if (length(cont_nested) > 1) "they are " else "it is ",
            "not part of the cell factor: a continuous nested variable enters ",
            "as a slope within the realized cells - `", cell_name, ":",
            cont_nested[1], "` - which is what nesting a quantity means. If ",
            "its values are labels rather than quantities, make it a factor ",
            "before declaring it.")
  counts <- table(data[[cell_name]])
  if (all(counts == 1))
    message("every realized cell contains a single observation: the model will be ",
            "saturated (zero residual degrees of freedom), so standard errors and ",
            "tests will be unavailable. Estimates remain valid as cell-mean ",
            "contrasts; inference needs trial-level data.")

  out <- structure(list(data = data, outcome = outcome, formula_in = formula,
                 parent = parent, crossed = crossed,
                 families = families, cat_families = cat_families,
                 cont_nested = cont_nested, covariates = covariates,
                 cov_by_cell = cov_by_cell, term_labels = labs,
                 cell_name = cell_name, cell_vars = cell_vars,
                 cells = cells, cell_levels = cell_levels,
                 cells_by_family = cells_by_family,
                 data_name = data_name, fit = fit, family = family,
                 random_original = random_original,
                 build = nestimand_build),
            class = "nesting_spec")
  ## Does the formula cross the structure fully, or does it restrict it? The
  ## count is taken from the very design the fit will use - term_span() is the
  ## width of that design - so what is reported here and what is fitted cannot
  ## disagree. They used to: this message counted one thing and nest_fit() built
  ## another, and a `+` between two nesting variables was quietly upgraded to a
  ## `*` while the user was told it had been fitted as written.
  span <- tryCatch(term_span(out, declared_terms(out)), error = function(e) NA_integer_)
  if (!is.na(span) && span < nrow(cells))
    message("the formula spans ", span, " of the ", nrow(cells), " realized ",
            "conditions, so it is fitted with one coefficient per effect it ",
            "names rather than one per condition. Effects the design cannot ",
            "inform are left out, so what remains is full rank and nothing is ",
            "held at zero. The coefficients are named for those columns: ",
            "nest_summary() reads them back as effects of your own variables, ",
            "and reduced_design() shows the columns themselves. (This counts ",
            "the mean structure alone; covariates and random terms multiply ",
            "those columns.)")
  out
}

## The ancestors of a nested variable, root first. Position in the family
## vector is not ancestry once a parent holds more than one child, so every
## caller that needs the strata above a variable asks here.
nest_ancestors <- function(spec, v) {
  p <- spec$parent
  out <- character(0)
  while (!is.null(p) && v %in% names(p)) { v <- unname(p[[v]]); out <- c(v, out) }
  out
}

## A family reads as `a > b > c` while each variable holds at most one child,
## and as its declarations otherwise, which stays unambiguous when it branches.
family_label <- function(x, fam) {
  mark <- function(v) if (v %in% x$cont_nested) paste0(v, " (continuous)") else v
  if (length(fam) == 1L) return(paste(mark(fam), "(crossed)"))
  kids <- vapply(fam, function(v) sum(x$parent == v), 1L)
  if (all(kids <= 1))
    return(paste(vapply(fam, mark, ""), collapse = " > "))
  paste(vapply(fam[-1], function(v)
    paste(unname(x$parent[[v]]), ">", mark(v)), ""), collapse = ", ")
}

print.nesting_spec <- function(x, ...) {
  for (i in seq_along(x$cat_families)) {
    cat("Nesting:", family_label(x, x$families[[i]]), "\n")
    cat("  realized cells:", nrow(x$cells_by_family[[i]]), "of",
        prod(vapply(x$cat_families[[i]],
                    function(v) length(unique(x$data[[v]])), 1L)),
        "in the full crossing\n")
  }
  cat("Cell factor `", x$cell_name, "`: ", length(x$cell_levels),
      " realized cells of ",
      prod(vapply(x$cell_vars, function(v) length(unique(x$data[[v]])), 1L)),
      " in the full crossing\n", sep = "")
  if (length(x$covariates)) cat("Covariates:", paste(x$covariates, collapse = ", "), "\n")
  md <- tryCatch(as.character(fitting_mode(x)), error = function(e) "cells")
  ## What is reported is always stated in the variables the user wrote. `cell`
  ## is an internal factor whose levels are the realized conditions: efficient
  ## to fit on and meaningless to read, so it belongs in the emitted code and
  ## nowhere a person is being told what the model is.
  tm <- tryCatch(declared_terms(x), error = function(e) character(0))
  cat("Structure fitted:", paste(x$outcome, "~", paste(tm, collapse = " + ")), "\n")
  if (identical(md, "cells"))
    cat("  (", nrow(x$cells), " of the ", nrow(x$cells), " realized cells - the ",
        "formula crosses the structure fully,\n   so every condition has its own ",
        "mean and the fit uses the `", x$cell_name, "` factor,\n   whose ",
        "coefficients are those means)\n", sep = "")
  else if (identical(md, "reduced"))
    cat("  (", term_span(x, tm), " of the ", nrow(x$cells),
        " realized cells - fitted as written, one column per effect the\n",
        "   design can inform. reduced_design() shows them)\n", sep = "")
  else
    cat("  (effects parameterization: requested by the engine or the priors)\n")
  invisible(x)
}

## Convert a variable's structurally undefined rows to an explicit sentinel.
## NA cannot be left as it is: R's model machinery deletes those rows casewise
## and silently, and the rows are the design rather than an accident. `where`
## marks the rows the structure says are undefined, so an NA outside them is
## genuine missing data and is refused rather than coded as a condition.
##
## This is not a step the user takes. Where the structure is declared, the
## declaration says which rows should be undefined; where it is inferred, the
## inference has just established it. Either way the check has an answer to
## check against, which is what asking the user to do it by hand never did.
set_sentinel <- function(data, var, where, sentinel = "none") {
  x <- data[[var]]
  if (any(is.na(x) & !where))
    stop(sum(is.na(x) & !where), " NA value(s) in `", var, "` fall outside the ",
         "rows the structure says are undefined. Those are missing data rather ",
         "than structure - R would delete them casewise and silently - and ",
         "coding them as a condition would analyse them as one. Resolve them, ",
         "or declare a structure that accounts for them.")
  if (any(!is.na(x) & where))
    stop(sum(!is.na(x) & where), " value(s) of `", var, "` occur in rows where ",
         "the structure says it is undefined. Either the declaration or the ",
         "data is wrong; resolve the inconsistency before fitting.")
  if (is.numeric(x)) {
    if (!is.numeric(sentinel)) sentinel <- 0
    x[where] <- sentinel
  } else {
    x <- as.character(x); x[where] <- as.character(sentinel)
    x <- factor(x, levels = c(as.character(sentinel),
                              sort(setdiff(unique(x), as.character(sentinel)))))
  }
  data[[var]] <- x
  data
}


## `||` is refused where a factor is involved, because it does not do what it is
## usually taken to do. The operator acts while the formula is parsed, before
## any model matrix exists, so it cannot tell a factor from a numeric variable:
## it splits the term and leaves the levels within each part correlated. The
## result is not a diagonal covariance, and with several factors it can carry
## more parameters than the single-bar term it was meant to simplify. On numeric
## terms it behaves as expected and is allowed through.
check_double_bar <- function(bar_labs, data) {
  for (b in bar_labs) {
    if (!grepl("||", b, fixed = TRUE)) next
    lhs <- trimws(strsplit(b, "||", fixed = TRUE)[[1]][1])
    vars <- all.vars(stats::as.formula(paste("~", lhs)))
    fac <- vars[vapply(vars, function(v)
      v %in% names(data) && (is.factor(data[[v]]) || is.character(data[[v]])), TRUE)]
    if (!length(fac)) next
    stop("`||` in `(", b, ")` involves the factor(s) `",
         paste(fac, collapse = "`, `"),
         "`, and does not give the diagonal covariance it is usually taken to ",
         "give. The operator acts while the formula is parsed, before any model ",
         "matrix exists, so it cannot tell a factor from a numeric variable: it ",
         "splits the term and leaves the levels within each part correlated. ",
         "For a diagonal covariance state it directly - `diag(...)` in lme4 ",
         "2.0-0 and later, or one indicator variable per level - or use ",
         "afex::mixed(expand_re = TRUE) or glmmTMB. For an unstructured ",
         "covariance over the realized cells, declare the full structure with a ",
         "single bar and let nestimand translate it.")
  }
  invisible(TRUE)
}
