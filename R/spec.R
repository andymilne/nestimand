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

nesting_spec <- function(data, formula, nests,
                         fit = c("lm", "glm", "lmer", "glmer", "clm", "clmm", "brm"),
                         family = NULL, random = NULL,
                         cell_name = "cell") {
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
  ne <- substitute(nests)
  ne_txt <- paste(deparse(ne), collapse = " ")
  ne_val <- tryCatch(eval(ne, parent.frame()), error = function(e) NULL)
  nests <- if (is.character(ne_val)) ne_val
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
  for (v in nest_vars) {
    if (!v %in% names(data)) stop("variable not in data: ", v)
    if (anyNA(data[[v]]))
      stop("`", v, "` contains NA. Code the undefined state as an explicit ",
           "sentinel (a factor level such as \"none\", or 0 for a numeric); ",
           "NA triggers silent casewise deletion in R's model machinery. ",
           "If the NA values mark structurally undefined rows, ",
           "apply_sentinel() converts them safely.")
  }
  for (ch in names(parent))
    if (is.numeric(data[[parent[[ch]]]]))
      stop("a continuous variable cannot nest another variable: ", parent[[ch]])
  if (cell_name %in% names(data))
    stop("`", cell_name, "` is already a column of the data. The realized-cell ",
         "factor is constructed under that name and would overwrite it; supply ",
         "cell_name = to choose another, or rename the column.")

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
            "cells, so it asks for less than the saturated structure - a `+` ",
            "where the design would admit a `*`. It is fitted as written: the ",
            "design the formula gives over the realized cells, with the columns ",
            "it cannot inform removed, which is full rank and needs nothing held ",
            "at zero. `~ 0 + cell` is used instead only when the formula does ",
            "cross the structure fully, where the two are the same model. ",
            "reduced_design() shows the columns and the effect each stands for. ",
            "(This counts the mean structure alone; covariates and random terms ",
            "multiply those columns.)")
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

apply_sentinel <- function(data, var, where = NULL, sentinel = "none") {
  ## Convert NA to a sentinel where the design says the variable is undefined.
  ## `where` marks the structurally undefined rows, e.g. data$chord_type == "aug".
  ## Supplying it is worthwhile: NA conflates two things this package keeps
  ## apart, and with `where` given, an NA outside those rows is genuine
  ## missingness and is refused, while a non-NA inside them is inconsistent
  ## coding and is likewise refused. Omitting it converts every NA and warns.
  x <- data[[var]]
  if (is.null(where)) {
    n <- sum(is.na(x))
    if (n)
      warning("converting all ", n, " NA value(s) in `", var, "` to the sentinel ",
              "without checking which are structural. Any that are genuine ",
              "missing data are now coded as an existing condition, and will be ",
              "analysed as one. Supplying `where` - a logical vector marking the ",
              "rows in which `", var, "` is undefined by design, such as ",
              "`data$parent == \"level\"` - restricts the conversion to those ",
              "rows and refuses the rest.", call. = FALSE)
    where <- is.na(x)
  } else {
    if (any(is.na(x) & !where))
      stop(sum(is.na(x) & !where), " NA value(s) in `", var,
           "` fall outside the declared undefined rows. Those are genuine ",
           "missing data, not structure; do not convert them to a sentinel.")
    if (any(!is.na(x) & where))
      stop(sum(!is.na(x) & where), " non-NA value(s) in `", var,
           "` occur inside the declared undefined rows; resolve the ",
           "inconsistent coding before applying a sentinel.")
  }
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
