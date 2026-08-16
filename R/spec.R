## nestimand: declaration surface -----------------------------------------
## The declaration surface is carried over from the chain-based prototype
## unchanged in its user-facing form (bare or quoted arguments, validations,
## teaching errors, saturation and outcome-type guards). What changed beneath
## it is the parameterization: the declared structure is now compiled into a
## realized-cell factor rather than into a chain of interaction terms.

nestimand_build <- "2026-08-15.1"

nesting_spec <- function(data, formula, nests,
                         fit = c("lm", "glm", "lmer", "glmer", "clm", "clmm", "brms"),
                         family = NULL, random = NULL,
                         cell_name = "cell") {
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
  ne <- substitute(nests)
  nests <- if (is.call(ne) && identical(ne[[1]], as.name("c")) &&
               !tryCatch(is.character(eval(ne, parent.frame())), error = function(e) FALSE))
    vapply(as.list(ne)[-1], function(x) paste(deparse(x), collapse = " "), "")
  else if (is.call(ne) && identical(ne[[1]], as.name("%in%")))
    paste(deparse(ne), collapse = " ")
  else eval(ne, parent.frame())
  parse1 <- function(s) {
    p <- trimws(strsplit(s, "%in%")[[1]]); stopifnot(length(p) == 2)
    stats::setNames(p[2], p[1])                       # child -> parent
  }
  parent <- do.call(c, lapply(nests, parse1))

  ## --- validation --------------------------------------------------------
  nest_vars <- unique(c(names(parent), parent))
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
  roots <- setdiff(parent, names(parent))
  chain_of <- function(root) {
    ch <- root
    repeat {
      kid <- names(parent)[parent == ch[length(ch)]]
      if (!length(kid)) break
      if (length(kid) > 1) stop("branching nests not supported: ", ch[length(ch)])
      ch <- c(ch, kid)
    }
    ch
  }
  families <- lapply(unique(roots), chain_of)

  ## --- outcome and fit guards (unchanged) --------------------------------
  outcome <- deparse(formula[[2]])
  if (fit %in% c("clm", "clmm") && !(outcome %in% names(data) && is.ordered(data[[outcome]])))
    stop("fit = '", fit, "' needs an ordered-factor outcome; `", outcome, "` is not one.")
  if (fit == "brms" && outcome %in% names(data) && is.ordered(data[[outcome]]) &&
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
  tt   <- stats::terms(formula)
  labs <- attr(tt, "term.labels")
  bar_labs <- c(labs[grepl("|", labs, fixed = TRUE)],
                if (!is.null(random))
                  grep("|", attr(stats::terms(stats::as.formula(
                    paste("~", random))), "term.labels"), value = TRUE, fixed = TRUE))
  labs <- labs[!grepl("|", labs, fixed = TRUE)]
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

  covariates <- setdiff(all.vars(formula[[3]]), c(unlist(families), "|"))
  covariates <- setdiff(covariates, unlist(lapply(bar_labs, function(b)
    all.vars(stats::as.formula(paste("~", gsub("|", "+", b, fixed = TRUE)))))))
  covariates <- unique(c(covariates, cont_nested))

  ## does the declared formula ask for covariate x structure interaction?
  cov_by_cell <- vapply(covariates, function(cv)
    any(vapply(strsplit(labs, ":"), function(vs)
      cv %in% vs && any(unlist(families) %in% vs), TRUE)), TRUE)
  names(cov_by_cell) <- covariates

  counts <- table(data[[cell_name]])
  if (all(counts == 1))
    message("every realized cell contains a single observation: the model will be ",
            "saturated (zero residual degrees of freedom), so standard errors and ",
            "tests will be unavailable. Estimates remain valid as cell-mean ",
            "contrasts; inference needs trial-level data.")

  structure(list(data = data, outcome = outcome, formula_in = formula,
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
}

print.nesting_spec <- function(x, ...) {
  for (i in seq_along(x$cat_families)) {
    cat("Nesting chain:", paste(x$cat_families[[i]], collapse = " > "), "\n")
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
  cat("Fitting formula:", paste(deparse(cell_formula(x)), collapse = " "), "\n")
  invisible(x)
}

apply_sentinel <- function(data, var, where, sentinel = "none") {
  ## Convert NA to a sentinel ONLY where the design says the variable is
  ## undefined. `where`: logical vector marking the structurally undefined
  ## rows (e.g. data$chord_type == "aug"). Both mismatches are errors:
  ## NA outside `where` is genuine missingness, not structure; non-NA
  ## inside `where` means the coding is already inconsistent.
  x <- data[[var]]
  if (any(is.na(x) & !where))
    stop(sum(is.na(x) & !where), " NA value(s) in `", var,
         "` fall outside the declared undefined rows. Those are genuine ",
         "missing data, not structure; do not convert them to a sentinel.")
  if (any(!is.na(x) & where))
    stop(sum(!is.na(x) & where), " non-NA value(s) in `", var,
         "` occur inside the declared undefined rows; resolve the ",
         "inconsistent coding before applying a sentinel.")
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
