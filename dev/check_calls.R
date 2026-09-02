## nestimand: does the package's own code call the package correctly? --------
## tests/test_brms.R runs only where Stan compiles, so it drifts against the
## signatures while the rest of the suite keeps pace: three separate faults in
## it - a spec passed positionally into `target`, a `weights` that was not an
## argument, contrast labels from a superseded direction convention - were each
## found by a user running it rather than by anything here.
##
## This is what can be checked without running anything: every call to a
## nestimand function in tests/ and dev/ is matched against that function's
## formals. A named argument the function does not have, or more positional
## arguments than it can take, is reported with its file and line. Argument
## *order* is beyond static analysis - `check_target()` catches the important
## case of that at runtime - so this narrows the gap rather than closing it.
##
##   Rscript dev/check_calls.R

for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
known <- ls(envir = globalenv())
known <- known[vapply(known, function(n) is.function(get(n)), TRUE)]

problems <- list()
note <- function(...) problems[[length(problems) + 1]] <<- paste0(...)

## Which of a file's own symbols hold a declaration and which hold a fit, taken
## from what they are assigned. Enough to tell a shifted call from a right one.
built_by <- function(ex, fname) {
  out <- character(0)
  rec <- function(e) {
    if (!is.call(e)) return(invisible())
    if (length(e) == 3L && as.character(e[[1]])[1] %in% c("<-", "=", "<<-") &&
        is.name(e[[2]]) && is.call(e[[3]]) && is.name(e[[3]][[1]]) &&
        as.character(e[[3]][[1]]) == fname)
      out <<- c(out, as.character(e[[2]]))
    for (i in seq_along(e)) tryCatch(rec(e[[i]]), error = function(x) NULL)
    invisible()
  }
  for (i in seq_along(ex)) rec(ex[[i]])
  unique(out)
}
specs <- character(0); fits <- character(0)

walk <- function(e, file, line) {
  if (!is.call(e)) return(invisible())
  ln <- attr(e, "srcref"); if (!is.null(ln)) line <- ln[1]
  fn <- e[[1]]
  if (is.name(fn) && as.character(fn) %in% known) {
    nm <- as.character(fn)
    fmls <- names(formals(get(nm)))
    dots <- match("...", fmls, nomatch = length(fmls) + 1L)
    args <- as.list(e)[-1]
    given <- names(args); if (is.null(given)) given <- rep("", length(args))
    ## R matches an argument name exactly, then by unique prefix
    matched <- vapply(given, function(g) {
      if (!nzchar(g)) return(NA_character_)
      if (g %in% fmls) return(g)
      hit <- fmls[startsWith(fmls, g)]
      if (length(hit) == 1L) hit else NA_character_
    }, "")
    unknown <- given[nzchar(given) & is.na(matched)]
    if (length(unknown) && dots > length(fmls))
      note(file, ":", line, "  ", nm, "() has no argument `",
           paste(unknown, collapse = "`, `"), "`")
    ## the positional arguments fill the formals the names did not take, in
    ## order and stopping at the dots. One landing on a formal already given by
    ## name is the mistake a shifted call makes, and R reports it only when the
    ## call runs - which for tests/test_brms.R means only where Stan compiles.
    free <- setdiff(fmls[seq_len(dots - 1L)], stats::na.omit(matched))
    npos <- sum(!nzchar(given))
    if (npos > length(free))
      note(file, ":", line, "  ", nm, "() has ", length(free),
           " argument(s) left for position after the named ones, and ", npos,
           " were passed positionally: `",
           paste(utils::tail(fmls[seq_len(dots - 1L)], 1), collapse = ""),
           "` onwards would be matched twice or overflow into `...`")
    ## Where a positional argument lands, when what is passed is recognisable.
    ## A file that builds a spec with nesting_spec() and a fit with nest_fit()
    ## says which of its own symbols are which, and a spec reaching any formal
    ## but `spec` - or a fit reaching any but `model` - is a shifted call. This
    ## is the fault that put a nesting_spec into `target` in four places, and it
    ## is not an error R can raise until the call runs. It also catches a
    ## declaration handed to nest_fit(), which now takes the arguments a
    ## declaration is made from rather than the declaration itself.
    pos <- args[!nzchar(given)]
    slot <- free[seq_along(pos)]
    for (k in seq_along(pos)) {
      if (!is.name(pos[[k]]) || is.na(slot[k])) next
      v <- as.character(pos[[k]])
      if (v %in% specs && !identical(slot[k], "spec"))
        note(file, ":", line, "  ", nm, "() would take `", v,
             "`, which this file builds with nesting_spec(), as `", slot[k],
             "`.")
      if (v %in% fits && "model" %in% fmls && !identical(slot[k], "model"))
        note(file, ":", line, "  ", nm, "() would take `", v,
             "`, which this file builds with nest_fit(), as `", slot[k], "`.")
    }
  }
  for (i in seq_along(e))
    tryCatch(walk(e[[i]], file, line), error = function(x) NULL)
  invisible()
}

files <- c(list.files("tests", pattern = "[.]R$", full.names = TRUE),
           list.files("dev", pattern = "[.]R$", full.names = TRUE))
files <- setdiff(files, "dev/check_calls.R")
for (f in files) {
  ex <- tryCatch(parse(f, keep.source = TRUE), error = function(e) {
    note(f, ":  does not parse: ", conditionMessage(e)); NULL })
  if (is.null(ex)) next
  specs <- built_by(ex, "nesting_spec")
  fits  <- built_by(ex, "nest_fit")
  for (i in seq_along(ex)) walk(ex[[i]], f, getSrcLocation(ex[i], "line")[1])
}

if (!length(problems)) cat("calls consistent with the signatures\n") else {
  cat(paste(unlist(problems), collapse = "\n"), "\n\n",
      length(problems), " call problem(s)\n", sep = "")
}
