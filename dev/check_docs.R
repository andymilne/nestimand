## Documentation audit: run before committing a signature change.
##   Rscript dev/check_docs.R
## Checks that every \usage matches the real signature, every argument has an
## \item, every export is documented, and no \link points nowhere.
for (f in list.files("R", full.names = TRUE)) source(f)
fail <- 0
say <- function(...) { cat(...); fail <<- fail + 1 }

al_of <- function(fl) sub("^\\\\alias\\{(.*)\\}$", "\\1",
                          grep("^\\\\alias\\{", readLines(fl), value = TRUE))
pages <- list.files("man", full.names = TRUE)
aliases <- unlist(lapply(pages, al_of))

## 1. usage matches the signature
for (fl in pages) {
  x <- readLines(fl)
  i <- grep("^\\\\usage\\{$", x); if (!length(i)) next
  j <- grep("^\\}$", x); j <- j[j > i][1]
  u <- paste(x[(i + 1):(j - 1)], collapse = " ")
  for (a in al_of(fl)) {
    if (!exists(a) || !is.function(get(a))) next
    if (!grepl(paste0("(^| )", gsub("\\.", "\\\\.", a), "\\("), u)) next
    for (p in setdiff(names(formals(get(a))), "..."))
      if (!grepl(paste0("[ (]", gsub("\\.", "\\\\.", p), "( |=|,|\\))"), u))
        say(basename(fl), ": `", a, "` usage omits `", p, "`\n")
  }
}
## 2a. no argument documented twice
for (fl in pages) {
  it <- sub("^\\s*\\\\item\\{([^}]*)\\}.*", "\\1",
            grep("\\\\item\\{", readLines(fl), value = TRUE))
  for (d in names(which(table(it) > 1)))
    say(basename(fl), ": `", d, "` is documented twice\n")
}
## 2. every argument has an \item
for (fl in pages) {
  x <- readLines(fl)
  it <- unlist(strsplit(sub("^\\s*\\\\item\\{([^}]*)\\}.*", "\\1",
                            grep("\\\\item\\{", x, value = TRUE)), ",\\s*"))
  fns <- Filter(function(a) exists(a) && is.function(get(a)), al_of(fl))
  need <- unique(unlist(lapply(fns, function(a)
    setdiff(names(formals(get(a))), "..."))))
  for (p in setdiff(need, it)) say(basename(fl), ": `", p, "` has no \\item\n")
}
## 3. exports documented, links resolve
ex <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", readLines("NAMESPACE"), value = TRUE))
for (e in setdiff(ex, aliases)) say("export `", e, "` is undocumented\n")
lk <- unique(unlist(lapply(pages, function(fl) {
  t <- paste(readLines(fl), collapse = " ")
  sub("\\\\link\\{(.*)\\}", "\\1", regmatches(t, gregexpr("\\\\link\\{[^}]*\\}", t))[[1]])
})))
for (l in setdiff(lk, aliases)) say("\\link{", l, "} points nowhere\n")
## 4. Rd validity
for (fl in pages) { r <- tools::checkRd(fl)
  if (length(r)) { say(basename(fl), ": ") ; print(r) } }

cat(if (fail) sprintf("\n%d documentation problem(s)\n", fail)
    else "documentation consistent\n")
