## Rewrite every \usage{} block from the functions' actual signatures.
## Run after changing any exported signature:  Rscript dev/sync_usage.R
for (f in list.files("R", full.names = TRUE)) source(f)

sig <- function(n) {
  a <- deparse(args(get(n)), width.cutoff = 60L)
  a <- paste(a[-length(a)], collapse = "\n")
  a <- sub("^function ", n, a)
  ## indent continuation lines under the opening parenthesis
  ln <- strsplit(a, "\n")[[1]]
  if (length(ln) > 1) {
    pad <- strrep(" ", nchar(n) + 1L)
    ln[-1] <- paste0(pad, trimws(ln[-1]))
  }
  paste(ln, collapse = "\n")
}

changed <- character(0)
for (fl in list.files("man", full.names = TRUE)) {
  x <- readLines(fl)
  i <- grep("^\\\\usage\\{$", x); j <- grep("^\\}$", x)
  if (!length(i)) next
  j <- j[j > i][1]
  al <- sub("^\\\\alias\\{(.*)\\}$", "\\1", grep("^\\\\alias\\{", x, value = TRUE))
  fns <- al[vapply(al, function(a) exists(a) && is.function(get(a)) &&
                     any(grepl(paste0("(^|[^A-Za-z0-9._])",
                                      gsub("\\.", "\\\\.", a), "\\("),
                               x[(i + 1):(j - 1)])), TRUE)]
  if (!length(fns)) next
  new <- c("\\usage{", unlist(lapply(fns, function(a) c(sig(a), ""))), "}")
  new <- new[-(length(new) - 1)]
  if (!identical(x[i:j], new)) {
    writeLines(c(x[seq_len(i - 1)], new, x[(j + 1):length(x)]), fl)
    changed <- c(changed, basename(fl))
  }
}
cat(if (length(changed)) paste("rewritten:", paste(changed, collapse = ", "))
    else "all usage sections already match", "\n")
