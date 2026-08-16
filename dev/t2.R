for (f in list.files("R", full.names = TRUE)) source(f)
source("dev/demo_data.R")
suppressPackageStartupMessages(library(marginaleffects))
sp <- nesting_spec(dat, response ~ chord_type * inversion + training,
                   nests = "inversion %in% chord_type")
m  <- lm(cell_formula(sp), data = sp$data)

est <- function(pol) {
  g <- counterfactual_grid(sp, sp$data, pol)
  e <- avg_predictions(m, newdata = g, by = "chord_type", wts = g$.w,
                       hypothesis = "pairwise")
  as.data.frame(e)
}
for (k in c("equal", "proportional", "hierarchical")) {
  p <- nest_policy(sp, "chord_type", k)
  e <- est(p)
  cat(sprintf("%-14s aug - maj = %.4f\n", k, e$estimate[e$term == "aug - maj"]))
}
print(nest_policy(sp, "chord_type", "equal"))
## supplied policy: convex combination check
p3 <- nest_policy(sp, "chord_type", c("0" = 0.5, "1" = 0.3, "2" = 0.2))
print(p3)
e3 <- est(p3); v3 <- e3$estimate[e3$term == "aug - maj"]
## single-version vertices
vert <- sapply(c("0","1","2"), function(iv) {
  p <- nest_policy(sp, "chord_type", "nominated", at = c(inversion = iv))
  e <- est(p); e$estimate[e$term == "aug - maj"]
})
cat("vertices aug-maj:", paste(sprintf("%.4f", vert), collapse = " "), "\n")
cat("supplied p = 0.5/0.3/0.2 :", sprintf("%.6f", v3),
    " convex combination:", sprintf("%.6f", sum(vert * c(.5,.3,.2))), "\n")
cat("bounds:", sprintf("%.3f to %.3f", min(vert), max(vert)),
    " equal-weight point:", sprintf("%.3f", est(nest_policy(sp,"chord_type","equal"))$estimate[3]), "\n")
