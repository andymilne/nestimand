## nestimand: one-screen diagnostic of the installed marginaleffects ---------
## Run from the package root:  source("tests/diagnose.R")
## Paste the whole output back; it identifies exactly what changed.
suppressPackageStartupMessages(library(marginaleffects))
local({
  want <- c("spec.R", "translate.R", "policy.R", "estimand.R", "fit.R",
            "latent.R", "priors.R")
  f <- list.files(".", pattern = "[.][Rr]$", recursive = TRUE, full.names = TRUE)
  for (p in f[basename(f) %in% want]) source(p)
})
set.seed(1); n <- 40
.c <- data.frame(chord_type = c(rep(c("dim","min","maj"), each = 3), "aug"),
                 inversion = c(rep(c("0","1","2"), 3), "none"),
                 mu = c(4.16,4.00,3.73, 4.29,4.20,3.71, 4.98,4.62,4.39, 3.92))
dat <- .c[rep(seq_len(nrow(.c)), each = n), ]
dat$participant <- factor(rep(seq_len(n), times = nrow(.c)))
.tr <- runif(n, 0, 10); dat$training <- .tr[as.integer(dat$participant)]
dat$response <- rnorm(nrow(dat), dat$mu + 0.12 * dat$training, 1.2)
dat$chord_type <- factor(dat$chord_type, levels = c("aug","dim","min","maj"))
dat$inversion <- factor(dat$inversion, levels = c("none","0","1","2"))

cat("R              ", R.version.string, "\n")
cat("marginaleffects", as.character(packageVersion("marginaleffects")), "\n")
cat("hypothesis form", if (mfx_formula_hypothesis()) "~pairwise" else '"pairwise"', "\n\n")

sp <- nesting_spec(dat, response ~ chord_type * inversion + training,
                   "inversion %in% chord_type")
m  <- lm(cell_formula(sp), data = sp$data)
pol <- nest_policy(sp, "chord_type", "equal")
g <- counterfactual_grid(sp, sp$data, pol)
e <- avg_predictions(m, newdata = g, by = "chord_type", wts = g$.w,
                     hypothesis = mfx_hypothesis("pairwise"))
d <- as.data.frame(e)
cat("--- columns returned ---\n"); print(names(d))
cat("\n--- first four rows, label and estimate columns ---\n")
lab <- mfx_term_column(d)
cat("label column located:", lab, "\n")
print(utils::head(d[, intersect(c(lab, "term", "hypothesis", "contrast",
                                  "estimate", "std.error"), names(d))], 4))
cat("\n--- expected: a row labelled 'aug - maj' with estimate -0.6779 ---\n")
cat("found:", if (!is.na(lab)) paste(d[[lab]], collapse = " | ") else "no label column", "\n")
