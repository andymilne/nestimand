for (f in list.files("R", full.names = TRUE)) source(f)
source("dev/demo_data.R")
suppressPackageStartupMessages(library(lme4))
set.seed(9)
## multi-trial data: 6 trials per participant per cell, so a rich RE structure
## is estimable in principle
dat6 <- do.call(rbind, lapply(1:6, function(i) {
  d <- dat; d$response <- d$response + rnorm(nrow(d), 0, 0.4); d }))
dat6$cell <- droplevels(interaction(dat6$chord_type, dat6$inversion, sep = "."))
cat("participants:", nlevels(dat6$participant), " rows:", nrow(dat6),
    " trials per participant per cell:", nrow(dat6)/nlevels(dat6$participant)/10, "\n\n")

## (1) the random slope the chain parameterization would need
X <- model.matrix(~ 0 + chord_type + chord_type:inversion, dat6)
cat("chain RE design columns:", ncol(X), " rank:", qr(X)$rank,
    " all-zero columns:", sum(colSums(X != 0) == 0), "\n")
r1 <- tryCatch({
  m <- lmer(response ~ 0 + cell + (0 + chord_type + chord_type:inversion | participant),
            data = dat6, control = lmerControl(check.nobs.vs.nRE = "ignore"))
  paste("fitted; RE variances:", ncol(VarCorr(m)$participant))
}, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 90)))
cat("chain random slope:", r1, "\n\n")

## (2) the cell random effect: unstructured over realized cells
m_cell <- lmer(response ~ 0 + cell + (0 + cell | participant), data = dat6,
               control = lmerControl(check.nobs.vs.nRE = "ignore"))
V <- VarCorr(m_cell)$participant
cat("cells unstructured: ", nrow(V), " variances + ", nrow(V)*(nrow(V)-1)/2,
    " correlations = ", nrow(V) + nrow(V)*(nrow(V)-1)/2, " parameters\n", sep = "")

## (3) the chain grouping submodel, which is what the chain form could express
m_chain <- lmer(response ~ 0 + cell + (1 | participant) +
                (1 | participant:chord_type) + (1 | participant:cell),
                data = dat6)
cat("chain grouping   : 3 variances, implying compound symmetry\n")
cr <- cov2cor(V)
cat("  cell-form correlations range:", paste(round(range(cr[upper.tri(cr)]), 3), collapse = " to "), "\n")
cat("  logLik  cells:", round(logLik(m_chain), 1), " ->", round(logLik(m_cell), 1),
    "  (df", attr(logLik(m_chain), "df"), "->", attr(logLik(m_cell), "df"), ")\n")
cat("  LRT p:", format.pval(anova(m_chain, m_cell)$`Pr(>Chisq)`[2], digits = 3), "\n")

## (4) diagonal-only, the middle option
m_diag <- lmer(response ~ 0 + cell + (0 + cell || participant), data = dat6,
               control = lmerControl(check.nobs.vs.nRE = "ignore"))
cat("cells diagonal   : 10 variances, no correlations; logLik",
    round(logLik(m_diag), 1), "\n")
