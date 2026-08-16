# nestimand

Estimands for partially nested designs, by translation into realized-cell space.

Declare the nesting structure once. The package fits the design in the
realized-cell parameterization, where estimation is unconditionally well posed —
full rank, no aliasing, no identification constraints — and translates grids,
labels, draws, priors, and results back into the original variable space in
which research questions, priors, and reports are stated.

The decisive reason is the random effects. Written as a chain of interaction
terms, a random slope over the nesting structure is rank-deficient by
construction: the columns for impossible combinations are identically zero for
every group, leaving dimensions of the covariance that no amount of data can
identify. No engine refuses such a structure — lme4 fits it and reports a
convergence warning, easily mistaken for sparse data — so the chain form is
confined in practice to nested grouping factors, which impose compound symmetry.
Placing the whole categorical structure on a single realized-cell factor makes
every random-effects structure R offers available and identified: unstructured,
diagonal, any brms covariance, or the grouping chain as a constrained submodel.

## Install

From GitHub:

```r
# install.packages("remotes")
remotes::install_github("<user>/nestimand")
```

Or from a local copy of the repository, with the working directory set to its
parent:

```r
install.packages("nestimand", repos = NULL, type = "source")
```

`marginaleffects` is needed for the estimand functions; `lme4`, `ordinal`, and
`brms` for the corresponding engines. None is required to install the package.

## Use

```r
library(nestimand)

sp <- nesting_spec(dat, response ~ chord_type * inversion + training,
                   "inversion %in% chord_type")
m  <- nest_fit(sp)
e  <- estimand(m, chord_type, policy = "equal")

e             # contrasts in the original variable names, with bounds
show_code(e)  # the code that produced them, for saving and adapting
```

Every across-boundary contrast requires a weighting policy — a distribution over
the versions of a compound condition — because no such contrast is defined
without one. The policy is an argument, never a silent default, and the printed
output carries the partial-identification bounds over all admissible policies
alongside the point estimate.

## Help

`?nestimand` gives an overview, with the functions grouped by how often they are
needed: the five that cover an ordinary analysis, those the analysis may call
for, those that report what the translation did, and those exported only so that
the code from `show_code()` will run.

## Tests

From the package root:

```r
source("tests/test_core.R")   # seconds, no sampling
source("tests/test_brms.R")   # needs a working brms
```

The brms checks compile Stan programs. Installing `cmdstanr` and setting
`options(brms.backend = "cmdstanr")` reduces compilation from minutes to
seconds; the test script detects and announces it automatically.

`NOTES.md` records what has been verified by execution and what has not.
