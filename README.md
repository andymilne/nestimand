# nestimand

Fitting and estimation for partially nested designs — designs containing a
factor whose levels exist only within certain levels of another. Chord
inversion, undefined for augmented triads. A dose that exists only in the
treatment arm. A follow-up question asked only of those who answered yes.

Give `nestimand` your formula and your data. It fits the model your formula
names over the conditions your design actually realizes — full rank, nothing
aliased, nothing held at zero — and translates grids, labels, draws, priors and
results back into the variables you wrote.

Random effects are where this matters most. A random slope written over the
original factors carries columns for combinations that do not exist: identically
zero for every group, so no amount of data can identify those dimensions of the
covariance. No engine refuses such a structure — lme4 fits it and reports a
convergence warning, easily mistaken for sparse data. Written over the columns
the design does identify, every parameter is estimable, and any covariance
structure R offers is available: unstructured, diagonal, any brms covariance.

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

Three verbs cover an analysis.

```r
library(nestimand)

m <- nest_fit("lm", response ~ chord_type * inversion + training, dat)
e <- nest_estimand(m, chord_type, policy = "equal")

e                # contrasts in the original variable names, with bounds
nest_summary(m)  # the coefficient table, read back as effects
show_code(e)     # the code that produced them, for saving and adapting
```

The nesting structure need not be stated: where a variable is undefined it has
no value, so a data frame carrying `NA` says what the structure is, and what was
read off it is reported in the syntax you would have written. Declare it
yourself when the data cannot settle it, or when you mean something else:

```r
m <- nest_fit("lm", response ~ chord_type * inversion, dat,
              nests = "inversion %in% chord_type")
```

Every across-boundary contrast requires a weighting policy — a distribution over
the versions of a compound condition — because no such contrast is defined
without one. The policy is an argument, never a silent default, and the printed
output carries the partial-identification bounds over all admissible policies
alongside the point estimate.

## Help

`?nestimand` gives an overview: the three verbs an ordinary analysis needs, how
the structure is declared or read, and what the model is fitted on. Everything
else the package exports is machinery, exported so that the code `show_code()`
emits will run.

## Tests

From the package root:

```r
source("tests/test_core.R")   # seconds, no sampling
source("tests/test_brms.R")   # needs a working brms
```

The brms checks compile Stan programs. Installing `cmdstanr` and setting
`options(brms.backend = "cmdstanr")` reduces compilation from minutes to
seconds; the test script detects and announces it automatically.

`dev/check_docs.R` and `dev/check_calls.R` audit the help pages and the
package's own calls against the signatures; run both before committing a
signature change. `NOTES.md` records what has been verified by execution and
what has not.
