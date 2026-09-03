\name{nestimand-package}
\alias{nestimand-package}
\alias{nestimand}
\docType{package}
\title{Estimands for Partially Nested Designs by Cell Translation}
\description{
A partially nested design contains a factor whose levels exist only within
certain levels of another: chord inversion, say, which is undefined for
augmented triads. Default modelling practice mishandles such designs quietly,
producing aliased coefficients with misleading labels and marginal summaries
that change when factor levels are permuted.

\pkg{nestimand} fits the model your formula names over the conditions your
design realizes, and translates grids, labels, draws, priors and results back
into the variables you wrote. Every contrast that crosses a structural boundary
requires an explicit weighting policy, because no such contrast is defined
without one.
}
\section{Three verbs}{
An analysis is three calls.

\tabular{ll}{
  \code{\link{nest_fit}} \tab fit the model, keeping the code that produced it \cr
  \code{\link{nest_estimand}} \tab contrasts in the original variable names, with bounds \cr
  \code{\link{nest_summary}} \tab the coefficient table, read back as effects \cr
}

\preformatted{m <- nest_fit("lm", rating ~ chord_type * inversion + training, dat)
e <- nest_estimand(m, chord_type, policy = "equal")
nest_summary(m)}

The same call gives the predictions the contrasts are formed from, under the same
policy, so the means and the differences reported cannot disagree:

\preformatted{nest_estimand(m, chord_type, policy = "equal", contrast = "none")}

\code{\link{show_code}} returns the script that produced any of them, for
saving, adapting and re-running.
}
\section{Declaring the structure, or leaving it to be read}{
The structure need not be stated. Where a variable is undefined it has no
value, and those rows are exactly the ones some other variable's levels pick
out, so \code{\link{nest_fit}} reads the structure off a data frame carrying
\code{NA}:

\preformatted{## inversion is NA wherever chord_type is "aug"
m <- nest_fit("lm", rating ~ chord_type * inversion, dat)
#> nesting structure read from the data: nests = "inversion \%in\% chord_type"}

What was inferred is always said, in the syntax you would have written, and an
\code{NA} pattern no combination of the other variables accounts for is refused
as missing data rather than read as structure. Declare it yourself when the data
cannot settle it, or when you mean something else:

\preformatted{m <- nest_fit("lm", rating ~ chord_type * inversion, dat,
              nests = "inversion \%in\% chord_type")}

\code{\link{nesting_spec}} makes the declaration without fitting, which is a way
of seeing the structure and the conditions it realizes before committing to a
long fit.
}
\section{Deeper and branching structures}{
A nested variable may itself hold nested variables, to any depth, and one parent
may hold several. Each link is one declaration, and links sharing a parent may
be gathered:

\preformatted{## inversion and X1 are both undefined for augmented triads;
## Z exists only within an inversion
nest_fit("lm", rating ~ chord_type * inversion * X1 * Z, dat,
         nests = c("inversion \%in\% chord_type", "X1 \%in\% chord_type",
                   "Z \%in\% inversion"))

## the two siblings, gathered - these four are the same declaration
c(inversion, X1) \%in\% chord_type
(inversion + X1) \%in\% chord_type
inversion + X1 \%in\% chord_type
c("inversion \%in\% chord_type", "X1 \%in\% chord_type")}

Only the nesting is declared. A categorical variable the formula crosses with
the declared structure is part of the design and joins the cell factor by
itself - \code{top_note} below needs no entry - so it can be the target of an
estimand like any other design variable:

\preformatted{nest_fit("lm", rating ~ chord_type * inversion * top_note, dat,
         nests = "inversion \%in\% chord_type")}

Three kinds of variable stay covariates: one entering additively, whose effect
is thereby declared common to every condition; a numeric one; and an ordered
factor, whose contrasts say it is meant as a quantity. Naming one in
\code{nests} on its own folds it into the cell factor deliberately.

The strata of a nested variable are its own ancestors, and a sibling is not one
of them: \code{nest_estimand(m, X1, by = "chord_type")} compares the levels of
\code{X1} within each chord type, while \code{nest_estimand(m, Z, by =
c("chord_type", "inversion"))} compares within each chord type and inversion. A
variable declared inside two parents is refused, as is \code{\%in\%} chained
within one entry - \code{Z \%in\% inversion \%in\% chord_type} is two
declarations, since \R reads it as \code{(Z \%in\% inversion) \%in\% chord_type}.
}
\section{What the model is fitted on}{
The model is fitted over the conditions the design realizes rather than over the
full crossing, so every coefficient is estimable and none is aliased. Where the
formula crosses the structure fully, the coefficients are the condition means
themselves; where it asks for less, there is one coefficient per effect the
formula names, and a random term is fitted on the effects its bar names.

This does not change how results are read. Everything \code{\link{nest_estimand}}
and \code{\link{nest_summary}} report is stated in the variables you wrote, and
\code{\link{nesting_spec}} prints the structure that will be fitted before you
fit it.
}
\section{Other exported functions}{
Everything else the package exports is machinery. Some of it is worth looking at
to see what was done - \code{\link{cell_formula}}, \code{\link{effect_basis}},
\code{\link{fitting_mode}}, \code{\link{reduced_design}}. The rest is exported
only because the code \code{\link{show_code}} emits has to run, and that code
calls the package rather than inlining cell tables and weight vectors as
literals. None of it is intended to be typed.
}
\references{
Milne, A. J. Partially nested variables: modelling and estimation in R.
Companion document, distributed with the package.
}
\author{Andrew J. Milne}
\keyword{package}
