\name{nest_summary}
\alias{nest_summary}
\alias{print.nestimand_summary}
\title{The fitted model in the original parameterization}
\description{
A cell fit carries one coefficient per realized condition, which is what makes
estimation well posed but is not what a report states. The identified effect
coefficients are an exact linear function of the cell means, so
\code{nest_summary} produces the coefficient table a reader expects without
refitting anything, and with the standard errors a chain fit would give.
}
\usage{
nest_summary(model, space = c("effects", "cells"), conf_level = 0.95, 
             random = FALSE, data = NULL, spec = NULL)
}
\arguments{
  \item{model}{a fitted model, ideally from \code{\link{nest_fit}}.}
  \item{spec}{a \code{nesting_spec} object, needed only for a model fitted outside \code{\link{nest_fit}}.}
  \item{space}{\code{"effects"} for the identified chain coefficients, \code{"cells"} for the realized-cell means.}
  \item{conf_level}{confidence level for the intervals.}
  \item{data}{data from which the cell grid is built; the spec's data by default.}
  \item{x}{an object to print.}
  \item{digits}{digits for the printed estimates and standard errors.}
  \item{random}{also read out the random-effects covariance, in the same parameterization, and print it beneath the fixed effects. See \code{\link{random_covariance}}.}
}
\details{
In effect space the \code{term} labels are those R derives from the design
matrix, and in a partially nested design they mislead: \code{chord_typemaj}
names a comparison at one inversion, not an average over them, and which one
depends on the pivot. Each row therefore carries a \code{meaning} column giving
the combination of realized conditions that the coefficient equals, read from
the design rather than from the label. In cell space the labels are the
conditions themselves and need no such gloss.

These are the model's own coefficients, reparameterized: the transformation is
exact and no refit or reference grid is involved.

A covariate that interacts with the conditions contributes one slope per cell.
That is a vector over the same space as the cell means, so it translates by the
same map: in effect space it becomes a reference-cell slope, named for the
covariate alone, plus the differences from it, which is what a chain fit
reports. A covariate entering additively is one slope common to every condition
and is left as fitted, as are ordinal thresholds. As with any intercept, a cell mean is the linear
predictor where the covariates are zero, so an uncentred covariate makes it an
extrapolation; centring covariates gives cell means at the covariate average
instead. Under a threshold family the reference
condition has no free coefficient, its latent mean being absorbed into the
thresholds; that row reads as exactly zero and is flagged.
}
\value{
A data frame of class \code{nestimand_summary}, printed with a header naming
the parameterization, the fitting function, and the formula that was fitted -
since the same table can come from any of the engines, and which one matters for
reading it. The engine is also in the \code{nestimand_fit} attribute with one row per coefficient,
giving the estimate, standard error, test statistic, p-value, interval, and the
combination of conditions the coefficient equals. The transformation matrix is
in the \code{nestimand_map} attribute.
}
\seealso{
\code{\link{nest_estimand}} for contrasts under a stated weighting policy,
\code{\link{effect_basis}} for the matrix relating the two spaces,
\code{\link{nest_fit}}.
}
\examples{
\dontrun{
m <- nest_fit(dat, response ~ chord_type * inversion + training,
              "inversion \%in\% chord_type")
nest_summary(m)                # effect coefficients, with their meanings
nest_summary(m, "cells")       # the realized-cell means
}
}
