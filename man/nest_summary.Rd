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
nest_summary(model, spec = NULL, space = c("effects", "cells"),
             conf_level = 0.95, data = NULL)
}
\arguments{
  \item{model}{a fitted model, ideally from \code{\link{nest_fit}}.}
  \item{spec}{a \code{nesting_spec} object. Optional: a model fitted by \code{\link{nest_fit}} carries its own declaration.}
  \item{space}{\code{"effects"} for the identified chain coefficients, \code{"cells"} for the realized-cell means.}
  \item{conf_level}{confidence level for the intervals.}
  \item{data}{data from which the cell grid is built; the spec's data by default.}
}
\details{
Coefficient names inherited from the design matrix are unreliable in a partially
nested design: \code{chord_typemaj} names a comparison at one inversion, not an
average over them, and which one depends on the pivot. Each row therefore
carries a \code{meaning} column giving the combination of realized conditions
that the coefficient equals, so no label has to be read in isolation.

Covariate slopes, and any covariate-by-condition interaction, are slopes rather
than means and are reported untranslated. Under a threshold family the reference
condition has no free coefficient, its latent mean being absorbed into the
thresholds; that row reads as exactly zero and is flagged.
}
\value{
A data frame of class \code{nestimand_summary} with one row per coefficient,
giving the estimate, standard error, test statistic, p-value, interval, and the
combination of conditions the coefficient equals. The transformation matrix is
in the \code{nestimand_map} attribute.
}
\seealso{
\code{\link{estimand}} for contrasts under a stated weighting policy,
\code{\link{effect_basis}} for the matrix relating the two spaces,
\code{\link{nest_fit}}.
}
\examples{
\dontrun{
sp <- nesting_spec(dat, response ~ chord_type * inversion + training,
                   "inversion \%in\% chord_type")
m  <- nest_fit(sp)
nest_summary(m)                # effect coefficients, with their meanings
nest_summary(m, space = "cells")   # the realized-cell means
}
}
