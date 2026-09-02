\name{nest_prior}
\alias{nest_prior}
\alias{print.nestimand_prior}
\alias{prior_audit}
\alias{prior_for_estimand}
\alias{prior_stanvars}
\alias{prior_statement}
\alias{prior_from_brms}
\alias{split_brms_prior}
\alias{parse_normal_prior}
\title{Priors stated in one space and fitted in another}
\description{
A prior is stated where it is thought - on named effects, or on cell means - and fitted
where estimation is well posed. For an elliptical prior the translation is exact and
needs no reparameterization: if the effects carry m ~ N(m0, D) then the cell means carry
mu ~ N(A m0, A D A'), and Stan accepts a correlated prior directly.
}
\usage{
nest_prior(spec, mean, sd = NULL, on = c("effects", "cells"),
           family = c("normal", "student_t"), df = 3, covariate_mean = 0,
           covariate_sd = NULL)

prior_audit(prior, contrasts = NULL) 

prior_for_estimand(prior, target, policy = "equal", at = NULL, contrast = "pairwise") 

prior_stanvars(prior) 

prior_statement(prior) 

prior_from_brms(spec, prior, covariate_mean = 0, covariate_sd = NULL)
}
\arguments{
  \item{spec}{a \code{nesting_spec} object.}
  \item{mean}{prior mean: one value for all, a full vector, or a named vector, optionally with a \code{.default} element.}
  \item{sd}{prior standard deviations, for independent priors.}
  \item{on}{the space in which the prior is stated.}
  \item{family}{\code{"normal"} or \code{"student_t"}.}
  \item{df}{degrees of freedom for the t family.}
  \item{covariate_mean, covariate_sd}{the prior on the covariate coefficients. A \code{brms} prior of class \code{"b"} covers every population-level coefficient, so a prior stated only on the cells would leave Stan with mismatched dimensions.}
  \item{prior}{a \code{nestimand_prior} object; for \code{prior_from_brms}, a \pkg{brms} prior written as it would be for \code{brm()}.}
  \item{txt}{a prior string, for \code{parse_normal_prior}.}
  \item{contrasts}{named list of contrast vectors over the realized cells.}
  \item{target}{name of the variable whose levels are contrasted, quoted or bare.}
  \item{policy}{a distribution over the versions of a compound condition: one of \code{"equal"}, \code{"proportional"}, \code{"hierarchical"}, \code{"nominated"}, or a named numeric vector of masses.}
  \item{at}{for \code{policy = "nominated"}, a named vector giving the version, e.g. \code{c(inversion = "0")}.}
  \item{contrast}{which comparisons to report.}
  \item{x}{an object to print.}
}
\section{Priors written in brms syntax}{
Priors may be stated once, in the usual \pkg{brms} way, and passed to
\code{\link{nest_fit}} as \code{prior}. The \code{class = "b"} entry describes
the mean structure: written in the original variables, as it would be for the
model the analyst has in mind, it is translated into cell space and reaches Stan
through \code{stanvars} rather than as a prior argument. Every other class -
the thresholds, the random-effects standard deviations and correlations, sigma -
is passed through untouched, and the two are combined into one call.

One combination is impossible: a translated prior places a multivariate normal
on the coefficients, and \code{sample_prior = "yes"} asks \pkg{brms} to draw one
prior sample per class into a scalar, which Stan rejects as an ill-typed
assignment. Draw from the prior in a separate run with
\code{sample_prior = "only"}, which works. Only elliptical priors translate
exactly, so
\code{normal} and \code{student_t} are accepted and others refused with the
alternative named; and because the translation is joint, the class must be
described in one statement rather than coefficient by coefficient.
}
\details{
Under a threshold family - the ordinal ones - the fitted coefficients are
contrasts of the cell means against the reference condition rather than the
means themselves, since the thresholds carry the overall level. The translation
is still linear and is applied in the same way, with one qualification: a prior
on the overall level is not transmitted, because no coefficient corresponds to
it. State that on the thresholds separately if it matters.
}
\value{
\code{nest_prior} an object of class \code{nestimand_prior}, holding the prior in
both spaces and over the full coefficient vector. \code{prior_audit} and
\code{prior_for_estimand} return data frames of implied means and standard deviations;
\code{prior_stanvars} and \code{prior_statement} return the \pkg{brms} objects.
}
\seealso{
\code{\link{nest_fit}}, \code{\link{nest_estimand}}
}
\keyword{internal}
