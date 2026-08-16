\name{latent_estimand}
\alias{latent_estimand}
\alias{latent_draws}
\alias{policy_contrast_matrix}
\title{Estimands as linear functionals of the coefficients}
\description{
On a linear scale - the latent scale of an ordinal model, the link scale of a
generalized linear model - a policy contrast is exactly c'b, where c is a weighted
difference of design-matrix rows formed on the counterfactual grid. Because the scale is
linear, averaging predictions over that grid and averaging its rows are the same
operation, so this is G-computation rather than an approximation to it.
}
\usage{
latent_estimand(model, spec, target, policy = "equal", at = NULL,
                contrast = "pairwise", data = spec$data, conf_level = 0.95)

latent_draws(model, spec, target, policy = "equal", at = NULL,
             contrast = "pairwise", data = spec$data)

policy_contrast_matrix(spec, target, policy, data = spec$data, model = NULL)
}
\arguments{
  \item{model}{a fitted model.}
  \item{spec}{a \code{nesting_spec} object.}
  \item{target}{name of the variable whose levels are contrasted, quoted or bare.}
  \item{policy}{a distribution over the versions of a compound condition: one of \code{"equal"}, \code{"proportional"}, \code{"counterfactual"}, \code{"hierarchical"}, \code{"nominated"}, or a named numeric vector of masses.}
  \item{at}{for \code{policy = "nominated"}, a named vector giving the version, e.g. \code{c(inversion = "0")}.}
  \item{contrast}{\code{"pairwise"}, \code{"reference"}, or \code{"sequential"}.}
  \item{data}{data over which the grid is formed.}
  \item{conf_level}{confidence level for the delta-method interval.}
}
\details{
This route never calls the prediction machinery, so it is immune to the per-category
expansion that ordinal fits provoke and to the guard rails that follow from it. The
contrast vectors it returns are also what the heterogeneity of an effect requires:
c'Sigma c over the random-effects covariance.
}
\value{
\code{latent_estimand} a data frame of contrasts with delta-method intervals, the
contrast vectors in its \code{nestimand_cvecs} attribute; \code{latent_draws} a data
frame of posterior draws, one column per contrast; \code{policy_contrast_matrix} the
weighted design-matrix rows, one per stratum.
}
\seealso{
\code{\link{estimand}}, \code{\link{nest_policy}}
}
