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
latent_estimand(model, target, policy = "equal", at = NULL, contrast = "pairwise", 
                data = NULL, conf_level = 0.95, spec = NULL, cells = NULL)

latent_draws(model, target, policy = "equal", at = NULL, contrast = "pairwise", 
             data = NULL, spec = NULL, cells = NULL)

policy_contrast_matrix(spec, target, policy, data = spec$data, model = NULL, 
                       cells = spec$cells)
}
\arguments{
  \item{model}{a fitted model.}
  \item{spec}{a \code{nesting_spec} object, needed only for a model fitted outside \code{\link{nest_fit}}.}
  \item{target}{name of the variable whose levels are contrasted, quoted or bare.}
  \item{policy}{a distribution over the versions of a compound condition: one of \code{"equal"}, \code{"proportional"}, \code{"hierarchical"}, \code{"nominated"}, or a named numeric vector of masses.}
  \item{at}{for \code{policy = "nominated"}, a named vector giving the version, e.g. \code{c(inversion = "0")}.}
  \item{contrast}{\code{"pairwise"}, \code{"reference"}, \code{"sequential"}, or \code{"interaction"}, the last taking two or more targets and giving a difference of differences among cells.}
  \item{data}{data over which the grid is formed.}
  \item{conf_level}{confidence level for the delta-method interval.}
  \item{cells}{the realized cells to work over. Chosen automatically: a marginal contrast of a nested variable is pooled only over the strata in which it varies, and an interaction over those in which every target varies.}
}
\details{
This route never calls the prediction machinery, so it is immune to the per-category
expansion that ordinal fits provoke and to the guard rails that follow from it. The
contrast vectors it returns are also what the heterogeneity of an effect requires:
c'Sigma c over the random-effects covariance.
}
\section{Frequentist and Bayesian summaries}{
The contrast is the same object either way - \code{C b}, with \code{C} one row
per comparison - and only the summary differs. For a frequentist fit the
coefficients and their covariance give a point estimate, a delta-method interval,
a test statistic and a p-value. For a \pkg{brms} fit the same \code{C} is
applied draw by draw, and the summary is the posterior of the contrast: its mean,
its standard deviation, a quantile interval, and \code{pd}, the probability of
direction - the larger of the posterior probabilities that the contrast is
positive and that it is negative, so at least a half. No test statistic or
p-value is reported there, a Wald test computed from a posterior covariance
having no Bayesian reading.
}
\value{
\code{latent_estimand} a data frame of contrasts: with delta-method intervals, the
contrast vectors in its \code{nestimand_cvecs} attribute; \code{latent_draws} a data
frame of posterior draws, one column per contrast; \code{policy_contrast_matrix} the
weighted design-matrix rows, one per stratum.
}
\seealso{
\code{\link{estimand}}, \code{\link{nest_policy}}
}
