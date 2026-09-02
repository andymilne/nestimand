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
                data = NULL, conf_level = 0.95, spec = NULL, cells = NULL,
                ndraws = NULL, weights = NULL, route = "g_computation", re_formula = NA)

latent_draws(model, target, policy = "equal", at = NULL, contrast = "pairwise", 
             data = NULL, spec = NULL, cells = NULL)

policy_contrast_matrix(spec, target, policy, data = spec$data, model = NULL, 
                       cells = spec$cells, weights = NULL, route = "g_computation")
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
  \item{ndraws}{thin the posterior of a \pkg{brms} fit to this many draws. The engine's own spelling, so the same argument serves on either route.}
  \item{weights}{unit weights, as a column name or one value per row. They multiply the row weights, standardizing the estimand to another population; the contrast is a weighted average of design rows either way.}
  \item{route}{over which rows the design is averaged, as in \code{\link{nest_estimand}}. The two coincide unless a design column is a nonlinear function of a covariate.}
  \item{re_formula}{\code{NA}, the default, evaluates the contrast at the average group. \code{NULL} averages over the groups in the data instead, by adding their draw-wise mean deviation before the contrast is taken; available for a \pkg{brms} fit, where the group-level parameters are in the draws.}
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
\code{marginaleffects::avg_predictions}, which \code{\link{nest_estimand}} uses on
the response scale, where this route uses the coefficients directly.

\code{\link{nest_estimand}}, \code{\link{nest_policy}}
}
