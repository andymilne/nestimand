\name{nest_policy}
\alias{nest_policy}
\alias{print.nestimand_policy}
\alias{counterfactual_grid}
\alias{policy_weights}
\alias{policy_vertices}
\alias{boundary_vars}
\title{Weighting policies for across-boundary contrasts}
\description{
A contrast that crosses the structural boundary compares a condition with a single
version against one with several, and is undefined until a distribution over those
versions is stated. \code{nest_policy} builds that distribution, and
\code{counterfactual_grid} attaches it to the G-computation grid as row weights. Each
policy names a stochastic intervention, so the estimand is a causal contrast between
two well-defined policies.
}
\usage{
nest_policy(spec, target, policy = "equal", at = NULL, data = spec$data, 
            cells = spec$cells)

counterfactual_grid(spec, data = spec$data, policy = NULL, cells = spec$cells) 

policy_weights(spec, grid, policy) 

policy_vertices(spec, target) 

boundary_vars(spec) 
}
\arguments{
  \item{spec}{a \code{nesting_spec} object.}
  \item{target}{name of the variable whose levels are contrasted, quoted or bare.}
  \item{policy}{a distribution over the versions of a compound condition: one of \code{"equal"}, \code{"proportional"}, \code{"hierarchical"}, \code{"nominated"}, or a named numeric vector of masses. How the versions are weighted; where the model is evaluated is the separate \code{route} argument of \code{\link{estimand}}.}
  \item{at}{for \code{policy = "nominated"}, a named vector giving the version, e.g. \code{c(inversion = "0")}.}
  \item{data}{the data over which empirical version frequencies are taken.}
  \item{grid}{a grid carrying the original factors.}
  \item{cells}{the realized cells to work over, when they have been restricted: a marginal contrast of a nested variable is pooled only over the strata in which it varies.}
  \item{x}{an object to print.}
}
\details{
A stratum with only one realized version has no choice to make, so the policy is
vacuous there and its single condition takes all the mass. A supplied distribution that
is silent on a realized version is refused rather than renormalized, since restricting
and renormalizing is a different intervention.
}
\value{
\code{nest_policy} an object of class \code{nestimand_policy}, holding the mass on
each realized version within each stratum; \code{counterfactual_grid} a data frame with
the weights in \code{.w} and the originating row in \code{.row};
\code{policy_vertices} the single-version policies whose contrast range is the
partial-identification region; \code{boundary_vars} the nested variables.
}
\seealso{
\code{\link{estimand}}, \code{\link{add_bounds}}, \code{\link{prior_for_estimand}}
}
