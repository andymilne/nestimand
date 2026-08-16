\name{estimand}
\alias{estimand}
\alias{print.nestimand_estimand}
\alias{add_bounds}
\alias{show_code}
\alias{show_code.default}
\alias{show_code.nestimand_estimand}
\title{Estimands in the original variable space}
\description{
\code{estimand} returns a labelled contrast table in the original variable names, under
a stated weighting policy, with partial-identification bounds over all admissible
policies and a reorder self-check. The analysis is assembled as text and then evaluated,
so \code{show_code} returns the code that ran rather than a reconstruction of it.
}
\usage{
estimand(model, target, policy = "equal", at = NULL, contrast = c("pairwise", 
         "reference", "sequential", "within"), route = c("g_computation",
         "observed", "cells"), weights = NULL, scale = c("response",
         "latent"), data = NULL, bounds = TRUE, self_check = TRUE,
         dry_run = FALSE, ..., spec = NULL, .env = parent.frame())

add_bounds(est, model, spec, target, contrast = "pairwise", scale = c("response", 
           "latent"), cells = spec$cells, ...)

show_code(x, ...) 
}
\arguments{
  \item{model}{a fitted model, ideally from \code{\link{nest_fit}}.}
  \item{spec}{a \code{nesting_spec} object, needed only for a model fitted by calling the engine directly: a fit from \code{\link{nest_fit}} carries its own.}
  \item{target}{name of the variable whose levels are contrasted, quoted or bare.}
  \item{policy}{a distribution over the versions of a compound condition: one of \code{"equal"}, \code{"proportional"}, \code{"hierarchical"}, \code{"nominated"}, or a named numeric vector of masses.}
  \item{at}{for \code{policy = "nominated"}, a named vector giving the version, e.g. \code{c(inversion = "0")}.}
  \item{route}{over which rows the model is evaluated. \code{"g_computation"} crosses every observed row with every realized cell, so each condition is averaged over the same covariate distribution, giving the population-averaged effect. \code{"observed"} uses the observed rows as they stand, reweighted so that the versions enter at the policy's proportions, predicting for no row that was never run. \code{"cells"} uses one row per realized cell with covariates at their means: the conditional effect at an average covariate value.}
  \item{weights}{unit weights, as a column name or a numeric vector of one weight per row. Where the policy weights the versions of a condition, these weight the rows, standardizing the estimand to a population other than the sample: survey weights, or post-stratification to a target distribution of the covariates. Not available on the \code{"cells"} route, which evaluates one row per condition and so has no units to weight.}
  \item{contrast}{which comparisons to form. \code{"within"} gives per-stratum contrasts of a nested variable, which cross no boundary and take no policy.}
  \item{scale}{\code{"response"} goes through the prediction machinery; \code{"latent"} computes the contrast as a linear functional of the coefficients, which is exact on a linear scale and necessary for ordinal fits.}
  \item{data}{data over which the estimand is averaged.}
  \item{bounds}{attach the partial-identification region.}
  \item{self_check}{re-run the estimand with nested factor levels permuted, and warn if it moves.}
  \item{dry_run}{return the code without computing anything. \code{show_code()} prints the code of a result that has already been computed, since the inner call is evaluated first; this is the way to see the code beforehand.}
  \item{...}{passed to the estimand engine and written into the saved code.}
  \item{.env}{environment in which the assembled code is evaluated.}
  \item{est}{an estimand table, for \code{add_bounds}.}
  \item{cells}{the realized cells the bounds are computed over, matching any restriction applied to the estimand itself.}
  \item{x}{an object carrying code, for \code{show_code}.}
}
\details{
The policy and the route answer different questions and are chosen separately.
In a linear model with covariates balanced across conditions the three routes
agree; they separate under imbalance, under covariate-by-condition interaction,
and in any nonlinear model. Where the covariate distribution differs across
conditions, \code{"observed"} compares subpopulations rather than the same units
under two interventions, and only \code{"g_computation"} supports a causal
reading.

Contrast direction is fixed by the package rather than inherited from the engine:
comparisons run in declared factor-level order, earlier level minus later, whatever the
installed version returns. Without this, the same analysis could report opposite-signed
effects on two machines.
}
\value{
\code{estimand} an estimand table of class \code{nestimand_estimand}, carrying the
policy, contrast, scale, bounds, self-check result, and code in its \code{nestimand}
attribute. \code{show_code} prints the code and returns it invisibly.
}
\seealso{
\code{\link{nest_fit}}, \code{\link{nest_policy}}, \code{\link{latent_estimand}}
}
\examples{
\dontrun{
sp <- nesting_spec(dat, response ~ chord_type * inversion + training,
                   "inversion %in% chord_type")
m  <- nest_fit(sp)
e  <- estimand(m, chord_type, policy = "equal")   # the fit carries its spec
e
show_code(e)
}
}
