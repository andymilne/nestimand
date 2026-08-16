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
         "reference", "sequential", "within", "interaction"), route = c("g_computation",
         "cells"), weights = NULL, scale = c("response", "latent"),
         data = NULL, bounds = TRUE, self_check = TRUE, dry_run = FALSE,
         ..., spec = NULL, .env = parent.frame())

add_bounds(est, model, spec, target, contrast = "pairwise", scale = c("response", 
           "latent"), cells = spec$cells, route = c("g_computation",
           "cells"), data = spec$data, ...)

show_code(x, ...) 
}
\arguments{
  \item{model}{a fitted model, ideally from \code{\link{nest_fit}}.}
  \item{spec}{a \code{nesting_spec} object, needed only for a model fitted by calling the engine directly: a fit from \code{\link{nest_fit}} carries its own.}
  \item{target}{the variable whose levels are contrasted, quoted or bare. Several may be named at once, and the formula operators carry their usual meaning: \code{c(chord_type, inversion)} computes each in turn, \code{chord_type:inversion} gives the interaction alone, and \code{chord_type * inversion} gives both variables and their interaction, as in a model formula.}
  \item{policy}{a distribution over the versions of a compound condition: one of \code{"equal"}, \code{"proportional"}, \code{"hierarchical"}, \code{"nominated"}, or a named numeric vector of masses.}
  \item{at}{for \code{policy = "nominated"}, a named vector giving the version, e.g. \code{c(inversion = "0")}.}
  \item{data}{data over which the bounds are computed, for \code{add_bounds}.}
  \item{route}{over which rows the model is evaluated. \code{"g_computation"} crosses every observed row with every realized cell, so each condition is averaged over the same covariate distribution, giving the population-averaged effect. \code{"cells"} uses one row per realized cell with covariates at their means: the conditional effect at an average covariate value.}
  \item{weights}{unit weights, as a column name or a numeric vector of one weight per row. Where the policy weights the versions of a condition, these weight the rows, standardizing the estimand to a population other than the sample: survey weights, or post-stratification to a target distribution of the covariates. Not available on the \code{"cells"} route, which evaluates one row per condition and so has no units to weight.}
  \item{contrast}{which comparisons to form. \code{"within"} gives per-stratum contrasts of a nested variable, and \code{"interaction"} a difference of differences over two or more targets; neither crosses a structural boundary, so neither takes a policy.}
  \item{scale}{\code{"response"} goes through the prediction machinery; \code{"latent"} computes the contrast as a linear functional of the coefficients, and uses no prediction at all. For an ordinal family the two differ in kind - the latent scale gives one number per contrast, the response scale one per outcome category - so the latent scale is the default there, and the response scale elsewhere. Where the response scale is used, the engine's own name for it is supplied automatically, since it differs between model classes and between releases.

Note that the latent scale is not the engine's \code{type = "link"}. For a cumulative model the linear predictor is one value per threshold, so \code{"link"} returns per-category output whose contrasts still carry the thresholds; the latent scale is the difference in the mean structure alone, thresholds cancelled.}
  \item{data}{data over which the estimand is averaged.}
  \item{bounds}{attach the partial-identification region. It is computed on the same route as the estimand, one evaluation per single-version policy, so its cost grows with the number of versions: on a design admitting many, \code{bounds = FALSE} is the way to skip it.}
  \item{self_check}{re-run the estimand with nested factor levels permuted, and warn if it moves. For a mixed or Bayesian fit the check runs on the fixed-effects shadow model - the same cell formula without the random terms - because order instability is a fixed-effects phenomenon, and refitting a large random structure is both expensive and liable to settle on a different optimum, which would report a failure having nothing to do with level order.}
  \item{dry_run}{return the code without computing anything. \code{show_code()} prints the code of a result that has already been computed, since the inner call is evaluated first; this is the way to see the code beforehand.}
  \item{...}{passed to the estimand engine and written into the saved code.}
  \item{.env}{environment in which the assembled code is evaluated.}
  \item{est}{an estimand table, for \code{add_bounds}.}
  \item{cells}{the realized cells the bounds are computed over, matching any restriction applied to the estimand itself.}
  \item{x}{an object carrying code, for \code{show_code}.}
  \item{digits}{decimal places for the printed table. The object keeps full precision; only the display is rounded.}
}
\section{Mixed models}{
An estimand here is population-level: it describes a typical group rather than
any sampled one, so the random effects are excluded from the predictions. That
is stated in the emitted code rather than left to a default, and in the spelling
the engine expects - \code{re.form = NA} for the frequentist engines,
\code{re_formula = NA} for \pkg{brms} - since each silently ignores the other's.
Supplying either yourself overrides it. Section 6.6 of the companion document
sets out the alternatives, which average over the sampled groups or integrate
over the fitted distribution instead.
}
\section{Policy and route}{
Two separate choices determine what an across-boundary contrast means, and both
are arguments rather than defaults, because no correct default exists.

\strong{The policy} says how the versions of a compound condition are weighted.
Where a nested variable is undefined, a level of the nesting variable is a
single condition; where it is defined, that level is several conditions at once.
A contrast between the two therefore compares one condition with a set of them,
and the set has no single value until its members are given weights. Each choice
of weights names a different intervention, and so a different estimand.

\tabular{ll}{
  \code{"equal"} \tab each realized version counts equally \cr
  \code{"proportional"} \tab versions counted at their frequency in the data \cr
  \code{"hierarchical"} \tab equal at each level of a nesting chain, not over the leaves; differs from \code{"equal"} only at depth two and beyond \cr
  \code{"nominated"} \tab all weight on one version, named in \code{at} \cr
  a named vector \tab weights supplied from outside the design, such as corpus frequencies \cr
}

\strong{The route} says over which rows the model is evaluated. This is not a
weighting but a question of which population the estimand describes.

\tabular{ll}{
  \code{"g_computation"} \tab every observed row crossed with every realized cell, so each condition is averaged over one common covariate distribution \cr
  \code{"cells"} \tab one row per realized cell, covariates at their means \cr
}

In a linear model the two agree. They separate in any nonlinear model, where
averaging predictions over a covariate distribution is not the same as
evaluating at its mean: \code{"g_computation"} gives the population-averaged
effect, \code{"cells"} the conditional effect at an average covariate value.

A third route, averaging the fitted values over each condition's own rows, was
offered and withdrawn. Under the cell parameterization the fixed part is
saturated over the cells, so that average returns the weighted raw cell means
exactly: the covariate adjustment cancels, because each condition's average
restores its own covariate mean. It summarizes the data rather than estimating
anything from the model, and \code{data =} with a subset covers the cases where
a summary of particular rows is wanted.

An interaction contrast is a different kind of quantity: every cell it uses
exists, so no boundary is crossed, no policy applies, and there are no bounds to
report. A policy supplied alongside one is dropped rather than refused, with a
message saying so, since \code{a * b} asks for marginal contrasts as well, and
there the policy does apply. It is computed only over the strata in which every target varies, so the
degenerate stratum drops out by construction.

\code{weights} is a third and independent choice: it weights the rows rather
than the versions, standardizing the estimand to a population other than the
sample.

Whatever the policy, the printed output carries the range of the estimand across
all admissible policies, so a reader can see what the choice was worth.
}
\details{
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
