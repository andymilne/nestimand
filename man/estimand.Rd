\name{estimand}
\alias{estimand}
\alias{print.nestimand_estimand}
\alias{add_bounds}
\alias{show_code}
\alias{show_code.default}
\alias{show_code.nestimand_estimand}
\alias{engine_types}
\alias{print.nestimand_types}
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
         "cells"), weights = NULL, type = NULL, subsample = NULL,
         data = NULL, bounds = TRUE, self_check = TRUE, dry_run = FALSE,
         ..., spec = NULL, .env = parent.frame())

add_bounds(est, model, spec, target, contrast = "pairwise", scale = c("response", 
           "latent"), cells = spec$cells, route = c("g_computation",
           "cells"), data = spec$data, ...)

show_code(x, ...) 

engine_types(model = NULL) 
}
\arguments{
  \item{model}{a fitted model, ideally from \code{\link{nest_fit}}.}
  \item{spec}{a \code{nesting_spec} object, needed only for a model fitted by calling the engine directly: a fit from \code{\link{nest_fit}} carries its own.}
  \item{target}{the variable whose levels are contrasted, quoted or bare. Several may be named at once, and the formula operators carry their usual meaning: \code{c(chord_type, inversion)} computes each in turn, \code{chord_type:inversion} gives the interaction alone, and \code{chord_type * inversion} gives both variables and their interaction, as in a model formula.}
  \item{policy}{a distribution over the versions of a compound condition: one of \code{"equal"}, \code{"proportional"}, \code{"hierarchical"}, \code{"nominated"}, or a named numeric vector of masses.}
  \item{at}{for \code{policy = "nominated"}, a named vector giving the version, e.g. \code{c(inversion = "0")}.}
  \item{data}{data over which the bounds are computed, for \code{add_bounds}.}
  \item{subsample}{average over a random sample of this many units rather than all of them, on the \code{"g_computation"} route. It estimates the same quantity, with Monte Carlo error in place of exactness, and is the way to keep a population-averaged estimand affordable on a nonlinear scale, where the averaging cannot be done before the model is evaluated. The rows drawn are written into the emitted code.}
  \item{route}{over which rows the model is evaluated. \code{"g_computation"} crosses every observed row with every realized cell, so each condition is averaged over the same covariate distribution, giving the population-averaged effect. \code{"cells"} uses one row per realized cell with covariates at their means: the conditional effect at an average covariate value.}
  \item{subsample}{average over a random sample of this many units rather than all of them, on the G-computation route. The same quantity, with Monte Carlo error in place of exactness; the seed is written into the emitted code, so the same rows come back. Which conditions the sampled rows happen to carry does not matter, since every row is crossed with every realized cell; and a proportional policy is counted from the whole data rather than from the sample, so the weights are unaffected.}
  \item{weights}{unit weights, as a column name or a numeric vector of one weight per row. Where the policy weights the versions of a condition, these weight the rows, standardizing the estimand to a population other than the sample: survey weights, or post-stratification to a target distribution of the covariates. Not available on the \code{"cells"} route, which evaluates one row per condition and so has no units to weight.}
  \item{contrast}{which comparisons to form. \code{"within"} gives per-stratum contrasts of a nested variable, and \code{"interaction"} a difference of differences over two or more targets; neither crosses a structural boundary, so neither takes a policy.}
  \item{type}{the quantity to contrast. \code{"eta"} is nestimand's own name for
the linear predictor - the mean structure, with thresholds and link excluded -
computed from the coefficients as c'b. It is exact, gives one number per
contrast, is available for every model class including those the prediction
machinery cannot handle, and costs one matrix product rather than a prediction
per row. Every other value is passed to \code{marginaleffects::avg_predictions}
verbatim, and means there exactly what that package means by it: \code{"link"}
is per-category on an ordinal \code{clm} fit and an error on an \code{lm},
because that is how \pkg{marginaleffects} defines it. The default is
\code{"response"}, or \code{"eta"} for an ordinal family. Where an engine
quantity would give the same numbers as \code{"eta"} more slowly, a note says
so; where the response scale is asked for on an ordinal fit, the engine's own
spelling for it is resolved against the model, since it differs between classes
and releases.}
  \item{route}{over which rows the model is evaluated. \code{"g_computation"} crosses every observed row with every realized cell, so each condition is averaged over the same covariate distribution, giving the population-averaged effect. \code{"cells"} uses one row per realized cell with covariates at their means: the conditional effect at an average covariate value.}
  \item{subsample}{average over a random sample of this many units rather than all of them, on the G-computation route. The same quantity, with Monte Carlo error in place of exactness; the seed is written into the emitted code, so the same rows come back. Which conditions the sampled rows happen to carry does not matter, since every row is crossed with every realized cell; and a proportional policy is counted from the whole data rather than from the sample, so the weights are unaffected.}
  \item{weights}{unit weights, as a column name or a numeric vector of one weight per row. Where the policy weights the versions of a condition, these weight the rows, standardizing the estimand to a population other than the sample: survey weights, or post-stratification to a target distribution of the covariates. Not available on the \code{"cells"} route, which evaluates one row per condition and so has no units to weight.}
  \item{contrast}{which comparisons to form. \code{"within"} gives per-stratum contrasts of a nested variable, and \code{"interaction"} a difference of differences over two or more targets; neither crosses a structural boundary, so neither takes a policy.}
  \item{type}{the quantity to contrast, in the engines' own vocabulary. \code{"link"} - equivalently \code{"latent"} or \code{"linear.predictor"} - is the linear predictor, computed from the coefficients as c'b: exact, one number per contrast, and independent of an engine's ordinal support. Where the link is the identity - \code{lm}, \code{lmer}, a gaussian family with the default link - the response scale and the linear predictor are the same quantity, so \code{type = "response"} is computed from the coefficients: the averaging is done on the design matrix rather than one prediction per row, which is the same answer for much less work. This is gated on the link rather than the family, since \code{gaussian(link = "log")} is not the identity, and set aside whenever the engine's own output is needed: a \code{hypothesis} or any other argument meant for it, a within-stratum contrast, or a model whose coefficients are not the cell design's - one fitted by hand in crossed or chain form. Otherwise anything but the linear predictor is passed to \code{marginaleffects::avg_predictions}: \code{"response"} for the expected outcome, \code{"prediction"} for the posterior predictive - which on an ordinal family draws category codes, so that contrasting them treats the rating scale as an interval one, and is noted when asked for - \code{"mean.class"} and the rest where a model offers them. The default is \code{"response"}, or \code{"link"} for an ordinal family, where the response scale gives one number per outcome category rather than one per contrast. For \code{"response"} on an ordinal fit the engine's own spelling is resolved against the model, since it differs between classes and releases. Some classes the prediction machinery does not handle at all - \code{clmm} among them - and there no quantity but the linear predictor is available, which is one reason it is computed from the coefficients rather than by prediction.}
  \item{data}{data over which the estimand is averaged.}
  \item{bounds}{attach the partial-identification region. It is computed on the same route as the estimand, one evaluation per single-version policy, so its cost grows with the number of versions: on a design admitting many, \code{bounds = FALSE} is the way to skip it.}
  \item{self_check}{re-run the estimand with nested factor levels permuted, and warn if it moves. For a mixed or Bayesian fit the check runs on the fixed-effects shadow model - the same cell formula without the random terms - because order instability is a fixed-effects phenomenon, and refitting a large random structure is both expensive and liable to settle on a different optimum, which would report a failure having nothing to do with level order.}
  \item{dry_run}{return the code without computing anything. \code{show_code()} prints the code of a result that has already been computed, since the inner call is evaluated first; this is the way to see the code beforehand.}
  \item{...}{passed to \code{marginaleffects::avg_predictions}, which computes the estimand, and written into the saved code: \code{conf_level}, \code{p_adjust}, \code{vcov}, \code{type}, \code{df} and the rest are supplied this way. See \code{?marginaleffects::avg_predictions} for what it accepts. \code{ndraws} thins the posterior of a \pkg{brms} fit and works on either route: passed to the prediction function, or applied to the draws directly. It is the other lever on a large problem, alongside \code{subsample}, and usually the cheaper, since a posterior has far more draws than the quantities they inform. A \code{hypothesis} supplied here replaces the comparison that \code{contrast} would have formed - the two cannot both be sent - and the direction convention goes with it, the labels then being whatever the engine returns. This is the way to ask for comparisons within each level of a grouped output, as an ordinal fit on the response scale produces: \code{hypothesis = ~ pairwise | group}. When \code{type} names the linear predictor the estimand is computed from the coefficients instead, and \code{conf_level} is the only one of these that applies.}
  \item{.env}{environment in which the assembled code is evaluated.}
  \item{est}{an estimand table, for \code{add_bounds}.}
  \item{cells}{the realized cells the bounds are computed over, matching any restriction applied to the estimand itself.}
  \item{x}{an object carrying code, for \code{show_code}.}
  \item{digits}{decimal places for the printed table. The object keeps full precision; only the display is rounded.}
  \item{scale}{internal, for \code{add_bounds}: whether the bounds are computed from the coefficients or through the prediction function. \code{\link{estimand}} sets it from \code{type}.}
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
\section{What each engine accepts}{
The \code{type} vocabulary is \pkg{marginaleffects}', and differs by model class
and between its releases. \code{\link{engine_types}} reports it for the version
installed; at the time of writing it was:

\tabular{ll}{
  any \tab \code{eta} - nestimand's own: the linear predictor, from the coefficients \cr
  \code{lm} \tab \code{response} \cr
  \code{glm} \tab \code{response}, \code{link}, \code{invlink(link)} \cr
  \code{lmerMod} \tab \code{response} \cr
  \code{glmerMod} \tab \code{response}, \code{link} \cr
  \code{clm} \tab \code{prob}, \code{cum.prob}, \code{linear.predictor} \cr
  \code{clmm} \tab none: the class is unsupported, so no quantity computed by prediction is available \cr
  \code{brmsfit} \tab \code{response}, \code{link}, \code{prediction}, \code{average} \cr
}

Note that \code{linear.predictor} on a \code{clm} is one value per outcome
category, and is not the contrast \code{eta} computes; and that a \code{clmm}
admits no type at all, so \code{eta} is the only quantity available there.
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

In a linear model the two agree, and on the linear predictor they agree exactly whatever the family, since averaging the design rows and averaging the predictions are then the same operation. They separate in any nonlinear model, where
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
