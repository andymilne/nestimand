\name{nest_estimand}
\alias{nest_estimand}
\alias{print.nestimand_estimand}
\alias{add_bounds}
\title{Estimands in the original variable space}
\description{
\code{nest_estimand} returns a labelled contrast table in the original variable names, under
a stated weighting policy, with partial-identification bounds over all admissible
policies and a reorder self-check. The analysis is assembled as text and then evaluated,
so \code{show_code} returns the code that ran rather than a reconstruction of it.
}
\usage{
nest_estimand(model, target, policy = "equal", at = NULL, contrast = c("pairwise", 
         "reference", "sequential"), route = c("g_computation",
         "cells"), weights = NULL, type = NULL, subsample = NULL,
         data = NULL, bounds = TRUE, self_check = TRUE, dry_run = FALSE,
         by = NULL, ..., .env = parent.frame(), .restrict = NULL)

add_bounds(est, model, spec, target, contrast = "pairwise", scale = c("response", 
           "latent"), cells = spec$cells, route = c("g_computation",
           "cells"), data = spec$data, ...)


}
\arguments{
  \item{model}{a fitted model, ideally from \code{\link{nest_fit}}.}
  \item{spec}{a \code{nesting_spec} object, needed only for a model fitted by calling the engine directly: a fit from \code{\link{nest_fit}} carries its own.}
  \item{target}{the variable whose levels are contrasted, quoted or bare. Several may be named at once, and the formula operators carry their usual meaning: \code{c(chord_type, inversion)} computes each in turn, \code{chord_type:inversion} gives the interaction alone, and \code{chord_type * inversion} gives both variables and their interaction, as in a model formula. The target is expanded as a formula expands its right-hand side, by the formula machinery itself, so any combination of operators and brackets means what it would mean there: \code{chord_type * inversion * top_note} returns the three variables, the three two-variable interactions and the three-variable one, and \code{chord_type * (inversion + top_note)} returns the three variables and the two interactions with \code{chord_type}, but not \code{inversion:top_note}. The results are named for the terms, in the order a formula writes them. An interaction naming only some of the design variables is formed on conditions averaged over the rest.}
  \item{by}{design variables whose levels the estimand is computed within, quoted or bare, as \code{by} groups an average in \pkg{marginaleffects}: the target's contrasts are formed inside each group, with the policy weighted over that group's conditions alone, and the groups are returned stacked and labelled by the grouping variables. It is the only way to group an estimand, and it groups by anything, so for a nested target it is its parent that is named. Each group's contrast is standardized to the whole sample, so the groups are comparable with each other. A group in which the target does not vary is left out and said. Both routes and every scale accept it.}
  \item{.restrict}{internal: the cells one group of a \code{by} call is computed over.}
  \item{policy}{a distribution over the versions of a compound condition: one of \code{"equal"}, \code{"proportional"}, \code{"hierarchical"}, \code{"nominated"}, or a named numeric vector of masses.}
  \item{at}{for \code{policy = "nominated"}, a named vector giving the version: one level for every variable nested below the target, so \code{c(inversion = "0")} at depth one and \code{c(inversion = "0", X1 = "a", Z = "z1")} where the target has three variables below it.}
  \item{data}{data over which the estimand is averaged; for \code{add_bounds}, the data the bounds are computed over.}
  \item{subsample}{average over a random sample of this many units rather than all of them, on the G-computation route. The same quantity, with Monte Carlo error in place of exactness; the seed is written into the emitted code, so the same rows come back. Which conditions the sampled rows happen to carry does not matter, since every row is crossed with every realized cell; and a proportional policy is counted from the whole data rather than from the sample, so the weights are unaffected.}
  \item{type}{the quantity to contrast. \code{"eta"} is nestimand's own name for
the linear predictor - the mean structure, with thresholds and link excluded -
computed from the coefficients as c'b: the group deviations held at zero, which
describes the average group rather than any sampled one. It is
exact, gives one number per contrast, is available for every model class including those the prediction
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
  \item{weights}{unit weights, as a column name or a numeric vector of one weight per row. Where the policy weights the versions of a condition, these weight the rows, standardizing the estimand to a population other than the sample: survey weights, or post-stratification to a target distribution of the covariates. Not available on the \code{"cells"} route, which evaluates one row per condition and so has no units to weight.}
  \item{contrast}{which comparisons to form among the target's levels: all pairs, each against the first, or each against the one before. That is all it decides. An interaction - a difference of differences - is asked for by the shape of the target, \code{nest_estimand(m, a:b)}, which crosses no structural boundary and so takes no policy; and the levels an estimand is computed within are named by \code{by}.}
  \item{bounds}{attach the partial-identification region. It is computed on the same route as the estimand, one evaluation per single-version policy, so its cost grows with the number of versions: on a design admitting many, \code{bounds = FALSE} is the way to skip it.}
  \item{self_check}{re-run the estimand with nested factor levels permuted, and warn if it moves. For a mixed or Bayesian fit the check runs on the fixed-effects shadow model - the same cell formula without the random terms - because order instability is a fixed-effects phenomenon, and refitting a large random structure is both expensive and liable to settle on a different optimum, which would report a failure having nothing to do with level order.}
  \item{dry_run}{return the code without computing anything. \code{show_code()} prints the code of a result that has already been computed, since the inner call is evaluated first; this is the way to see the code beforehand.}
  \item{...}{passed to \code{marginaleffects::avg_predictions}, which computes the estimand, and written into the saved code: \code{conf_level}, \code{p_adjust}, \code{vcov}, \code{type}, \code{df} and the rest are supplied this way. See \code{?marginaleffects::avg_predictions} for what it accepts. \code{ndraws} thins the posterior of a \pkg{brms} fit and works on either route: passed to the prediction function, or applied to the draws directly. It is the other lever on a large problem, alongside \code{subsample}, and usually the cheaper, since a posterior has far more draws than the quantities they inform. A \code{hypothesis} supplied here replaces the comparison that \code{contrast} would have formed - the two cannot both be sent - and the direction convention goes with it, the labels then being whatever the engine returns. This is the way to ask for comparisons within each level of a grouped output, as an ordinal fit on the response scale produces: \code{hypothesis = ~ pairwise | group}. When \code{type} names the linear predictor the estimand is computed from the coefficients instead, and \code{conf_level} is the only one of these that applies.}
  \item{.env}{environment in which the assembled code is evaluated.}
  \item{x}{an estimand, for printing.}
  \item{est}{an estimand table, for \code{add_bounds}.}
  \item{cells}{the realized cells the bounds are computed over, matching any restriction applied to the estimand itself.}
  \item{digits}{decimal places for the printed table. The object keeps full precision; only the display is rounded.}
  \item{scale}{internal, for \code{add_bounds}: whether the bounds are computed from the coefficients or through the prediction function. \code{\link{nest_estimand}} sets it from \code{type}.}
}
\section{Mixed models}{
Whether the group-level effects enter the predictions is the engine's default to
make: \code{re.form} and \code{re_formula} are passed through untouched, like
every other engine argument. What the package does is say which default is in
force when none is given, since it settles whether the estimand describes a
typical group or the ones in the data - \pkg{brms} and \pkg{lme4} both include
them unless told otherwise. \code{NA} gives the typical group. The spelling
differs between engines and neither accepts the other's, so the note uses the
one that fits the model. \code{"eta"} evaluates the contrast with the group deviations at zero - the
average group, as a coefficient is read with the other terms at zero - which
corresponds to \code{NA} rather than \code{NULL}. On a linear scale the two give
the same estimate, since the deviations average to zero; under a link they do
not, and that is when the choice matters.

The spread is where the choice bites, and that depends on the engine. On a
\pkg{brms} fit the draws carry the group-level parameters, so \code{NULL} gives
a wider posterior than \code{NA}. On a frequentist mixed fit
\pkg{marginaleffects} forms the interval from the fixed-effect covariance
either way, so the sampled-group spread is not available on any route there.
\code{"eta"} is the average-group quantity by default. On a \pkg{brms} fit
\code{re_formula = NULL} is accepted there too: the group-level draws are added
before the contrast is taken, so the estimate is unchanged and the posterior
widens, at the speed of the linear route. On any other engine there is nothing
to average over and it is refused rather than ignored. A \pkg{brms} fit carries the group-level parameters in its
draws, so there the uncertainty does follow the choice. Section 6.6 of the
companion document sets out what the three answers mean.
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
  \code{"hierarchical"} \tab equal at each node of the declared structure, not over the leaves. Each variable's split is made conditional on its own ancestors, so variables sharing a parent are independent choices whose probabilities multiply, and the weighting does not depend on the order the declarations were written. Over a set of siblings it therefore coincides with \code{"equal"}; it separates from it where one variable is nested inside another, since the number of levels to choose among then depends on the parent's value \cr
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

On the linear predictor the two routes give the same matrix whenever every
column of the design is a linear function of the variables held in the grid -
the ordinary case, covering a covariate entering additively or interacted with
the conditions. Averaging those columns over the observed rows and evaluating
them at the covariate means then coincide exactly, so \code{route} costs nothing
there and \code{"cells"} is the cheaper way to the same number. They separate
only where a design column is a nonlinear function of a covariate, since the
average of a function is not the function of the average.

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
\code{nest_estimand} an estimand table of class \code{nestimand_estimand}, carrying the
policy, contrast, scale, bounds, self-check result, and code in its \code{nestimand}
attribute; \code{\link{show_code}} prints it.
}
\seealso{
\code{\link{nest_fit}}, \code{\link{nest_policy}}, \code{\link{latent_estimand}}
}
\examples{
\dontrun{
m <- nest_fit("lm", response ~ chord_type * inversion + training, dat,
              nests = "inversion %in% chord_type")
e <- nest_estimand(m, chord_type, policy = "equal")    # the fit carries its declaration
e
show_code(e)

## A branching, three-deep structure. The target decides which strata the
## contrast is pooled over: X1 within each chord type, Z within each chord
## type and inversion.
m2 <- nest_fit("lm", response ~ chord_type * inversion * X1 * Z, dat,
               nests = c("inversion %in% chord_type", "X1 %in% chord_type",
                         "Z %in% inversion"))
nest_estimand(m2, X1, by = "chord_type")               # one contrast per chord type
nest_estimand(m2, top_note, by = chord_type)           # the same, strata named, for a
                                                  # variable that has no ancestors
nest_estimand(m2, top_note, by = c(chord_type, inversion))
nest_estimand(m2, Z,  by = c("chord_type", "inversion"))  # per chord type and inversion
nest_estimand(m2, Z,  policy = "equal")                # pooled, and so needs a policy
nest_estimand(m2, chord_type, policy = "nominated",
         at = c(inversion = "0", X1 = "a", Z = "z1"))
}
}
