\name{nest_fit}
\alias{nest_fit}
\alias{random_terms}
\alias{random_slope_rank}
\alias{has_thresholds}
\alias{print.nestimand_code}
\title{Fit a declared design, keeping the call}
\description{
\code{nest_fit} assembles the model call as text and evaluates it, so the returned fit
carries the code that produced it and \code{\link{show_code}} can show the whole
pipeline. Arguments in \code{...} reach the engine unaltered.
\code{random_terms} translates a declared random structure.
}
\usage{
nest_fit(spec, formula = NULL, nests = NULL, mode = NULL,
         random_structure = c("cells", "reduced", "chain", "as_declared"),
         data = NULL, dry_run = FALSE, ..., fit = NULL, family = NULL,
         random = NULL, cell_name = "cell", priors = NULL,
         prior_space = c("effects", "cells"), .env = parent.frame())

random_terms(spec, structure = c("cells", "chain", "chain_slope", 
             "as_declared"))

random_slope_rank(spec, data = spec$data) 

has_thresholds(spec) 
}
\arguments{
  \item{spec}{a \code{\link{nesting_spec}}, or the data to build one from. Given data, \code{formula} is required and the remaining declaration arguments are passed to \code{nesting_spec}, so a model can be declared and fitted in one call; the declaration is then part of the emitted code. Given a \code{nesting_spec}, the declaration arguments must be left out, since it already carries them.}
  \item{formula}{the model formula, when the declaration is being made here.}
  \item{nests}{the nesting structure, when the declaration is being made here. May be omitted, in which case it is read from the data - see \code{\link{nesting_spec}}.}
  \item{fit}{the engine, when the declaration is being made here. It belongs to the declaration, so passing it alongside a \code{nesting_spec} is refused rather than allowed to disagree with it.}
  \item{family}{the family, when the declaration is being made here.}
  \item{random}{the random structure, when the declaration is being made here and it is not written into the formula with a bar.}
  \item{cell_name}{the name of the realized-cell factor, when the declaration is being made here.}
  \item{mode}{the parameterization; chosen by \code{\link{fitting_mode}} when \code{NULL}, which is what it should normally be left as. A declaration that crosses the structure fully is fitted on the cell factor, whose coefficients are the cell means; one that asks for less is fitted on the reduced design - one column per identified effect, so the model is the one the formula names, no coefficient is aliased, and nothing has to be held at zero.}
  \item{random_structure}{how a declared random structure is translated. Chosen to match the mean structure when not given: a declaration that crosses the structure fully gets \code{"cells"}, one that asks for less gets \code{"reduced"}. \code{"reduced"} puts the random effects on the same columns the mean structure is fitted on, so they are deviations in exactly the effects declared and the two sides have the same dimension; a random term naming something the mean structure does not contain is refused, since that is a different model and not what this option expresses. \code{"cells"} gives an unstructured covariance over realized cells, \code{"chain"} for the grouping-factor submodel - one rung per ancestor path, so a chain gives \code{(1|p) + (1|p:a) + (1|p:a:b)} and siblings enter crossed, \code{(1|p:a:b) + (1|p:a:c) + (1|p:a:b:c)}, neither of them made the coarser division - \code{"chain_slope"} for the chain counterpart of the cell covariance, \code{"as_declared"} to pass it through unchanged.}
  \item{data}{data to fit; the spec's data by default.}
  \item{priors}{the package's own translated prior, as \code{\link{nest_prior}} returns. Called without a declaration it returns one that is dimensioned here, against the declaration this call builds or was given, so a prior can be written in the same call as the fit. It governs the mean structure; \pkg{brms}'s own \code{prior} remains the route to everything else - the random-effects standard deviations and correlations, the thresholds, sigma - and the two are combined into one call. A \code{class = "b"} entry in \code{prior} is translated into cell space rather than passed through, since those coefficients are the realized conditions rather than the effects of the original variables; supplying one alongside \code{priors} is refused as two priors for the same structure. It sits after \code{...} and must be named in full, so that \pkg{brms}'s own \code{prior} argument reaches the engine rather than being partial-matched to it.}
  \item{dry_run}{return the code without fitting, which matters when the fit is a sampling run.}
  \item{...}{passed to the fitting function named by \code{fit} in the declaration - \code{lm}, \code{glm}, \code{lme4::lmer}, \code{lme4::glmer}, \code{ordinal::clm}, \code{ordinal::clmm} or \code{brms::brm} - and written into the saved code: this is how \code{prior}, \code{sample_prior}, \code{save_pars}, \code{init}, \code{threads}, \code{file}, \code{chains}, \code{iter}, \code{seed} and the rest reach \pkg{brms}, and their counterparts reach the other engines.}
  \item{.env}{environment in which the assembled call is evaluated.}
  \item{structure}{as \code{random_structure}.}
  \item{x}{an object to print.}
  \item{prior_space}{the space a \code{class = "b"} prior in \code{prior} is written in: \code{"effects"}, the original variables, in which case it is translated to the fitted coefficients, or \code{"cells"}, in which case it is left as written. Named in full, since \pkg{brms}'s own \code{prior} would otherwise partial-match it.}
}
\details{
Only terms crossing the structural boundary are translated: a random intercept, a slope
on a fully crossed variable, or a hand-written grouping structure is returned exactly as
declared. A structure varying with a nested variable but not with its parent is refused
rather than enlarged, and the alternatives are named. A covariance-structure wrapper such
as \code{diag()} survives translation.

The \code{"chain"} rungs are ordered by depth and then by name, so a design emits the
same formula however its declarations were written; the \code{"cells"} form does not
depend on the order in any case.

\preformatted{sp <- nesting_spec(dat, rating ~ chord_type * inversion * X1 +
                          (chord_type * inversion * X1 | participant),
                   c(inversion, X1) \%in\% chord_type, fit = "lmer")
random_terms(sp, "cells")
## (0 + cell | participant)
random_terms(sp, "chain")
## (1 | participant) + (1 | participant:chord_type) +
## (1 | participant:chord_type:X1) + (1 | participant:chord_type:inversion) +
## (1 | participant:chord_type:X1:inversion)}
}
\value{
\code{nest_fit} the fitted model, classed \code{nestimand_fit} behind the
engine's own classes - so that packages dispatching or validating on
\code{class(model)[1]} see what they expect - with the code in its \code{nestimand_code} attribute and the
declaration in \code{nestimand_spec}, so that \code{\link{nest_estimand}} and
\code{\link{nest_summary}} can be called on the model alone. An
\code{update} method carries the declaration across a refit, so a simplified
model remains usable without the declaration being supplied again; whether the
updated model still corresponds to it is checked where that matters. The extra
class changes no behaviour: the engine's own \code{predict}, \code{summary},
\code{vcov}, and \pkg{marginaleffects} methods apply as before. S4 fits are not reclassified; a \code{brmsfit} is an ordinary
S3 object and is, but since \pkg{brms} supplies its own \code{update} method
that one takes precedence, so an updated brms fit needs \code{spec} supplied, or an object of class \code{nestimand_code} when \code{dry_run = TRUE}.
\code{random_terms} a character string; \code{random_slope_rank} a list reporting the
columns, rank, and identically zero columns of a structural random slope.
}
\seealso{
\code{\link{nesting_spec}}, \code{\link{nest_estimand}}, \code{\link{nest_summary}}
}
