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
nest_fit(spec, mode = NULL, random_structure = c("cells", "chain", 
         "as_declared"), data = NULL, engine = "marginaleffects",
         dry_run = FALSE, ..., priors = NULL, prior_space = c("effects",
         "cells"), .env = parent.frame())

random_terms(spec, structure = c("cells", "chain", "chain_slope", 
             "as_declared"))

random_slope_rank(spec, data = spec$data) 

has_thresholds(spec) 
}
\arguments{
  \item{spec}{a \code{nesting_spec} object.}
  \item{mode}{the parameterization; chosen by \code{\link{fitting_mode}} when \code{NULL}.}
  \item{random_structure}{how a declared random structure is translated: \code{"cells"} for an unstructured covariance over realized cells, \code{"chain"} for the grouping-factor submodel, \code{"chain_slope"} for the chain counterpart of the cell covariance, \code{"as_declared"} to pass it through unchanged.}
  \item{data}{data to fit; the spec's data by default.}
  \item{engine}{the estimand engine the fit is intended for.}
  \item{priors}{a \code{\link{nest_prior}} or \code{\link{chain_priors}} object - a nestimand prior, translated or derived. It governs the mean structure; \pkg{brms}'s own \code{prior} remains the route to everything else - the random-effects standard deviations and correlations, the thresholds, sigma - and the two are combined into one call. A \code{class = "b"} entry in \code{prior} draws a note, since those coefficients are the realized conditions rather than the effects of the original variables, and supplying one alongside \code{priors} is refused as two priors for the same structure. It sits after \code{...} and must be named in full, so that \pkg{brms}'s own \code{prior} argument reaches the engine rather than being partial-matched to it.}
  \item{dry_run}{return the code without fitting, which matters when the fit is a sampling run.}
  \item{...}{passed to the engine and written into the saved code: this is how \code{prior}, \code{sample_prior}, \code{save_pars}, \code{init}, \code{threads}, \code{file}, \code{chains}, \code{iter}, \code{seed} and the rest reach \pkg{brms}, and their counterparts reach the other engines.}
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
}
\value{
\code{nest_fit} the fitted model, classed \code{nestimand_fit} behind the
engine's own classes - so that packages dispatching or validating on
\code{class(model)[1]} see what they expect - with the code in its \code{nestimand_code} attribute and the
declaration in \code{nestimand_spec}, so that \code{\link{estimand}} and
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
\code{\link{nesting_spec}}, \code{\link{estimand}}, \code{\link{chain_priors}}
}
