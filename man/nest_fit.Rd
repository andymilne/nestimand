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
nest_fit(spec, mode = NULL,
         random_structure = c("cells", "chain", "chain_slope", "as_declared"),
         data = NULL, engine = "marginaleffects", priors = NULL,
         dry_run = FALSE, ..., .env = parent.frame())

random_terms(spec, structure = c("cells", "chain", "chain_slope", "as_declared"))
random_slope_rank(spec, data = spec$data)
has_thresholds(spec)
}
\arguments{
  \item{spec}{a \code{nesting_spec} object.}
  \item{mode}{the parameterization; chosen by \code{\link{fitting_mode}} when \code{NULL}.}
  \item{random_structure}{how a declared random structure is translated: \code{"cells"} for an unstructured covariance over realized cells, \code{"chain"} for the grouping-factor submodel, \code{"chain_slope"} for the chain counterpart of the cell covariance, \code{"as_declared"} to pass it through unchanged.}
  \item{data}{data to fit; the spec's data by default.}
  \item{engine}{the estimand engine the fit is intended for.}
  \item{priors}{a \code{\link{nest_prior}} or \code{\link{chain_priors}} object.}
  \item{dry_run}{return the code without fitting, which matters when the fit is a sampling run.}
  \item{...}{passed to the engine and written into the saved code.}
  \item{.env}{environment in which the assembled call is evaluated.}
  \item{structure}{as \code{random_structure}.}
}
\details{
Only terms crossing the structural boundary are translated: a random intercept, a slope
on a fully crossed variable, or a hand-written grouping structure is returned exactly as
declared. A structure varying with a nested variable but not with its parent is refused
rather than enlarged, and the alternatives are named. A covariance-structure wrapper such
as \code{diag()} survives translation.
}
\value{
\code{nest_fit} the fitted model, classed \code{nestimand_fit} ahead of the
engine's own class, with the code in its \code{nestimand_code} attribute and the
declaration in \code{nestimand_spec}, so that \code{\link{estimand}} and
\code{\link{nest_summary}} can be called on the model alone. An
\code{update} method carries the declaration across a refit, so a simplified
model remains usable without the declaration being supplied again; whether the
updated model still corresponds to it is checked where that matters. The extra
class changes no behaviour: the engine's own \code{predict}, \code{summary},
\code{vcov}, and \pkg{marginaleffects} methods apply as before. S4 fits, brms
among them, are not reclassified, so an updated brms fit needs \code{spec}
supplied, or an object of class \code{nestimand_code} when \code{dry_run = TRUE}.
\code{random_terms} a character string; \code{random_slope_rank} a list reporting the
columns, rank, and identically zero columns of a structural random slope.
}
\seealso{
\code{\link{nesting_spec}}, \code{\link{estimand}}, \code{\link{chain_priors}}
}
