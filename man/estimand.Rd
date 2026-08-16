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
estimand(model, spec, target, policy = "equal", at = NULL,
         contrast = c("pairwise", "reference", "sequential", "within"),
         scale = c("response", "latent"), data = NULL,
         bounds = TRUE, self_check = TRUE, ..., .env = parent.frame())

add_bounds(est, model, spec, target, contrast = "pairwise",
           scale = c("response", "latent"), ...)
show_code(x, ...)
}
\arguments{
  \item{model}{a fitted model, ideally from \code{\link{nest_fit}}.}
  \item{spec}{a \code{nesting_spec} object.}
  \item{target}{name of the variable whose levels are contrasted, quoted or bare.}
  \item{policy}{a distribution over the versions of a compound condition: one of \code{"equal"}, \code{"proportional"}, \code{"counterfactual"}, \code{"hierarchical"}, \code{"nominated"}, or a named numeric vector of masses.}
  \item{at}{for \code{policy = "nominated"}, a named vector giving the version, e.g. \code{c(inversion = "0")}.}
  \item{contrast}{which comparisons to form. \code{"within"} gives per-stratum contrasts of a nested variable, which cross no boundary and take no policy.}
  \item{scale}{\code{"response"} goes through the prediction machinery; \code{"latent"} computes the contrast as a linear functional of the coefficients, which is exact on a linear scale and necessary for ordinal fits.}
  \item{data}{data over which the estimand is averaged.}
  \item{bounds}{attach the partial-identification region.}
  \item{self_check}{re-run the estimand with nested factor levels permuted, and warn if it moves.}
  \item{...}{passed to the estimand engine and written into the saved code.}
  \item{.env}{environment in which the assembled code is evaluated.}
  \item{est}{an estimand table, for \code{add_bounds}.}
  \item{x}{an object carrying code, for \code{show_code}.}
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
e  <- estimand(m, sp, chord_type, policy = "equal")
e
show_code(e)
}
}
