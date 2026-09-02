\name{random_covariance}
\alias{random_covariance}
\alias{random_heterogeneity}
\alias{print.nestimand_random}
\title{The random-effects covariance, in the original terms}
\description{
A cell random effect carries one deviation per realized condition, so its
covariance is already in the design's own terms - but labelled by the fitted
factor rather than by the conditions, and printed by the engines in a form that
wraps or truncates once the structure is large. \code{random_covariance} reads it
out in full, in the original labels, and translates it into effect space where
that is meaningful. \code{random_heterogeneity} gives the spread of a named
effect across units.
}
\usage{
random_covariance(model, spec = NULL, space = c("cells", "effects")) 

random_heterogeneity(model, target, policy = "equal", at = NULL, contrast = "pairwise", 
                     spec = NULL)
}
\arguments{
  \item{model}{a mixed fit, from \code{\link{nest_fit}} or otherwise.}
  \item{spec}{a \code{nesting_spec}, needed only if the fit does not carry one.}
  \item{space}{\code{"cells"} keeps one dimension per realized condition; \code{"effects"} translates the covariance to the identified effects, as A^-1 Sigma A^-T.}
  \item{target}{the variable whose contrasts are of interest.}
  \item{policy}{the weighting policy defining those contrasts.}
  \item{at}{the nominated version, for \code{policy = "nominated"}.}
  \item{contrast}{which comparisons to form.}
  \item{x}{an object to print.}
  \item{digits}{digits for the printed standard deviations and correlations.}
}
\details{
Under a grouping-factor structure, heterogeneity in an effect is a named
variance component that can be read from the summary. Under a cell random
effect there is no such component: the covariance describes variation across
conditions, and the spread of an \emph{effect} is the quadratic form c'Sigma c,
where c is the same weight vector that defines the corresponding fixed-effect
contrast. It therefore depends on the weighting policy, and changes with every
contrast reported. \code{random_heterogeneity} builds c from the same policy the
estimand uses, so the two cannot disagree.

A random term that is not a covariance over the realized conditions - a
grouping-factor variance, or a slope on the nesting variable - has nothing to
translate, and is printed as it stands. The \code{translated} attribute of each
matrix records which case it was, for code that needs to know.
}
\value{
\code{random_covariance} a list of covariance matrices, one per grouping factor,
of class \code{nestimand_random}. \code{random_heterogeneity} a data frame of
standard deviations, one row per grouping factor and contrast.
}
\seealso{
\code{\link{nest_summary}}, which shows the covariance alongside the fixed
effects when \code{random = TRUE}; \code{\link{latent_estimand}} for the
contrast vectors themselves.
}
\keyword{internal}
