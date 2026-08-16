\name{nesting_spec}
\alias{nesting_spec}
\alias{print.nesting_spec}
\alias{apply_sentinel}
\title{Declare a partially nested design}
\description{
\code{nesting_spec} validates a declared design, records its nesting families, and
constructs the factor of realized cells that the fitting functions use.
\code{apply_sentinel} converts structurally undefined \code{NA} values to an explicit
sentinel level, refusing genuine missingness.
}
\usage{
nesting_spec(data, formula, nests, fit = c("lm", "glm", "lmer", 
             "glmer", "clm", "clmm", "brms"), family = NULL, random = NULL,
             cell_name = "cell")

apply_sentinel(data, var, where = NULL, sentinel = "none") 
}
\arguments{
  \item{data}{a data frame containing the outcome, the nesting variables, and any covariates.}
  \item{formula}{the desired model, written freely in ordinary R syntax; crossed terms are welcome and are transformed rather than fitted as given. Random-effects terms may be written here in \code{lme4}/\code{brms} syntax.}
  \item{nests}{the nesting links, quoted or bare: \code{"inversion \%in\% chord_type"}, or \code{c(inv \%in\% chord, doub \%in\% inv)} for a chain.}
  \item{fit}{the modelling engine.}
  \item{family}{family for the engines that need one, quoted or as a bare call, e.g. \code{cumulative("probit")}.}
  \item{random}{random-effects terms, if not written into \code{formula}.}
  \item{cell_name}{name for the constructed cell factor. An existing column of that name is refused rather than overwritten.}
  \item{var}{for \code{apply_sentinel}, the variable to convert.}
  \item{where}{logical vector marking the rows in which \code{var} is structurally undefined. Optional: when omitted, every \code{NA} is converted and a warning is given, since any that were genuine missing data are then coded as an existing condition. When supplied, an \code{NA} outside those rows is refused as genuine missingness, and a non-\code{NA} inside them is refused as inconsistent coding.}
  \item{sentinel}{the sentinel value; numeric variables default to \code{0}.}
  \item{x}{an object to print.}
}
\details{
Validation is deliberately strict. An \code{NA} in a nesting variable is refused, since
R's model machinery would delete those rows silently; a continuous variable in any
non-leaf position is refused; an ordered-factor outcome sent to an engine that cannot
model one is refused; and \code{||} on a factor is refused, because it does not give
the diagonal covariance it is usually taken to give.
}
\value{
\code{nesting_spec} returns an object of class \code{nesting_spec}: a list whose
\code{data} element is the supplied data with the cell factor added, together with the
nesting families, the realized-cell table, the covariates, and the declared engine.
\code{apply_sentinel} returns the modified data frame.
}
\seealso{
\code{\link{nest_fit}}, \code{\link{estimand}}, \code{\link{cell_formula}}
}
\examples{
\dontrun{
sp <- nesting_spec(dat, response ~ chord_type * inversion + training,
                   "inversion %in% chord_type")
sp
}
}
