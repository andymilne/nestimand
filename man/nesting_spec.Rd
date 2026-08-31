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
             "glmer", "clm", "clmm", "brm"), family = NULL, random = NULL,
             cell_name = "cell")

apply_sentinel(data, var, where = NULL, sentinel = "none") 
}
\arguments{
  \item{data}{a data frame containing the outcome, the nesting variables, and any covariates.}
  \item{formula}{the desired model, written freely in ordinary R syntax; crossed terms are welcome and are transformed rather than fitted as given. Random-effects terms may be written here in \code{lme4}/\code{brms} syntax.}
  \item{nests}{the categorical design - which variable sits inside which, and which sit inside nothing. Quoted or bare: \code{"inversion \%in\% chord_type"}, or \code{c(inv \%in\% chord, doub \%in\% inv)} for a chain. A parent may hold several nested variables, declared one per entry or together as \code{c(inversion, X1) \%in\% chord_type}; \code{inversion + X1 \%in\% chord_type} and \code{(inversion + X1) \%in\% chord_type} read the same way, the brackets being optional since the left side is a list of names rather than a model term. Anything that is not a list of names - an interaction, say - is refused. Each variable is nested in exactly one parent, and \code{\%in\%} is not chained within one entry - \code{Z \%in\% inversion \%in\% chord_type} is two entries, since R would read it as \code{(Z \%in\% inversion) \%in\% chord_type}. Only the nesting has to be declared: a categorical variable the formula crosses with the declared structure - \code{chord_type * inversion * top_note} - is part of the design and joins the cell factor without being named here. An entry with no \code{\%in\%} folds one in deliberately, which is how an ordered factor or a variable entering additively is made a design variable rather than a covariate.}
  \item{fit}{the modelling engine, named for its fitting function: \code{"lm"}, \code{"glm"}, \code{"lmer"}, \code{"glmer"}, \code{"clm"}, \code{"clmm"}, \code{"brm"}. \code{"brms"} is accepted for the last, being the package's name rather than the function's.}
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

A declared family is a tree: one parent may hold several nested variables, and each
of those may hold further ones, to any depth. The strata of a nested variable are its
own ancestors - a sibling is not one of them - so \code{inversion} and \code{X1}
both inside \code{chord_type} are compared within chord type, not within
chord type crossed with each other.

The cell factor is the realized categorical design, and it is formed from the
declared nesting variables together with every categorical variable the formula
crosses with them. Folding those in changes no fitted model - \code{cell + cell:x}
and a cell factor over the enlarged design span the same space - but it drops the
columns of any combination the design does not realize, which the interaction form
would carry as aliased coefficients, and it makes the variable a legal target of an
estimand. Three kinds of variable are left as covariates: one entering additively,
whose effect is thereby declared common to every condition; a numeric one; and an
ordered factor, whose contrasts say it is meant as a quantity. Naming any of them in
\code{nests} folds it in.

A declared nesting variable that is \emph{numeric} is not part of the cell factor: a
continuous nested variable enters as a slope within the realized cells,
\code{cell:x}, which is what nesting a quantity means. \code{nesting_spec} says so
when it happens, and the printed structure marks the variable, since a variable whose
values are labels is easily left numeric by accident - convert it to a factor first if
that is what it is.

A formula that spans fewer dimensions than there are realized cells - \code{+}
between two variables the design admits an interaction between - asks for a
restricted mean structure, and is fitted as written. The cell factor is saturated by
construction and cannot express a restriction, so such a declaration is fitted in the
effects parameterization instead, with each nested variable given its parent:
\code{chord_type * (inversion + top_note)} is fitted as \code{chord_type +
top_note + chord_type:inversion + chord_type:top_note}, since \code{inversion} has no
effect outside the strata in which it varies. \code{fitting_mode()} reports which
parameterization a declaration leads to, and the printed spec shows the formula that
will be fitted. On \pkg{brms} a restricted structure needs
\code{\link{chain_priors}} to hold its uninformative columns at zero, as the chain
form always does.
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

## Two variables inside one parent, and one inside another. The declarations
## are given one per entry, bare or quoted.
sp2 <- nesting_spec(dat, response ~ chord_type * inversion * X1 * Z,
                    c("inversion %in% chord_type", "X1 %in% chord_type",
                      "Z %in% inversion"))
sp2                                   # chord_type > inversion, chord_type > X1,
                                      # inversion > Z

## Links sharing a parent may be gathered instead. These three declare the same
## structure, and `+` is read as a list of children, not as a model term.
nesting_spec(dat, response ~ chord_type * inversion * X1,
             c(inversion, X1) %in% chord_type)
nesting_spec(dat, response ~ chord_type * inversion * X1,
             inversion + X1 %in% chord_type)
nesting_spec(dat, response ~ chord_type * inversion * X1,
             c("inversion %in% chord_type", "X1 %in% chord_type"))
}
}
