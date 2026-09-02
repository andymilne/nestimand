\name{nesting_spec}
\alias{nesting_spec}
\alias{infer_nests}
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
nesting_spec(data, formula, nests = NULL, fit = c("lm", "glm", "lmer", 
             "glmer", "clm", "clmm", "brm"), family = NULL, random = NULL,
             cell_name = "cell")

apply_sentinel(data, var, where = NULL, sentinel = "none") 

infer_nests(data, vars, max_parents = 2) 
}
\arguments{
  \item{data}{a data frame containing the outcome, the nesting variables, and any covariates.}
  \item{formula}{the desired model, written freely in ordinary R syntax; crossed terms are welcome and are transformed rather than fitted as given. Random-effects terms may be written here in \code{lme4}/\code{brms} syntax.}
  \item{vars}{for \code{infer_nests}, the design variables among which to look for structure.}
  \item{max_parents}{for \code{infer_nests}, the largest number of variables a nesting is allowed to be inside.}
  \item{nests}{the categorical design - which variable sits inside which, and which sit inside nothing. May be omitted, in which case it is read off the data: a variable is taken to be nested within the smallest set of design variables such that, within every combination of their levels, it is either always \code{NA} or never \code{NA}. What is inferred is reported in the syntax you would have written, and the undefined rows are given a sentinel automatically, so a data frame carrying \code{NA} where a variable does not apply needs no declaration at all. Where the gaps are not structural - an \code{NA} pattern no combination of the other variables picks out - it stops rather than inventing a structure, which is the distinction \code{apply_sentinel} exists to draw; and where two variables explain the same gaps equally well it says so and asks for a declaration, that being something the data cannot settle. Inference sees the sample rather than the design, so gaps that are regular but not structural would be read as a nesting: that is why it is announced. Quoted or bare: \code{"inversion \%in\% chord_type"}, or \code{c(inv \%in\% chord, doub \%in\% inv)} for a chain. A parent may hold several nested variables, declared one per entry or together as \code{c(inversion, X1) \%in\% chord_type}; \code{inversion + X1 \%in\% chord_type} and \code{(inversion + X1) \%in\% chord_type} read the same way, the brackets being optional since the left side is a list of names rather than a model term. Anything that is not a list of names - an interaction, say - is refused. Each variable is nested in exactly one parent, and \code{\%in\%} is not chained within one entry - \code{Z \%in\% inversion \%in\% chord_type} is two entries, since R would read it as \code{(Z \%in\% inversion) \%in\% chord_type}. Only the nesting has to be declared: a categorical variable the formula crosses with the declared structure - \code{chord_type * inversion * top_note} - is part of the design and joins the cell factor without being named here. An entry with no \code{\%in\%} folds one in deliberately, which is how an ordered factor or a variable entering additively is made a design variable rather than a covariate.}
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

The formula says what to fit and the \code{nests} argument says which combinations
exist. Between them they determine the model, and nothing is added to it: a \code{+}
where the design would admit a \code{*} is a restriction the user stated, and it is
fitted as the restriction it is. \code{response ~ chord_type + inversion} fits one
inversion effect shared by every chord type, not a separate one for each.

What is fitted is the design that formula gives over the realized cells, with the
columns that design cannot inform removed - identically zero first, for conditions
that do not exist, then linearly redundant. What remains is full rank, so no
coefficient is aliased and nothing has to be held at zero on any engine, \pkg{brms}
included. \code{~ 0 + cell} is used when, and only when, the formula crosses the
structure fully, where the two are the same model in different coordinates.
\code{\link{reduced_design}} shows the columns and the effect each stands for, and
the printed spec shows the structure that will be fitted, counted in realized cells.

The level a nested variable takes where it is undefined - the sentinel - is not a
value it takes. Rows carrying it contribute nothing to any column that variable helps
to build, and its contrasts run among the levels that exist. So in
\code{chord_type + inversion} the inversion coefficients read as inversion effects -
second inversion versus root position, shared across chord types - rather than as
comparisons against whichever stratum the sentinel stands for.

The random side follows the same rule: a random structure that asks for less than the
saturated one is fitted on the same columns rather than enlarged to an unstructured
covariance over every realized condition. See \code{random_structure} in
\code{\link{nest_fit}}.
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
