\name{show_code}
\alias{show_code}
\alias{show_code.default}
\alias{show_code.nestimand_estimand}
\alias{engine_types}
\alias{print.nestimand_types}
\title{The code that produced a result}
\description{
Every fit and every estimand carries the script that produced it. \code{show_code}
prints that script and returns it invisibly, so a result can be saved, adapted,
checked by a reader, or handed to a collaborator who does not use this package
beyond installing it.

The code is what ran, not a reconstruction of it: \code{\link{nest_fit}} and
\code{\link{nest_estimand}} assemble the call as text and evaluate that text, so
the two cannot drift apart. It calls \pkg{nestimand} rather than inlining cell
tables and weight vectors as literals, on the ground that a reader of the code
can be expected to have the package installed - which is why the package exports
more than the three functions an analysis uses.

\code{engine_types} reports which \code{type} values the installed
\pkg{marginaleffects} accepts for a given model, since the vocabulary differs by
model class and between releases.
}
\usage{
show_code(x, ...) 

engine_types(model = NULL) 
}
\arguments{
  \item{x}{a fit or an estimand, or anything else carrying a \code{nestimand_code} attribute.}
  \item{model}{a fitted model, whose class decides which types are available. \code{NULL} lists the types the package itself knows about.}
  \item{...}{unused, for method consistency.}
}
\value{
\code{show_code} returns the lines of code invisibly, having printed them.
\code{engine_types} returns a character vector, classed for printing.
}
\seealso{
\code{\link{nest_fit}}, \code{\link{nest_estimand}}, \code{\link{nest_summary}}
}
\examples{
\dontrun{
m <- nest_fit(dat, rating ~ chord_type * inversion + training)
e <- nest_estimand(m, chord_type, policy = "equal")

show_code(m)   # the fit, including the declaration it was built from
show_code(e)   # the policy, the grid, and the contrast

## the code without computing anything first
nest_estimand(m, chord_type, policy = "equal", dry_run = TRUE)
}
}
