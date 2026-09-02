\name{cell_formula}
\alias{cell_formula}
\alias{chain_terms}
\alias{reduced_design}
\alias{with_reduced}
\alias{uses_reduced}
\alias{add_cells}
\alias{cell_grid}
\alias{effect_basis}
\alias{translation}
\alias{print.nestimand_translation}
\alias{fitting_mode}
\title{The translation between original variables and realized cells}
\description{
These functions own the maps between the space in which questions are stated and the
space in which estimation is well posed. \code{cell_formula} gives the fitting formula
in either parameterization; \code{add_cells} recomputes the cell factor on a grid;
\code{effect_basis} returns the matrix A relating cell means to identified effects; and
\code{fitting_mode} states which parameterization is chosen and why;
\code{reduced_design} gives the design a restricted declaration implies, one row per
realized cell and one column per identified effect; \code{with_reduced} carries
those columns alongside the data, looked up by cell; and \code{uses_reduced} says
whether the fit this declaration produces needs them, which the random side can
require even when the mean structure does not.
}
\usage{
cell_formula(spec, mode = c("cells", "reduced"), intercept = NULL) 

chain_terms(spec) 

reduced_design(spec) 

with_reduced(spec, data = spec$data) 

uses_reduced(spec) 

add_cells(spec, newdata, unrealized = c("drop", "error")) 

cell_grid(spec, data = spec$data, covariates = c("mean", "keep")) 

effect_basis(spec) 

translation(spec) 

fitting_mode(spec) 
}
\arguments{
  \item{spec}{a \code{nesting_spec} object.}
  \item{mode}{\code{"cells"} places the whole categorical structure on one realized-cell factor, which is saturated by construction; \code{"reduced"} is what a declaration asking for less than that is fitted as - the same design, one column per identified effect of the declared structure, with the columns the data cannot inform removed, so nothing is aliased and nothing has to be held at zero.}
  \item{intercept}{whether the coding carries an intercept. Chosen from the family when \code{NULL}: ordinal families already carry one in their thresholds.}
  \item{newdata}{a grid or data frame in the original variable space.}
  \item{unrealized}{what to do with rows naming combinations that do not exist: drop them, or stop.}
  \item{data}{data from which covariate means are taken.}
  \item{covariates}{hold covariates at their means, or keep the supplied values.}
  \item{x}{an object to print.}
}
\value{
\code{cell_formula} a formula; \code{chain_terms} a character vector of terms;
\code{reduced_design} a matrix with one row per realized cell, whose columns are
named \code{dm_} plus the design-matrix column they came from, with \code{:}
written \code{.} so that the name is ordinary in a formula, and whose
\code{effect_names} attribute names the effect each column stands for - which is
what \code{\link{nest_summary}} reports;
\code{with_reduced} the data frame with those columns appended;
\code{uses_reduced} \code{TRUE} or \code{FALSE};
\code{add_cells} the grid with the cell factor added and the number of dropped rows in
the \code{nestimand_dropped} attribute; \code{cell_grid} a data frame with one row per
realized cell; \code{effect_basis} a square full-rank matrix; \code{translation} an
object holding A and its inverse; \code{fitting_mode} a string with the reason in its
\code{reason} attribute.
}
\seealso{
\code{\link{nesting_spec}}, \code{\link{nest_fit}}, \code{\link{nest_estimand}}
}
\keyword{internal}
