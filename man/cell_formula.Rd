\name{cell_formula}
\alias{cell_formula}
\alias{chain_terms}
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
\code{fitting_mode} states which parameterization is chosen and why.
}
\usage{
cell_formula(spec, mode = c("cells", "effects"), intercept = NULL)
chain_terms(spec)
add_cells(spec, newdata, unrealized = c("drop", "error"))
cell_grid(spec, data = spec$data, covariates = c("mean", "keep"))
effect_basis(spec)
translation(spec)
fitting_mode(spec, engine = "marginaleffects", priors = NULL)
}
\arguments{
  \item{spec}{a \code{nesting_spec} object.}
  \item{mode}{\code{"cells"} places the whole categorical structure on one realized-cell factor; \code{"effects"} uses the ancestry-closed chain form, whose coefficients read as effects.}
  \item{intercept}{whether the coding carries an intercept. Chosen from the family when \code{NULL}: ordinal families already carry one in their thresholds.}
  \item{newdata}{a grid or data frame in the original variable space.}
  \item{unrealized}{what to do with rows naming combinations that do not exist: drop them, or stop.}
  \item{data}{data from which covariate means are taken.}
  \item{covariates}{hold covariates at their means, or keep the supplied values.}
  \item{engine}{the estimand engine the fit is intended for; \code{"emmeans"} requires the chain form.}
  \item{priors}{declared priors, which may pin the parameterization.}
}
\value{
\code{cell_formula} a formula; \code{chain_terms} a character vector of terms;
\code{add_cells} the grid with the cell factor added and the number of dropped rows in
the \code{nestimand_dropped} attribute; \code{cell_grid} a data frame with one row per
realized cell; \code{effect_basis} a square full-rank matrix; \code{translation} an
object holding A and its inverse; \code{fitting_mode} a string with the reason in its
\code{reason} attribute.
}
\seealso{
\code{\link{nesting_spec}}, \code{\link{nest_fit}}, \code{\link{estimand}}
}
