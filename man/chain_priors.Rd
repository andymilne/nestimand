\name{chain_priors}
\alias{chain_priors}
\alias{print.nestimand_chain_priors}
\alias{show_code.nestimand_chain_priors}
\alias{chain_prior_object}
\alias{zero_columns}
\alias{chain_random_zeros}
\title{Declarations for the chain parameterization in brms}
\description{
The chain parameterization keeps the original variables as predictors, so its
coefficients read as effects and \pkg{emmeans} sees the nesting. The design then carries
columns the data cannot inform: those for conditions that do not exist, and a further set
reconstructable from the others. \code{chain_priors} derives which, from the design
matrix rather than by assumption, and distinguishes the two kinds.
}
\usage{
chain_priors(spec, regularize = "normal(0, 5)", data = NULL) 

chain_prior_object(x) 

zero_columns(formula, data, drop_intercept = TRUE, spec_for_ref = NULL) 

chain_random_zeros(spec, data = NULL) 
}
\arguments{
  \item{spec}{a \code{nesting_spec} object.}
  \item{regularize}{prior for the remaining coefficients, or \code{NULL} for none. Identification comes from the declarations; this does only its ordinary regularizing work.}
  \item{data}{data from which the design matrix is built.}
  \item{x}{a \code{nestimand_chain_priors} object.}
  \item{formula}{a model formula, for \code{zero_columns}.}
  \item{drop_intercept}{drop the intercept column before examining the design. The intercept is never a candidate for constraining, but it is retained for the rank analysis, which would otherwise understate the number of constraints.}
  \item{spec_for_ref}{a \code{nesting_spec}, used to take the reference condition from the declared level order rather than from a pivot.}
}
\details{
A structural zero says a condition does not exist. An identification constraint is a
coding choice, like a reference level, which leaves every estimand unchanged. The fixed
and random parts are built from the same design matrix, so the same columns are held in
both.
}
\value{
\code{chain_priors} an object of class \code{nestimand_chain_priors}, whose
\code{table} names each held coefficient and its kind and whose \code{code} is the block
as \pkg{brms} code. \code{chain_prior_object} the \pkg{brms} prior;
\code{zero_columns} a list of the structural, identification, and retained columns.
}
\seealso{
\code{\link{nest_prior}}, \code{\link{nest_fit}}, \code{\link{fitting_mode}}
}
