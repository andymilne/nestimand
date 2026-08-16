\name{nestimand-package}
\alias{nestimand-package}
\alias{nestimand}
\docType{package}
\title{Estimands for Partially Nested Designs by Cell Translation}
\description{
A partially nested design contains a factor whose levels exist only within
certain levels of another factor: chord inversion, say, which is undefined for
augmented triads. Default modelling practice mishandles such designs quietly,
producing aliased coefficients with misleading labels and marginal summaries
that change when factor levels are permuted.

\pkg{nestimand} declares the design once, fits it in the realized-cell
parameterization, where estimation is unconditionally well posed, and
translates grids, labels, draws, priors, and results back into the original
variable space in which research questions and reports are stated. Every
contrast that crosses the structural boundary requires an explicit weighting
policy, because no such contrast is defined without one.
}
\section{Declaring the design}{
\code{\link{nesting_spec}} is the entry point: it validates the declaration,
constructs the realized-cell factor, and refuses the coding errors that
otherwise pass silently. \code{\link{apply_sentinel}} converts structurally
undefined \code{NA} values to an explicit sentinel level.
}
\section{Fitting}{
\code{\link{nest_fit}} dispatches to the declared engine and attaches the call
it ran, so \code{\link{show_code}} returns the whole pipeline.
\code{\link{cell_formula}} gives either parameterization,
\code{\link{random_terms}} translates a declared random structure, and
\code{\link{fitting_mode}} states which parameterization is chosen and why.
}
\section{Estimands}{
\code{\link{estimand}} returns a labelled contrast table in the original
variable names, with partial-identification bounds
(\code{\link{add_bounds}}) and a reorder self-check.
\code{\link{nest_policy}} builds the weighting policy;
\code{\link{counterfactual_grid}} forms the G-computation grid;
\code{\link{latent_estimand}} and \code{\link{latent_draws}} compute the same
quantity as a linear functional of the cell coefficients, which is exact on any
linear scale and necessary for ordinal fits.
}
\section{Priors}{
\code{\link{nest_prior}} translates an elliptical prior between effect space
and cell space; \code{\link{prior_audit}} and
\code{\link{prior_for_estimand}} report what a stated prior implies.
\code{\link{chain_priors}} derives the \code{constant(0)} declarations that the
chain parameterization requires in \pkg{brms}.
}
\section{Function index}{
\tabular{ll}{
  \code{\link{nesting_spec}} \tab declare the design \cr
  \code{\link{apply_sentinel}} \tab convert structural \code{NA} to a sentinel \cr
  \code{\link{cell_formula}} \tab fitting formula, cells or effects \cr
  \code{\link{chain_terms}} \tab the ancestry-closed chain terms \cr
  \code{\link{add_cells}} \tab recompute the cell factor on a grid \cr
  \code{\link{cell_grid}} \tab one row per realized cell \cr
  \code{\link{effect_basis}} \tab the matrix A relating cells to effects \cr
  \code{\link{translation}} \tab A and its inverse, with names \cr
  \code{\link{fitting_mode}} \tab which parameterization, and why \cr
  \code{\link{boundary_vars}} \tab the variables that cross the boundary \cr
  \code{\link{nest_policy}} \tab a distribution over versions \cr
  \code{\link{counterfactual_grid}} \tab the G-computation grid \cr
  \code{\link{policy_weights}} \tab policy weights for a grid \cr
  \code{\link{policy_vertices}} \tab the single-version policies \cr
  \code{\link{nest_fit}} \tab fit, with the call attached \cr
  \code{\link{random_terms}} \tab translate a random structure \cr
  \code{\link{random_slope_rank}} \tab rank of a structural random slope \cr
  \code{\link{has_thresholds}} \tab does the family carry thresholds \cr
  \code{\link{estimand}} \tab contrasts in the original variables \cr
  \code{\link{add_bounds}} \tab partial-identification bounds \cr
  \code{\link{show_code}} \tab the code that produced a result \cr
  \code{\link{latent_estimand}} \tab contrast as a linear functional \cr
  \code{\link{latent_draws}} \tab the same, draw by draw \cr
  \code{\link{policy_contrast_matrix}} \tab stratum rows of the design \cr
  \code{\link{nest_prior}} \tab translate a prior between spaces \cr
  \code{\link{prior_audit}} \tab what a prior implies elsewhere \cr
  \code{\link{prior_for_estimand}} \tab what a prior implies for a contrast \cr
  \code{\link{prior_stanvars}} \tab the brms data blocks \cr
  \code{\link{prior_statement}} \tab the brms prior statement \cr
  \code{\link{chain_priors}} \tab constant(0) declarations \cr
  \code{\link{chain_prior_object}} \tab the same, as a brms object \cr
  \code{\link{zero_columns}} \tab uninformative design columns \cr
  \code{\link{chain_random_zeros}} \tab the same, per grouping factor \cr
  \code{\link{mfx_canonical}} \tab normalize engine contrast labels \cr
  \code{\link{mfx_hypothesis}} \tab the accepted hypothesis spelling \cr
}
}
\references{
Milne, A. J. Partially nested variables: modelling and estimation in R.
Companion document, distributed with the package.
}
\author{Andrew J. Milne}
\keyword{package}
