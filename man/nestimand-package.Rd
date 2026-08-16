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
\section{The everyday workflow}{
Five functions cover an ordinary analysis.

\tabular{ll}{
  \code{\link{nesting_spec}} \tab declare the design once: which factor nests within which \cr
  \code{\link{nest_fit}} \tab fit it, keeping the call that was run \cr
  \code{\link{estimand}} \tab contrasts in the original variable names, with bounds \cr
  \code{\link{show_code}} \tab the code that produced a result, for saving and adapting \cr
  \code{\link{nest_summary}} \tab the coefficient table, in the original parameterization \cr
  \code{\link{apply_sentinel}} \tab convert structurally undefined \code{NA} to a sentinel \cr
}

\preformatted{sp <- nesting_spec(dat, response ~ chord_type * inversion + training,
                   "inversion \%in\% chord_type")
m  <- nest_fit(sp)
e  <- estimand(m, chord_type, policy = "equal")
show_code(e)}
}
\section{As the analysis requires}{
\tabular{ll}{
  \code{\link{latent_estimand}} \tab the contrast as a linear functional: exact, and necessary for ordinal fits \cr
  \code{\link{latent_draws}} \tab the same, draw by draw, from a \pkg{brms} posterior \cr
  \code{\link{nest_policy}} \tab build and inspect a weighting policy before using it \cr
  \code{\link{nest_prior}} \tab state a prior in effect space or cell space \cr
  \code{\link{prior_audit}} \tab what a stated prior implies in the other space \cr
  \code{\link{prior_for_estimand}} \tab what it implies for a named contrast \cr
  \code{\link{chain_priors}} \tab the \code{constant(0)} declarations chain-mode \pkg{brms} requires \cr
  \code{\link{random_terms}} \tab how a declared random structure is translated \cr
  \code{\link{random_slope_rank}} \tab the rank deficiency of a structural random slope \cr
}
}
\section{Inspecting the translation}{
Not needed for an analysis, but useful for seeing what the package did.

\tabular{ll}{
  \code{\link{cell_formula}} \tab the fitting formula, in either parameterization \cr
  \code{\link{effect_basis}} \tab the matrix A relating cell means to identified effects \cr
  \code{\link{translation}} \tab A and its inverse, with the cell and effect names \cr
  \code{\link{fitting_mode}} \tab which parameterization is chosen, and why \cr
  \code{\link{chain_terms}} \tab the ancestry-closed chain terms \cr
  \code{\link{boundary_vars}} \tab the variables that cross the structural boundary \cr
  \code{\link{zero_columns}} \tab which design columns the data cannot inform \cr
  \code{\link{cell_grid}} \tab one row per realized cell \cr
  \code{\link{has_thresholds}} \tab whether the declared family carries thresholds \cr
}
}
\section{Called by emitted code}{
The code returned by \code{\link{show_code}} has to run, so the functions it calls
are exported. They are not intended to be typed.

\tabular{ll}{
  \code{\link{add_cells}} \tab recompute the cell factor on a grid \cr
  \code{\link{counterfactual_grid}} \tab the G-computation grid, with policy weights \cr
  \code{\link{policy_weights}} \tab the weights a policy places on a grid \cr
  \code{\link{policy_vertices}} \tab the single-version policies bounding the estimand \cr
  \code{\link{policy_contrast_matrix}} \tab weighted design-matrix rows, one per stratum \cr
  \code{\link{add_bounds}} \tab attach the partial-identification region \cr
  \code{\link{mfx_canonical}} \tab fix contrast labels and direction across engine versions \cr
  \code{\link{mfx_term_column}} \tab locate the contrast label column \cr
  \code{\link{mfx_hypothesis}} \tab the \code{hypothesis} spelling the installed version accepts \cr
  \code{\link{mfx_formula_hypothesis}} \tab whether that version takes a formula \cr
  \code{\link{chain_prior_object}} \tab the declarations as a \pkg{brms} prior \cr
  \code{\link{chain_random_zeros}} \tab the uninformative columns, per grouping factor \cr
  \code{\link{prior_stanvars}} \tab the \pkg{brms} data blocks for a translated prior \cr
  \code{\link{prior_statement}} \tab the \pkg{brms} prior statement \cr
}
}
\references{
Milne, A. J. Partially nested variables: modelling and estimation in R.
Companion document, distributed with the package.
}
\author{Andrew J. Milne}
\keyword{package}
