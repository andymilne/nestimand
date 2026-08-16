\name{mfx_canonical}
\alias{mfx_canonical}
\alias{mfx_term_column}
\alias{mfx_hypothesis}
\alias{mfx_formula_hypothesis}
\title{Engine label and argument conventions}
\description{
Presentation conventions have moved between \pkg{marginaleffects} versions: the
\code{hypothesis} argument changed from a string to a formula, the contrast label moved
from \code{term} to \code{hypothesis}, labels gained parentheses, and the direction of
the difference reversed. These functions normalize all four, so that the same analysis
reports the same signs and labels on any installed version.
}
\usage{
mfx_canonical(d, levs = NULL)
mfx_term_column(d)
mfx_hypothesis(contrast = "pairwise", version = NULL)
mfx_formula_hypothesis(version)
}
\arguments{
  \item{d}{an estimand table.}
  \item{levs}{the declared factor levels, which fix the contrast direction: earlier level minus later.}
  \item{contrast}{the comparison type.}
  \item{version}{package version to test; the installed one by default.}
}
\value{
\code{mfx_canonical} the table with a \code{term} column in declared level order,
estimates and intervals negated where the engine ran the difference the other way, and
the number of reversed rows in its \code{nestimand_flipped} attribute.
\code{mfx_hypothesis} the spelling the installed version accepts.
}
\seealso{
\code{\link{estimand}}
}
