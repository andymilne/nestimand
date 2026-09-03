\name{nest_fit}
\alias{nest_fit}
\alias{random_terms}
\alias{has_thresholds}
\alias{print.nestimand_code}
\title{Fit a partially nested design}
\description{
\code{nest_fit} takes the engine, a formula and the data, and fits the model the
formula names. Where one variable is undefined at some levels of another - inversion
for an augmented triad, a dose that exists only in the treatment arm - the conditions
that do not exist are removed from the design rather than estimated, and nothing is
left aliased or held at zero.

\preformatted{m <- nest_fit("lm", rating ~ chord_type * inversion + training, dat)}

The structure is read off the data, from the rows where a variable is \code{NA}. Say
it yourself with \code{nests} when the data does not carry those gaps, or when the
structure is not the one the data suggests.

What comes back is the engine's own model object, with the declaration attached to it -
an \code{lm} is an \code{lm}, a \code{brmsfit} a \code{brmsfit} - so every method the
engine and \pkg{marginaleffects} provide works on it unchanged. Printing it therefore
shows the engine's own coefficient names, which are those of the columns the model was
fitted on. \code{\link{nest_summary}} is how those coefficients are read back as
effects of your own variables.

The call is assembled as text and evaluated, so the returned fit carries the code that
produced it and \code{\link{show_code}} can show the whole pipeline - the code that
ran, not a reconstruction of it. Arguments in \code{...} reach the engine unaltered.
}
\usage{
nest_fit(fit, formula, data, family = NULL, nests = NULL,
         random = NULL, dry_run = FALSE, ..., .env = parent.frame())

random_terms(spec)

has_thresholds(spec)
}
\arguments{
  \item{fit}{the engine, named for its fitting function: \code{"lm"}, \code{"glm"}, \code{"lmer"}, \code{"glmer"}, \code{"clm"}, \code{"clmm"} or \code{"brm"}.}
  \item{formula}{the model formula, in your own variables. Random terms may be written into it with a bar.}
  \item{data}{the data frame.}
  \item{family}{the family, where the engine takes one, e.g. \code{"binomial"} or \code{cumulative("probit")}.}
  \item{nests}{the nesting structure, written \code{inversion \%in\% chord_type} - quoted or not, one entry or several. Omit it and the structure is read from the data: a variable is taken to be nested where it is \code{NA} for whole strata of another variable and never otherwise. What was read is reported, in the syntax you would have written, so it can be checked.}
  \item{random}{the random structure, if you would rather give it separately than write bars into the formula: \code{random = "(1 | participant)"}.}
  \item{dry_run}{return the code without fitting, which matters when the fit is a sampling run.}
  \item{...}{passed to the fitting function and written into the saved code: this is how \code{prior}, \code{sample_prior}, \code{save_pars}, \code{init}, \code{threads}, \code{file}, \code{chains}, \code{iter}, \code{seed} and the rest reach \pkg{brms}, and their counterparts reach the other engines. A \pkg{brms} \code{prior} for the random-effects standard deviations and correlations, the thresholds or sigma is handed over untouched; one of \code{class = "b"} describes the mean structure, so it is translated to the coefficients the model is fitted on and reaches Stan through stanvars.}
  \item{.env}{environment in which the assembled call is evaluated.}
  \item{spec}{a \code{\link{nesting_spec}}.}
  \item{x}{an object to print.}
}
\details{
On the random side, a term is translated only if it reaches across the structure. A
random intercept, a slope on a variable that every unit meets at every level, or a
hand-written grouping structure such as \code{(1 | participant:chord_type)} comes back
exactly as declared. A slope that does reach across is written on the columns the
declared structure identifies: over the original factors it would carry columns for
conditions that do not exist, which no amount of data can inform. A covariance-structure
wrapper such as \code{diag()} survives translation.

A random term naming something the model formula does not contain is refused.

\preformatted{sp <- nesting_spec(dat, rating ~ chord_type * inversion +
                     (chord_type * inversion | participant),
                   inversion \%in\% chord_type, fit = "lmer")
random_terms(sp)
## (1 + dm_chord_typedim + ... | participant)}
}
\value{
\code{nest_fit} the fitted model, with the code in its \code{nestimand_code}
attribute and the declaration in \code{nestimand_spec}, so that
\code{\link{nest_estimand}} and \code{\link{nest_summary}} can be called on the model
alone. \code{update} carries the declaration across a refit, so a simplified model
remains usable. The class \code{nestimand_fit} is added behind the engine's own, where
it changes no behaviour: \code{predict}, \code{summary}, \code{vcov} and the
\pkg{marginaleffects} methods apply as before. S4 fits are not reclassified. With
\code{dry_run = TRUE} an object of class \code{nestimand_code} instead.

\code{random_terms} a character string: the random structure as it will be written to
the engine.
}
\seealso{
\code{\link{nest_estimand}}, \code{\link{nest_summary}}, \code{\link{nesting_spec}}
}
