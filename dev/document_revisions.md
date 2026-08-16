# Revisions to `partially_nested_variables.md`

*Four edits, in document order: an amendment to the front-page Purpose note; a new
closing subsection for Section 1; a rewritten Section 6; and one checklist item.
All numbers cited are executed — R 4.3.3, lme4 1.1.35.1, on a
six-trials-per-participant-per-cell version of the demonstration data (2400 rows,
40 participants). Generating scripts: `dev/re_evidence.R` and `dev/re_closure.R`
in the package repository.*

---

## Edit 1 — Purpose note, front page

The current note describes the document's chain parameterization as accepting
"some residual awkwardness (most visibly for random effects, Section 6) as the
price of parameters a reader can interpret directly". That understates the
matter, and Section 1 below now says why. Replace that clause with:

> …and develops interpretable, hand-codeable remedies – the chain
> parameterization, with `emmeans` and `marginaleffects` recipes – accepting one
> substantive limit as the price of parameters a reader can interpret directly:
> the chain form cannot express every random-effects structure the design admits
> (Section 1, "A second fact about the design", and Section 6).

---

## Edit 2 — new subsection, at the end of Section 1

*Placed after "The unavoidable choice", as the second thing a reader should know
about such a design.*

### A second fact about the design

The choice above concerns fixed effects: which conditions exist, and how to
average across them. A second consequence of the same ten-cell structure
concerns random effects, and it is not a choice at all but a limit on what can be
written down.

A participant-level random effect is, at bottom, a set of personal deviations —
one number for each condition a participant can encounter. This design has ten
such conditions, so a participant's deviations form a vector of ten numbers, and
the general description of how participants vary is an arbitrary 10 × 10
covariance matrix. Every particular random-effects structure is a *constraint* on
that matrix: a random intercept alone says all ten deviations are the same
number; the additive grouping structure of Section 6 says cells sharing a chord
type are equally correlated; an unstructured covariance imposes nothing at all.

The question that decides the parameterization is therefore a question about
coverage. **Can every one of those structures be written in the chain form?**

It cannot. The chain form's random basis — a slope on `chord_type` together with
a slope on `chord_type:inversion` — has sixteen columns, one for each cell of the
full crossing, of which six do not exist. Those six columns are identically zero
for every participant, so the sixteen columns span only ten dimensions, and a
covariance placed on them has six directions along which the data are silent.
No sample size alters this: the column for "augmented chords in first inversion"
is zero because no participant ever met that condition, not because too few did.
To keep the random side identified, the chain form must retreat to additive
grouping terms, and those impose one variance per level with equal correlations
within a stratum. That is a proper subset of the 10 × 10 matrices. Structures
outside it — an unstructured covariance, a diagonal one, most of what `brms`
offers — have no chain-form equivalent.

Put the same ten conditions on a *single* factor, one level per realized cell,
and the limit disappears. A deviation over realized cells simply *is* a vector of
ten numbers, so `(0 + cell | participant)` with an unrestricted covariance is the
general case by construction, and every other structure is a constrained
submodel of it. Nothing can be specified that cannot be represented. The cell
parameterization is closed under arbitrary random-effects structures; the chain
parameterization is not.

This is the strongest argument for the companion package, and it is worth being
plain about why. The fixed-effect difficulties of Sections 3 and 4 are real but
surmountable by hand: a careful analyst who follows the recipes there will reach
correct answers. The closure failure is not surmountable by care, because it is
not a matter of discipline but of what the parameterization can express — and the
remedy, working in cell space, exacts a translation burden that is mechanical for
software and unreasonable by hand. A worked instance appears in Section 6.3:
under the grouping structure, "how much do participants differ in the chord-type
effect?" is answered by reading off a named variance component, whereas under the
cell parameterization the same question becomes a quadratic form c′Σc, where c is
the *same* ten-element weight vector that defines the corresponding fixed-effect
contrast — and therefore depends on the weighting choice of Section 4.2, differs
for every contrast reported, and must be rebuilt whenever the question changes.
That is precisely the kind of bookkeeping a package should carry and a person
should not.

---

## Edit 3 — Section 6, rewritten

## 6. Random effects in the partially nested design

The structural zeros do not disappear in a mixed model. They constrain the random
side too, and more severely than the fixed side, because the constraint is one of
expressiveness rather than of care: the chain parameterization cannot state every
structure the design admits, and neither the formula syntax nor the fitting
engines give a clear signal when an inexpressible one is attempted.

### 6.1 The closure question

Section 1 posed it: a participant's deviations form a vector of one number per
realized condition, so the general participant-level structure is an arbitrary
10 × 10 covariance, and every named structure is a constraint on it. Can the
chain form express them all?

The chain form's random basis has sixteen columns of rank ten, six columns being
identically zero. An unconstrained covariance on that basis has six unidentified
directions; identifiability is regained only by retreating to additive grouping
terms,

```r
(1 | participant) + (1 | participant:chord_type) + (1 | participant:cell)
```

whose implied covariance has one variance per level and equal correlation between
any two cells sharing a stratum — compound symmetry. That is a proper subset of
the admissible matrices, so the answer is no.

The cell form answers yes, trivially. `cell` is a complete basis for the realized
categorical structure, so `(0 + cell | participant)` with unrestricted Σ *is* the
general case, and diagonal, compound-symmetric, grouping-chain, and every `brms`
covariance structure are constrained submodels of it. On complete data it is an
exact reparameterization of the crossed random-slope form, verified to 10⁻⁹ in
log-likelihood; in a partially nested design it differs only by omitting the
combinations that do not exist.

One qualification keeps this from over-reaching. Terms on fully crossed variables
need no translation: every participant hears every chord type, so
`(1 + chord_type | participant)` is identified as written and passes through
unchanged. Translation is required only for terms spanning the structural
boundary — and such a term corresponds in cell space to a *reduced-rank*
covariance over the ten cells, which is representable, being a constrained 10 × 10
matrix. That is the closure property doing its work.

### 6.2 Why the restriction goes unnoticed

An unidentified fixed effect announces itself: R prints `NA` in the coefficient
table, and Section 3.2 traced what follows. The random side has no such signal.

Fitted on 2400 rows — six trials per participant per cell, generous by any
standard — `lmer` returns the sixteen-dimensional participant covariance, a
Hessian with twenty-one negative eigenvalues, and a convergence warning. Nothing
in that output states that the structure is not identifiable. A convergence
warning is the most familiar message in mixed modelling and is routinely read as
a sign of sparse data, met by simplifying the random structure until the warning
stops. The simplification usually succeeds, which is exactly the trouble: the
analyst concludes the data were thin, when the model in fact asked a question the
design cannot answer and would still be unable to answer with ten times as much
data.

### 6.3 What the cell parameterization restores, and what it costs

The gain is not that the richest structure is the right one. On the demonstration
data it is not: the unstructured covariance costs fifty-two additional parameters
for twenty-nine log-likelihood units, which a likelihood-ratio test does not
distinguish from chance (p = 0.221). The fitted correlations range from −0.319 to
0.453, but with forty participants that spread is consistent with sampling noise
around a common value, and compound symmetry is a defensible description of these
data.

The gain is that the question became askable. Under the chain form, "does this
design require an unstructured covariance, or will compound symmetry serve?"
cannot be posed, because only one of the two candidates can be fitted. Under the
cell form both can be fitted, compared, and reported. A decision formerly settled
by the parameterization returns to the data and to the analyst.

The cost is interpretive, and it falls on reporting rather than on fitting. Under
the grouping structure, participant heterogeneity in the chord-type effect is a
named variance component, read straight off the summary. Under the cell
parameterization there is no such component: Σ describes variation across cells,
and the heterogeneity of any *effect* is the quadratic form c′Σc, where c is the
same weight vector that defines the corresponding fixed-effect contrast. For the
equal-weight augmented-minus-major contrast, c places +1 on the single augmented
cell and −1/3 on each of the three major cells,

```
+1.00  0  0  0  0  0  0  −0.33  −0.33  −0.33
```

giving a participant standard deviation of 1.380 on that effect, against 0.909 on
the equal-weight diminished-minus-major effect. The quantities are exactly the
ones a report needs, and in `brms` each carries a full posterior. But c changes
with every contrast reported, and depends on the weighting policy chosen in
Section 4.2, so the vectors must be rebuilt whenever the question changes. Doing
this by hand for a table of contrasts is error-prone in a way that is difficult to
detect, since a mistaken c yields a plausible number rather than an error. It is
mechanical for software, and it is the reason the companion package treats the
random structure as one more translation: `nestimand`'s `random_terms()` renders
a declared structure as an unstructured cell covariance, as the grouping chain, or
unchanged, and reports the rank deficiency of any structural slope passed through
unaltered rather than leaving a convergence warning to be misread.

### 6.4 Retained material

*The following are kept from the current Section 6, following the above: the
sentinel-slope demonstration — `(1 + inversion | participant)` asks for a
participant-level slope on the sentinel level, fitted singular, AIC 2862.9
against 2842.1 for chord-type slopes and 2865.3 for a bare intercept; the
`participant:cell` grouping factor and its 300 levels; the single-trial
requirement, where `lmer` stops outright rather than returning a singular fit; and
the closing note that the fixed-effect summaries of Section 4 apply unchanged
under any random structure, with `re.form = NA` for frequentist fits and
`re_formula = NA` for `brms`. The prediction-mode paragraph — conditional,
finite-population, and integrated, with the Louis shrinkage result — is unchanged
and follows.*

*One sentence of the current text requires amendment. It states that "for a
partially nested design, the full covariance matrix was never going to be well
identified anyway". Over the realized cells it is identified — the fit above
estimates all fifty-five parameters — and the sentence should be narrowed to the
crossed slope form, where the claim is exactly right and is the substance of
Section 6.1.*

---

## Edit 4 — Section 7 checklist

The current mixed-model item reads:

> - [ ] For mixed models: does any random slope span the structural boundary?
>   Slopes on the nesting factor are fine; slopes on the partially nested factor
>   are not (Section 6).

Replace with:

> - [ ] For mixed models: does any random slope span the structural boundary?
>   Slopes on the nesting factor are fine; slopes on the partially nested factor
>   are not (Section 6.1). A convergence warning is not evidence of sparse data
>   here — check the rank of the random-effects design before simplifying
>   (`nestimand::random_slope_rank()`).
> - [ ] Is the random structure actually the one intended, or the one the chain
>   parameterization permits? Compound symmetry is the only structure the chain
>   form can identify; if an unstructured or diagonal covariance is wanted, the
>   cell parameterization is required (Section 6.1).
