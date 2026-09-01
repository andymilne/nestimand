# nestimand: implementation status

*Build 2026-08-15.11. Two kinds of claim are distinguished below: things that were
run in the development container and observed to work, and things that are
written but could not be run here and need running on a machine with a working
brms. Environment: R 4.3.3, marginaleffects 0.18.0, emmeans 1.10.0,
brms 2.20.4, lme4 1.1.35.1, ordinal 2023.12.4, posterior 1.5.0.*

## Why cells: the random-effects argument

*Corrected 2026-08-15.11. An earlier version of this note claimed that the chain
form cannot express an unstructured covariance at all. That is false, and the
test is recorded below: deleting the six uninformative columns and placing an
unrestricted covariance on the ten that remain reproduces the cell fit to
optimizer tolerance. The accurate claim is that the ten columns cannot be
produced by any formula, so reaching them means building them by hand, or
silencing six coefficients by prior and reading around an over-sized covariance
afterwards.*

The cell parameterization is often introduced through the fixed effects — full
rank, no aliasing, no identification constraints — but the decisive reason is the
random effects, and it is a reason of kind rather than degree. The chain form can
express only a subset of the random-effects structures R offers; the cell form
can express all of them.

**Under the chain form, a random slope over the nesting structure is
rank-deficient by construction.** The columns for impossible combinations are
identically zero for every group, so the covariance has more dimensions than any
amount of data can identify. On the demonstration data the structural random
slope has 16 columns of rank 10, three of them identically zero: six dimensions
of the covariance are unidentified, and no sample size changes that.

**No engine refuses it.** Fitting that structure in lme4 on 2400 rows with six
trials per participant per cell produces a 16-dimensional covariance, a
degenerate Hessian with 21 negative eigenvalues, and a convergence warning — the
kind of warning routinely attributed to sparse data and worked around by
simplifying the model. The structure is not sparse-data-limited; it is
unidentified in principle. `nestimand` now reports this explicitly, with the
column count, the rank, and the remedy, whenever such a structure is passed
through unchanged.

**So the chain form is confined *in formula terms* to nested grouping factors** —
`(1|p) + (1|p:chord) + (1|p:chord:inversion)` — which impose compound symmetry:
one variance per level, and equal correlations among cells sharing a stratum.

A correction to an earlier draft of this note, which claimed the demonstration
data exhibit a correlation pattern compound symmetry cannot represent. The fitted
correlations do range from −0.319 to 0.453, but the unstructured covariance costs
fifty-two parameters for twenty-nine log-likelihood units, which a
likelihood-ratio test does not distinguish from chance (p = 0.221). With forty
participants that spread is consistent with sampling noise around a common value.
The argument does not need the stronger claim and should not make it: the point
is that the chain form cannot express the alternatives, so the question of which
structure the data support cannot be asked at all. Under the cell form both
candidates can be fitted, compared, and reported, and on these data compound
symmetry is a defensible answer.

**The cell form places the whole categorical structure on one factor**, where the
entire menu is available and every dimension identified: unstructured
(`(0 + cell | p)`, 10 variances and 45 correlations), diagonal
(`(0 + cell || p)`), any brms covariance structure, or the grouping chain itself
as a constrained submodel via `random_structure = "chain"`. Choosing among them
becomes a modelling decision about the data rather than a limitation imposed by
the parameterization — which is the whole point of the translation architecture,
arrived at from the random-effects side.

## What is implemented

**`R/spec.R` – declaration surface.** Carried over from the prototype with its
user-facing form intact: bare or quoted `nests`, the NA refusal with the
`apply_sentinel` remedy, the continuous-parent refusal, the outcome-type and
saturation guards, `apply_sentinel` itself. What changed beneath it is that the
declared structure now compiles into a realized-cell factor rather than into a
chain of interaction terms. The cell factor is written into the returned data
under `cell_name` (default `"cell"`); a name collision with an existing column
is refused rather than overwritten.

**`R/translate.R` – the translation module.** Replaces `model_formula()`.
`cell_formula(spec, mode)` emits either parameterization; `chain_terms()`
retains the ancestry-closed chain form for the effect basis and the emmeans
engine; `add_cells()` recomputes the cell factor from the crossed original
factors, dropping unrealized combinations by default and refusing them on
request; `cell_grid()` builds the realized grid in original space;
`effect_basis()` returns A, QR-pivot-selected and rank-checked;
`translation()` bundles A and its inverse with the cell and effect names;
`fitting_mode()` states the basis choice and its reason.

**`R/policy.R` – the weighting policy.** One object: a distribution p over the
versions of a compound condition. `nest_policy()` builds p from an alias
(`equal`, `proportional` / `counterfactual`, `hierarchical`, `nominated`) or from
a supplied distribution; `counterfactual_grid()` forms the G-computation grid
with p attached as row weights; `policy_vertices()` enumerates the single-version
policies whose contrast range is the partial-identification region.

**`R/estimand.R` – the estimand function and its code view.** `estimand()` takes a
fitted model and a declared spec and returns a labelled estimand table in the
original variable names, with a policy, partial-identification bounds, a reorder
self-check, and provenance attached. The code view is guaranteed rather than
maintained: the function assembles the analysis as text and then evaluates that
text, so `show_code()` returns the object that was run, and cannot drift from it.
Non-core arguments in `...` are deparsed into the same text, so they reach the
destination function unaltered and appear in the saved code, which is what makes
the code re-usable and adaptable. `contrast = "within"` emits per-stratum
contrasts, which cross no boundary and so take no policy.

**`R/fit.R` – the fitting side.** `nest_fit()` dispatches to the declared engine,
assembling the call as text and evaluating it, so the fitted model carries the
call that made it; `estimand()` finds that code and joins the two into a single
pipeline in the code view. `...` reaches the engine unaltered, which is how
`chains`, `iter`, `seed`, and their kind are supplied. `dry_run = TRUE` returns
the code without fitting, which matters when the fit is a sampling run.
`random_terms()` translates a declared random structure: `cells` gives the
unstructured covariance over realized cells, `chain` the progressive grouping
submodel, `as_declared` passes the structure through untouched. A clean
intercept term is left alone.

**`R/latent.R` – estimands as linear functionals.** On a linear scale – the
latent scale of an ordinal model, the link scale of a generalized linear model –
a policy contrast is exactly c′b, where c is a weighted difference of
design-matrix rows formed on the counterfactual grid. Because the scale is
linear, averaging predictions over that grid and averaging its design-matrix rows
are the same operation, so this is G-computation rather than an approximation to
it. `latent_estimand()` returns the contrast with a delta-method interval;
`latent_draws()` applies the same c draw by draw to a `brms` posterior, which is
the Bayesian counterpart and the entry point for the derived-draws output.
`estimand(..., scale = "latent")` routes through it, bounds included.

## Run and observed to work

- The cell parameterization is full rank with no aliased coefficients, where the
  chain form aliases six; the two are the same fit (equal log-likelihood).
- The effect basis is square and full rank at depth one (10 × 10) and depth two
  (13 × 13); μ = A m reproduces the cell means, and A⁻¹ recovers the effect
  coefficients from them, to 10⁻¹⁰.
- `by =` on the non-predictor original factors returns estimands in original
  labels from a cell fit: aug − maj = −0.6779, matching the prototype.
- `equal`, `proportional`, `counterfactual`, and `hierarchical` all return
  −0.6779 at depth one; `hierarchical` returns −0.6283 at depth two, matching the
  prototype; `equal` and `hierarchical` differ at depth two, as the tree-versus-
  leaves distinction requires.
- A supplied p = (0.5, 0.3, 0.2) returns exactly the convex combination of the
  three single-version contrasts (agreement to 10⁻¹⁰).
- The partial-identification region on the demonstration data is −1.048 to
  −0.245, with the equal-weight point at −0.678.
- The estimand is unchanged under permutation of every nested factor's levels.
- Grids: unrealized combinations dropped (six of sixteen) or refused; a grid
  missing a nesting variable refused; recomputed cells agree with the fitted
  factor.
- The code returned by `show_code()` re-runs in a fresh environment and
  reproduces the estimand it came from.
- A non-core argument (`conf_level = 0.9`) reaches `marginaleffects`, changes the
  interval, and appears in the saved code.
- `contrast = "within"` returns three contrasts in each of the three strata in
  which inversion varies, with no sentinel contrast, and passes the reorder check.

- `nest_fit()` produces a full-rank cell fit whose code re-states the
  parameterization and the reason for it; `engine = "emmeans"` selects the
  effect basis instead, as designed.
- A declared random structure `(chord_type * inversion + training | participant)`
  translates to `(0 + cell + training | participant)` under `cells` and to the
  progressive grouping chain under `chain`.
- Ordinal engines: `clm` fits the threshold-aware coding at 15 parameters and no
  longer warns; the brms cumulative family is threaded into the emitted call;
  non-core arguments reach brms and appear in the saved code.

- The linear map agrees with the prediction route to 5 × 10⁻¹⁴ on estimates and
  1.6 × 10⁻⁹ on standard errors, and reproduces the convex-combination identity.
- On an ordinal `clm` fit the latent route returns one number per contrast, with
  bounds and a passing reorder check, where the prediction route refuses.

- A prior of independent cell means with sd 1.5 implies effect sds of 1.5√2 with
  non-zero covariances, and translating back recovers the stated sds exactly.
- The implied prior on the equal-weight aug − maj contrast is 1.5√(1 + 1/3), and
  on the within-family dim − maj contrast 1.5√(2/3), both to 10⁻¹⁰.
- A prior silent on a parameter is refused rather than quietly widened.

105 checks, all passing (`tests/test_core.R`).

## Findings from the fitting side

**Ordinal engines needed a different intercept convention.** A zero-intercept
cell coding is one parameter too many where thresholds already carry the
intercept: `clm()` warns and reinstates the intercept silently, and a brms
ordinal model would sample a weakly identified parameterization. The cell factor
under default contrasts is full rank either way – dummy coding of a single factor
cannot produce structural zeros – so the coding switches to `~ cell` for `clm`,
`clmm`, and the ordinal brms families, and no estimand changes.

**The clmm-by-marginaleffects support boundary is real**, as the design
document's checklist anticipated. On ordinal fits marginaleffects returns
per-category output even on the latent scale, and version 0.18.0 refuses the
pairwise comparison as too large. `estimand()` therefore requires the scale to be
stated for ordinal fits rather than assuming one, and offers the remedy that the
translation layer makes available: `scale = "latent"`, which computes the
contrast as a linear functional of the cell coefficients. That route is exact,
costs one matrix product, is immune to the per-category expansion, and – since it
never calls the prediction machinery – is not subject to the engine's guard
rails at all. It is now the recommended route for ordinal fits, and it doubles as
the mechanism the derived-draws output will need.

**`R/priors.R` – prior translation.** A prior is stated where it is thought – on
named effects, or on cell means – and fitted where estimation is well posed. For
an elliptical prior the translation is exact and needs no reparameterization: if
the effects carry m ~ N(m0, D) then the cell means carry mu ~ N(A m0, A D A'),
and Stan accepts a correlated prior directly, so the whole thing reduces to two
data blocks and one prior statement. `prior_audit()` reports what a stated prior
implies in the other space; `prior_for_estimand()` reports what it implies for a
named contrast, which is the form a report actually states.

## Engine version differences, found by running on a second machine

The development container has marginaleffects 0.18.0; a target machine running
0.32.0 exposed three changes at once, all of which the package now absorbs.

1. **The `hypothesis` argument.** The string form (`"pairwise"`) was accepted up
   to 0.18.x and the formula form (`~pairwise`) is required from 0.19.0. The
   package emits whichever the installed version accepts, in the computation and
   in the code view alike.
2. **The label column.** The contrast label lived in `term` and now lives in
   `hypothesis`. The label is located by content rather than by name, and a
   `term` column is guaranteed on every estimand returned.
3. **The direction of the contrast, which reversed.** 0.18.0 returns
   `aug - maj = -0.6779`; 0.32.0 returns `(maj) - (aug) = +0.6779`. This one is
   not cosmetic: inheriting the engine's convention would mean the same analysis
   reporting opposite-signed effects on two machines, with nothing to alert
   either user. The package therefore fixes the convention itself – contrasts
   run in declared factor-level order, earlier level minus later – negating the
   estimate, the statistic, and the swapped confidence bounds where the engine
   ran the difference the other way. The number of reversed rows is recorded on
   the result.

The general lesson is recorded here because it will recur: the estimand layer
must not inherit presentation conventions from whichever engine version happens
to be installed. Anything a user would report is normalized by the package.

## Verified on the target machine (R 4.6.1, marginaleffects 0.32.0, cmdstanr)

`tests/test_brms.R`, eight checks, all passing.

- The cell parameterization samples cleanly in brms, eleven population-level
  coefficients, and the fitted call travels with the model.
- **The draw-wise translation agrees with the linear map exactly.** Posterior
  mean −0.6788 with sd 0.2076, against a delta-method point estimate of −0.6788
  with se 0.2076. The same functional computed two ways.
- **The translated prior samples and reproduces its audit table.** Drawing from
  mu ~ N(A m, A D A') with no data gives sd 1.644 on the across-boundary
  contrast against 1.732 stated, and 1.150 on a within-family contrast against
  1.225, each centred where stated. This is the empirical confirmation that the
  prior translation is correct.
- The latent route is unaffected by the presence of random effects.

Two findings came out of that run. brms names cell coefficients with the
punctuation intact (`b_cellaug.none`), where the draw-matching rule had stripped
it; matching is now exact with a punctuation-insensitive fallback. And a brms
prior of `class = "b"` covers every population-level coefficient, covariates
included, so a prior stated only on the cells left Stan with mismatched
dimensions — reported merely as chains finishing unexpectedly, with nothing to
indicate the cause. The prior now spans the full vector, and the dimension is
checked before sampling begins.

### The prior scale, checked and cleared

A first run of the prior check, at 1000 draws, returned sampled standard
deviations about five per cent below the stated ones for both contrasts tested,
in the same direction — enough to be worth ruling out, since a systematic loss of
scale in A D A' would invalidate the prior route. Repeating at 10,000 draws
returned a ratio of 0.986 on the across-boundary contrast and ratios of
0.969, 1.001, 0.986, 0.997, 0.980, 1.000 across all six, scattering either side
of one with no pattern by contrast type. A discrepancy that shrinks as the draw
count rises is Monte Carlo noise; a systematic error would have persisted. The
residual spread slightly exceeds the naive Monte Carlo error of 1/sqrt(2n),
which indicates an effective sample size below the nominal draw count; the check
now reports that figure where the `posterior` package is available.

The prior check runs 10,000 draws and a five per cent tolerance for this reason:
comparing a sampled standard deviation against a stated one is a low-precision
comparison, and a tolerance loose enough to accommodate noise at 1000 draws would
be too loose to detect a real error.

## Not run here: needs a machine that can compile Stan

`tests/test_brms.R` collects four checks that require sampling. The container
compiles Stan too slowly to finish one within a session, so the brms path is
verified here only as far as the emitted call, the family threading, the prior
objects, and the argument passthrough. The most consequential of the four is the
second: `latent_draws()` assumes brms names cell coefficients `b_` followed by
the cell label with non-alphanumeric characters removed, and cell labels contain
full stops, so that assumption is the one most likely to be wrong.

## Branching families and depth (added 2026-08-31, run and observed to work)

A declared family is now a tree rather than a chain: one parent may hold several
nested variables. `chain_of()` became `family_of()`, a depth-first walk, and the
family vector is no longer ancestry — `nest_ancestors()` is. Four places had been
reading a position in that vector as ancestry, and each is now asked: the sentinel
search, `degenerate_strata()`, the `within` strata, and `spec_nests()`. The
declaration accepts `c(inversion, X1) %in% chord_type` and the
`inversion + X1 %in% chord_type` sugar (carried as text, since `%in%` binds tighter
than `+`, so R's parse tree does not match the reading); a variable declared inside
two parents, and `%in%` chained within one entry, are refused with the remedy named.

Two bugs the one-level chain never exposed, both fixed and covered by
`tests/test_core.R`: `degenerate_strata()` reports a composite stratum key, and
three call sites restricted the cell table by its first variable only, so at depth
three the pooled estimand of the deepest variable matched no cells and failed inside
`combn()`; and `add_bounds()` and `estimand_values()` called `latent_estimand()`
without `cells`, which errored at depth three and, at depth two, silently computed
the bounds on the unrestricted cell table.

Run here on a four-variable design — `chord_type > {inversion, X1}`,
`inversion > Z` — with `lm`: 19 to 41 realized cells, effect basis square and full
rank, no aliased coefficients, every target's estimand and its emitted script
agreeing, and the reorder self-check passing at depth three. The whole suite is 487
checks, 0 failures.

`random_terms(structure = "chain")` builds its rungs from ancestor paths rather than
from prefixes of the family vector. The two coincide on a chain - `p:a`, `p:a:b`,
`p:a:b:c`, unchanged - but where a parent holds two children the prefix form made one
of them the coarser division and the other the finer one, which the design does not
say and which the order of the declaration then decided: declaring `inversion` first
gave `p:chord:inversion` as the intermediate rung, declaring `X1` first gave
`p:chord:X1`. Siblings now enter crossed, one variance each, with the whole structure
as the finest rung:

    (1|p) + (1|p:chord_type) + (1|p:chord_type:inversion) + (1|p:chord_type:X1)
          + (1|p:chord_type:X1:inversion)

One variance component more than the ladder it replaces, and no longer a function of
how the declaration was typed - the rungs are ordered by depth and then by name, so
the same design emits the same formula however it was written. Run here on an `lmer`
fit: five components, one per rung, and the `cells` form is unchanged.

### Hierarchical weighting under branching

`hierarchical_weights()` split the version *label* position by position, so with two
variables under one parent it read one of them as dividing the other, and the weights
depended on which was declared first: for versions `0.a, 0.b, 1.a, 2.a` it gave
(1/6, 1/6, 1/3, 1/3) reading inversion first and (1/6, 1/2, 1/6, 1/6) reading X1
first. The fix is the same one as elsewhere - condition on a variable's own ancestors
rather than on its position - so the function now works from the cell table:

    w(cell) proportional to prod over v of 1 / n(v | v's ancestors in that cell)

which is the generative reading of the policy, one uniform draw per node of the
declared structure. Siblings are then independent choices whose probabilities
multiply; the weights are renormalized over the realized versions, as everywhere
else. Consequences, all checked: on a chain the numbers are unchanged (the depth-one
anchor -0.6779 and the depth-two anchor -0.6283 both hold); over a set of siblings
hierarchical coincides with equal, which is right, since with independent draws and a
complete crossing every leaf has the same probability; and it still separates from
equal wherever one variable is nested inside another. Working from the cell table
also removes a latent fragility - the old version re-split dot-joined labels, which a
level containing a full stop would have broken.

### A random slope cost a variable its fixed effect (fixed 2026-08-31)

`nesting_spec()` took the covariates from `all.vars(formula[[3]])` and then
subtracted every variable named in a bar. The intention was to drop the grouping
factors; the effect was to drop the slope variables too. Declaring

    response ~ chord_type * inversion * top * GMSI + (chord_type * inversion * top | id)

with `inversion %in% chord_type` therefore fitted `~ cell + cell:GMSI` - `top` gone
from the mean structure entirely, while keeping a random slope on it, and nothing
said. The covariates are now read from the fixed term labels, which already exclude
the bars, so a bar can no longer remove anything from the fixed side.

The same declaration exposed a second reduction, on the random side.
`random_terms()` drops the structural terms of a bar as subsumed by the cell factor,
which is right for `chord_type:inversion` but not for `chord_type:top`: `cell` alone
carries no slope that varies by condition. A covariate declared crossed with the
structure now keeps the crossing - `(0 + cell + cell:top | id)` rather than
`(0 + cell + top | id)` - the same distinction `cov_by_cell` draws on the fixed side.
One left additive in the bar stays additive. The compound-symmetric submodel
(`structure = "chain"`) crosses nothing, by design.

### Two more from the same session (fixed 2026-08-31)

`estimand(m, a * b * c)` refused. `a * b * c` parses as `(a * b) * c`, and the
operands were read off the top call alone, so `a * b` stood as though it were a
variable name. They are now gathered through the nesting.

`nest_summary()` left a *factor* covariate crossed with the cells untranslated, and
called every column of it a threshold. The columns are named `cell<k>:x<level>` -
the variable plus a level suffix - and the block matcher required the name to end at
the variable, so it found none of them; they fell to the leftover block, whose label
was `ifelse(name %in% covariates, "common slope", "threshold")`. Each level is now
its own block, matched allowing the suffix, longest covariate name first so a name
that is a prefix of another does not swallow its columns. The leftover label no
longer guesses: only `Intercept[k]` and `a|b` are called thresholds, and anything
else says it was left as fitted. A wrong label is worse than an unhelpful one - these
rows were being read as ordinal cut points.

### A design variable that is nested in nothing (added 2026-08-31)

The declaration surface had no way to say "this factor is part of the categorical
design and is crossed with the rest of it". Such a variable was left out of `nests`,
which made it a covariate: it entered the fit as `cell:x`, which is the right mean
structure, but it was not in the cell factor and so could not be the target of an
estimand - there were no conditions to weight over. The user's case was
`inversion_top_notes`, defined at every level of chord type and inversion.

Two changes. A categorical variable the formula crosses with the declared structure
is now folded into the cell factor automatically - nothing to declare, since the
formula already said the conditions are the combinations of all three. Folding
changes no fitted model: `cell + cell:x` and a cell factor over the enlarged design
span the same space, verified by comparing projections. Where a combination is
unrealized it is strictly better - on a ragged case the interaction form gave 20
columns of rank 17, three of them aliased, against 17 clean cells. Left as
covariates: a variable entering additively, whose effect is thereby declared common
to every condition; a numeric one; and an ordered factor, whose contrasts say it is
meant as a quantity.

An entry with no `%in%` also names such a variable explicitly - `c("inversion %in%
chord_type", "inversion_top_notes")`, or unquoted - which is how an additive or
ordered one is folded in deliberately. It becomes a one-variable family, joins
`cell_vars`, and everything downstream follows without further change - 20 realized
cells rather than 10, a square full-rank effect basis, `estimand(m, top)` with
policies and bounds, and the reorder check passing once `spec_nests()` was taught to
restate the crossed declarations (it rebuilt the spec from the parent map alone,
which would have dropped them). A numeric variable is refused there, since a numeric
variable crossed with the structure is exactly what a covariate already is.

The estimand error for an undeclared target now says this rather than only refusing:
it names the variable, says it is entering as a covariate, and prints the `nests`
vector that would make it a target.

### Interaction contrasts beyond two variables (fixed 2026-08-31)

`interaction_matrix()` read `lv[[1]]` and `lv[[2]]` and looked its corners up in a key
built from *all* the variables it was given. With three, nothing matched, and the
failure came out as "no two levels of `chord_type` share two levels of `inversion`" -
false of the design, and pointing at the wrong thing. `estimand(m, a * b * c)` reaches
this, since `*` returns each variable on its own and their interaction.

`*` also returned only the whole-set interaction, where a formula crosses: it now
returns every interaction among the targets, so `a * b * c` gives seven results - the
three variables, the three pairs, and the triple - in the order a formula writes them.

`interaction_matrix()` now takes two levels of every variable and forms the product of
the simple contrasts over the 2^n corners they define, which is a difference of
differences at n = 2 and a difference of those at n = 3. The corner set must be realized, which is
what excludes the augmented stratum here: nine three-way contrasts among dim, min and
maj. Two-variable output is unchanged, in values, labels and order. A design with no
realized corner set now says so and says which variables it could not span; more than
500 contrasts is refused with the alternatives named, since the count grows as the
product of the pair counts.

### An interaction over some of the design variables averaged nothing

Reached as soon as the cell factor held a variable the interaction did not name -
which the automatic fold made common. `latent_contrast_matrix()` built its comparison
matrix from the raw cell table, so when several cells shared one combination of the
targets, `match()` took the first and the contrast was formed at one arbitrary level
of the remaining variables. On a demonstration design with a large three-way effect
the interaction of chord type and inversion came out 0.515 where the answer is 3.278,
and the prediction route agreed with it, since on an identity link it takes the same
shortcut through the coefficients. The rows are now averaged within each combination
of the targets, as the prediction route does by averaging predictions, and the two
routes agree with the value computed by hand from the cell means. Two more
`deg$vars[1]` stratum keys were fixed in the same function.

### The sentinel reached an interaction contrast

An interaction is computed only over the strata in which every one of its variables
varies; the code took the restriction of `target[length(target)]` alone, on the
reasoning that in a chain the deepest target's restriction implies the others'. That
stopped being true as soon as a target could sit outside the chain. With
`inversion:inversion_top_notes` the last target is the crossed variable, whose
restriction is empty, so nothing was restricted and the augmented stratum stayed in:
the output offered `(none - 0) x (t2 - t1)`, comparing the absence of an inversion
against root position, which is the comparison the package exists to refuse.

`degenerate_strata_multi()` intersects the restrictions of all the targets and
returns them keyed over the union of their ancestors, so the existing plumbing - the
emitted `subset()`, the restriction note - is unchanged. `inversion:top` now gives
three contrasts and no sentinel; `chord_type:top` keeps the augmented stratum, which
is right, since both variables vary there.

### A p-value on a Bayesian coefficient table

`nest_summary()` had no Bayesian branch. On a brms fit it took `fixef()` and
`vcov()`, applied the linear map, and reported `2 * pnorm(-|est / se|)` - a Wald test
against a normal approximation to the posterior, which answers a question about
repeated sampling that the model was not fitted to ask, computed from a Gaussian
stand-in for a posterior whose draws were sitting right there. The estimand side had
known this since it was written; the coefficient table had simply never been given
the same treatment.

It now maps the draws through the same matrix and summarizes each row's posterior:
mean, standard deviation, quantile interval, and probability of direction. No
p-value, no test statistic. A row the parameterization holds at zero - the intercept
under an ordinal family - gets `pd = NA` rather than a meaningless 0. The print method
shows `pd` where a frequentist summary shows `p.value`. Frequentist fits are
untouched, which the suite checks.

`pd` is monotone in the Wald p-value for a symmetric posterior, so this does not add
information so much as stop mislabelling it. A ROPE percentage or a Bayes factor
would answer a different question, and both need a decision about what counts as
negligible on the latent scale; neither is implemented.

### A bracketed random slope lost its grouping factor

The declared random terms were split out of `random_original` with a regex that
matched a parenthesised group, `\\([^()]*\\)`. That cannot match a bar whose left
side is itself bracketed - `(chord_type * (inversion + top) | participant)`, which is
how a slope over a set of variables gets written - because the character class stops
at the inner parenthesis. The regex matched the inner group instead, so the term
arrived as `inversion + top` with no bar at all, `grp` became NA, and the partial
structure check refused a term the user had not written: "the random term
`(inversion + inversion_top_notes | NA)` varies with ...".

`bar_terms_of()` reads them from the parsed expression instead - flattening `+`,
peeling `(` and any covariance wrapper, then taking the two sides of `|` or `||`.
Used by `random_terms()`, `chain_random_zeros()` and `grouping_vars()`, all three of
which had their own variant of the same regex. Checked on a bracketed left side, two
terms bracketed differently, a `diag()` wrapper, and a nested grouping factor
`school:class`; the translation is identical however the left side is bracketed.

### The effect basis followed the formula, which the fit does not

`chain_terms()` built the effect basis from the declared term labels. The cell
parameterization fits `~ 0 + cell` whatever the formula says - `nesting_spec()`
already messaged as much - so a formula that restricts the structure left the basis
smaller than the thing it is meant to reparameterize:
`chord_type * (inversion + inversion_top_notes)` spans 14 of 20 realized cells, and
`effect_basis()` refused a perfectly well formed spec with "not square and full rank
(20 cells, 14 identified effects)". Reached through the prior translation, so a
declared `class = "b"` prior on a restricted formula could not be fitted at all.

The basis is now the saturated one over the realized cells, independent of the
formula: every subset of the cell variables that is closed under ancestry, which is
the same list the fully crossed formula used to produce. What the declared structure
still decides is which terms may appear - a variable is named only alongside its own
ancestors, so no term is a marginal effect the design cannot support. A chain still
gives its prefix ladder, and more than twelve categorical design variables is refused
rather than enumerated.

### A restricted mean structure is now fitted as written

Until now the declared formula did not decide the mean structure: cells mode fits
`~ 0 + cell`, which is saturated, so `chord_type * (inversion + top)` - a model
without the three-way interaction, and a reasonable thing to want when the data are
thin - was silently enlarged to the saturated one and only a message said so.

Two functions now do what one did. `chain_terms()` is the translation basis, the
saturated ancestry-closed set over the realized cells; `declared_terms()` is the
structure the formula asks for, closed under ancestry so that a nested variable is
always given its parent - `inversion` alone is not a term this design admits and
becomes `chord_type:inversion`. `cell_formula(mode = "effects")` fits the second, and
`fitting_mode()` chooses it whenever the declared structure spans fewer dimensions
than there are cells, since the cell factor cannot express a restriction.

Checked on a design of 20 realized cells with a declaration spanning 14: the fitted
values and residual degrees of freedom match a hand-written restricted `lm`, the 12
columns the chain form carries that the data cannot inform are dropped as aliased
rather than fitted, estimands agree between the prediction and coefficient routes,
and `nest_summary()` reports the coefficients as fitted rather than refusing to
translate a fit that needs no translation. `chain_priors()` derives 6 structural
zeros and 6 identification constraints for the same design, which is what brms needs
to fit it.

### The constant(0) block is now derived rather than requested

A restricted structure is fitted in the effects parameterization, which carries
coefficients the data cannot inform: columns for conditions that do not exist, and one
redundant column per stratum. `lm` drops them as aliased; brms cannot, and samples an
improper posterior along them. The package knew this and said so - "chain_priors(spec)
derives the constant(0) block; pass it as priors =" - which leaves the default fit
wrong and the correct one an opt-in.

`nest_fit()` now derives and applies the block itself on brms in effects mode, and
reports what it did: how many coefficients, how many are structural zeros (conditions
the design does not realize) and how many identification constraints (a coding choice,
like a reference level, which leaves every estimand unchanged). A `class = "b"` prior
supplied by the user becomes the prior on the remaining coefficients rather than being
displaced, and appears once in the call; supplying `priors` replaces the block
entirely. The emitted code reads `chain_prior_object(chain_priors(sp))`, so it stands
on its own rather than referring to an object only the fitting environment held - and
it carries the arguments the block was derived with, which the first version did not.
`nest_fit()` runs the code it emits, so re-deriving with the default regularizer put a
second `class = "b"` prior beside the user's and brms refused the whole call as
duplicated. It now emits `chain_priors(sp, regularize = NULL)` where the user has
supplied one. Checked by evaluating the emitted expression and handing the result to
`brms::validate_prior()`, which accepts it, with and without a prior of the user's.

### Two counts of two different things, neither saying which

`nesting_spec()` reported "spans 14 of the 20 realized cells" and `nest_fit()` then
reported 24 coefficients held at zero, and nothing said the two were counting
different things. They are consistent: 14 of 20 is the rank of the declared *mean
structure* over the cell means, while 24 is a count of *design columns* - 12 in the
mean structure and 12 in the random structure for one grouping factor, each structural
term being crossed with the covariates it interacts with, so one condition that does
not exist can carry several columns. On this design the fixed 12 are 6 structural
zeros and 6 identification constraints, and the random side repeats them.

Both messages now say what they count. The fit message breaks the total down by where
the coefficients sit - the mean structure, then each grouping factor by name - and
says outright that it counts coefficients rather than cells; the spec message says it
counts cell means, and that the covariates and random terms multiply them.

### Why these kept surfacing, and what was done about it

Four faults in a row had one shape: a design matrix built one way and coefficients
fitted another. The prior block re-derived with different arguments than the object it
stood for; an interaction matrix built over the raw cell table against predictions
aggregated over targets; design rows built in the cell parameterization for a fit in
the effects one; and finally design rows built from the data as it stands for a fit
whose data had been releveled.

Two things let them through. The effects parameterization was until now a niche path,
exercised almost entirely by dry runs - the code was checked as text rather than run -
while the cell path had 500 executed checks. And every test used demonstration data in
which the sentinel is already the first factor level, so the relevel `nest_fit()`
applies for an effects fit changed nothing, and a grid coded without it agreed by
accident. Real data ordered `0, 1, 2, none` does not.

The remedy is an invariant rather than four more instances. `tests/test_core.R` now
runs, for both positions of the sentinel in the level order: a restricted declaration
(fitted as effects) and a saturated one (fitted as cells and, forced, as effects);
every estimand route on each; the requirement that the two parameterizations of one
model give the same numbers; and the requirement that a design built for a fit carries
no column that fit lacks. Ten checks, and any of the four faults above would have
failed one of them.

### The target expression, expanded by R rather than by hand

The same shape again, in the one place the invariant added above does not reach: the
target was read by a hand-rolled walk over the parse tree. It handled `a * b`, then
`a * b * c` once taught to recurse, and then failed on `a * (b + c)` - the brackets a
reader writes for exactly the structure this package is about - passing
`(inversion + inversion_top_notes)` on as though it were a variable name.

R expands a right-hand side already. The target is now handed to `terms()` and the
results are its term labels: `a * (b + c)` gives a, b, c, a:b, a:c and no b:c, which
is what a formula means by it, and `a:b` alone stays a single interaction. Six checks
compare the returned names against `terms()` over five bracketings, which is the
invariant rather than another instance.

### A grid over one stratum had no contrasts to code

Asking for an estimand within one level of something - the top-note effect within
each chord type - means building a design over a slice of the cells, and a factor with
one level in that slice has no contrasts, so `model.matrix()` refused. The columns
have to be the ones the fit has whatever the grid contains, so `design_rows()` now
gives every factor the levels the fit saw rather than the levels the slice happens to
hold. Same fault line as the relevel: a design built from what is in front of it
rather than from what the fit was built from. Two more checks on the invariant, over
both sentinel positions and both parameterizations.

With that, a per-stratum estimand can be assembled from exported functions:

    per <- do.call(rbind, lapply(levels(dat$chord_type), function(k) {
      cs  <- sp$cells[sp$cells$chord_type == k, , drop = FALSE]
      pol <- nest_policy(sp, "inversion_top_notes", "proportional", cells = cs)
      e   <- latent_estimand(m, "inversion_top_notes", pol, spec = sp, cells = cs)
      cbind(chord_type = k, as.data.frame(e))
    }))

`estimand()` now has `by` for it, in the sense marginaleffects gives the word: the
target's contrasts are formed inside each group, with the policy weighted over that
group's conditions alone, and the groups come back stacked and labelled. It could not
simply be passed through - `estimand()` already gives `by = target` to
`avg_predictions()`, and the coefficient route makes no such call - so it is
implemented as a restriction of the cells, one group at a time, over the same
machinery that restricts them when a nested target does not vary everywhere. Both
routes, every scale, several grouping variables, bounds per group, and one runnable
script per group from `show_code()`.

`contrast = "within"` is the same estimand with the strata taken from a nested
target's ancestors, and the two agree where both apply - checked, since they are
computed by different paths. What `by` adds is naming the strata, which is what a
variable crossed with the structure needs, having no ancestors to take them from. A
group in which the target does not vary is dropped and said; a covariate cannot group
an estimand, having no conditions; nor can the target group itself.

### An unknown argument was reported as somebody else's

Anything in `...` that the coefficient route could not use was described as "an
argument of marginaleffects::avg_predictions", whether or not it was one.
`method = "within"` - a reasonable guess at how to ask for per-stratum contrasts -
drew a message sending the reader to the documentation of a function that has no such
argument. The two cases are now told apart: an argument the prediction function does
have is named as one and explained; an argument neither has says so, with a
suggestion where the name is close to a real one (`contrst` draws `contrast`, while
`method` draws nothing rather than a spurious guess). Either way the message ends by
naming the vocabulary that does exist - `contrast`, including `contrast = "within"`,
and `by`.

## Not yet implemented

The prior translation and audit
functions, the prediction modes, the four typed outputs, the random-effects cell
structure, and the emmeans effect-basis path.

## Decisions raised by the implementation

1. **`hierarchical` is absent from the design document's policy table** but is
   retained here, since it is a distinct p – uniform at each level of the chain
   rather than uniform over the leaves – and the depth-two prototype anchor
   (−0.6283) is stated in terms of it. The table should either gain the row or
   the alias should be retired with its anchor restated.
2. **Incomplete coverage of a supplied p is refused rather than renormalized.**
   Restricting p to the named versions and renormalizing is a different
   intervention, and choosing it silently would violate the no-silent-defaults
   rule. A degenerate stratum, having no choice to make, takes its single
   realized version at mass one.
3. **Settled: functions, plus a code view.** The package computes through
   functions, and the code that a function ran is retrievable with `show_code()`
   for saving, adapting, and re-use. Emitted code calls `nestimand` rather than
   inlining cell tables and weight vectors as literals, on the ground that a
   reader of the code can be expected to have the package installed. The reorder
   self-check accordingly runs inside `estimand()` rather than being pasted into
   a script, and is reported in the printed output.
4. **The reorder check refits the model**, and so is skipped with a stated note
   for `brms` fits, where a refit doubles sampling time. Under the cell
   parameterization order instability is impossible by construction, so the check
   now verifies the translation layer rather than the fit.

## One parameterization: the reduced design

A declaration that asks for less than the saturated structure could not be
written on the cell factor - the factor is saturated by construction - so it was
fitted on the original factors, in the effects parameterization. That form
carries columns the data cannot inform: `lm` drops them as aliased, and brms
cannot, so every such fit needed a derived block of `constant(0)` priors. Two
parameterizations, two ways of building a design, and each of the last few faults
was the same shape - a design matrix built one way against coefficients fitted
another.

The restriction can be carried per cell after all. The declared structure's
design over the realized cells, with the identically-zero columns and then the
redundant ones removed, is full rank by construction: same fitted values, same
residual degrees of freedom, nothing aliased, nothing held at zero, on any
engine. It is the same kind of object the cell factor is - one row per realized
condition, computed once and looked up by cell - and the cell factor is its
saturated special case. So there is now one construction rather than two:

- `reduced_design(spec)` is the matrix, its `effect_names` attribute naming the
  effect each column stands for; `with_reduced()` carries the columns alongside
  the data, and `counterfactual_grid()` recomputes them for every grid, since a
  column copied across with a row would still describe the cell the row came
  from.
- Column names are syntactic (`z_` prefix, `:` written `.`) so that every engine
  carries them unaltered; `nest_summary()` reports the effect, not the column.
- A covariate crossed with the structure gets one slope per column *and* one for
  the intercept the design does not carry - the reference condition. Without that
  term the crossing is one slope short of the structure it is crossed with.
- `mode = "effects"` remains, but nothing reaches it by default: only
  `engine = "emmeans"`, which is formula-driven, and priors stated
  coordinate-wise in effect space. The `constant(0)` block belongs to that form
  alone.

The claim is checked rather than assumed: for both positions of the sentinel in
the level order, the reduced fit must drop nothing, must have the same fitted
values and residual degrees of freedom as the effects fit of the same
declaration, and must give the same estimands on both routes.

## `target` is checked where every route to a policy passes

A `nesting_spec` passed in the `target` slot - these functions take the model
first and the target second, since a fit from `nest_fit()` carries its own
declaration - failed several frames later inside `versions_of()`, on
`if (target %in% f)`, an error naming neither argument. The check now sits in
`versions_of()` itself, which every route to a policy runs through, and
`latent_draws()` runs it before the engine check so that a misplaced argument is
reported as one rather than as whatever engine happens to be in hand.

Two defects surfaced beside it, both in `latent_draws()` and both of the same
shape - a second implementation of something `latent_estimand()` already did:

- it passed `weights = weights` to `policy_contrast_matrix()` without having a
  `weights` argument of its own, so any path that forced the promise failed on a
  missing object. `weights` and `route` are now arguments, as they are on
  `latent_estimand()`.
- it left `cells` unresolved and passed `NULL` on, where `latent_estimand()`
  dropped the strata in which the target does not vary. `estimand_cells()` now
  makes that choice for both, so the two routes are over the same cells.

## The brms test file had gone stale against its own package

`tests/test_brms.R` called `latent_draws(mb, spb, "chord_type", "equal")` and
three others in the same shape: the argument order these functions had before
the fit began carrying its own declaration. Nothing caught it, because the file
cannot run anywhere Stan cannot compile, and the rest of the suite had long since
been updated. The calls are now `f(model, "chord_type", "equal", spec = spb)`,
which is what the current signature takes.

This is the standing gap rather than a fault in the code: `test_core.R` runs
everywhere and is the reason the frequentist paths hold, while `test_brms.R`
runs only on a machine with a working Stan toolchain, so it drifts silently
between runs. It is worth running after any change to the latent route.

## `dev/check_calls.R`: what can be checked without running anything

Three separate faults in `tests/test_brms.R` were found by a user running it
rather than by anything here - a spec passed positionally into `target` in five
places, a `weights` that was not an argument, and contrast labels
(`aug - maj`) from a superseded direction convention. The file only runs where
Stan compiles, so it drifts against the signatures while the rest of the suite
keeps pace.

`dev/check_calls.R` closes part of that gap statically. It parses every file in
`tests/` and `dev/`, matches each call to a nestimand function against that
function's formals, and reports a named argument the function does not have, a
call with more positional arguments than free formals, and - the useful one -
a positional argument that would land in the wrong slot. It knows which of a
file's own symbols hold a declaration and which hold a fit, because it can see
what they are assigned, so a `nesting_spec` reaching any formal but `spec`, or a
fit reaching any but `model`, is reported with its file and line. That is what
caught `estimand(mb, spb, chord_type, ...)`, which R would not have complained
about at all: it would have taken `chord_type` as `at` and failed elsewhere.

Argument *order* in general is beyond static analysis, and stale string
literals - a contrast label from an old convention - are beyond it entirely.
`check_target()` catches the important case of the first at runtime. Run
`check_calls.R` beside `check_docs.R` after any signature change.
