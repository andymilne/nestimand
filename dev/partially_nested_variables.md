# Partially Nested Variables: Modelling and Estimation in R

**Purpose.** This document and its companion software divide one problem between them. The document serves *manual* analysis: it explains why partially nested variables break default modelling practice and develops interpretable, hand-codeable remedies – the chain parameterization, with `emmeans` and `marginaleffects` recipes – accepting one substantive limit as the price of parameters a reader can interpret directly: the chain form cannot express every random-effects structure the design admits (Section 1, "A second fact about the design", and Section 6). The companion package `nestimand` serves *automated* analysis and inverts the trade: it fits the fully robust but uninterpretable realized-cell parameterization, and its principal purpose is translation – from the original, intuitive, research-question-pertinent variable and parameter space into the cell space, where estimation is unconditionally well posed; inference is performed there; and results are translated back into the original space in which questions, priors, and reports are stated. The Appendix documents an earlier chain-based prototype of `nestimand`; the cell-translation design is specified in the accompanying `nestimand_design.md`.

---

## 1. What a partially nested variable is

Some designs contain a factor whose levels exist only within certain levels of another factor. The emptiness of the remaining cells may be a matter of logic – the combination is meaningless, so no observation *could* arise there – or a matter of protocol – the combination is perfectly coherent but was, deliberately, never implemented. The canonical music example is of the first kind:

- **Chord type**: major, minor, diminished, augmented
- **Inversion**: root position (0), first inversion (1), second inversion (2)

Inversion is meaningless for augmented triads: the augmented triad is symmetric under inversion, so there is no augmented-chord-in-second-inversion condition to run. Other logically empty cases include bowing technique (strings only), key or mode (musical stimuli only), and pitch accuracy (sung responses only). The by-design kind is just as common: a manipulation applied in only some arms of an experiment, a control condition that crosses with no treatment sub-factor, or – as in Section 5.1 – a factor varied within one stratum although it could have been varied within others. For everything in this document the two kinds are equivalent, because what matters statistically is that the empty cells are empty *with certainty, as a fixed feature of the design*; the kinds differ only in interpretation, since an unimplemented cell has an in-principle population value that the design cannot identify, whereas an impossible cell has none. Both are distinct from *missing data* – values that exist but were not observed – for which none of what follows is a substitute. The inversion factor is **partially nested** within chord type: its levels exist only inside some levels of the nesting factor. Each level of the nesting factor defines a **stratum** – a slice of the design within which the nested factor either varies or is undefined. Here each chord type is a stratum: the major stratum contains three inversion conditions, while the augmented stratum contains a single condition with no inversion at all. A comparison among inversions stays *within* a stratum; a comparison between augmented and any other chord type crosses the **structural boundary** – out of the strata in which inversion is defined, into the one in which it is not – and only these across-boundary comparisons force the weighting decision below. Nothing in the definition requires the nested variable to be categorical, or the nesting to stop at one level: a *continuous* variable can be defined only within some strata (a vibrato rate that exists only for stimuli with vibrato), and a nested variable can itself contain further nested variables, to any depth. Sections 1–4 develop the one-level categorical case; Section 5 shows that the deeper and continuous cases follow the same rules. This is a statement about *fixed effects* – which conditions exist – and is distinct from the "nesting" of grouping factors (trials within participants) familiar from multilevel modelling.

**A note on terminology.** The labels in the literature do not line up neatly, and the relevant work is scattered across several of them. The clinical-trials and psychotherapy-methods literature uses *partially nested design* for a trial in which clustering (therapy groups, shared clinicians) exists in one arm only. That is not a competing meaning so much as the same idea applied to a different kind of variable: there, a *random grouping factor* exists only within some levels of another factor; here, a *fixed factor* does. This document in effect generalizes the term from the random side to the fixed side – and onward to continuous variables and nesting within nesting (Section 5), with the random side rejoining in Section 6. In both usages the structure itself is a design property – which cells exist, and where clustering occurs, are facts no analysis can alter – while everything done about it is a modelling decision; the trials literature accordingly uses the term for both, speaking of partially nested *models* as readily as designs. The fixed-factor structure also travels under other names: it is an *incomplete factorial design* in the clinical-trials sense (some treatment combinations deliberately or necessarily unimplemented – the classic case being a control condition that cannot be crossed with the treatment factors), a factorial with *structurally empty cells* or *structural zeros* in the linear-models and categorical-data literatures, and *nested fixed factors* in the `emmeans` documentation, which is the operational label, since `nesting = "inversion %in% chord_type"` is how the structure is declared in code. Strictly, it is not classical nesting either – in a classically nested design each level of the nested factor occurs within exactly one level of the nesting factor, whereas here inversion is *crossed* with chord type over three of its levels and undefined for the fourth. When searching the literature, "incomplete factorial" and "empty cells" find the relevant work; "partially nested" mostly finds the clustering design. (On forums the problem also circulates informally as "nested variables" or "nested predictors".) This document keeps *nesting* as its umbrella term: it is the most intuitive of the candidates, it is the word the software uses, it is continuous with the trials-literature usage rather than in conflict with it, and it is the only candidate that extends gracefully to continuous nested variables and to nesting within nesting.

The following simulated dataset is used throughout (40 observations per cell, a participant identifier, and a continuous covariate `training`):

```r
library(emmeans)
library(marginaleffects)

set.seed(1); n <- 40
cells <- data.frame(
  chord_type = c(rep(c("dim", "min", "maj"), each = 3), "aug"),
  inversion  = c(rep(c("0", "1", "2"), times = 3), "none"),
  mu         = c(4.16, 4.00, 3.73,   # dim, inversions 0/1/2
                 4.29, 4.20, 3.71,   # min
                 4.98, 4.62, 4.39,   # maj
                 3.92))              # aug
dat <- cells[rep(seq_len(nrow(cells)), each = n), ]
dat$participant <- factor(rep(seq_len(n), times = nrow(cells)))
tr <- runif(n, 0, 10)
dat$training <- tr[as.integer(dat$participant)]
dat$response <- rnorm(nrow(dat), dat$mu + 0.12 * dat$training, 1.2)
dat$chord_type <- factor(dat$chord_type, levels = c("aug", "dim", "min", "maj"))
dat$inversion  <- factor(dat$inversion,  levels = c("none", "0", "1", "2"))
```

Note that the undefined level is coded as an explicit sentinel, `"none"`, not as `NA` – Section 3 explains why. The first thing to do with any such design is to count the cells:

```r
table(dat$chord_type, dat$inversion)
```

```
      none  0  1  2
  aug   40  0  0  0
  dim    0 40 40 40
  min    0 40 40 40
  maj    0 40 40 40
```

**Ten non-empty cells, not sixteen.** The saturated mean structure has exactly ten parameters, and any analysis that implicitly requires more is asking for something the data cannot supply.

### The unavoidable choice

Suppose the aim is to compare augmented with major chords. Augmented is *one* condition; major is *three*, and their cell means differ (in this simulation: 5.73, 5.42, and 4.93 for inversions 0, 1, and 2, against 4.68 for augmented). Before the comparison can be made, "the major chord mean" must be defined as a single number:

- the **unweighted average** of the three inversions;
- the average **weighted by presentation frequency**;
- major **at one nominated inversion** only.

None of these is wrong; they are different questions. The choice is made in every analysis – the only issue is whether it is made deliberately (Section 4) or left, unexamined, to the software (Sections 2 and 3).

### A second fact about the design

The choice above concerns fixed effects. The same ten-cell structure has a second consequence, for *random* effects, and this one is not a choice but an obstruction.

A mixed model treats some factor as a sample from a population – participants here, but equally stimuli, classrooms, or sessions. Its units may respond differently from condition to condition, and how much freedom that variation is given is a modelling decision, written in the formula:

```r
(1 | participant)                                 # compound symmetry: one offset
                                                  #   per unit
(1 | participant) + (1 | participant:chord_type)  # compound symmetry: offsets
                                                  #   also varying by chord type
(chord_type | participant)                        # unstructured over chord type
(chord_type * inversion | participant)            # unstructured over all ten
                                                  #   conditions: each its own
                                                  #   variance, each pair its own
                                                  #   correlation
```

These four formulas request four different *covariance structures* – the constraints placed on the matrix of variances and correlations among a unit's deviations. The first three are ordinary requests and are estimated without difficulty. The fourth is not, and the reason has nothing to do with sample size: it asks for a deviation on every combination of chord type and inversion, including the six that do not exist, and a deviation for a condition that never occurred is not estimable from any amount of data. Rewriting it in the chain form that Section 4 recommends for the fixed effects – `(chord_type + chord_type:inversion | participant)` – does not help. Both forms produce sixteen columns of rank ten; they differ only in how the deficiency is arranged. `lme4` fits the request and warns only of convergence trouble, which reads like ordinary sparse data. `brms` can hold the six absent conditions at zero, through a prior specification naming each in turn. The fit is then correct, but the absent conditions remain in the model: on the demonstration design, 776 parameters are sampled where 455 suffice, and the summary reports 120 correlations of which 45 mean anything. Section 6.2 sets out what that costs. Neither engine offers a way to leave them out.

**The obstruction lies in the notation, not in the experiment.** The ten conditions support every one of these covariance structures, including the least constrained; what fails is the attempt to write it in terms of two factors, six of whose combinations were never presented. Indexing the deviations by the ten conditions themselves removes the difficulty: same model, same fit, every parameter estimable. Section 6 gives the formula and the parameter count for each of the three standard covariance structures – unstructured, diagonal, and compound symmetry. Of these, only compound symmetry is directly available in a partially nested design. The other two are reached either by rewriting the random term and supplying priors that hold the absent conditions at zero, which `brms` alone permits, or by reparameterizing the model onto the ten conditions; Section 6 sets out the grounds on which `nestimand` chooses.

This is the strongest motivation for the companion package. The fixed-effect difficulties of Sections 3 and 4 are surmountable by hand. This one is not a matter of care but of what the notation reaches, and the translation it demands is mechanical for software and unreasonable for a person.

---

## 2. Why splitting the analysis is deficient

The usual first instinct is two models:

```r
# Model A: inversion effects – augmented chords excluded
nonaug <- droplevels(subset(dat, chord_type != "aug"))
mA <- lm(response ~ chord_type + chord_type:inversion + training, data = nonaug)

# Model B: chord-type effects – inversion omitted
mB <- lm(response ~ chord_type + training, data = dat)
```

Both fit cleanly, and if inversion is genuinely the only question of interest, Model A alone is defensible. The deficiencies appear when Model B is used for chord type:

**Model B removes the weighting choice from the analyst's control.** A natural objection: if inference proceeds through predictions – `avg_comparisons` or `emmeans` – rather than raw coefficients, does the weighting problem not evaporate? It does not, and the reason is instructive. For a linear model, prediction-based summaries of Model B return exactly its coefficient. On an unbalanced variant of the data (root position presented 60 times, the others 20 each), the coefficient, `avg_comparisons`, and `emmeans` all give maj − aug = 0.6864, precisely the count-weighted average of the major cell means, where the equal-weight answer is 0.5347. The weighting was fixed when the model was fitted: omitting `inversion` marginalizes over it with observed-frequency weights, and prediction-based tools can only reweight over variables the model contains. The usual lever is disconnected: `emmeans(mB, ~ chord_type, weights = "equal")` still returns the count-weighted mean (4.9615, not the equal-weighted 4.8098), because there is no inversion dimension in the reference grid for `weights` to act on. With the full model, by contrast, the same choice is made at the marginalization step (Section 4.2), where it remains open to be made, changed, and stated.

**Model B inflates the residual error.** Systematic variation between inversions is unexplained, so it lands in the residual. To keep this comparison strictly like-for-like, compare Model B without the covariate (`response ~ chord_type`) against a model that differs in exactly one respect – it adds the inversion structure, `response ~ chord_type + chord_type:inversion` – so that any difference is attributable to modelling inversion and to nothing else. The residual SD is 1.246 against 1.216 on the balanced data, and 1.288 against 1.252 on the unbalanced variant. For the identical estimand – maj − aug under observed-frequency weighting – both models return the identical point estimate (0.6779 balanced; 0.6864 unbalanced), but the two-factor model with the smaller standard error (0.2221 against 0.2274 balanced; 0.2341 against 0.2409 unbalanced). The gains are modest here because the simulated inversion effects are modest relative to the noise; they scale with the amount of within-stratum variation Model B leaves unexplained.

**Shared quantities are estimated twice.** The covariate effect is fitted once per model – 0.1157 in the full model against 0.1160 in Model A – and the two are quietly allowed to disagree. In a mixed model the same applies to each participant's random effect, so there is no coherent single estimate of between-participant variance; in an ordinal model each fit has its own thresholds, so the two latent scales cannot be compared at all.

**What a single model provides.** Taken together: Model B answers one particular weighting question, unchangeably and less precisely than a two-factor model answers the same question, and the split as a whole duplicates estimates that should be shared. The alternative is a single model with the same two-factor structure as Model A, but fitted to all 400 rows – augmented chords included – rather than to the 360 non-augmented rows. A reader who has attempted this may protest that it is impossible: with inversion coded as `NA` for augmented chords, R silently deletes those 40 rows, and the "full" model is Model A in disguise. It is possible – it requires coding the undefined level as an explicit sentinel, as in the data setup above, for reasons Section 3.1 explains. Once fitted, including the augmented data changes none of the within-stratum inversion contrasts: the maj inversion-0 − inversion-1 contrast is 0.3074, with matching standard errors, whether the augmented rows are in the model or not. (Model B plays no part in this comparison – having no inversion terms, it can estimate no inversion contrast at all.) The single model thus matches the defensible half of the split exactly, replaces the indefensible half, and does so while sharing one residual variance, one covariate estimate, and one set of participant effects, and while making the cross-boundary weighting available to choose and state. That two-factor model, introduced in Section 3.2, is the vehicle for everything in Section 4. One caution before proceeding: *fitting it is not the remedy*. Nothing in this section licenses running the model and reading its output – its coefficient table is actively misleading (Section 3.2), and so are the default marginal summaries built on it (Sections 3.3–3.5). The model is the right object; the rest of this document is about the right summaries of it.

---

## 3. Intuitive solutions that do not work as expected

### 3.1 Coding the undefined level as `NA`

```r
dat$inv_na <- dat$inversion
dat$inv_na[dat$chord_type == "aug"] <- NA
m_na <- lm(response ~ chord_type * inv_na + training, data = dat)
nobs(m_na)   # 360, not 400
```

R's default `na.action = na.omit` performs casewise deletion: the moment `inv_na` enters the formula, every augmented row is silently removed: `nobs` drops to 360 and `aug` is absent from the model frame – with **no warning of any kind**, because the remaining design is rectangular and full rank. The estimated "effect of chord type" simply omits one of its levels. The blame lies upstream of any modelling package: the deletion happens in base R's model-preparation machinery – `stats::model.frame`, applied with the global `na.action` option, whose default is `na.omit` – which `lm`, `lmer`, and `clm` all call. The rows are gone before estimation begins: `nrow(model.frame(response ~ chord_type * inv_na, dat))` is already 360. The one partial exception is `brms`, which prepares its data through its own validation: it excludes the same rows, but at least emits a generic warning ("Rows containing NAs were excluded from the model") – easy to wave through, since the missingness looks intentional. No choice of package offers protection, because the problem sits at the data level – and so does the fix: an explicit sentinel level (`"none"`, `"undefined"`) in place of `NA`.

### 3.2 Fitting the two-factor model anyway

```r
m <- lm(response ~ chord_type + chord_type:inversion + training, data = dat)
```

The formula is written in *chain form* – `chord_type + chord_type:inversion` rather than the crossed `chord_type * inversion` – for reasons Section 4 sets out; everything in this section afflicts both forms. The formula requests 16 mean-structure parameters where only 10 are estimable. (The chain form is no leaner here than the crossed form: with no inversion main effects present, its interaction expands over every chord type, reference included, giving 1 + 3 + 12 = 16 columns – the same total the crossed form reaches as 1 + 3 + 3 + 9. The two allocate the sixteen differently, which is Section 4.5's subject, but they request the same number and share the same rank of 10.) The fit does not stop: the pivoted decomposition inside the fitting routine detects the rank of 10, estimates the ten coefficients its pivot selects, and marks the six beyond the rank `NA`, with one quiet note about singularities. The other engines fit anyway too, but announce it differently: `glm()` and `clm()` likewise mark the six coefficients `NA`; `lmer()` and `glmer()` drop the columns outright, leaving a coefficient vector six entries shorter and a one-line message; `brms` retains all sixteen (Section 4.5). Those six are **aliased** – each of their design-matrix columns is a linear combination of columns already in the model, so the data cannot attribute an effect to them separately, and R's pivoting drops them. (The term comes from experimental design, where effects that share a design-matrix column are said to be aliased with one another; R's `alias()` function will list them.) The six fall into two kinds: three are the augmented-interaction columns, which are all-zero because those cells do not exist, and three are one interaction column per non-augmented chord type, chosen by the pivot – a partition Section 4.5 returns to. The counts generalize: a design in which several strata lack the nested variable contributes one all-zero triple per degenerate stratum, and one pivot choice per stratum that has it. The fit itself is correct – but **the labels are not**. Comparing against the raw cell means:

| Coefficient | Label implies | Actually estimates |
|---|---|---|
| `(Intercept)` | grand mean or baseline | aug cell mean (4.6799) |
| `chord_typemaj` | main effect of maj | maj **at inversion 2** − aug (0.2452) |

Which cell gets absorbed depends on which columns R's pivoting happens to drop – determined by the ordering of the factor levels: reordering the inversion levels changes `chord_typemaj` from 0.2452 to 0.7405, the same label now naming maj at inversion 1 minus aug. A significant `chord_typemaj` is a statement about one inversion, not about major chords, and nothing in the output says so.

### 3.3 The reorder test

Refit with the levels of the nested factor permuted. The fit must be unchanged; anything that changes was never identified.

```r
dat$inv_b <- factor(dat$inversion, levels = c("none", "2", "0", "1"))
m_b <- lm(response ~ chord_type + chord_type:inv_b + training, data = dat)
logLik(m); logLik(m_b)   # identical: -626.114 in both
```

Now ask each model for the "average effect of chord type" with `marginaleffects` defaults:

```r
avg_comparisons(m,   variables = list(chord_type = "pairwise"))
avg_comparisons(m_b, variables = list(chord_type = "pairwise"))
```

The results:

| Contrast | Default order | Reordered |
|---|---|---|
| dim − aug | −0.084 (p = .697) | −0.037 (p = .863) |
| min − aug | −0.079 (p = .712) | −0.004 (p = .984) |
| maj − aug | 0.635 (p = .003) | 0.684 (p = .001) |
| min − dim | 0.004 (p = .978) | 0.033 (p = .830) |

Same model, same likelihood, same predictions at every real cell – and every estimate moves with the order of the factor labels: maj − aug by eight per cent of its value, min − dim by a factor of eight. The shifts are modest only by the accident of this section's parameterization; under the fit-identical crossed form `chord_type * inversion` – the formula many analysts type first – a second contamination route opens (Section 3.4) and the same table becomes disastrous: dim − aug flips sign, from −0.473 (p = .063) to +0.019 (p = .940), and maj − aug moves from null to significant (0.245, p = .350, to 0.740, p = .005). Order-dependence of any size is disqualifying, because nothing about the data changed.

The failure is not universal across tools. The same request put to `emmeans` – `pairs(emmeans(m, ~ chord_type))` – returns identical estimates under both orderings (aug − maj is −0.6779 either way, agreeing to 10⁻¹⁵), because it detects the nesting from the chain formula and computes only quantities the realized cells support; under the crossed form, where detection fails, it refuses outright (`nonEst`) rather than producing a movable number. Why one tool is immune and the other is not is the subject of the next subsection.

### 3.4 Why this happens

`avg_comparisons()` works by **G-computation**: it copies the dataset, sets `chord_type` to each level in turn, predicts, and averages the differences. Setting `chord_type = "dim"` on the 40 augmented rows requests a prediction for *diminished at the sentinel level* – a cell that does not exist. The model obliges using whichever aliased column pivoting dropped, and the answer moves when the aliasing moves. Under the chain form this is the *only* contamination route: the reverse substitution is harmless, because an augmented row's prediction involves no inversion terms at all – the fitted augmented value is 4.6631 at every inversion – which is why the shifts in the table above stay below 0.05 (40 contaminated rows in 400, with partial cancellation between the two sides of each within-boundary contrast). Under the crossed form the reverse substitution is contaminated too – `chord_type = "aug"` on a non-augmented row requests the augmented chord at a real inversion, itself unrealized – and with 360 further contaminated predictions per contrast the shifts grow tenfold, producing the sign flip quoted above. This asymmetry between the two fit-identical parameterizations is one of the reasons the chain form is this document's standard (Section 4). The general lesson: **G-computation over a rectangular grid assumes the grid is realizable.** With structural zeros it is not, and every marginal-effects tool inherits that assumption by default. Any method that never requests a prediction for an unrealized cell is safe; any method that does is not.

This is the mechanism behind the contrast in behaviour just seen. `emmeans` never walks onto the rectangle at all under the chain form: it detects the nesting from the formula (printing "A nesting structure was detected: `inversion %in% chord_type`") and computes only realized-cell quantities, so even an undeclared `emmeans(m, ~ chord_type)` returns correct averages, identical under every factor ordering. Under the crossed form, where detection fails, it falls back to refusing – `nonEst` for anything that depends on an unrealized cell. `avg_comparisons()` complies under either form: it returns a plausible-looking number that changes with the factor ordering, exactly as in the table above. (Under the hood, its predictions come from `predict()`, which does warn – "prediction from rank-deficient fit" – but the warning never reaches the `marginaleffects` output.) Correct by default, or failing safe, in one package; a quiet wrong answer in the other. This shapes the next section: `emmeans` needs at most a documenting declaration (Section 4.1), whereas `marginaleffects` needs its prediction grid *rebuilt by hand* to exclude the impossible rows (Section 4.2).

### 3.5 Three more things that look like fixes but are not

All three are demonstrated on the dataset above.

**Restricting `newdata` to real rows does nothing.** `avg_comparisons(m, ..., newdata = unique(dat[, c("chord_type", "inversion", "training")]))` returns min − dim = 0.0043, identical to the unrestricted call and equally order-dependent (0.004 under one ordering, 0.033 under the other). `newdata` restricts only the *starting* rows; the level substitution then recreates every impossible cell.

**Releveling does nothing.** Making a real inversion the reference (`relevel(dat$inversion, ref = "0")`) merely changes which columns drop: the same contrast now returns 0.0298 – a third different answer, not a stable one.

**`by = "chord_type"` mixes stable and unstable rows.** `avg_comparisons(m, variables = list(inversion = "pairwise"), by = "chord_type")` returns, for major chords, the genuine within-stratum contrasts (stable under reordering) *interleaved with* contrasts against the sentinel such as `mean(0) − mean(none)`, which returned 0.803 under one ordering and 0.307 under the other – each being a real contrast wearing a wrong label. Nothing in the table distinguishes the two kinds of row.

---

## 4. Correct strategies

One principle covers everything: **estimates come only from realized cells, and the weighting is made explicit.** In practice that means either declaring the nesting to `emmeans`, or restricting the prediction grid in `marginaleffects` so that no impossible cell is ever requested. (A third route – collapsing the design into a single ten-level cell factor – also removes the rank deficiency, but the factors then no longer exist as named variables, so every contrast must be hand-indexed and the substantive meaning of the design is lost; that route is not pursued here.)

A word on parameterization before the recipes, because this document standardizes on one for manual analysis: the **chain form**, `chord_type + chord_type:inversion` (extended as `+ chord_type:inversion:doubling` and so on at depth). The crossed form `chord_type * inversion` fits identically – same likelihood, same fitted values, and every estimand in this section agrees between the two to machine precision – but the chain form has four advantages the equivalence does not touch: each structurally impossible cell receives a dedicated coefficient of its own, which Section 4.5 requires (under the crossed form no such coefficients exist at all – Section 4.5 explains the mechanics – so `brms` rejects the priors); `emmeans` detects the nesting from it without a declaration; the augmented predictions are constant across inversions, closing one of Section 3.4's two contamination routes; and `lm`'s own pivot then drops exactly the three structural and three conventional coefficients that Section 4.5 distinguishes. Its one drawback is unfamiliarity – chain formulas are awkward to write by hand at depth – a burden the Appendix's generator removes by writing them, to any depth, covariates included. Reference levels, reassuringly, need no careful choosing for *validity*: the safe recipes of this section, the augmented-prediction constancy, and `emmeans`' auto-detection are all unaffected by the choice, even when a real inversion, rather than the sentinel, is made the nested factor's reference. What placing the sentinel first secures is *legibility*: the structural zeros are then exactly the all-zero columns named `stratum:level`, one triple per degenerate stratum, which is the form Section 4.5's declaration assumes; with a real level as reference the same information smears across a different set of zero and linearly dependent columns, and the prior recipe no longer applies verbatim. This is why the data setup in Section 1 lists the sentinel first. Keep the model from Section 3.2 throughout:

```r
m <- lm(response ~ chord_type + chord_type:inversion + training, data = dat)
```


**An alternative parameterization: cell means.** Fitting the realized cells directly – `response ~ 0 + cell + training`, with `cell` the factor of realized chord-type–inversion combinations – removes the pathologies of Section 3 at their source. The design matrix is full rank by construction, no coefficient is ever aliased, no identification constraints (and, in `brms`, no `constant(0)` priors) are required, and reordering can change nothing because no order-dependent pivot exists; requests that would mislabel across the structural boundary fail loudly rather than silently, because the nested variable is no longer a predictor. These are considerable advantages. The cost is manual: every research question, prior, and label must then be translated by hand between the cell space and the variables in which the questions are actually posed, and formula-driven tools – notably `emmeans` – no longer see those variables. This guide therefore proceeds with the chain form, which keeps parameters interpretable and both estimation packages usable for manual analysis, while noting that the cell parameterization is the natural foundation for software. The principal purpose of `nestimand`, in its redesigned form, is precisely that translation: from the original, intuitive, research-pertinent variable and parameter space into the realized-cell space, where estimation is unconditionally well posed; inference is performed there; and the results are translated back into the original space in which questions, priors, and reports are stated.

### 4.1 Within-stratum contrasts

```r
en <- emmeans(m, ~ inversion | chord_type)   # the chain formula itself
pairs(en)                                    # declares the nesting
```

The chain formula *is* the declaration: `emmeans` reads `inversion %in% chord_type` out of the term structure (printing a note to that effect), so no `nesting` argument is needed – and none is needed at greater depth or for continuous leaves either, where the whole chain is detected the same way (Section 5). The call produces only contrasts among inversions that exist within each chord type – no sentinel rows at all – and every estimate is identical under factor reordering (for major: 0.307, 0.803, 0.495, matching the raw cell differences).

The `marginaleffects` analogue is *not* the `by = "chord_type"` call – Section 3.5 showed that its output mixes genuine contrasts with sentinel ones, and restricting `newdata` does not remove them, because the sentinel level is reintroduced by the level substitution itself, not by the rows supplied. Instead, stay with predictions on observed rows, one stratum at a time:

```r
avg_predictions(m, newdata = subset(dat, chord_type == "maj"),
                by = "inversion", hypothesis = ~ pairwise)
```

No substitution occurs, so nothing can extrapolate: this returns exactly the `emmeans` values above, with matching standard errors, and is stable under factor reordering. (When comparing its output across refits, match contrasts by label, not by row position – pairwise tables are ordered by factor-level order, so reordering permutes the rows and flips some signs, in both packages.)

An explicit `nesting = "inversion %in% chord_type"` argument reproduces the identical grid and results, and adding it anyway serves as in-code documentation of the design; it becomes necessary only under the crossed parameterization, where detection fails and everything returns `nonEst`.

**Pooled contrasts of the nested variable.** The within-stratum contrasts above have a pooled counterpart – "the inversion effect, averaged over chord types" – and it requires one rule: **compute it only over the strata in which the nested variable varies.** Degenerate strata must be excluded, because there the variable holds only its sentinel level, and any contrast against that level is a comparison *between strata* wearing the nested variable's label – on the demonstration data, a naive `by = "inversion"` over all rows reports "none − 0" = −0.500, which is nothing about inversions: it is the augmented chords versus the pooled root-position triads. Restricting to the strata where inversion varies removes such rows from the request rather than the output:

```r
tri <- subset(dat, chord_type != "aug")
avg_predictions(m, newdata = tri, by = "inversion", hypothesis = ~ pairwise)
```

giving 0 − 1 = 0.163, 0 − 2 = 0.734, 1 − 2 = 0.571 – observed-frequency-pooled over the three triad types, stable under factor reordering to 10⁻¹⁵, and free of sentinel rows by construction. The exclusion generalizes structurally, with no reference to the sentinel's label: a degenerate stratum is any combination of the nesting variables in which the nested variable takes a single level. `emmeans` offers no pooled marginal of a nested variable – in its nesting semantics the levels are meaningful only within their strata, which is a defensible position; the per-stratum contrasts above are its answer – so this estimand is a `marginaleffects` construction, and the strata pooled over should be stated when it is reported.

### 4.2 Across-boundary contrasts: an explicit choice of weights

Each weighting scheme is available in both packages.

**Equal weights** – each inversion counts equally within its chord type. In `emmeans`:

```r
pairs(emmeans(en, ~ chord_type))
```

In `marginaleffects`, the same estimand comes from averaging predictions over a grid containing one row per realized cell:

```r
g_eq <- datagrid(model = m, chord_type = levels(dat$chord_type),
                 inversion = levels(dat$inversion))
g_eq <- subset(g_eq, (chord_type == "aug" & inversion == "none") |
                     (chord_type != "aug" & inversion != "none"))
avg_predictions(m, newdata = g_eq, by = "chord_type", hypothesis = ~ pairwise)
```

**Proportional weights** – inversions weighted by their observed frequency. In `emmeans`:

```r
pairs(emmeans(en, ~ chord_type, weights = "proportional"))
```

In `marginaleffects`, no grid is needed at all – average the model's predictions over the observed rows themselves:

```r
avg_predictions(m, newdata = dat, by = "chord_type", hypothesis = ~ pairwise)
```

All four calls evaluate the model only at cells that exist, and all are stable under factor reordering. On the balanced data every one returns aug − maj = −0.6779, equal to the hand-computed value; on the unbalanced variant the two weightings separate exactly as they should – equal −0.5347, proportional −0.6864 – with the two packages agreeing on each. (One caveat when covariates are in the model: the `emmeans` calls and the cell grid evaluate covariates at their means, whereas the observed-rows call carries each row's own covariate values. In a linear model with covariates balanced across conditions these coincide; the distinction becomes real in the situations the next route is built for.)

**G-computation over the counterfactual grid.** The recipes above only ever evaluate the model where data exist. The counterfactual route instead asks, for every observed row, what the model predicts with that row's chord type set to each level in turn – build the substituted grid first, *then* delete the impossible combinations the substitution creates, then average:

```r
cf <- datagrid(model = m, newdata = dat,
               chord_type = levels(dat$chord_type),
               grid_type  = "counterfactual")
cf <- subset(cf, (chord_type == "aug" & inversion == "none") |
                 (chord_type != "aug" & inversion != "none"))
avg_predictions(m, newdata = cf, by = "chord_type", hypothesis = ~ pairwise)
```

This is the safe counterpart of the failed `newdata` fix in Section 3.5: because the grid is built before subsetting, the impossible pairings created by the substitution can be removed *after* they arise. This gives aug − maj = −0.6779, and all six contrasts agree between the two factor orderings to within 3.6 × 10⁻¹⁵ – machine precision.

Three properties define what this route does and does not add:

1. **It adds no new weighting.** Each chord type's inversions still enter at their observed frequencies, so its point estimate reproduces the proportional-weight answer above exactly (−0.6864 on the unbalanced data). The contribution of the additional construction lies in the next two properties.
2. **The observed covariate distribution is preserved.** The grid rows keep their real `training` values, so every substituted chord-type level is averaged over the same covariate distribution: the mean of `training` is 5.1452 for all four levels, identical to the sample mean. The contrast is standardized to a common covariate population by construction, rather than evaluated at a single covariate value.
3. **This matters beyond the linear model.** For a linear model with no covariate-by-factor interaction, averaging over the covariates and evaluating at their mean coincide. For anything nonlinear they do not: in a logistic version of this analysis, G-computation over the observed covariates gives an aug − maj risk difference of −0.175, where `emmeans` evaluated at the mean of `training` gives −0.186. The G-computation number is the population-averaged (marginal) effect; the `emmeans` number is a conditional effect at one covariate value. Both are legitimate – but they are different estimands, and only one of them respects the covariate distribution actually observed.

One subtlety: in the recipe above, the augmented predictions are averaged over the augmented rows' covariate values and the other chord types' predictions over their own rows'. In a within-subjects design where every participant contributes to every condition this is exactly balanced (as the identical `training` means above show). If covariates differ across conditions, the fixed-population variant evaluates all four chord types over one and the same set of rows:

```r
cfv <- datagrid(model = m, newdata = subset(dat, chord_type != "aug"),
                chord_type = levels(dat$chord_type),
                grid_type  = "counterfactual")
cfv$inversion[cfv$chord_type == "aug"] <- "none"
avg_predictions(m, newdata = cfv, by = "chord_type", hypothesis = ~ pairwise)
```

On these data it gives the identical result (−0.6779), identical covariate means for all four levels, and is stable under reordering.

**Covariates may themselves interact with the nested structure.** A training effect that differs by chord type enters the covariate down the same chain as the means:

```r
m_mod <- lm(response ~ chord_type + chord_type:inversion +
              training + chord_type:training, data = dat)
```

with a further term, `+ chord_type:inversion:training`, if the slope also varies by inversion. The covariate chain mirrors the mean structure, so the impossible covariate slopes are all-zero columns just as the impossible means are, and the `brms` declaration of Section 4.5 extends to them mechanically. Nothing about the recipes changes. On a variant of the demonstration data in which the training slope differs by chord type, both packages return the identical averaged contrast, stable under both factor orderings (aug − maj = −1.6468): `avg_predictions` on the restricted grid exactly as above, or `pairs(emmeans(m_mod, ~ chord_type, weights = "proportional"))`. The moderation itself is reported by evaluating the nominated-cell comparison of Section 4.4 at chosen covariate values:

```r
g2 <- datagrid(model = m_mod, chord_type = c("aug", "maj"),
               inversion = c("none", "0"), training = 2)
g2 <- subset(g2, (chord_type == "aug" & inversion == "none") |
                 (chord_type == "maj" & inversion == "0"))
avg_predictions(m_mod, newdata = g2, by = "chord_type", hypothesis = ~ pairwise)

eg <- emmeans(m_mod, ~ chord_type * inversion,
              at = list(training = 2), nesting = NULL)
pairs(eg[keep])   # keep: the two nominated rows, found as in Section 4.4
```

Both return maj-at-inversion-0 − aug = 1.0688 at training = 2 and 2.0481 at training = 8.


### 4.3 The causal status of an across-boundary contrast

A causal contrast requires each of its arms to name a well-defined intervention applied to a common set of units (Rubin, 1980). Across a structural boundary, one arm fails that requirement. "Set the chord to augmented" is complete: exactly one such stimulus exists. "Set the chord to major" is not: it names three distinct stimuli, one per inversion. The quantity E[Y(major)] therefore has no value until a distribution over those versions is specified. This is the problem of *multiple versions of treatment*, and with it a failure of the consistency assumption in its usual form (Cole and Frangakis, 2009; VanderWeele, 2009; VanderWeele and Hernán, 2013; Hernán, 2016).

The deficiency is definitional rather than statistical, which has a practical consequence: covariate adjustment does not repair it. Adjustment addresses a different threat, namely that the units observed in each stratum may differ in composition, and randomization, adjustment, and G-computation handle that threat completely. Neither they nor any other estimation technique can supply a meaning for an intervention that names three stimuli at once.

**The problem is universal; its unsolvability here is not.** Versions do not stop at inversion: timbre, loudness, register, and duration are also glossed over by "set the chord to augmented", and at a fine enough grain every treatment in every science is compound (Hernán, 2016, makes the point with water, whose versions include drinking and drowning; Holland, 1986, draws the general moral). Causal claims are therefore always relative to a declared grain of description, and consistency is required only up to *treatment-variation irrelevance*: the undeclared versions must either not affect the outcome or be drawn from a common distribution across the arms being compared (VanderWeele, 2009; VanderWeele and Hernán, 2013). Ordinary stimulus design supplies that condition in either of two ways. A version may be held fixed – one timbre, register, and duration used for every chord in the experiment – in which case its distribution is degenerate but identical across the arms; the contrast is causal, and its scope is the effect of chord type *at that timbre, register, and duration*, with no licence to generalize to any others. Alternatively the versions may be varied deliberately: every chord type presented in the same set of timbres, registers, and durations, so that the version distribution is non-degenerate and, again, common to every arm. The estimand is then the effect of chord type averaged over that stimulus distribution, and it generalizes as far as the distribution sampled – the substance of the stimulus-sampling tradition (Clark, 1973), whose statistical corollary is that the sampled versions belong in the model as random effects (Section 6) if the intervals are to reflect the variation across them. What breaks the condition is neither constancy nor variation but *asymmetry*: were augmented chords rendered on one instrument and major chords on another, timbre would differ systematically between the arms and the contrast would confound chord type with timbre. What distinguishes a partially nested variable is that the balancing move is unavailable for it, and the two kinds of emptiness distinguished in Section 1 make it unavailable in two different senses. Where the cells are empty by logic, balance is impossible *in principle*: the augmented stratum has no inversions to equate, and no experiment could supply them, so E[Y(augmented, inversion = 1)] has no referent at all. Where the cells are empty by protocol, the corresponding potential outcomes exist and could have been realized by a differently designed experiment; balance is impossible only *in this design*, which leaves them unidentified. The practical requirement is the same in both cases – the distribution over versions must be declared explicitly, since it does not cancel – but what the declaration acknowledges differs: in the first case that a common distribution does not exist, in the second that this design cannot inform one. The second case admits, in principle, a modelling assumption bridging the unimplemented cells; Section 3 shows what happens when such an assumption is made silently and by accident rather than deliberately and in the open. That declaration is the subject of the remainder of this section.

Set beside one another, the balanced versions and the boundary-crossing one differ in a single respect: whether the distribution cancels. Every estimand is defined relative to distributions over the versions it does not name, and where those distributions are common to both arms – timbre, register, duration – they cancel in the contrast and may be left implicit; re-weighting them to some other target distribution remains meaningful, and is what standardization to a population does, but it is optional, because the contrast is causal without it. Where the versions cannot be made common, whether in principle or by this design, the distribution no longer cancels: it enters the estimand itself, so specifying it is no longer optional. The construction below is therefore not a special device for boundary crossings but the general account, applied at the point where its choices stop being invisible.

**Contrasts of policies.** Let *p* be a distribution over the versions of a compound condition – over the inversions of a major triad. Then

E[Y(major, *I* ~ *p*)] = Σ_k *p*_k E[Y(major, inversion = k)]

is a well-defined causal quantity: the expected outcome under a *stochastic intervention*, which assigns the chord type and draws its version from *p* (Muñoz and van der Laan, 2012; Díaz and van der Laan, 2013; Kennedy, 2019). A contrast between two such policies – augmented, which has a single version, against major under *p* – is a causal effect in the ordinary sense, estimated by the G-computation recipe of Section 4.2 with weights *p* applied to the counterfactual grid (Robins, 1986; Hernán and Robins, 2020). The weighting is thus not an obstacle to causal interpretation but the specification of what is being intervened upon, and it is on that understanding that the reporting requirement of Section 4.2 rests.

**The weighting schemes are policies.** Each scheme in Section 4.2 is this construction at a particular *p*: equal weighting takes *p* uniform over the realized versions; proportional weighting takes *p* at the empirical version frequencies; a nominated level (Section 4.4) takes *p* at a vertex of the simplex. A distribution may equally come from outside the design – from corpus frequencies of inversions in the repertoire of interest – which standardizes the estimand to a target population and answers an applied question directly (Pearl and Bareinboim, 2014; Westreich et al., 2017). For the demonstration data, *p* = (0.5, 0.3, 0.2) over inversions 0, 1, and 2 gives aug − maj = −0.7952, the convex combination 0.5(−1.0479) + 0.3(−0.7405) + 0.2(−0.2452) of the three single-version contrasts.

**Bounds.** Every mixture estimand is a convex combination of the single-version contrasts, so the interval those contrasts span is the range of the causal effect across all admissible policies – a partial-identification region in the sense of Manski (1990). One additional line of output therefore establishes what no choice of weights can alter:

> Augmented triads are rated between 0.25 and 1.05 below major triads, depending on the inversion in which major is realized; under equal weighting over inversions, 0.68.

The point summary states the policy chosen; the interval states the range available to any other choice. Where a single-version question is acceptable, Section 4.4 is cleaner still, both arms then naming atomic interventions. Describing the stimuli instead by properties defined for every stimulus, such as roughness or harmonicity, removes the structural boundary from the exposure space altogether, at the price of answering a question about features rather than about categories; that route lies outside this document.

### 4.4 Comparison at a nominated level

Often the clearest report is not an average at all but "augmented versus major *in root position*". Two recipes deliver it, and the distinction between them matters. The first uses observed rows only, with no level substitution anywhere:

```r
sub0 <- subset(dat, (chord_type == "aug" & inversion == "none") |
                     inversion == "0")
avg_predictions(m, newdata = sub0, by = "chord_type", hypothesis = ~ pairwise)
```

This is a comparison between two *subpopulations* – the augmented trials and the root-position trials – and is therefore vulnerable to any covariate imbalance between them. The second evaluates every unit under both conditions, which is the G-computation form and the one that supports a causal reading (Section 4.2):

```r
g0 <- datagrid(model = m, newdata = dat,
               chord_type = levels(dat$chord_type),
               inversion = levels(dat$inversion), grid_type = "counterfactual")
g0 <- g0[do.call(paste, g0[, c("chord_type", "inversion")]) %in%
         do.call(paste, dat[, c("chord_type", "inversion")]), ]   # realized cells
g0 <- g0[g0$inversion %in% c("none", "0"), ]                      # nominated version
avg_predictions(m, newdata = g0, by = "chord_type", hypothesis = ~ pairwise)
```

Both return maj − aug = 1.0479 here, under either factor ordering, equal to the cell-mean difference – they coincide because this design is balanced and the covariate is distributed alike across conditions. That coincidence is a property of these data, not of the recipes: under imbalance, or with a covariate related to condition, only the second is a contrast of the same units under two interventions. Contrast this with the substitution-based version in Section 3.5, which produced two different wrong answers.

The `emmeans` analogue subsets the cell grid to the nominated cells before taking pairs:

```r
eg <- emmeans(m, ~ chord_type * inversion,   # the 16-cell grid (a grid
              nesting = NULL)                # specification, not a model
s  <- summary(eg)                            # formula); nesting = NULL keeps
keep <- which((s$chord_type == "aug" & s$inversion == "none") |   # it rectangular
              (s$chord_type != "aug" & s$inversion == "0"))
pairs(eg[keep])
```

This returns the identical estimates and standard errors, stable under reordering, and it fails safe: if the subset accidentally includes an unrealized cell, every contrast touching that cell comes back `NA` rather than wrong. (A substitution-based version – `avg_comparisons` with `newdata` restricted to the inversion-0 rows – also happens to be safe under the chain form, because the augmented predictions are constant across inversions; the observed-rows recipe is preferred because its validity does not depend on the parameterization.) The `nesting = NULL` argument is essential here, not decoration: because the chain formula triggers auto-detection, omitting it yields a *nested* grid, and in emmeans 1.10.0 subsetting a nested grid returns mislabelled contrasts – a silent wrong answer (a numeric contrast vector on a nested grid is rejected outright). Rectangular-grid subsetting, by contrast, is safe and fails safe. Whichever level is nominated should be stated in the text: a specific, defensible comparison beats an ambiguous general one.

### 4.5 Bayesian models in `brms`

Bayesian estimation changes the fitting machinery but not the estimand problem, and it adds one genuine capability: the structural zeros can be *declared* rather than silently discovered. Two facts to get right:

**The chain form is what makes the declaration possible.** Both parameterizations have sixteen columns; they allocate them differently. The crossed form allocates three to inversion main effects and nine to interactions (non-reference chord types × non-reference inversions only); the chain form allocates none to main effects and twelve to interactions – over *every* chord type, reference included – and the three extra columns are precisely the augmented cells' own, all-zero because those cells hold no data. The difference is between shared and dedicated parameters. Under the crossed form an impossible cell's prediction is assembled from parameters that also describe real cells – aug-at-inversion-0 is `(Intercept) + inversion0`, the same `inversion0` that serves every real inversion-0 cell – so no coefficient's sole job is the impossible cell, there is nothing for `constant(0)` to attach to, and `brms` stops with an error ("The following priors do not correspond to any model parameter"). The same entanglement is why augmented predictions are contaminated under the crossed form but constant under the chain form (Section 3.4). Running `get_prior()` and reading the coefficient names before writing priors is a prudent precaution either way.

**Fixing the impossible coefficients alone is necessary but not sufficient.** The design matrix has 16 mean-structure columns and rank 10. `constant(0)` priors on the three augmented-interaction coefficients declare the impossible cells, but three redundant directions remain – one per non-augmented chord type, because those chords never occur at the sentinel level – and under flat priors these leave the posterior improper. Two resolutions:

*Fully identify with three more constants.* Fix one interaction coefficient per non-augmented chord type to zero as well – here the `:inversion0` terms:

```r
library(brms)
pri <- c(prior(constant(0), class = "b", coef = "chord_typeaug:inversion0"),
         prior(constant(0), class = "b", coef = "chord_typeaug:inversion1"),
         prior(constant(0), class = "b", coef = "chord_typeaug:inversion2"),
         prior(constant(0), class = "b", coef = "chord_typedim:inversion0"),
         prior(constant(0), class = "b", coef = "chord_typemin:inversion0"),
         prior(constant(0), class = "b", coef = "chord_typemaj:inversion0"))
```

Sixteen columns minus six constants leaves ten free mean parameters for ten cells: the model is exactly identified, and no prior is left doing structural work. It samples cleanly even with `brms`' default flat priors on the remaining coefficients (maximum R̂ ≤ 1.005, no divergent transitions), and the restricted-grid contrast below returns aug − maj = −0.667. The three additional zeros are different in kind from the first three, and the difference belongs in a methods section: the augmented zeros declare that cells *do not exist*, whereas the others are an *identification constraint* – a coding choice, like picking a reference level, one per non-augmented chord type. Which coefficient is constrained is arbitrary and changes nothing that matters: constraining the `:inversion2` terms instead alters every Section 4 estimand only within Monte Carlo error (below 0.007 with 2,000 draws here), because under flat priors the induced posterior on the ten cell means is identical either way.

In practice, combine the constants with ordinary prior hygiene – keep all six `constant(0)` declarations and add a weakly informative prior on the remaining coefficients:

```r
pri <- c(prior(normal(0, 5), class = "b"), pri)   # regularization + the six constants
```

Identification comes from the constants; the prior does only its ordinary regularizing work. This is the specification to use; the flat-prior fit above is best read as the demonstration that the constants alone identify the model. It samples cleanly (maximum R̂ ≤ 1.003, no divergent transitions), returns aug − maj = −0.673 [95% CrI −1.12, −0.25], and – unlike the regularization-only route below – its coefficient posteriors are data-driven throughout: the posterior SD of `b_chord_typedim` is 0.262, on the scale of the frequentist standard error, rather than the 2.496 reported below.

*Or keep only the structural zeros and regularize.* Where choosing a constraint is unappealing, weakly informative priors on the remaining coefficients make the posterior proper instead:

```r
pri <- c(prior(normal(0, 5), class = "b"),
         prior(constant(0), class = "b", coef = "chord_typeaug:inversion0"),
         prior(constant(0), class = "b", coef = "chord_typeaug:inversion1"),
         prior(constant(0), class = "b", coef = "chord_typeaug:inversion2"))

mb <- brm(response ~ chord_type + chord_type:inversion + training,
          data = dat, prior = pri, chains = 2, iter = 2000, seed = 1)
```

The model samples cleanly (maximum R̂ = 1.004, bulk ESS ≥ 680, no divergent transitions), the three declared coefficients are exactly zero in every draw, and – the practically important point – **the same estimand machinery applies unchanged**: the restricted counterfactual grid from Section 4.2, passed to `avg_predictions(mb, ...)`, returns aug − maj = −0.669 [95% CrI −1.10, −0.24], closely matching the frequentist −0.6779. The self-documentation is the principal benefit: a reader of the model code can see which parameters were declared absent and why.

This model's predictions and estimands are essentially identical to the fully constrained model's – aug − maj = −0.669 [−1.10, −0.24] here against −0.673 [−1.12, −0.25] there, a difference within Monte Carlo error – because the inflated coefficient uncertainties cancel in every identified combination: the correlations of ±0.99 mean that the flat-direction variance each coefficient carries is removed the moment the coefficients are summed. The clearest single demonstration: the posterior SD of `b_dim + b_dim:inversion0` (which is dim-at-inversion-0 − aug) is 0.267 in this model, from components with SDs near 2.5, and the same quantity appears in the fully constrained model as the lone coefficient `b_dim`, with SD 0.262. The drawback of this route is therefore confined to the coefficient table, but there it is severe: along the three unresolved directions the posteriors are pure prior read-back – SDs of about 2.49, almost exactly the 2.5 that a `normal(0, 5)` prior implies along a flat direction. Report estimands from this model, never its coefficients.

This whole apparatus is the Bayesian face of Section 3.2's aliasing. The redundancy is a property of the design matrix, not of any particular coefficient – which six coefficients `lm` reports as `NA` was always an artefact of its order-dependent pivoting. `brms` drops nothing, so the same six-dimensional redundancy sits in the posterior until it is resolved. The full set of `constant(0)` declarations amounts to choosing one's own pivot – explicitly, and immune to factor ordering – and the regularization route is a decision not to choose, with the consequences just described.

---

## 5. Deeper nesting and nested continuous variables

The machinery of Sections 3–4 was built for one categorical nested factor, but nothing in it depends on that. Nesting can recurse: a nested categorical variable has levels, and further variables can be nested within those. A nested *continuous* variable, by contrast, terminates a chain – having no levels, nothing can nest within it – so it is the leaf case, treated second.

### 5.1 Nesting within nesting

Let a third factor enter: for root-position triads only, the chord's bass note is doubled at the `octave` or the `fifth` (`doubling`, sentinel `"none"` elsewhere). This is nesting of the by-design kind from Section 1 – doubled first-inversion triads are perfectly coherent, merely unimplemented – and, as promised there, nothing in the analysis is different for it. The nesting is now a chain – doubling within inversion within chord type – and the realized design has 13 cells of the 48 in the full crossing. The rules do not change; three points deserve attention.

**The formula generalizes as a chain, and the reorder test permutes every nested factor.**

```r
m <- lm(y ~ chord + chord:inv + chord:inv:doub, data = dd)
```

This fits with 35 of its 48 coefficients aliased. Permuting the levels of `inv` *and* `doub` leaves the likelihood identical, and the naive `avg_comparisons` chord contrast again moves (0.706 to 0.749) while every recipe below is stable to 9 × 10⁻¹⁶.

**Auto-detection scales with the chain.** With the chain formula, `emmeans` prints "A nesting structure was detected: `inv %in% chord, doub %in% (chord*inv)`" and returns correct nested answers even undeclared – though declaring the chain, `nesting = c("inv %in% chord", "doub %in% inv")`, still documents the design and is recommended.

**"Equal weights" becomes ambiguous at depth, and the intended sense needs stating.** With two doubling cells inside inversion 0, there are now two different equal-weight estimands for a chord-type mean: *hierarchical* – equal at each level of the tree, so inversions count equally and doublings count equally within inversion 0 – and *flat* – all realized cells counted equally. They differ (aug − maj: −0.6283 hierarchical, −0.7574 flat on the demonstration data). `emmeans`' nested averaging gives the hierarchical one; `marginaleffects` gives every option from one grid of the 13 realized cells:

```r
cellgrid   <- unique(dd[, c("chord", "inv", "doub")])   # the 13 realized cells
cellgrid$w <- ifelse(cellgrid$doub == "none", 2, 1)     # tree weights: inversion 0
                                                        # splits its share between
                                                        # its two doubling children
avg_predictions(m, newdata = cellgrid, by = "chord",
                hypothesis = ~ pairwise)                # flat:         -0.7574
avg_predictions(m, newdata = cellgrid, by = "chord",
                wts = cellgrid$w, hypothesis = ~ pairwise)  # hierarchical: -0.6283
```

The hierarchical result matches `emmeans` to the standard error and is stable under reordering of both nested factors. (Weights are normalized within each `by` group, so the weight on the augmented chord's single row is immaterial; in general a cell's hierarchical weight is the product of 1/(number of siblings) down its branch.) Observed-frequency weighting is unambiguous and is what the observed-rows and restricted-counterfactual recipes deliver (−0.7574 here, coinciding with flat only because the cells are balanced); the restricted grid is built by keeping the counterfactual rows whose (chord, inv, doub) triple occurs in the data. Comparisons at a nominated cell work unchanged at any depth: predictions on the observed aug rows and the maj/inversion-0/octave rows return their raw (adjusted) cell difference, 1.0797, stable under every reordering.

### 5.2 Continuous nested variables

Suppose stimuli either have vibrato or do not, and for those that do, a continuous vibrato `rate` (4–8 Hz) is manipulated. `rate` is a nested continuous variable: undefined – not zero, not missing – for the no-vibrato stimuli. (This is the structure of "years since quitting" for ex-smokers, or "husband's income" for the married.) The demonstrations below use 150 no-vibrato and 250 vibrato trials, with a true vibrato effect of 0.3 + 0.25 × rate.

**Every failure mode of Section 3 recurs here, one-for-one.** Coding `rate` as `NA` for the undefined rows deletes *all* of them – and because that removes an entire level of `vibrato`, `lm` now fails with the baffling error "contrasts can be applied only to factors with 2 or more levels", which blames the factor coding rather than the deletion. The sentinel here is **0**. With zero-coding, the crossed formula `y ~ vibrato * rate` is rank deficient – `rate` and `vibrato:rate` are identical columns – and the surviving coefficient labelled `rate` (0.2009) is the vibrato-condition slope wearing a main-effect label. The nested-style formula is the right one:

```r
m <- lm(y ~ vibrato + vibrato:rate, data = d)
```

This carries one all-zero column (`vibratonone:rate`), aliased and reported `NA` – and, exactly as in Section 4.5, that coefficient exists by name in `brms` under this formula and is declarable with `constant(0)`. One more lying label: `vibratopresent` (0.703) is the fitted difference *at rate 0*, an extrapolation far below the observed 4–8 Hz range. Centring the rate within the defined rows (`rate_c = rate − mean(rate among vibrato trials)`, kept at 0 elsewhere) leaves the fit unchanged and makes the coefficient the difference at the mean defined rate (1.9313).

**The defaults fail in the familiar, asymmetric ways.** The reference grid reduces `rate` to its pooled mean – **3.82 Hz, an average over the structural zeros**, a rate that is neither the sentinel nor anything realized. `emmeans` then refuses (`nonEst`, because its grid places that rate on the no-vibrato row, touching the aliased column); `avg_comparisons()` complies, returning 1.4707 – the effect at 3.82 Hz, an estimand nobody asked for.

**The clean recipes are the familiar ones.** Averaged over the observed rates (the proportional analogue): predictions on observed rows,

```r
avg_predictions(m, newdata = d, by = "vibrato", hypothesis = ~ pairwise)
```

giving 1.9313 – the value the centred coefficient reports directly. At a nominated rate (the nominated-level analogue): a two-row grid via `datagrid(vibrato = ..., rate = 6)` gives 1.9086, and the `emmeans` cell-grid subset from Section 4.4 – `emmeans(m, ~ vibrato * rate, at = list(rate = c(0, 6)))`, keeping the (none, 0) and (present, 6) rows – returns the identical estimate and standard error, correctly marking the impossible (none, 6) cell `nonEst`. The prediction for the no-vibrato condition is invariant to whatever rate sits on its grid row, because under the nested-style formula rate enters only through the vibrato interaction. There is no natural "equal weights" here – a continuous variable has no levels to weight equally – so the estimand menu is: at a nominated rate, or averaged over a stated rate distribution (the observed one being the usual choice).

---

## 6. Random effects in the partially nested design

The structural zeros constrain the random side as well as the fixed one, and less visibly, because the formula syntax accepts requests the design cannot support.

A *unit* is a level of the grouping factor – a participant, a stimulus, a classroom – and a *condition* is one of the ten realized cells. A unit's random effect is its set of deviations, one per condition it can meet.

### 6.1 Three standard covariance structures

Ten conditions give a unit ten deviations, whose joint distribution is a 10 × 10 covariance matrix. The constraints placed on that matrix are its *covariance structure*. Three are standard.

**Unstructured.** Every condition has its own variance and every pair its own correlation: 10 + 45 = 55 parameters.

```r
(0 + cell | participant)
```

**Diagonal.** Variances free, all correlations fixed at zero: 10 parameters.

The double bar is the usual notation for this, and **it is not recommended with factors.** Because `||` acts while the formula is parsed, before any model matrix exists, it cannot tell a factor from a numeric variable: it splits the term literally, leaving the levels within each part still correlated (Bates, Mächler, Bolker and Walker, 2015; Singmann and Kellen, 2019). `(0 + cell || participant)` therefore estimates all 55 parameters, in `lme4` and `brms` alike. What it does produce with several factors – separate blocks, correlated within and uncorrelated between – is seldom the intended structure, and can carry more parameters than the single-bar model it was meant to simplify.

From `lme4` 2.0-0 a diagonal covariance is stated directly, and this is the clearest form:

```r
diag(0 + cell | participant)
```

On lme4 2.0-6 that gives the intended 10 parameters, as does `(0 + cell || participant)` once `options(lme4.doublevert.default = "diag_special")` is set. The package default remains `split`, however, so the behaviour described above is what an unmodified installation still does, and the documentation recommends stating `diag()` explicitly rather than relying on `||`. On earlier versions the same structure requires one indicator variable per condition,

```r
(0 + d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8 + d9 + d10 || participant)
```

which works on every version. `afex::mixed(..., expand_re = TRUE)` and `glmmTMB` have long produced the intended structure directly (Singmann and Kellen, 2019).

**Compound symmetry.** A variance for each level of the nesting, the deviation at each level shared by every condition beneath it: 3 parameters here. Written with nested grouping factors:

```r
(1 | participant) + (1 | participant:chord_type) + (1 | participant:cell)
```

The name describes what the implied covariance looks like: all conditions share one variance, and any two sharing a chord type share one correlation. Both restrictions come with the notation; neither is estimated.

Compound symmetry and diagonal are each special cases of unstructured, so any of them may be compared against it by likelihood-ratio test wherever both fits are available.

### 6.2 The reach of the original variables

A random-effects term written in chord type and inversion produces sixteen columns, one per combination in the full crossing. Six carry no information: three because those conditions do not exist, three because they can be reconstructed from the others. This is the same redundancy that leaves six fixed-effect coefficients unestimable in Section 3.2 – both parts are built from the same design matrix, so one examination settles both.

Nested grouping factors are unaffected: only realized combinations become levels of `participant:cell`, so nothing empty ever arises. Compound symmetry is therefore available in the original variables, from any engine.

Unstructured and diagonal are a different matter. Both can be had in the original variables – delete the six uninformative columns and place the covariance on the ten that remain, which fits identically to `(0 + cell | participant)`, agreeing on the demonstration data to optimizer tolerance (−1932.99 against −1932.97, both with 55 parameters). What no formula produces is those ten columns: any formula in chord type and inversion regenerates all sixteen. In `lme4` they must be built as ten variables by hand – which is the cell factor under other names.

`brms` offers an alternative. A prior specification can hold a coefficient at zero,

```r
set_prior("constant(0)", class = "sd", group = "participant",
          coef = "chord_typeaug:inversion1")
```

and six such declarations silence the uninformative columns. This is what the prototype of the Appendix does, and the fit is correct. But silencing is not deletion, and four costs follow.

*Sampling.* The unit deviations are stored at sixteen per unit rather than ten, and the correlation matrix at sixteen dimensions rather than ten. On the demonstration design that is 776 sampled parameters against 455, so each iteration does roughly seventy per cent more work for the same inferences.

*Diagnostics.* Of the 120 correlations, 75 are invisible to the likelihood and simply reproduce their prior. They are not ill-behaved – a parameter sampling from a proper prior mixes perfectly well – but they must be excluded before any summary of convergence over the model as a whole means anything.

*Output.* Every posterior summary lists all 120, and a reader has no way to tell from the table which 45 carry information.

*Priors.* This is the one that affects the answers. A prior placed on a sixteen-dimensional correlation matrix does not say the same thing about the ten real conditions as the same prior placed on a ten-dimensional one: the implied prior on the smaller matrix takes a different shape, so a specification intended as weakly informative in one parameterization is not the same specification in the other. Matching them requires adjusting the prior's parameter, not merely restating it.

Terms that do not cross the boundary need none of this. Every unit meets every chord type, so `(1 + chord_type | participant)` is identified as written.

### 6.3 Which parameterization `nestimand` uses

The two ways of writing the model are equivalent, so the choice is practical, and the emitted code states it.

**Cells by default.** Nothing to derive, nothing held at zero, no coefficient named, no choice of which columns to constrain. Ten standard deviations and forty-five correlations, all meaningful, in place of sixteen and 120 of which most are not; diagnostics read as they stand; fewer parameters to sample. Every frequentist engine is served on the same footing as `brms`.

**The chain form when something requires it.** `emmeans` works from the formula and needs the original factors as predictors. A prior stated independently on each effect, or with bounded support, is tied to the coordinates it is written in. And a user may want coefficients that read as effects: with the constraints in place, `b_chord_typemaj` is major-at-root-position minus augmented, where a cell coefficient is a cell mean. `chain_priors()` derives the declarations from the design matrix, covering the fixed and random parts with the same six columns, and distinguishes structural zeros from identification constraints – the first says a condition does not exist, the second is a coding choice like a reference level.

Fixed-effect estimands agree between the two to Monte Carlo error.

**The random structure remains the user's.** Only terms crossing the boundary are translated; `(1 | participant)`, `(chord_type | participant)`, and a hand-written grouping structure are returned unchanged. A structure varying with inversion but not chord type is refused rather than enlarged, since over the realized conditions it is a covariance of reduced rank that no formula states, and the alternatives are named in the error.

**Identification is only half the problem.** A fully identified model still returns a number for a condition that does not exist – on the demonstration data, diminished chords at the sentinel level come back as 4.2854 – and a grid built by crossing the factors still counts the augmented chord four times. Asked for a chord-type contrast over such a grid, the identified model gives −0.5697 where the answer is −0.6779. Restricting the grid and choosing the weighting are questions of what is being estimated, not of whether it can be, and Section 4 settles them for every parameterization alike.

### 6.4 The cost of each parameterization

An unidentified fixed effect prints `NA`. The random side gives no such signal: on 2400 rows – six trials per unit per condition – `lmer` returns the sixteen-dimensional covariance, a Hessian with twenty-one negative eigenvalues, and a convergence warning. That warning is the most familiar message in mixed modelling and is usually read as sparse data, met by simplifying until it stops. The simplification succeeds, and the analyst concludes the data were thin rather than that the request was unanswerable.

The least constrained structure is not always the right one. On the demonstration data it is not: unstructured costs fifty-two parameters more than compound symmetry for twenty-nine log-likelihood units, which a likelihood-ratio test does not separate from chance (p = 0.221). The fitted correlations range from −0.319 to 0.453, but with forty units that spread is consistent with noise around a common value.

The cell parameterization makes the comparison routine rather than laborious. In the original variables, unstructured needs ten hand-built columns, or six declarations and an over-sized covariance to read around; in cell coordinates it is one formula. A comparison available in principle but tedious in practice tends not to be made.

The cost falls on reporting. Under compound symmetry, heterogeneity in the chord-type effect is a named variance component. Under the cell parameterization it is the quadratic form c′Σc, where Σ is the covariance over conditions and c is the same weight vector that defines the corresponding fixed-effect contrast. For the equal-weight augmented-minus-major contrast, c places +1 on the augmented condition and −1/3 on each major one, giving a unit standard deviation of 1.380 on that effect against 0.909 on diminished-minus-major. These are the quantities a report needs, and in `brms` each carries a posterior. But c changes with every contrast and every weighting policy, and a mistaken c yields a plausible number rather than an error; `nestimand` builds it from the same policy that defines the contrast.

### 6.5 Fitting the identified structures

The demonstrations use a version of the dataset with three trials per unit per condition and genuine variance at unit, unit × chord-type, and unit × condition level.

**Slopes on the nesting factor are safe.** Every unit meets every chord type, so `(1 + chord_type | participant)` is identified and fits cleanly: 10 parameters, no singularity.

**Slopes on the partially nested factor are not.** `(1 + inversion | participant)` asks for a slope on the sentinel level, treating "being an augmented chord" as a kind of inversion, identified only from the augmented trials. On the demonstration data it fitted singular, with a worse AIC (2862.9) than the chord-type slopes (2842.1) and barely better than a bare random intercept (2865.3). A singular fit is easily misread as over-parameterization rather than structural impossibility.

**Compound symmetry**, the structure available in the original variables:

```r
dd$cell <- droplevels(interaction(dd$chord_type, dd$inversion, sep = "_"))

lmer(response ~ chord_type + chord_type:inversion +
       (1 | participant) +
       (1 | participant:chord_type) +
       (1 | participant:cell),
     data = dd)
```

`cell` holds only the ten realized combinations, so `participant:cell` has 30 × 10 = 300 levels and the impossible ones never arise. Three variance parameters; on the demonstration data it fitted without singularity and gave the best AIC of the four (2827.6). Note that `cell` enters as a *grouping* variable only – the fixed effects remain the two-factor chain form, so this is not the cell parameterization.

`(1 | participant:cell)` needs more than one observation per unit per condition. On this document's main dataset, with one trial each, `lmer` stops outright: "number of levels of each grouping factor must be < number of observations". With single trials, drop that term.

**The fixed-effect summaries of Section 4 apply unchanged.** The `emmeans` equal-weight chord-type contrast and the `marginaleffects` realized-cell grid – both with the random effects set to zero, `re.form = NA` for frequentist fits and `re_formula = NA` for `brms`, so that each describes the population rather than a particular unit – return the identical value, aug − maj = −0.9036, stable under factor reordering. The reorder test applies regardless of the random structure.

### 6.6 Which population does a prediction describe?

A further distinction becomes consequential once a model carries random effects, and it is a choice of estimand rather than of algorithm. A population-level prediction may set the random effects to zero, giving the effect for a typical group; it may average over the groups actually sampled, using their estimated deviations; or it may integrate over the fitted random-effect distribution, giving the population-averaged effect. In a linear model with balanced data the three coincide, which is why the choice can go unnoticed. In nonlinear models – logistic, ordinal, and the rest – they diverge, and the second and third answer different questions: an average over the groups actually sampled, and an average over the population from which they were drawn. How the second is computed matters. Point estimates of group deviations – best linear unbiased predictors, or posterior means – are shrunk towards zero, so their ensemble is under-dispersed relative to the parameters it represents: in a hierarchical model with a shrinkage factor of 0.39, a true group standard deviation of 0.883 appears as 0.555 across the point estimates, and as 0.896 across single posterior draws, one per group (Louis, 1984; Shen and Louis, 1998). Averaging a nonlinear quantity over the shrunken ensemble therefore distorts it – in a probit illustration, a risk difference of 0.304 with the random effects at zero becomes 0.229 averaged over point estimates against 0.225 integrated over the fitted distribution – whereas averaging draw-wise, each posterior draw using its own group deviations, carries the correct spread and leaves only the intended finite-population difference. A Bayesian fit supports this directly and also propagates uncertainty in the variance component; a frequentist fit, having only the shrunken point estimates, cannot reproduce the same quantity. In `marginaleffects` the first is the `re.form = NA` setting of Section 6.5, the second follows from supplying observed rows or a grid crossed with the grouping factor, and the third requires new group levels drawn from the fitted distribution, which `brms` supports directly. The redesigned `nestimand` makes this an explicit argument rather than a default, on the same principle as the weighting declaration of Section 4.2.

## 7. Checklist

Where the Appendix's generator automates an item, its function is named in parentheses.

Before modelling:

- [ ] `table(factorA, factorB)` – how many cells actually exist? (`nesting_spec()` reports the realized-cell count for every declared family.)
- [ ] Is the undefined level an explicit sentinel (`"none"`), never `NA` – and listed first among the factor's levels (Section 4)? (`apply_sentinel()` performs the conversion, refusing genuine missingness; `nesting_spec()` refuses `NA` outright.)
- [ ] Is the model formula in chain form – `g1 + g1:g2 + …`, never crossed (`g1 * g2`)? Section 4 explains why. (`model_formula()` derives it from a freely written formula.)
- [ ] For a nested *continuous* variable: sentinel is 0, no free-standing main effect (`g + g:x`, never `g * x`), and check what value the reference grid assigned it before trusting any covariate-adjusted contrast (Section 5.2). (The same two functions apply; a continuous variable is refused in any non-leaf position.)

After fitting:

- [ ] Compare `nobs(model)` with `nrow(data)` – did casewise deletion remove a condition?
- [ ] **Run the reorder test** (Section 3.3): permute the levels of *every* nested factor, refit, and recompute the reported estimates. Identical `logLik` is necessary but not the test – the fit is always unchanged; it is the *estimates* that must not move, and any that do were never identified. (Scripts from `emit_estimand()` end with this check built in, at the estimate level.)
- [ ] For `brms`: do the structurally impossible coefficients carry `constant(0)` priors – or deliberate regularization – derived from the actual design matrix rather than assumed (Section 4.5)? (`brms_priors()` derives the full block, distinguishing structural from conventional constraints.)
- [ ] For mixed models: does any random slope span the structural boundary? Slopes on the nesting factor are fine; slopes on the partially nested factor are not (Section 6.1). A convergence warning is not evidence of sparse data here – check the rank of the random-effects design before simplifying (`nestimand::random_slope_rank()`).
- [ ] Is the random structure the one intended, or the one a formula in the original variables happens to reach? Compound symmetry, through nested grouping factors, is what such a formula gives directly; richer covariances require either the cell parameterization or, in `brms`, a block of `constant(0)` declarations (Sections 6.1 and 6.2).

Before reporting:

- [ ] Within-stratum contrasts: from `emmeans` (the chain formula itself declares the nesting), or from predictions restricted to realized cells – never from default `avg_comparisons` output containing sentinel rows. (`emit_estimand(weights = "within")` emits both forms.)
- [ ] Across-boundary contrasts: which weighting – equal, proportional, G-computation over observed rows, or a nominated level – and is it stated in the methods? For nominated-cell comparisons via the `emmeans` cell grid, remember `nesting = NULL` (Section 4.4). (`emit_estimand()` makes the weighting an explicit argument and emits the matching calls.)
- [ ] For nonlinear models: population-averaged (G-computation) or conditional at covariate values, and does the text say which?

### Reporting

A methods section does not need the machinery in this document – reviewers need the estimand, not the pivoting story. Three sentences suffice, whichever package produced the numbers: that the design is partially nested and which factor is undefined where; that contrasts were computed over the realized cells only; and which weighting was used for any comparison crossing the structural boundary (equal, observed-frequency, averaged over the observed data, or at a nominated level). Where covariates entered a nonlinear model, a fourth states whether effects are population-averaged or conditional at covariate values.

---


## References

Bates, D., Mächler, M., Bolker, B., and Walker, S. (2015). Fitting linear mixed-effects models using lme4. *Journal of Statistical Software*, 67(1), 1–48.

Bates, D., Kliegl, R., Vasishth, S., and Baayen, H. (2015). Parsimonious mixed models. arXiv:1506.04967.

Clark, H. H. (1973). The language-as-fixed-effect fallacy: a critique of language statistics in psychological research. *Journal of Verbal Learning and Verbal Behavior*, 12(4), 335–359.

Cole, S. R., and Frangakis, C. E. (2009). The consistency statement in causal inference: a definition or an assumption? *Epidemiology*, 20(1), 3–5.

Díaz, I., and van der Laan, M. J. (2013). Assessing the causal effect of policies: an example using stochastic interventions. *The International Journal of Biostatistics*, 9(2), 161–174.

Hernán, M. A. (2016). Does water kill? A call for less casual causal inferences. *Annals of Epidemiology*, 26(10), 674–680.

Hernán, M. A., and Robins, J. M. (2020). *Causal Inference: What If*. Boca Raton: Chapman & Hall/CRC.

Holland, P. W. (1986). Statistics and causal inference. *Journal of the American Statistical Association*, 81(396), 945–960.

Kennedy, E. H. (2019). Nonparametric causal effects based on incremental propensity score interventions. *Journal of the American Statistical Association*, 114(526), 645–656.

Louis, T. A. (1984). Estimating a population of parameter values using Bayes and empirical Bayes methods. *Journal of the American Statistical Association*, 79(386), 393–398.

Manski, C. F. (1990). Nonparametric bounds on treatment effects. *American Economic Review*, 80(2), 319–323.

Muñoz, I. D., and van der Laan, M. J. (2012). Population intervention causal effects based on stochastic interventions. *Biometrics*, 68(2), 541–549.

Pearl, J., and Bareinboim, E. (2014). External validity: from do-calculus to transportability across populations. *Statistical Science*, 29(4), 579–595.

Robins, J. M. (1986). A new approach to causal inference in mortality studies with a sustained exposure period. *Mathematical Modelling*, 7(9–12), 1393–1512.

Rubin, D. B. (1980). Comment on "Randomization analysis of experimental data: the Fisher randomization test". *Journal of the American Statistical Association*, 75(371), 591–593.

Shen, W., and Louis, T. A. (1998). Triple-goal estimates in two-stage hierarchical models. *Journal of the Royal Statistical Society: Series B*, 60(2), 455–471.

Singmann, H., and Kellen, D. (2019). An introduction to mixed models for experimental psychology. In D. H. Spieler and E. Schumacher (eds), *New Methods in Cognitive Psychology*, 4–31. New York: Routledge.

VanderWeele, T. J. (2009). Concerning the consistency assumption in causal inference. *Epidemiology*, 20(6), 880–883.

VanderWeele, T. J., and Hernán, M. A. (2013). Causal inference under multiple versions of treatment. *Journal of Causal Inference*, 1(1), 1–20.

Westreich, D., Edwards, J. K., Lesko, C. R., Stuart, E., and Cole, S. R. (2017). Transportability of trial results using inverse odds of sampling weights. *American Journal of Epidemiology*, 186(8), 1010–1014.

## Appendix: the `nestimand` prototype

*Development note: this Appendix documents the earlier, chain-based version of `nestimand` – the prototype verified against the examples in this document. The current design – the cell-translation architecture stated in the Purpose note and motivated in Sections 4 and 6 – is specified in the accompanying design document (`nestimand_design.md`); the user-facing declarations below are intended to survive the transition.*

`nestimand` is a proof-of-concept companion to this document (a single R source file, not yet a package). Its design principle is to *generate scripts rather than wrap computation*: every function returns either a validated description of the design or plain R code that can be read, edited, and pasted into a supplement. Each generated analysis script ends with an embedded reorder test (Section 3.3), so its own validity check travels with it. All behaviour described below was verified by executing the functions against the examples in this document.

**`nesting_spec(data, formula, nests, fit = "lm", family = NULL, random = NULL)`**
Declares the analysis. `formula` is the *desired* model, written freely in ordinary R syntax – crossed terms welcome, e.g. `Y ~ X1 * X2 * (Z1 + Z2) + X3 * X4 + Z4`; it is transformed, not fitted as given. `nests` is a character vector of `"child %in% parent"` links, and several independent nesting families may be declared at once (`c("X2 %in% X1", "X4 %in% X3")`). Random-effects terms may be written directly in the formula in the usual `lme4`/`brms` syntax (or supplied via `random`); either way they are transformed by the same principle as the fixed part – covariate slopes are retained, while factor structure over nested variables becomes a chain of grouping terms over realized cells (`(chord_type * inversion | id)` becomes `(1 | id) + (1 | id:chord_type) + (1 | id:cell)`, with the `cell` factor constructed in the emitted script), since slopes over nested variables are unidentified across the structural boundary for any amount of data; the grouping chain models the same per-group cell deviations under an additive variance-components covariance – a constrained, always-identified *submodel* of the slope form's unstructured covariance, not an equivalent reparameterization (cells sharing a stratum are constrained to equal correlation; an identified unstructured alternative over realized cells is `(0 + cell | group)`, at a much larger parameter cost) (Section 6); the emitted script records the original and transformed structures side by side. A `fit` argument selects the modelling engine – `"lm"` (default), `"glm"`, `"lmer"`, `"glmer"`, `"clm"`, `"clmm"`, or `"brms"` – with `family` and `random` companions where the engine needs them (`family` accepts a string or a bare family call, e.g. `family = cumulative("probit")`); under an ordinal `brms` family the emitted `marginaleffects` calls contrast on the latent scale (`type = "link"`, one number per contrast), with the per-category probability alternative and its computational guardrail noted in the script; emitted scripts then load the right package, append the random terms, add `re.form = NA` to `marginaleffects` calls for mixed models (Section 6), embed the `brms_priors()` block for `brms` (with the reorder tolerance relaxed for Monte Carlo error), and the outcome-type validation becomes fit-aware (`"clm"` requires the ordered factor that `"lm"` refuses). Everything in the formula that belongs to no nesting family is treated as a covariate, and covariate moderation is expressed in the formula itself (`X1 * X2 * Z1` yields the covariate chain automatically). Validation is deliberately strict: any `NA` in a nesting variable is refused with a message explaining the sentinel convention (Section 3.1) – `apply_sentinel()` below performs the conversion safely – and a continuous variable in any non-leaf position is refused (Section 5). The `print` method reports each family's chain and realized-cell count, the covariates, and the transformed formula.

**`model_formula(spec)`**
Returns the safe, chain-form transformation of the declared formula, obtained by one rule applied to every term: *close it under nesting ancestry* – any term containing a nested variable gains that variable's ancestors – then collapse duplicates and complete each used chain's prefixes. Under this rule `X1 * X2` becomes `X1 + X1:X2`, `vibrato * rate` becomes `vibrato + vibrato:rate`, and `chord_type * inversion * training` becomes the full covariate chain of Section 4.2; the chain form of Section 4 is thus derived rather than imposed. The result is fit- and estimand-identical to the declared crossed formula wherever both are estimable, and it is the parameterization the document itself uses throughout (the note opening Section 4 gives the reasons).

**`brms_priors(spec)`**
Returns, as text, a complete `brms` prior block for the chain model, derived mechanically from the design matrix: all-zero columns become `constant(0)` priors annotated as *structural zeros* (cells that do not exist), the residual rank deficiency identified by QR pivoting yields further `constant(0)` priors annotated as *identification constraints* (a coding choice; estimands are unaffected by which coefficients are chosen), and a weakly informative `normal(0, 5)` prior covers the remaining coefficients – the recommended specification of Section 4.5, with its two kinds of zero distinguished in comments.

**`emit_estimand(spec, target, weights = "proportional", engine = c("marginaleffects", "emmeans"), at = NULL)`**
Returns, as text, a self-contained analysis script for contrasts of `target` levels – a single variable or a vector of them (`c(chord_type, inversion)`), in which case the model is fitted once, each target's estimand is computed in its own `local()` block, `est` is a named list, and the reorder check verifies every element – under a stated weighting: `"within"` (per-stratum contrasts of a nested `target` – the Section 4.1 quantity – emitted for either engine), `"proportional"` (predictions over observed rows), `"equal"` (flat over the realized-cell grid), `"hierarchical"` (the same grid with tree weights – equal at each level of the chain), `"counterfactual"` (G-computation over the restricted counterfactual grid of Section 4.2, preserving the observed covariate distribution; `marginaleffects` only, since the corresponding `emmeans` machinery fails on rank-deficient fits), or `"nominated"` (observed rows restricted by the `at` condition). Marginal weightings of a *nested* target are `marginaleffects`-only, and are computed over the strata in which the target varies: degenerate strata are excluded from the emitted grid (detected structurally, as parent combinations with a single target level), so no across-boundary sentinel contrasts appear in the output. `emmeans` has no pooled marginal for a nested variable and the engine refuses with a pointer to `"within"`. The `emmeans` engine emits the nested-declaration equivalents for the hierarchical and proportional cases; flat-equal and nominated cells are emitted for `marginaleffects`. The estimand vocabulary of Section 4 is thus the function's argument list, so the weighting disclosure required for reporting (Section 7) is forced at the point of code generation. Every emitted script fits the model, computes the estimand into `est`, and closes with a reorder test that permutes the non-reference levels of every nested factor and requires the estimates themselves to be unchanged. A `self_check` argument scales this to the fitting cost: `"full"` (the default), `"fast"` (for `brms`, the verification refit runs one short chain – gross instability is large, so reduced sampling suffices to detect it; for `lmer`, `glmer`, and `clmm`, the check runs instead on the fixed-effects *shadow model* – the same chain formula fitted by `lm`, `glm`, or `clm` without the random terms – since order-instability is a fixed-effects phenomenon and the random structure is orthogonal to it, Section 6), or `"none"`, which warns at emission and leaves a comment block in the script recording the omission, so the analysis record shows the check was waived rather than passed. The returned object is text – printing it shows the script with a banner saying so – and `run_estimand(script)` executes it in one call, first echoing the script it runs (`quiet = TRUE` suppresses the echo), so nothing executes unseen. Emitted scripts carry every consequential default on their face: the `brms` MCMC controls are written out explicitly rather than inherited, and the weakly informative prior is comment-flagged for scale review.

**`apply_sentinel(data, var, where, sentinel = "none")`**
Converts `NA` values in `var` to a sentinel, but only in the rows the `where` predicate declares structurally undefined (e.g. `data$chord_type == "aug"`); numeric variables default to a sentinel of 0. The predicate is mandatory by design, because `NA` conflates two things this document insists on separating: an `NA` outside the declared rows is genuine missingness and is refused rather than converted, and a non-`NA` value inside them is inconsistent coding and is likewise refused. For factors, the sentinel is placed first in the level order of the returned variable.

---

*Every numerical claim in this document was verified by executing the code under R 4.3.3 with emmeans 1.10.0, marginaleffects 0.18.0, lme4 1.1-35.1, ordinal 2023.12-4, and brms 2.20.4 (August 2026). The double-bar and `diag()` behaviour of Section 6.1 was additionally verified under R 4.6.1 with lme4 2.0-6. Code blocks use the current `marginaleffects` hypothesis syntax (`hypothesis = ~ pairwise`); the verification environment's release used the equivalent string syntax of its day. The core behaviours – the `NA` deletion, the aliasing, the reorder instability, and the Section 4 recipes – were additionally confirmed for ordinal models fitted with `clm()`, where `emmeans` and `marginaleffects` agree to machine precision whenever the estimand is defined identically (predicted category probabilities with covariates at their means). Package behaviour may change between versions: the reorder test should be rerun on the data and software at hand rather than these results being trusted.*
