# Using This Code

The scripts are built around fish passing a dam, but the models are general. If you mark animals and detect them again, count animals with imperfect detection, or estimate a probability that changes over time and between groups, at least one of these scripts fits your problem.

## What each model is good for

**01 First Bayesian model.** Any probability estimated from successes out of trials: survival, detection, occupancy, nest success. It's the template for reading any JAGS output: check Rhat and trace plots, then summarize the draws.

**02 Trap efficiency.** Any study where you release a known number of marked animals and count how many come back: trap efficiency trials, mark-resight, survey recapture. Use it to see how sampling noise moves the estimate, why a small number of recaptures needs the +1 adjustment, and why 1 / (average detection) and the average of 1 / detection differ.

**03 Logit GLM.** Any yes/no outcome with covariates: survival to a site, presence at a plot, nest success, a fish taking one route or another. The pattern (Bernoulli trial, logit link, linear predictor) is the building block for everything after it.

**04 State-space model.** Repeated counts over time where detection is below 1: daily trap catches, weekly redd counts, annual surveys. Two lessons carry over:

- Counts alone estimate an index (N × p), not N. To get N you need separate data on detection, such as efficiency trials, double observers, or marked animals.
- Missing days or years can be left as `NA`. The model estimates them from the days around them.

**05 Six-cell multistate model.** Any system where animals take one of several hidden routes or states, and you see them at some detection sites and not others: dam passage, river junctions, breeding states, habitat use. The route is the hidden state; the detection histories are the data.

**06 Hierarchical model.** Many small groups that are alike but not identical: weeks, sites, streams, years. Thin groups borrow strength from the rest, and well-sampled groups mostly speak for themselves.

**07 Day and night.** What to do when one state splits further (day and night here), and why it matters when groups borrow from each other. If two groups share the same underlying rates but differ in timing or habitat use, borrowing a pooled rate carries that difference across; borrowing the split rates doesn't. Equal detection day and night is a simplifying assumption here. A single group cannot freely identify both rates, but sharing guidance across groups with different timing can add identifying information. Check identifiability and precision for the actual model.

**08 Bootstrap.** Any estimate built by expanding a count by an estimated rate (catch / efficiency, counts / detection). If the rate is uncertain, the interval has to carry that uncertainty too, or it will be far too narrow.

**09 Recover the truth.** Not an optional exercise: a template for every model you adapt. Before trusting a modified model on real data, simulate data from it with known values and confirm it gets them back. The script's header lists the four things to swap in.

## Adapting a script to your data

1. **Run the script unchanged first.** Make sure JAGS works and you get the printed results and plot.
2. **Find the data the model needs.** Each JAGS call takes a `jags.data` list. That list is the only thing your real data have to match. For example:
   - `04`: a vector of daily counts `C` (with `NA` for missed days), plus trial releases `M` and recaptures `R` with the day each trial ran.
   - `05`: a vector of six cell counts `n` and their total `N_seen`.
   - `06`: an array `n[group, time, cell]`, plus covariates and a stratum index.
   - `07`: a matrix `n11[group, cell]` with 11 cells in a fixed order (day cells, night cells, then c5 pooled).
3. **Format your data into that shape.** This is where most analyses go wrong, not the math. `05_six_cell.R` section 2 starts from one row per tagged fish, assigns each fish to exactly one cell, and tabulates. With groups and weeks, `table()` on group, week, and cell builds the array; the comments in `06` and `07` show the layout. Check unique tags, nonmissing binary flags, and compatibility with the assumed routes before assigning cells. Contradictory histories need investigation, not arbitrary priority. Check that every detected fish lands in exactly one cell and that the counts add up.
4. **Replace the simulation section with your data.** Keep the simulation code: you'll need it for step 6.
5. **Change the model to match your system.** Redraw the detection tree for your sites, write one `pi_raw` line per detectable history, and keep the conditioning lines (`p_seen` and `pi_obs`) as they are.
6. **Simulate from your new model and refit it** using `09` as the template, before using real data. If it can't recover known values, it won't recover unknown ones. Do this every time you change the model, not once.

## Where to get real data

- **PIT tag detections:** PTAGIS, https://www.ptagis.org/. Detection records come one row per tag and site; step 3 above turns them into cell counts.
- **River flow and spill:** Fish Passage Center historical flow and spill query, https://www.fpc.org/smolt/smolt_queries/Q_river_historicflowspill_dataquery.php. Covariates like the spill effect in `03` and `06`.

## Gotchas that bite real studies

- **Marked animals that aren't a fair sample.** Animals tagged after the outcome you're estimating (fish tagged inside the bypass, after they were already guided) make the rate look far better than it is, with a tight interval that misses the truth. Check where and when every tag was applied, and drop tags whose release site guarantees a route (`05`, section 7).
- **Incomplete transport flags.** Check transport records and distinguish passage time from loading time. In `05`, moving unflagged barged fish from c6 to c1 leaves closed-form GE unchanged because those fish remain guided, while transport and bypass recovery change. That example cannot demonstrate the production claim that GE moved substantially. Establish that effect with the actual before/after data and model.

- **A parameter the data can't separate.** Without the right detection site, two parameters trade off along a ridge. The sampler still runs and Rhat still looks fine (`05`, section 8). Ask what part of the data identifies each parameter before you fit.
- **Expanding by an uncertain rate.** 1 / rate is curved, so expanding by the average rate and averaging the expansions give different answers, and neither is guaranteed to be unbiased (`02` and `08`). Carry the rate's uncertainty into the interval, and check the point estimate against simulated truth.
- **An inverse with no finite mean.** With no guided fish, the wild-only flat-GE prior in `07` yields an infinite posterior mean of `1 / GE`. A finite MCMC average hides this. Lesson `07` uses posterior medians for all three methods and reports zero-guided weeks. Changing the prior requires scientific justification, not just a desire for a finite mean.
- **A conditional expansion labeled as abundance.** `05` holds `N_seen` fixed and propagates detection-probability uncertainty through `N_seen / p_seen`. Its expansion interval is not a full abundance posterior.
- **An interval that ignores an input.** A bootstrap that holds an estimated rate fixed gives an interval that's far too narrow (`08`).
- **Borrowing a mixed quantity.** If groups differ in timing or habitat use, a pooled rate mixes the thing you want with the mix itself. Borrow the components, not the blend (`07`).
- **Monitoring the wrong parameter.** When offsets are centered to sum to zero, their common mean cancels and the data say nothing about it; its posterior is just its prior. Rhat on the raw offsets can look poor for that reason alone, even when everything the model uses has converged. Monitor the identified quantity, the centered offsets (`06`).
- **Dropping fits that fail diagnostics.** In a simulation study, the fits that fail are often the hard data sets. Retry them with longer chains, keep comparisons paired, and report the failure rate next to results that describe only the fits that passed (`06`, `07`, `09`).

## Checks to run every time

- **Convergence.** Use `fit_diagnostics()` from `helpers.R`: rank-normalized split Rhat <= 1.01 and bulk/tail effective sample sizes >= 400 are initial screens. Missing diagnostics fail the check. Inspect trace plots and Monte Carlo error for the quantities reported, including inverse-GE summaries. Run longer chains, improve parameterization, or simplify when needed. See the [modern diagnostic guidance](https://mc-stan.org/docs/reference-manual/analysis.html).
- **Prior vs posterior.** If a parameter's posterior looks just like its prior, the data aren't informing it. In a Bayesian fit an unidentified parameter still gets a posterior (BPA 7.9).
- **Coverage.** Measure how often intervals contain the truth under the conditions you simulate. Nominal 95% Bayesian intervals do not guarantee 95% coverage at every fixed parameter value. Report Monte Carlo uncertainty and failures; 50 repeats give only a coarse check.
- **Use the draws, not the summary.** When an answer is a nonlinear function of parameters (1 / p, a product, a ratio), compute it for every posterior draw and then summarize. That carries the uncertainty through. It doesn't by itself remove bias: the posterior mean of a curved function can sit away from the truth, and plugging in posterior means gives a different answer again (scripts 02 and 08, Jensen's inequality). Check against simulated truth.

## Common extensions

- Covariates on any probability: add terms to its logit linear predictor, as in `03` and `06`.
- Random effects for time or site: replace a fixed parameter with draws from a shared distribution, as in `06`.
- More detection sites: add branches to the tree and cells to the multinomial.
- Combining data sets: add another likelihood line to the same model. The parameters they share are estimated from both.

## Questions

David Smith, Idaho Department of Fish and Game, david.smith@idfg.idaho.gov
