# Answers to the Try This Exercises

Try each exercise before reading its answer. Lessons 05, 06, 07, and 09 now use low downstream recovery throughout. Their old numerical answers came from a different, high-recovery scenario and have been removed. The answers below describe what to investigate; new fitted values and coverage results await the full lesson runs.

Numerical examples retained for 01--04 and 08 are historical runs of unchanged simulation settings, not newly verified outputs under the current diagnostic checks. Inspect diagnostics before interpreting any numerical comparison. JAGS results can vary by software version. Lesson 07a has its own previously checked closed-form outputs and no JAGS fits.

## 01 First Bayesian Model

### 1. Ten marked fish, two recaptured

The 95% credible interval for p widens from 0.084–0.286 (width 0.20) to 0.060–0.518 (width 0.46), about 2.3 times wider. The posterior mean moves from 0.173 to 0.25 even though R / M is 0.2: with only 10 fish, the flat prior pulls the mean toward 0.5. The posterior mean is (R + 1) / (M + 2).

### 2. The Beta(3, 17) prior with ten fish

With 10 fish, the informative prior changes the answer a lot: the mean drops from 0.25 to 0.167 and the interval width from 0.46 to 0.26. With 50 fish it moved the mean only from 0.173 to 0.157. A prior matters most when data are thin.

### 3. Doubling the catch

Expanded passage doubles: the median goes from 3,565 to 7,131 and the interval from 2,102–7,340 to 4,204–14,679. The efficiency posterior doesn't change (mean 0.172), because the catch carries no information about efficiency in this model. Only the marked fish do.

## 02 Lincoln-Petersen

### 1. Fifteen marked fish

854 of 10,000 repeats recapture no marked fish (1 with 50 marked fish), so the simple expansion is infinite in those. Among repeats with recaptures it runs 21% high. The +1 version averages 3,702, 7.5% low, matching the formula's 3,703, and its middle 95% widens to 1,539–9,808.

### 2. Efficiency of 0.03

The +1 correction no longer lands near truth: it averages 3,134, 22% low. The formula N(1 − (1 − p)^(M + 1)) predicts 3,154, so the bias is expected. Its source is (1 − p)^(M + 1), the chance of recapturing no marked fish, which is now large: 21% of repeats have zero recaptures. A different formula won't fix that; more marked fish will.

## 03 Bernoulli Logit

### 1. Fifty fish

In this run the interval for b1 is −1.54 to −0.17, so it still excludes zero, but it is about 2.8 times wider than with 400 fish (−0.96 to −0.47). With a different seed it could include zero. The fitted curve's band is about 3 times wider.

### 2. No spill effect

The b1 estimate is −0.03 with interval −0.25 to 0.19, which includes zero. The fitted curve is nearly flat (0.29 to 0.30 across the spill range). The band stays about as wide as before, because the model is still uncertain about the level and the slope.

## 04 State-Space Screw Trap

### 1. A 25-day outage

Inside the outage, the 95% interval for daily passage is on average 2.7 times the true passage, compared with 1.1 times on days the trap fished (1.4 vs 1.0 with the 5-day outage). All 25 outage-day intervals contain the truth. The season-total interval widens from 5,190–7,145 to 3,843–10,300. Both models print a diagnostics warning: passage in the middle of the outage (days 37 to 48) and the process SD have bulk effective sample sizes under 400, so run longer chains before reporting.

### 2. Ten marked fish per trial

Efficiency gets less precise (average 95% width 0.18 vs 0.13), and so does passage: the season-total interval widens from 5,190–7,145 to 4,746–7,289. Model 1 doesn't change, because it never uses the trials. Model 2 fails diagnostics on 69 quantities (lowest tail effective sample size 91), so it needs longer chains before its results can be trusted.

## 05 Six-Cell Multistate Model

Default GOJ recovery is about 3.25% for nonguided fish and 9.86% for bypassed fish returned to river. These probabilities combine survival and downstream detection. They are not survival estimates alone.

### 1. Five thousand fish tagged in the facility

Adding fish selected because they already entered the bypass makes the sample disproportionately guided. Compare GE and its interval with and without those fish, then check whether either interval covers the simulated truth. A narrow interval need not reveal selection bias. Report numbers only after checking the revised fits.

### 2. A Beta(1, 9) prior on psi without GOJ

Without GOJ, the data identify a combination of GE and detection, not each probability separately. A prior change on GE can move the posterior for detection even when detection's own prior is unchanged. Inspect the joint ridge and diagnostics. Sampler agreement does not establish identifiability.

### 3. Six hundred tagged fish

Low recovery makes c2, c4, and c5 especially sparse. Compare interval widths relative to truth and count these cells; which parameter is weakest depends on the realized sample. Section 7 still adds 1,500 facility fish, so the selection problem becomes large relative to the remaining upstream sample. Do not reuse the old high-recovery interval widths or GE estimates.

## 06 Hierarchical Model

### 1. Better-sampled wild weeks

Increasing the wild sample usually reduces dependence on shared structure, but the amount of identifying information depends on downstream overlaps as well as total tags. Compare wild GE error, interval width, coverage, and diagnostic failures. Do not assume a narrow hierarchical interval is well calibrated.

### 2. More week-to-week variation

Larger departures from the shared pattern can make shrinkage more costly. Examine whether the model estimates greater variation and whether unusual wild weeks retain interval coverage. One season is not enough to establish comparative performance; optional repeats provide that check.

### 3. Allowing either spill-slope direction

`restrict.spill <- TRUE` imposes a nonpositive coefficient before the data are seen. With `FALSE`, the prior permits either direction; the likelihood and simulated data are unchanged. Compare the slope, wild GE, and diagnostics. In the extension with true `beta <- 0.4`, the restricted model cannot recover the positive coefficient because its prior excludes it. Other model terms may absorb the mismatch, so inspect their behavior too.

## 07 Day and Night

### 1. Equal night shares

The true pooled GE becomes equal for wild and hatchery. This removes the systematic timing mismatch. Finite-sample error remains, and estimating separate day/night rates can cost precision. Compare repeated, paired results rather than treating one simulated week as a guaranteed improvement.

### 2. No night effect

When day and night GE both equal 0.02, night share cannot change pooled GE. Splitting is unnecessary for timing correction in this scenario. Low recovery and sparse wild histories can still leave large uncertainty and prior influence.

### 3. Twenty wild fish

Look for zero guided fish and zero downstream overlaps, not just total detections. Compare uncertainty, prior sensitivity, and diagnostic failures. The posterior median of inverse GE stays defined when the posterior mean does not, but a finite median does not restore information absent from the sample.

## 07a Closed-Form Day/Night Comparison

This example already uses low nonguided recovery (123/3,785, about 3.25%). Its simulation and plot are unchanged by the low-recovery revision to the JAGS lessons.

The previously checked seed-2025 run retains 1,997 of 2,000 complete seasons. Pooled borrowing has median season error +33.3% (middle 90%: +26.8% to +40.6%); split borrowing has +1.3% (-5.8% to +9.1%). Three wild week-15 samples lack GRS--GOJ overlap. These are conditional simulation-error summaries, not posterior intervals.

Equal timing removes the imposed pooled timing mismatch. Reducing wild sample size increases the chance of failed estimation, so inspect exclusions as well as successful-only summaries. If wild and hatchery also differ in within-period GE, borrowing the split components need not remove bias.

## 08 Bootstrap

### 1. One hundred PIT fish

The fixed-GE interval covers the truth in only 14% of weeks (36% with 600 PIT fish), because it ignores the now much larger GE uncertainty. The GE-draw interval still covers about 95% (94.7%), but its median width grows from 26,772 to 69,649. The mean of the GE-draw replicates now runs 61% high and the median 13.5% high: with an uncertain GE, Jensen's inequality bites hard. Report the median, or carry the whole distribution.

### 2. Two hundred replicates

Coverage barely moves (fixed 35.7%, draws 96.0%, vs 36.0% and 95.3%), and median widths are nearly the same. The coverage estimate from 300 weeks has its own Monte Carlo noise of a couple of points, which is larger than the change from fewer replicates. The one-week interval endpoints shift by up to about 900 fish, partly from fewer replicates and partly from different random draws, so use more replicates for any single interval you report.

## 09 Recover the Truth

All repeats use nonguided GOJ recovery drawn from 1--6% and bypass recovery from 4--15%. These are chosen low-recovery simulation ranges around the talk's setting, not estimated uncertainty intervals. The signed logit contrast is derived from the two drawn probabilities. Recovery priors remain broad; the simulation ranges are not imposed as posterior bounds.

### 1. Five hundred tagged fish

Expect fewer downstream overlaps and less information, but do not assume any particular error ranking or coverage percentage. Compare bias, absolute error, interval width, coverage, and fit failures with the 3,000-fish default. Check whether failures concentrate at the lowest recovery or smallest guided counts. Coverage among successful fits alone may overstate the method's overall reliability.

### 2. Two hundred repeats

More repeats reduce Monte Carlo uncertainty in estimated performance, not uncertainty within each individual data set. For 95% coverage, the rough binomial sampling band is about ±6 percentage points at 50 independent successful repeats and ±3 at 200. Use the number of successful repeats, and report failures alongside conditional coverage. New numerical bias and coverage results await the revised model runs.
