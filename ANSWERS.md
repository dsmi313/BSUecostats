# Answers to the Try This Exercises

Try each exercise before reading its answer. Results were checked on 4 November 2026 with R 4.3.3, JAGS 4.3.2, jagsUI 1.6.3, and posterior 1.5.0 (posterior 1.7.0 for the two longer checks run afterward, `04_small_long` and `06_repeats_long`). The default lessons, every exercise, and the optional repeated studies were executed in isolated folders. Numerical JAGS results can vary across software versions. Unless labeled optional, the numbers use `run.repeats <- FALSE`. Enabling repeated studies consumes random draws, so examples later in the same script may use different simulated data.

Diagnostics use rank-normalized split Rhat ≤ 1.01 and bulk/tail ESS ≥ 400, with finite median Monte Carlo error. Passing this screen does not establish identifiability, lack of bias, or adequate model fit. Failed repeated fits remain in the failure counts; performance summaries describe successful fits. See [validation results](validation/results/) and the [rerun instructions](validation/README.md).

## 01 First Bayesian Model

### 1. Ten marked fish, two recaptured

The exact Beta posterior’s 95% credible interval for p widens from 0.084–0.286 (width 0.20) to 0.060–0.518 (width 0.46), about 2.3 times wider. The posterior mean moves from 0.173 to 0.25 even though R / M is 0.2: with only 10 fish, the flat prior pulls the mean toward 0.5. The posterior mean is (R + 1) / (M + 2).

### 2. The Beta(3, 17) prior with ten fish

With 10 fish, the informative prior changes the answer a lot: the mean drops from 0.25 to 0.167 and the interval width from 0.46 to 0.26. With 50 fish it moved the mean only from 0.173 to 0.157. A prior matters most when data are thin.

### 3. Doubling the catch

Expanded passage doubles: the median goes from 3,565 to 7,131 and the interval from 2,102–7,340 to 4,204–14,679. The efficiency posterior doesn't change (mean 0.172), because the catch carries no information about efficiency in this model. Only the marked fish do.

## 02 Lincoln-Petersen

### 1. Fifteen marked fish

854 of 10,000 repeats recapture no marked fish (1 with 50 marked fish), so the simple expansion is infinite in those. Among repeats with recaptures it runs 21% high. The +1 version averages 3,702, 7.5% low, matching the formula's 3,703, and its middle 95% widens to 1,539–9,808.

### 2. Efficiency of 0.03

The +1 correction no longer lands near truth: it averages 3,134, 22% low. The formula N(1 − (1 − p)^(M + 1)) predicts 3,154, so the bias is expected. The bias factor is (1 − p)^(M + 1). The probability of zero recaptures among the M actual releases is (1 − p)^M; 21.2% of these repeats have zero recaptures. The correction cannot replace missing information; increasing the marked sample helps.

## 03 Bernoulli Logit

### 1. Fifty fish

In this run the interval for b1 is −1.54 to −0.17, so it still excludes zero, but it is about 2.8 times wider than with 400 fish (−0.96 to −0.47). With a different seed it could include zero. The fitted curve's uncertainty band also widens.

### 2. No spill effect

The b1 estimate is −0.03 with interval −0.25 to 0.19, which includes zero. The fitted curve is nearly flat around a guidance probability of 0.29. Uncertainty remains in both its level and its slope.

## 04 State-Space Screw Trap

### 1. A 25-day outage

The original 60,000-iteration runs warn about insufficient bulk ESS in both models. The longer check uses `ni <- 240000; nb <- 120000; nt <- 30; nc <- 3`; both fits then pass. For the model using trials, the season-total 95% interval widens from the default 5,190–7,145 to 3,781–11,270. During the outage, average daily interval width is 2.82 times true passage, versus 1.06 on fishing days; all 25 outage-day intervals contain truth. The corresponding default ratios are about 1.42 and 1.00. The model fills the gap through its process assumptions, with much greater uncertainty where observations are missing.

### 2. Ten marked fish per trial

With the script's 60,000 iterations, the model using trials fails diagnostics on 69 quantities (lowest tail ESS 91), so its numbers can't be trusted. The longer check uses `ni <- 480000; nb <- 240000; nt <- 30; nc <- 3`; both fits then pass (largest Rhat 1.002; smallest bulk and tail ESS 1,323 and 1,356). Efficiency gets less precise: its average daily 95% interval width is 0.18, versus 0.13 with 30 marked fish. Passage follows: the season total is 5,943 (95% CrI 4,764–7,450) against truth 6,284, an interval about 37% wider than the default 5,190–7,145. The counts-only model doesn't change (1,235; 1,150–1,332), because it never uses the trials and still estimates N × p, not N. Fewer marked fish cost precision, and the sampler needs far longer chains to explore the less constrained daily efficiency.

## 05 Six-Cell Multistate Model

Default GOJ recovery is 3.25% for nonguided fish and 9.86% for bypassed fish returned to river. These combine survival and downstream detection. In the verified default fit, GE is 0.094 (95% CrI 0.078–0.111), against truth 0.100. Longer default chains now allow all four fits, including the facility-tag and no-GOJ examples, to pass the diagnostic screen.

### 1. Five thousand fish tagged in the facility

Estimated GE rises to 0.494 (0.447–0.536), compared with 0.094 (0.078–0.111) using upstream tags alone. Even adding the default 1,500 facility tags moves it to 0.268 (0.232–0.303). Neither mixed-sample interval contains truth 0.100. The fits pass sampling diagnostics while answering a question distorted by route-based selection.

### 2. A Beta(1, 9) prior on psi without GOJ

With flat priors, GE is 0.142 (0.030–0.214) and detection p is 0.626 (0.114–0.985). Changing only GE's prior to Beta(1, 9) moves them to 0.121 (0.021–0.209) and 0.524 (0.078–0.970). Both fits pass after longer sampling, but the data still identify only `psi / (psi + (1 - psi) * p)`. A prior on one component changes inference for the other along this ridge; sampler agreement cannot resolve that confounding.

### 3. Six hundred tagged fish

Observed c2, c4, and c5 are 2, 3, and 10. Bypass recovery has the widest 95% interval relative to its truth: 0.010–0.197, width about 1.90 times truth 0.099. GE is next: its 0.044–0.160 interval is about 1.16 times truth 0.100. Nonguided recovery spans 0.008–0.044, about 1.10 times truth 0.0325. The original short fit failed Rhat; the revised 24,000-iteration fit passes. Section 7 still adds 1,500 facility tags, pushing GE to 0.662 (0.457–0.793). That section illustrates selection bias, not recovery of the upstream-sample truth.

## 06 Hierarchical Model

In the default season, hierarchical sharing reduces wild weekly mean absolute GE error from 0.102 to 0.049 and average 95% interval width from 0.338 to 0.197. These are summaries of 16 weeks in one simulated season, not independent simulation-based coverage estimates. All four default fits pass diagnostics. The posterior predictive comparison gives 0.456 for matching data and 0.000 for deliberately rising detection fitted with constant detection, correctly flagging the latter mismatch.

### 1. Better-sampled wild weeks

Raising the upper tag limit to 1,500 improves both models. Separate versus hierarchical mean absolute error becomes 0.053 versus 0.032; mean interval width becomes 0.209 versus 0.119. More tags reduce errors and prior dependence, but shared structure still helps in this realized season. All fits pass. The exercise changes random draws and does not create a controlled experiment with exactly the same fish plus extra observations.

### 2. More week-to-week variation

The sigma_psi = 0.8 exercise requires longer sampling: its original hierarchical fit fails Rhat for p (1.013). With 40,000 iterations and 20,000 burn-in, all four fits pass. Separate versus hierarchical mean absolute error is 0.107 versus 0.072, and mean interval width is 0.321 versus 0.288. Estimated weekly SD is 1.037 (95% CrI 0.736–1.454), up from 0.367 in the default season. The model allows more variation and shrinks less; 14 of 16 hierarchical intervals contain truth in this one season, which is not a calibrated coverage study.

### 3. Allowing either spill-slope direction

With negative truth −0.8, the unrestricted fit estimates −0.798 (95% CrI −0.981 to −0.622), close to the restricted default's −0.800 (−0.974 to −0.633).

With positive truth +0.4, the restricted fit is forced against its prior boundary: beta = −0.036 (−0.132 to −0.001). Its estimated weekly SD rises to 0.654, absorbing some mismatch. The unrestricted fit estimates beta = +0.456 (+0.288 to +0.626), with weekly SD 0.353. All fits pass sampling diagnostics. A converged sampler cannot recover a coefficient excluded by its prior.

### Optional: twenty simulated seasons

All 20 paired seasons pass after the allowed retries. Within these seasons:

| Wild tags per week | Model | Mean absolute GE error | Mean 95% interval width | Coverage |
| --- | --- | ---: | ---: | ---: |
| 8–20 | Separate weeks | 0.114 | 0.436 | 86.2% |
| 8–20 | Hierarchical | 0.030 | 0.155 | 95.4% |
| 21–50 | Separate weeks | 0.058 | 0.257 | 91.9% |
| 21–50 | Hierarchical | 0.026 | 0.133 | 93.7% |
| 51–150 | Separate weeks | 0.036 | 0.160 | 89.3% |
| 51–150 | Hierarchical | 0.022 | 0.110 | 94.3% |

Sharing helps most in the thinnest weeks in this simulation. Weeks within a season are related, and only 20 seasons were simulated; do not treat every week as an independent replicate when quantifying Monte Carlo uncertainty.

The additional rising-detection predictive-check fit *after* the optional loop fails Rhat for p (1.013) at the script's 10,000 iterations. That is a different fit from the 20 paired comparisons and the default predictive check. Rerunning with 40,000 iterations and 20,000 burn-in for the two predictive-check fits only (validation case `06_repeats_long`) reproduces the same 20 seasons and table above, and all four direct fits then pass (largest Rhat 1.007). The matching season gives Bayesian p 0.445 (0.456 from the shorter default fit to the same data; the difference is Monte Carlo variation). The rising-detection season gives 0.000: its observed discrepancy (95% CrI 73–83) sits well above the simulated discrepancy (38–57), correctly flagging the mismatch. The original warning is retained in the validation records.

## 07 Day and Night

The default wild-only, pooled, and split fits pass diagnostics. True wild GE is 0.069. Posterior means are 0.122, 0.047, and 0.083; abundance errors based on posterior medians of inverse GE are −40.1%, +48.2%, and −16.0%. A single noisy week need not favor the method that performs best across repeats.

### 1. Equal night shares

Wild and hatchery now have the same true pooled GE, about 0.044. Pooled borrowing estimates 0.047 (95% CrI 0.038–0.056), with −5.4% abundance error. Split borrowing estimates 0.065 (0.041–0.098), with −29.9% error. All fits pass. The systematic timing mismatch disappears, and splitting costs precision in this realized week.

### 2. No night effect

Both true day/night GE are 0.020, so timing cannot change pooled GE. Pooled and split posterior means are 0.0248 and 0.0243; abundance errors are −18.7% and −17.2%. Both intervals contain truth, and all fits pass. Splitting offers no timing correction when there is no period effect.

### 3. Twenty wild fish

Wild-only inference is very uncertain and prior-sensitive. The original pooled fit fails Rhat for p (1.013). With 24,000 iterations and 8,000 burn-in, all seven fits pass. The first wild-only GE mean is 0.298 (95% CrI 0.075–0.600), versus truth 0.069; the prior has substantial influence with so little information. In the separate sample used by the prior table, changing Beta(1, 1) to Beta(1, 9) changes mean GE from 0.058 to 0.039 and the median-based abundance error from +67.8% to +149.8%. A stronger prior need not improve an individual estimate. The separate prior-sensitivity exercise also uses a new simulated sample, so do not confuse its GE estimate with the first week's estimate. Compare the same sample across priors, and inspect both zero guided fish and downstream overlaps.

### Optional: one hundred simulated weeks

All 100 paired weeks pass after the allowed retries; no simulated wild week has zero guided fish in this run.

| Method | Median abundance error | Middle 90% of errors across weeks |
| --- | ---: | ---: |
| Wild only | −31.1% | −59.3% to +39.2% |
| Borrow pooled GE | +50.5% | +26.1% to +78.7% |
| Borrow day/night GE, then blend | −2.4% | −32.5% to +51.5% |

These are sampling-error distributions across repeated data sets, not posterior intervals from one fit. Sharing day/night components removes most systematic timing error under the simulation's shared-component assumptions, but sparse wild histories still leave substantial variation.

## 07a Closed-Form Day/Night Comparison

This example uses low nonguided recovery (123/3,785, about 3.25%). Its default simulation reproduces the slide figure. The optional component-mismatch exercise now has an explicit argument that changes both wild data generation and wild truth.

The verified seed-2025 run retains 1,997 of 2,000 complete seasons. Pooled borrowing has median season error +33.3% (middle 90%: +26.8% to +40.6%); split borrowing has +1.3% (-5.8% to +9.1%). Three wild week-15 samples lack GRS--GOJ overlap. These are conditional simulation-error summaries, not posterior intervals.

### 1. Equal night shares

1,996 of 2,000 complete seasons remain estimable. Median season error is +0.6% for pooled borrowing (middle 90%: −4.8% to +6.7%) and +1.5% for split borrowing (−6.2% to +9.7%). Equal timing removes the imposed pooled mismatch; estimating separate rates can cost precision.

### 2. One-tenth as many wild fish

Only 285 of 2,000 complete seasons remain estimable. Among those, pooled median error is +33.4% (+27.0% to +41.4%) and split median error is +5.4% (−10.4% to +27.0%). These successful-only results omit 85.8% of seasons. Their exclusions are a central result, not evidence that splitting works reliably with such sparse data.

### 3. Wild day/night GE 20% below hatchery GE

Run `run_borrow_simulation(wild.ge.multiplier = 0.8)`. The new truth is 0.8 times the original wild GE. Splitting can correct the timing mixture while still borrowing the wrong components; it cannot remove that biological mismatch.

With 1,995 successful seasons, pooled median error is +6.6% (middle 90%: +1.5% to +12.4%); split median error is −19.1% (−24.7% to −12.9%). Pooled borrowing happens to have offsetting timing and component errors here. Its smaller error does not validate its assumptions.

## 08 Bootstrap

### 1. One hundred PIT fish

The fixed-GE interval covers the truth in only 14.3% of weeks (36% with 600 PIT fish), because it ignores the now much larger GE uncertainty. The GE-draw interval still covers about 95% (94.7%), but its median width grows from 26,772 to 69,649. The mean of the GE-draw replicates now runs 61% high and the median 13.5% high: with an uncertain GE, Jensen's inequality bites hard. The median is less affected by the extreme upper tail, but is not guaranteed unbiased. Report the uncertainty interval with any point summary.

### 2. Two hundred replicates

Coverage barely moves (fixed 35.7%, draws 96.0%, vs 36.0% and 95.3%), and median widths are nearly the same. The coverage estimate from 300 weeks has its own Monte Carlo noise of a couple of points, which is larger than the change from fewer replicates. The one-week interval endpoints shift by up to about 900 fish, partly from fewer replicates and partly from different random draws, so use more replicates for any single interval you report.

## 09 Recover the Truth

All repeats use nonguided GOJ recovery drawn from 1--6% and bypass recovery from 4--15%. These are chosen low-recovery simulation ranges around the talk's setting, not estimated uncertainty intervals. The signed logit contrast is derived from the two drawn probabilities. Recovery priors remain broad; the simulation ranges are not imposed as posterior bounds.

### Default: three thousand tagged fish, fifty repeats

48 of 50 fits pass after the allowed retry; failures are 4%. Replicates 30 and 44 fail Rhat for GE and p, with values around 1.013. Neither has a zero downstream overlap cell. Sparse overlap and diagnostic failure are related concerns, but zero counts do not explain every failure.

| Parameter | Mean error | Mean absolute error | 95% interval coverage among 48 passed fits |
| --- | ---: | ---: | ---: |
| GE, psi | +0.0035 | 0.0126 | 95.8% |
| GRS detection, p | +0.0048 | 0.0336 | 97.9% |
| Transport | −0.0015 | 0.0202 | 91.7% |
| Nonguided recovery | +0.0010 | 0.0042 | 95.8% |
| Bypass recovery | −0.0041 | 0.0164 | 95.8% |

These are absolute probability errors. With 48 retained independent data sets, a rough Monte Carlo band around 95% coverage is ±6 percentage points. The failed fits are saved for inspection; they are not silently called accurate or counted as covered.

### 1. Five hundred tagged fish

All 50 fits pass after the allowed retries, but that does not mean they contain more information than the 3,000-tag fits. Four samples have zero c2 counts. Mean 95% interval widths roughly double: GE 0.155 versus 0.072, p 0.442 versus 0.214, and bypass recovery 0.191 versus 0.084. Mean errors also grow, especially for p. Failure counts need not increase monotonically in a finite run with broad, prior-influenced posteriors.

| Parameter | Mean error | Mean absolute error | 95% interval coverage among 50 passed fits |
| --- | ---: | ---: | ---: |
| GE, psi | +0.0202 | 0.0355 | 94% |
| GRS detection, p | +0.0638 | 0.0985 | 96% |
| Transport | −0.0063 | 0.0344 | 100% |
| Nonguided recovery | +0.0065 | 0.0105 | 94% |
| Bypass recovery | +0.0081 | 0.0354 | 96% |

Both zero overlaps and prior sensitivity deserve attention even when every sampler passes its screen.

### 2. Two hundred repeats

More repeats reduce Monte Carlo uncertainty in estimated performance, not uncertainty within each data set. A rough band around 95% coverage is ±3 points at 200 successful independent repeats. Use the actual retained count, and report failures alongside conditional coverage.

With 200 repeats, 194 fits pass after retries and six fail (3%). There are no zero c2/c4/c5 cells. Among passed fits:

| Parameter | Mean error | Mean absolute error | 95% interval coverage |
| --- | ---: | ---: | ---: |
| GE, psi | +0.0034 | 0.0133 | 94.8% |
| GRS detection, p | +0.0149 | 0.0403 | 94.8% |
| Transport | −0.0037 | 0.0202 | 93.3% |
| Nonguided recovery | +0.0013 | 0.0044 | 93.3% |
| Bypass recovery | −0.0029 | 0.0195 | 94.8% |

The error ranking changes slightly as more data sets are simulated. None of these conditional coverage percentages differs convincingly from 95% at this simulation size. That does not establish calibration in other parameter ranges or repair the six failed fits.
