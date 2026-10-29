# Answers to the Try This Exercises

Try each exercise before reading its answer. Every number below comes from one run of the script with that single change, using the script's own seed. JAGS results can differ slightly on another computer or package version, so expect the same pattern, not identical digits. Several exercises run one simulated week or season, so some differences are partly chance. Where that matters, the answer says so.

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

### 1. Five thousand fish tagged in the facility

Estimated GE is 0.506 (95% CrI 0.493–0.519) against a truth of 0.10, compared with 0.277 with 1,500 facility fish. The interval is narrow and nowhere near the truth, so it gives no warning. More wrongly tagged fish make the model more confidently wrong.

### 2. A Beta(1, 9) prior on psi without GOJ

p's posterior moves from a mean of 0.647 (0.161–0.983) to 0.532 (0.086–0.970), even though its own prior and the data are unchanged. Without GOJ the data pin down only a combination of psi and p, so a prior on one shifts the other along the ridge. psi's mean moves from 0.146 to 0.122. The without-GOJ fit's Rhat is 1.025, another hint that something is off.

### 3. Six hundred tagged fish

phi_B has the widest interval relative to its true value (0.245–0.608, width 0.9 times truth, vs 0.3 times with 6,000 fish). psi is next (0.55 times). The sparse histories are c1 (15 fish) and c2 (11 fish), the bypassed fish returned to the river, and they are what informs phi_B. Section 7 still adds 1,500 facility fish, now more than twice the tagged sample, so its GE jumps to 0.749.

## 06 Hierarchical Model

### 1. Better-sampled wild weeks

Not every week is thick: wild fish tagged now range from 15 to 1,255, with 12 of 16 weeks above 50. The hierarchical advantage shrinks: mean absolute error for wild GE is 0.033 on its own vs 0.026 hierarchical (0.057 vs 0.031 before), and mean interval width is 0.161 vs 0.116 (0.258 vs 0.173). In weeks above 50 fish the errors are 0.032 vs 0.023. The four remaining thin weeks still get much narrower intervals from the hierarchical model (0.30 vs 0.16).

### 2. More week-to-week variation

The hierarchical model estimates sigma_psi at 1.08, so it shrinks weeks less toward their stratum. Its advantage shrinks: mean absolute error is 0.058 on its own vs 0.045 hierarchical, and interval width is 0.255 vs 0.227. Hierarchical coverage is 15 of 16 weeks in this one season. Borrowing helps less when weeks really are different.

## 07 Day and Night

### 1. Equal night shares

True wild and hatchery GE are now both 0.044. Borrowing pooled GE no longer has the timing problem: its wild abundance error falls from +36% to −12%. The remaining error comes from the shared estimate itself (0.050 this week), not from timing. Borrowing split GE does worse this week (−35%) because wild estimates its own night share from 59 fish. With equal night shares, splitting costs precision and buys nothing. One week is noisy; set `run.repeats <- TRUE` to see the average.

### 2. No night effect

With day and night GE equal (both 0.02), night share doesn't change pooled GE, so true wild and hatchery GE are the same. Both borrowing methods land close: 0.023 pooled (−12% abundance error) and 0.021 split (−6%). Wild alone estimates 0.072 (−70%), pulled up by its flat prior with few guided fish.

### 3. Twenty wild fish

Only 16 wild fish are detected. Wild alone gives GE 0.219 with an interval of 0.067–0.434 and a −67% abundance error. Borrowing split narrows that to 0.061–0.204 (−43%). Borrowing pooled is unchanged at 0.050 (+37%), because hatchery fish dominate it. This week happened to have guided wild fish; with 20 fish, some weeks won't, and `run.repeats <- TRUE` counts them. All fits passed diagnostics.

## 08 Bootstrap

### 1. One hundred PIT fish

The fixed-GE interval covers the truth in only 14% of weeks (36% with 600 PIT fish), because it ignores the now much larger GE uncertainty. The GE-draw interval still covers about 95% (94.7%), but its median width grows from 26,772 to 69,649. The mean of the GE-draw replicates now runs 61% high and the median 13.5% high: with an uncertain GE, Jensen's inequality bites hard. Report the median, or carry the whole distribution.

### 2. Two hundred replicates

Coverage barely moves (fixed 35.7%, draws 96.0%, vs 36.0% and 95.3%), and median widths are nearly the same. The coverage estimate from 300 weeks has its own Monte Carlo noise of a couple of points, which is larger than the change from fewer replicates. The one-week interval endpoints shift by up to about 900 fish, partly from fewer replicates and partly from different random draws, so use more replicates for any single interval you report.

## 09 Recover the Truth

### 1. Five hundred tagged fish

Every parameter gets less precise. Mean absolute error grows 2.4 to 3.2 times (psi 0.006 to 0.015, phi_B 0.026 to 0.068). Coverage stays between 94% and 98%, within binomial noise of 95%, and no fits failed diagnostics. With fewer fish the model is less precise but still honest about it.

### 2. Two hundred repeats

The binomial noise band around 95% shrinks from about ±6 to about ±3 points, and coverage lands between 94.5% and 97%. The ranking of mean absolute errors stays the same (psi smallest, then phi_S, p, trans, and phi_B largest). With more repeats, phi_B's small upward bias from the 50-repeat run (+0.014) mostly goes away (+0.002), so it was largely chance.
