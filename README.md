# How Many Fish Are There? Mark-Recapture with JAGS

Companion code for the BSU ecological statistics group talk on IDFG mark-recapture, from screw traps to Lower Granite Dam. **[View the slides](https://docs.google.com/presentation/d/1f0dWv8g_EwbsDDwEKLQE8IRd6JRxob2HT59992U6Y-Y/edit?slide=id.p1#slide=id.p1)**

Most scripts simulate data from known values, fit a model, and check the answer against the truth. **These are teaching examples, not the models IDFG uses for its estimates.** Written with help from Claude.

> **Start here.** Run [00_setup.R](00_setup.R) first, then:
>
> - **New to Bayesian models?** Run [01](01_first_bayesian_model.R), [02](02_lincoln_petersen.R), then [05](05_six_cell.R).
> - **Here for the talk's main results?** Run [05](05_six_cell.R), [06](06_hierarchical.R), then [07a](07a_day_night_closed_form.R).
> - **Have 20 minutes?** Take the [guided route](#a-20-minute-route) below.
> - **Want the standard capture-recapture models?** Run [02a](02a_closed_population.R) (population size), [04a](04a_cjs_survival.R) (survival), and [04b](04b_jolly_seber.R) (both), taken from Kéry and Schaub.

## A 20-minute route

After setup, run [05_six_cell.R](05_six_cell.R) top to bottom once (under a minute), then read these sections in order. Each stop has one thing to notice and one question to answer before moving on.

1. **05, sections 0 and 3b** (4 min). Notice GE fall from 0.12 to 0.11 to 0.05 as the row correction and then GOJ are added. *Why does GOJ cut GE by more than half when the row correction barely moves it?*
2. **05, sections 3 and 3c** (5 min). Notice that c4 and c5 differ only in whether GRS read the fish. *Explain why c5 informs detection, and which assumption in 3c that depends on.*
3. **05, section 5** (3 min). Notice the JAGS estimate of GE, 0.094 (0.078–0.111), against truth 0.10. *Which parameter has the widest interval relative to its true value?*
4. **05, section 8** (3 min). Notice the plot: without GOJ, GE and detection trade off. *What data would separate them?*
5. **[07a_day_night_closed_form.R](07a_day_night_closed_form.R), sections 0 and 4** (5 min). Run it (about 5 seconds). Notice GE of 1% by day and 26% at night, then the +33% error from borrowing hatchery's pooled GE against +1% when split. *Why does borrowing pooled GE overestimate wild abundance when both rears share day and night GE?*

## Setup

1. Install R (4.0 or later), then [JAGS 4.x](https://sourceforge.net/projects/mcmc-jags/). JAGS is a separate program. On Windows, install it before the R packages, and restart R afterward.
2. Download this repository (**Code → Download ZIP**) and extract it.
3. In R, set the working directory to the extracted folder, the one containing `00_setup.R`, and run `source("00_setup.R")`. It installs jagsUI, ggplot2, and posterior, checks that R can find JAGS, and fits a tiny test model.

Keep that folder as the working directory for every script. The JAGS lessons start with `source("helpers.R")`, which loads the shared diagnostic checks; if R says it cannot open `helpers.R`, the working directory is wrong. Run one script at a time and read the comments as you go. Keep the [glossary](GLOSSARY.md) open for the fish terms (GE, GRJ, spill) and the statistics terms (prior, Rhat, logit). Every JAGS fit is screened for convergence, and the script warns if one fails.

## Scripts

Times and printed numbers come from running each script unchanged on 4 November 2026. Times varied up to 2× between machines. Each "plot" link is a reference copy of one of the script's plots, so you can check yours.

### Basics: no Bayesian background needed

- [01_first_bayesian_model.R](01_first_bayesian_model.R): prior, likelihood, and posterior for one efficiency trial, by hand and in JAGS; how to read JAGS output. *Under 5 seconds.* Prints p 0.173 (95% CrI 0.084–0.286). [Plot](figures/01_first_bayesian_model.png)
- [02_lincoln_petersen.R](02_lincoln_petersen.R): the talk's Petersen estimate (plus Bailey's modification), trap-efficiency expansion, and Jensen's inequality. No JAGS. *Under 5 seconds.* Prints Petersen 10,000 and Bailey 9,533, and a +1-adjusted mean of 3,984 against truth 4,000. [Plot](figures/02_lincoln_petersen.png)
- [03_bernoulli_logit.R](03_bernoulli_logit.R): a logit GLM, with spill shifting guidance. *Under 5 seconds.* Prints spill slope −0.71 (−0.96 to −0.47), truth −0.8. [Plot](figures/03_bernoulli_logit.png)

### Capture-recapture models

- [02a_closed_population.R](02a_closed_population.R): population size from capture histories with data augmentation, models M0, Mt, Mb, and Mh (Kéry and Schaub, ch. 6). *About 3–4 minutes.* Prints N 100 (92–111), truth 100; with heterogeneous catchability, M0 gives 89 (86–93) and Mh 104 (90–145). [Plot](figures/02a_closed_population.png)
- [04_state_space_screwtrap.R](04_state_space_screwtrap.R): a state-space model for daily screw trap passage. Counts alone estimate N × p; efficiency trials recover N. *About 2–3 minutes.* Prints season total 6,060 (5,190–7,145), truth 6,284. [Plot](figures/04_state_space_screwtrap.png)
- [04a_cjs_survival.R](04a_cjs_survival.R): survival from recaptures with the Cormack-Jolly-Seber model, the m-array, identifiability, and a fit check (Kéry and Schaub, ch. 7). *About 30–60 seconds.* Prints survival 0.652 (0.577–0.731), truth 0.65; check p-values 0.486 for the right model and 0.005 for the wrong one; then fits real data on Leisler's bats. [Plot](figures/04a_cjs_survival.png)
- [04b_jolly_seber.R](04b_jolly_seber.R): survival, recruitment, and population size at each occasion with the Jolly-Seber model (Kéry and Schaub, ch. 10). *About 3–8 minutes, the slowest script.* Prints superpopulation 403 (366–451), truth 400. [Plot](figures/04b_jolly_seber.png)
- [04c_mark_recovery.R](04c_mark_recovery.R): survival from dead recoveries, like rings reported by the public or tags returned by anglers, then real data on Swiss red kites with two priors (Kéry and Schaub, ch. 8). *Under 10 seconds.* Prints red kite adult survival 0.84 (0.78–0.90); juvenile survival stays wide under either prior. [Plot](figures/04c_mark_recovery.png)
- [05_six_cell.R](05_six_cell.R): the multistate model at Lower Granite: six detection cells, conditioning on detection, and two ways it goes wrong. *10–35 seconds.* Prints the talk's GE steps 0.12 → 0.11 → 0.05, then GE 0.094 (0.078–0.111), truth 0.10. [Plot](figures/05_six_cell.png)

### Graduate level

- [06_hierarchical.R](06_hierarchical.R): borrowing strength across weeks and rear types, with a posterior predictive check. *15–30 seconds.* Prints weekly GE error 0.102 fitting each week alone vs 0.049 hierarchical. [Plot](figures/06_hierarchical.png)
- [07_day_night.R](07_day_night.R): splitting cells by day and night, and why wild should borrow hatchery's split GE, not its pooled GE. *Under 10 seconds.* Prints truth 0.069; wild alone 0.122, borrow pooled 0.047, borrow split 0.083. [Plot](figures/07_day_night.png)
- [07a_day_night_closed_form.R](07a_day_night_closed_form.R): the talk's season day/night GE (1% by day, 26% at night) and its 2,000-season day/night simulation. No JAGS. *About 5 seconds.* Prints median season error +33.3% pooled vs +1.3% split, and saves [the slide figure](figures/day_night_borrowing.png).
- [08_bootstrap.R](08_bootstrap.R): bootstrap intervals with GE fixed vs GE drawn. No JAGS. *About 10 seconds.* Prints 95% coverage 36% vs 95.3%. [Plot](figures/08_bootstrap.png)
- [09_recover_truth.R](09_recover_truth.R): a simulation study of bias and coverage over 50 data sets. Use it as a template before trusting an adapted model on real data. *30–45 seconds.* Prints 48 of 50 fits passing, GE coverage 95.8%. [Plot](figures/09_recover_truth.png)

## Exercises

Each script ends with **Try this** exercises. Predict the result, change one setting, rerun, and check the diagnostics before interpreting it. [ANSWERS.md](ANSWERS.md) gives checked answers.

## More

- [USING_THIS_CODE.md](USING_THIS_CODE.md): adapting a model to your own data, formatting detection histories, and gotchas.
- [TALK_NOTES.md](TALK_NOTES.md): which script matches each slide, and the settings behind the dam examples and the day/night figure.
- [validation/](validation/README.md): rerun every lesson and exercise, quick checks, and how fits are screened.
- Real data: PIT tag detections from [PTAGIS](https://www.ptagis.org/); flow and spill from the [Fish Passage Center](https://www.fpc.org/smolt/smolt_queries/Q_river_historicflowspill_dataquery.php).
- Reading: Kéry, M., and M. Schaub. 2012. *Bayesian Population Analysis using WinBUGS: A Hierarchical Perspective.* Academic Press. The authors' [code for JAGS and jagsUI](https://www.vogelwarte.ch/wp-content/uploads/2025/01/BPA-with-jagsUI-Kopie.txt). Each script ends with a "Where BPA teaches this" section.
