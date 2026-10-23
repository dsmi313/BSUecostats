# How Many Fish Are There? Mark-Recapture with JAGS

**First, download and install JAGS:** https://sourceforge.net/projects/mcmc-jags/ (JAGS is a separate program from R; install it before the R packages, then restart R.)

Companion code for the BSU ecological statistics group talk on IDFG mark-recapture, from screw traps to Lower Granite Dam.

**[View the presentation slides](https://docs.google.com/presentation/d/177eQCKUWrDXGobg0EtnxK9MQ_RDyjjAWeHS4S07U_pA/edit?slide=id.p1#slide=id.p1).**

Most scripts simulate data from known values, fit a model, and check the answer against the truth. Start at 01 and work up: each step adds one idea.

New to Bayesian statistics? Start with [01_first_bayesian_model.R](01_first_bayesian_model.R), and keep [GLOSSARY.md](GLOSSARY.md) open: it explains the fish terms (GE, GRJ, PIT tag, spill) and the statistics terms (prior, posterior, Rhat, logit) in plain language.

**These are teaching examples, not the models IDFG uses for its estimates.** They are simplified so each idea stands on its own. The concepts carry over to the production models. Written with help from Claude. I've tested every script and checked every reference.

## Requirements

- R (4.0 or later)
- JAGS 4.x, installed separately from R: https://sourceforge.net/projects/mcmc-jags/
- The jagsUI and ggplot2 packages: `install.packages(c("jagsUI", "ggplot2"))`

On Windows, install JAGS before jagsUI. If R can't find JAGS, restart R after installing it.

Run [00_setup.R](00_setup.R) first. It installs any missing packages, checks that R can find JAGS, fits a tiny test model, and prints your R, JAGS, and package versions. Save that output with your results so the environment can be rebuilt later.

## Quick start

1. Download this repository using **Code → Download ZIP**, then extract it.
2. Open R or RStudio and set the working directory to the extracted folder.
3. Run the setup check:

   ```r
   source("00_setup.R")
   ```

4. Begin with [01_first_bayesian_model.R](01_first_bayesian_model.R). Run one numbered lesson at a time, reading the comments and inspecting the plots.
5. Keep the [glossary](GLOSSARY.md) handy; see [Using This Code](USING_THIS_CODE.md) to adapt an example.

The numbered R scripts stay together so they can be run directly from this folder.

## Levels

- **Start here, no Bayesian background needed:** `01`, `02`, `03`. The basics: prior, likelihood, posterior, MCMC, and the logit.
- **Next steps:** `04`, `05`. State-space and multistate models, conditioning on detection, and identifiability.
- **Graduate level:** `06`, `07`, `08`, `09`. Hierarchical models, borrowing strength, propagating uncertainty, and simulation studies.

## Scripts

| Script | Idea | Talk section |
|---|---|---|
| [00_setup.R](00_setup.R) | Installs packages, checks JAGS, runs a test model, and records versions. | |
| [01_first_bayesian_model.R](01_first_bayesian_model.R) | Prior, likelihood, and posterior for one efficiency trial, solved exactly by hand, then the same model in JAGS. How to read JAGS output, trace plots, Rhat, and credible intervals. | Background for the whole talk |
| [02_lincoln_petersen.R](02_lincoln_petersen.R) | Trap efficiency from an efficiency trial, catch expanded by efficiency, a +1 adjustment that removes small-sample bias, and Jensen's inequality. Plain R, no JAGS. | Screw traps; Jensen's inequality |
| [03_bernoulli_logit.R](03_bernoulli_logit.R) | A logit GLM in JAGS: each fish is a Bernoulli trial, and spill shifts guidance on the logit scale. | Why model GE on the logit scale |
| [04_state_space_screwtrap.R](04_state_space_screwtrap.R) | A state-space model for daily screw trap passage. Counts alone estimate N × p; adding efficiency trials recovers N. Days with the trap out are estimated anyway. | Screw traps (extension) |
| [05_six_cell.R](05_six_cell.R) | The multistate model at Lower Granite: six detection cells, conditioning on being detected, and a multinomial likelihood. Then two ways it goes wrong: fish tagged after the guidance decision inflate GE, and without a downstream site p and GE can't be separated. | Steps 3 and 4; Gotcha: fish tagged after the decision |
| [06_hierarchical.R](06_hierarchical.R) | Borrowing strength across weeks and rear types. Compares each wild week on its own to the hierarchical fit. | Borrow strength; the hierarchical GE regression |
| [07_day_night.R](07_day_night.R) | Day and night: five cells split by passage time, c5 stays pooled, 11 cells in all. Shows why wild should borrow hatchery's day and night GE, not its pooled GE, when the two pass at night at different rates. | Splitting the cells by day and night; what splitting buys |
| [08_bootstrap.R](08_bootstrap.R) | From GE to abundance: a bootstrap with GE held fixed vs one where every replicate carries its own GE draw. Coverage of each over 300 simulated weeks. Plain R, no JAGS. | Where SCRAPI's interval comes from; Jensen's inequality |
| [09_recover_truth.R](09_recover_truth.R) | A small simulation study: bias and 95% interval coverage over 50 simulated data sets. | Does it recover GE? |

Every lesson script runs top to bottom and prints its results and a plot. Scripts that use JAGS write their models to `.jags` files in the working directory. Run times are a few seconds to a few minutes; `04`, `06`, `07`, and `09` are the slowest. In `06` and `07`, the single example and its plot run by default. The repeated simulations are off: set `run.repeats <- TRUE` at the top of either script to add them.

## Real data to try

The scripts use simulated data, but the same kinds of data are public:

- **PIT tag detections:** PTAGIS, https://www.ptagis.org/
- **River flow and spill at the dams:** Fish Passage Center historical flow and spill query, https://www.fpc.org/smolt/smolt_queries/Q_river_historicflowspill_dataquery.php

[USING_THIS_CODE.md](USING_THIS_CODE.md) explains how to turn fish-level detection records into the cell counts the models need.

## What to read alongside

Kéry, M., and M. Schaub. 2012. *Bayesian Population Analysis using WinBUGS: A Hierarchical Perspective.* Academic Press. The book's code was written for WinBUGS; the authors also provide a translation that runs in JAGS through jagsUI: [BPA code with jagsUI](https://www.vogelwarte.ch/wp-content/uploads/2025/01/BPA-with-jagsUI-Kopie.txt).

Each script ends with a `Where BPA teaches this` section pointing to the matching book sections.

- Chapter 1: the binomial observation process, and why simulate (scripts 02, 09)
- Chapter 2: Bayes rule, priors, posteriors, MCMC, and convergence (script 01)
- Chapter 3: GLMs, link functions, and the logit (scripts 03, 06)
- Chapter 4: random effects (script 06)
- Chapter 5: state-space models for counts (script 04)
- Chapter 6: closed-population size from capture-recapture (script 02)
- Chapter 7: time and group effects, identifiability, and the multinomial likelihood (scripts 05, 06, 09)
- Chapter 9: multistate models (scripts 05, 07)

See [USING_THIS_CODE.md](USING_THIS_CODE.md) for how to adapt these models to your own data, including how to format raw detection histories into the arrays the models need, and a list of gotchas that bite real studies.

**Before using any adapted model on real data, run it through [09_recover_truth.R](09_recover_truth.R).** It's a template, not just an exercise: if a model can't recover known values from simulated data, it won't recover unknown ones from yours.
