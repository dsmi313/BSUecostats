# Handoff for Claude Code: BSUecostats validation

Dave asked to push completed results immediately because he is running out of usage. This handoff records the remaining work without calling unfinished sampling validated.

## Completed

- Executed all default lessons 00–09, including 07a, every Try This exercise, the optional 20-season hierarchical study, and the optional 100-week day/night study.
- Updated ANSWERS.md with checked numerical results, failure rates, and limits. All four focused tests passed; all R files parse.
- Increased lesson 05 sampling; all four fits pass in each final default/exercise run. Added explicit checks for facility-tagging and no-GOJ fits. Longer sampling does not resolve selection bias or structural confounding.
- Added parameter-level diagnostic reasons and full failed-fit saving to 09, with counts/truth for every replicate. Default: 48/50 pass; 500-tag exercise: 50/50; 200-repeat exercise: 194/200. Never relax diagnostics or quietly omit failed fits.
- Added an explicit wild.ge.multiplier argument for 07a's component-mismatch exercise. Default slide results remain +33.3% pooled versus +1.3% split median season error, among 1,997/2,000 valid seasons.
- Fixed 06's predictive-check axes so both panels show the diagonal. Default predictive fits pass, with rounded Bayesian p-values 0.456 and 0.000.
- Added validation/run_all.py, run_lesson.R, summarize_runs.R, and committed validation CSVs. See validation/README.md.

## Remaining

Large direct-fit CSVs (`fit_diagnostics.csv` and `posterior_summary.csv`) remain in the local workspace; upload them in the follow-up. The answer key, repeated-study summaries, failure tables, source changes, and this handoff were prioritized for the immediate push. Local commit 6e9fc60 contains all 23 files.

1. **Long screw-trap exercise still running at handoff.** The original 10-marked-fish exercise completed but model 2 failed diagnostics on 69 quantities (minimum tail ESS about 91). Its numerical outputs must remain labeled inadequate until the longer fit is checked. The longer run uses 480,000 iterations, 240,000 burn-in, thinning 30, three chains. On the same workspace, inspect:
   `/tmp/bsu-validation/long/04_small_long/run.log` and `result.rds`.
   Background exec session was 5271. Do not assume this process or scratch files survive into another environment. If unavailable, rerun from the repo root:
   ```sh
   python validation/run_all.py --output /tmp/BSUecostats-validation --workers 1 --cases 04_small_long
   ```
   When done, inspect **both** fits and all parameter diagnostics. Update ANSWERS.md section 04.2 using only adequately sampled results; if still failing, report that honestly. Add the completed case to validation/results with summarize_runs.R. Do not remove existing validation records by collecting only a fresh partial directory; preserve or combine the existing results.
2. **A later PPC fit in 06_repeats warns.** All 20 paired comparison seasons passed after retries. The separate rising-detection PPC fit after the loop fails Rhat for p (1.013). Its printed p-value is not validated. The default PPC is a different, passing run. ANSWERS.md and validation records explicitly distinguish these. Investigate/retry separately if desired; do not invalidate the 20 passing paired comparisons or claim this later fit passed.
3. Final reread of ANSWERS.md after the last screw-trap result, then push a follow-up commit. No production model is to be released. These are teaching examples for MS ecology students; Dave represents the agency and deliberately keeps production code private.

## Current environment and data

R 4.3.3; JAGS 4.3.2; jagsUI 1.6.3; rjags 4-18; ggplot2 3.4.4; posterior 1.5.0.

Current repo: `/workspace/scratch/3f4bd71fdb0a/BSUecostats`.
Run folders: `/tmp/bsu-validation/runs`, `reruns`, `final05`, `extra`, `extra07`, `long`, `final06`. Later directories replace earlier cases when collecting results. Most full console logs and plots are outside Git; committed CSVs retain numerical evidence. The 09 local RDS files also contain failed full fits and the simulated inputs.

The root README's ZIP instructions are Dave's requested workflow: extract the ZIP contents directly into the current R working directory. Preserve those instructions.

## Important interpretation

Use rank-normalized split Rhat <= 1.01, bulk and tail ESS >= 400, and finite median MCSE. Passing is a screen, not proof of identifiability or calibration. Repeated-study coverage is conditional on passed fits; report exclusions. Within-season weeks are related, not independent simulation replicates. Enabling optional loops consumes RNG draws and can change later simulated examples. The answer key labels default versus optional results.
