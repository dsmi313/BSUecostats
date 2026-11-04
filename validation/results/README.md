# Validation results: 4 November 2026

Executed 43 cases: all 11 default scripts (00–09, including 07a), every Try This exercise, both optional repeated studies, and three completed longer-chain exercise checks. A fourth longer check, `04_small_long`, is still running at handoff and is not included as a completed result. See `../../CLAUDE_HANDOFF.md`. All completed without execution errors. The four focused tests also passed. R 4.3.3, JAGS 4.3.2, jagsUI 1.6.3, rjags 4-18, ggplot2 3.4.4, and posterior 1.5.0 were used; see `session.txt`.

All 20 directly fitted default examples passed the diagnostic screen. Lesson 09 separately retained 48 of 50 simulated fits, with two failures reported. Completion and sampling adequacy are different outcomes.

## What changed after running the code

- Lesson 05 uses longer chains and checks the mixed-tagging and no-GOJ fits explicitly. All four final fits pass in each default/exercise case. This improves sampling; it does not remove selection bias or structural confounding.
- Lesson 06's predictive-check panels now use identical axes and keep the diagonal reference visible. The matching-data and deliberately misspecified examples retain Bayesian p-values 0.456 and 0.000 in the default run.
- Lesson 07a adds a component-mismatch exercise with explicit wild GE scaling. Its default slide simulation still gives +33.3% pooled versus +1.3% split median season error, conditional on 1,997 successful seasons.
- Lesson 09 saves failed-fit parameters, counts, truth, diagnostics, and full failed objects locally. Checked failure tables from this validation are included here.
- The answer key replaces stale and provisional values with the executed results and their limits.

## Difficult exercises

| Original case | Fits failing the diagnostic screen | Longer check |
| --- | ---: | --- |
| `04_outage` | 2 | `04_outage_long`: 0 failures |
| `04_small` | 1 | `04_small_long` still running at handoff; original warning retained. |
| `06_repeats` (later predictive check) | 1 | Not resolved in this optional run; its 20 paired study fits all passed. Default predictive check passed separately. |
| `06_variable` | 1 | `06_variable_long`: 0 failures |
| `07_small` | 1 | `07_small_long`: 0 failures |

The original warnings remain in the records. Longer-run results are labeled in `ANSWERS.md`; thresholds were not relaxed. Intervals from an inadequately sampled fit should not be treated as checked numerical answers.

## Repeated recovery studies

| Case | Attempted | Passed | Failed |
| --- | ---: | ---: | ---: |
| `09_default` | 50 | 48 | 2 |
| `09_more` | 200 | 194 | 6 |
| `09_small` | 50 | 50 | 0 |

Coverage and errors in `recovery_summary.csv` describe passed fits. They do not establish unconditional calibration when failed fits are excluded. The 20-season hierarchical comparison and 100-week day/night comparison retain paired outcomes in `repeat_status.csv`; their summaries are in the answer key. Repeated weeks within a season share model structure and should not be treated as independent simulation replicates for Monte Carlo uncertainty.

## Files

- `run_status.csv`: execution status, warnings, direct-fit failure counts, and MD5 hashes of executed lesson/helper files. Exercise variants have different source hashes by design; the recorded 06 reruns that changed only plotting also have different hashes.
- `fit_diagnostics.csv`: rank-normalized split Rhat, bulk/tail ESS, and median MCSE for direct fits, including inverse GE where relevant.
- `posterior_summary.csv`: direct-fit posterior summaries, with a pass flag. Includes original flagged fits for transparency.
- `repeat_status.csv`: per-replicate diagnostic decisions for optional studies. The 06 individual pass flags are unavailable; its paired flag is retained.
- `recovery_summary.csv`: conditional bias, mean absolute error, interval width, coverage, and retained counts for 09.
- `failed_recovery_diagnostics.csv`: parameter-level reasons for final 09 failures after retries.
- `session.txt`: software environment.

Use the [runner instructions](../README.md) to regenerate logs, plots, numerical results, and the full 09 failed-fit objects outside the repository. Repeated studies consume random draws, so later examples can change when optional loops are enabled. The checked numerical answers distinguish default and optional settings. Version differences can also change Monte Carlo output.
