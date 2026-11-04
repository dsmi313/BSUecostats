# Running the lesson and answer checks

Run `00_setup.R` first. From the repository root:

```sh
Rscript validation/run_all.R --output /tmp/BSUecostats-validation --workers 2
Rscript validation/summarize_runs.R validation/results /tmp/BSUecostats-validation
```

On Windows, replace the output path with a folder on your computer. The runner copies each lesson and `helpers.R` to a separate case folder, applies one exercise change, and runs the whole lesson. It includes the optional 20-season and 100-week studies. These checks take substantially longer than the focused tests. Each JAGS lesson may use three cores; lower `--workers` if needed.

Use `--cases 05_default 09_small` to run selected cases. Each case saves a console log, plots, numerical summaries, diagnostics, and software versions. The R runner uses only base R and the bundled parallel package, supports Windows and Unix, and checkpoints execution status to manifest.csv and manifest.rds after each batch of cases. Direct-fit records are checkpointed after each fit, so a later interruption need not erase earlier diagnostic evidence. Lesson 09 also saves its own failed-fit records. Keep run folders outside the repository. To combine a full run and later reruns, give both directories to `summarize_runs.R`; the later directory takes precedence for duplicate case names.

A zero exit status means execution completed. It does **not** mean every fit passed diagnostics. Check `run_status.csv`, `fit_diagnostics.csv`, `repeat_status.csv`, and the failure counts in `recovery_summary.csv`. Repeated-study coverage is conditional on the fits that passed. The answer key explains these limits. Source MD5 hashes identify each executed lesson variant and helper file.

The runner also includes explicit longer-chain checks for the difficult 04 outage/trial-size, 06 variation, 06 repeated-study predictive check, and 07 sparse-wild exercises. They are separate cases, so the original warnings remain visible. In `repeat_status.csv`, lesson 06 reports the paired pass decision it retains; individual pass flags are unavailable there and left blank. Lesson 07 retains both individual and paired pass decisions. `posterior_summary.csv` contains the direct fits’ estimates and intervals, including flagged fits; interpret its `passed` column.

## How the lessons check fits

`helpers.R` screens each JAGS fit with rank-normalized split Rhat (<= 1.01), bulk and tail effective sample sizes (>= 400), and finite median Monte Carlo error, and the lessons warn when a fit fails. These are screens, not proof of convergence: inspect trace plots too. Passing them does not establish identifiability, lack of bias, or adequate fit.

In the repeated simulation studies (`06` and `07` with `run.repeats <- TRUE`, and `09`), a fit that errors or fails the screen is retried once and counted if it still fails. Comparisons remain paired. Lesson 07 compares posterior-median abundance expansions because the wild-only inverse-GE mean can be infinite when no guided fish are observed. Lesson 09 prints its failure rate and saves replicate statuses and failed-fit diagnostics in `diagnostics/`, which Git ignores.

## Quick checks

After setup, run these from the repository root. They take seconds and do not replace running the lessons.

- `Rscript tests/validation_runner.R`: R-only runner scheduling, exercise patches, logs, and failure handling (no JAGS needed).
- `Rscript tests/audit_checks.R`: diagnostic failures, detection-history validation, and the finite median used for inverse-GE expansion.
- `Rscript tests/talk_examples.R`: the slide arithmetic (base R only).
- `Rscript tests/day_night_closed_form.R`: the closed-form day/night example (base R only).
- `Rscript tests/low_recovery.R`: recovery ranges and route contrasts, without fitting JAGS models.
