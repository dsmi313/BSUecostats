# Running the lesson and answer checks

Run `00_setup.R` first. From the repository root, with Python 3.8 or later available:

```sh
python validation/run_all.py --output /tmp/BSUecostats-validation --workers 2
Rscript validation/summarize_runs.R validation/results /tmp/BSUecostats-validation
```

On Windows, replace the output path with a folder on your computer. The runner copies each lesson and `helpers.R` to a separate case folder, applies one exercise change, and runs the whole lesson. It includes the optional 20-season and 100-week studies. These checks take substantially longer than the focused tests. Each JAGS lesson may use three cores; lower `--workers` if needed.

Use `--cases 05_default 09_small` to run selected cases. Each case saves a console log, plots, numerical summaries, diagnostics, and software versions. Direct-fit records are checkpointed after each fit, so a later interruption need not erase earlier diagnostic evidence. Lesson 09 also saves its own failed-fit records. Keep run folders outside the repository. To combine a full run and later reruns, give both directories to `summarize_runs.R`; the later directory takes precedence for duplicate case names.

A zero exit status means execution completed. It does **not** mean every fit passed diagnostics. Check `run_status.csv`, `fit_diagnostics.csv`, `repeat_status.csv`, and the failure counts in `recovery_summary.csv`. Repeated-study coverage is conditional on the fits that passed. The answer key explains these limits. Source MD5 hashes identify each executed lesson variant and helper file.

The runner also includes explicit longer-chain checks for the difficult 04 outage/trial-size, 06 variation, 06 repeated-study predictive check, and 07 sparse-wild exercises. They are separate cases, so the original warnings remain visible. In `repeat_status.csv`, lesson 06 reports the paired pass decision it retains; individual pass flags are unavailable there and left blank. Lesson 07 retains both individual and paired pass decisions. `posterior_summary.csv` contains the direct fits’ estimates and intervals, including flagged fits; interpret its `passed` column.
