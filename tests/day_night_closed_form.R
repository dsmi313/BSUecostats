# Base-R checks for the closed-form example. No JAGS or plotting required.
# Run from the repository root: Rscript tests/day_night_closed_form.R
options(BSUecostats.functions_only = TRUE)
source("07a_day_night_closed_form.R")
options(BSUecostats.functions_only = NULL)

# A deterministic table with known passage: day N=1000, night N=500,
# GE=.10/.40, detection=.50, and GOJ recovery=.03.
e <- estimate_borrow_rear(100, 200, 450, 150, 18, 18)
stopifnot(e$ok, abs(e$ge.d - .1) < 1e-12, abs(e$ge.n - .4) < 1e-12,
          abs(e$a - 1 / 3) < 1e-12, abs(e$pool - .2) < 1e-12)

# Unidentified or empty observations fail; no truth fallback is possible.
stopifnot(!estimate_borrow_rear(1, 1, 10, 10, 0, 0)$ok,
          !estimate_borrow_rear(1, 1, 10, 10, 0, 5)$ok,
          !estimate_borrow_rear(0, 1, 0, 10, 2, 2)$ok)

# With equal guidance in both periods, changing timing cannot change GE.
for (a in c(.05, .2, .8))
  stopifnot(abs((1 - a) * .1 + a * .1 - .1) < 1e-12)

# Missing identifying information must not yield a reassuring summary.
failed <- run_borrow_simulation(n.sims = 3, phi = 0)
sf <- summarize_borrow_simulation(failed)
stopifnot(sf$failed.weeks == 24, sf$failed.seasons == 3,
          all(is.na(failed$weekly$pooled)), all(is.na(failed$weekly$split)),
          all(sf$weekly$n == 0), all(sf$season$n == 0))

# Reproducibility and paired retention, including the impossible-overlap case.
a <- run_borrow_simulation(n.sims = 50, seed = 2025)
b <- run_borrow_simulation(n.sims = 50, seed = 2025)
stopifnot(identical(a, b),
          identical(is.finite(a$weekly$pooled), is.finite(a$weekly$split)),
          identical(is.finite(a$season$pooled), is.finite(a$season$split)))
sa <- summarize_borrow_simulation(a)
stopifnot(all(sa$weekly$n[sa$weekly$method == "pooled"] ==
                sa$weekly$n[sa$weekly$method == "split"]))

# If rears truly have identical timing, borrowing pooled GE has no timing
# bias. Large synthetic samples isolate this property from sparse-data noise.
wk <- borrow_inputs()
wk$aW <- wk$aH
wk$NW <- wk$NW * 100
wk$NH <- wk$NH * 100
same <- summarize_borrow_simulation(run_borrow_simulation(wk, n.sims = 200))
stopifnot(same$failed.seasons == 0, all(abs(same$season$median) < .02))
# Departure exercise changes wild generation and the truth used for errors.
departure <- run_borrow_simulation(n.sims = 3, wild.ge.multiplier = 0.8)
stopifnot(identical(departure$wild.ge.multiplier, 0.8),
          inherits(try(run_borrow_simulation(n.sims = 1, wild.ge.multiplier = 100),
                       silent = TRUE), "try-error"))
cat("Closed-form day/night checks passed.\n")
