# Does the model recover the truth? A small simulation study
# (Run 00_setup.R once first: it installs packages and checks JAGS.)
# (Run time: about 30 to 45 seconds. Times vary up to 2x between computers.)
#
# "Truth is known, so we know what to expect from the analysis"
# (Kery and Schaub 2012, Section 1.5). Simulate data from known values, fit
# the model, and ask two questions over many repeats:
#   Bias: on average, does the estimate land on the truth?
#   Coverage: do 95% credible intervals contain the truth about 95% of the time?
#
# This uses the six-cell model from 05_six_cell.R. Each repeat draws its own
# true parameters, so the check isn't tuned to one lucky setting.
# The entire study uses low GOJ recovery: nonguided 1--6%, bypassed 4--15%.
# These are survival times downstream detection, not survival alone.
# Low overlap counts are part of the problem, not removed in data generation.
#
# Use this as a template, not just an exercise. Any time you change a model
# (new species, new sites, new covariates), check it here before trusting it
# on real data:
#   1. Replace the JAGS model in section 1 with yours.
#   2. Replace the simulation inside the loop with one that generates data
#      from your model, in the same layout as your real data.
#   3. Choose true values that cover the range you expect in the field.
#   4. Run, and check bias, coverage, and how many fits failed to converge.
# If your model can't recover known values, it won't recover unknown ones.
# Expect a few minutes. Requires JAGS and the jagsUI package.

library(jagsUI)
source("helpers.R")
library(ggplot2)

set.seed(7)

# Low recovery is the default for every simulated data set, not an optional
# stress test. Draw both route probabilities directly, then derive the signed
# logit contrast used by the fitted model. This permits bypass recovery to be
# higher than nonguided recovery, as in the talk's MY2025 example.
draw_low_recovery_truth <- function() {
  phi_S <- runif(1, 0.01, 0.06)
  phi_B <- runif(1, 0.04, 0.15)
  c(phi_S = phi_S, phi_B = phi_B, delta = qlogis(phi_S) - qlogis(phi_B))
}

# 1. The model (same as 05_six_cell.R) ----

cat(file = "six_cell.jags", "
model {
psi ~ dbeta(1, 1)
p ~ dbeta(1, 1)
trans ~ dbeta(1, 1)
phi_S ~ dbeta(1, 1)
# Reminder: JAGS's dnorm(mean, precision) uses precision = 1 / SD^2, not SD.
delta ~ dnorm(0, 0.1)
phi_B <- ilogit(logit(phi_S) - delta)

pi_raw[1] <- psi * (1 - trans) * (1 - phi_B)
pi_raw[2] <- psi * (1 - trans) * phi_B
pi_raw[3] <- (1 - psi) * p * (1 - phi_S)
pi_raw[4] <- (1 - psi) * p * phi_S
pi_raw[5] <- (1 - psi) * (1 - p) * phi_S
pi_raw[6] <- psi * trans

p_seen <- sum(pi_raw[1:6])
for (cc in 1:6) { pi_obs[cc] <- pi_raw[cc] / p_seen }
n[1:6] ~ dmulti(pi_obs[1:6], N_seen)
}
")

# 2. Simulate, fit, record ----

n.reps <- 50
n.tagged <- 3000
pars <- c("psi", "p", "trans", "phi_S", "phi_B")
results <- data.frame()
rhat.all <- numeric()
passed <- logical(n.reps)
overlap.counts <- matrix(0L, n.reps, 3, dimnames = list(NULL, c("c2", "c4", "c5")))
counts.all <- matrix(0L, n.reps, 6, dimnames = list(NULL, paste0("c", 1:6)))
truth.all <- matrix(NA_real_, n.reps, length(pars) + 1L,
                   dimnames = list(NULL, paste0("truth.", c(pars, "delta"))))
retried <- logical(n.reps)
diagnostics.all <- vector("list", n.reps)
failed.fits <- list()
failed.diagnostics <- data.frame(rep = integer(), variable = character(), rhat = numeric(),
                                  ess_bulk = numeric(), ess_tail = numeric(), mcse_median = numeric(),
                                  reason = character(), stringsAsFactors = FALSE)

for (i in 1:n.reps) {
  truth <- c(psi = runif(1, 0.03, 0.30), p = runif(1, 0.25, 0.60), trans = runif(1, 0.3, 0.8),
             draw_low_recovery_truth())
  pr <- with(as.list(truth), c(psi * (1 - trans) * (1 - phi_B), psi * (1 - trans) * phi_B,
                                (1 - psi) * p * (1 - phi_S), (1 - psi) * p * phi_S,
                                (1 - psi) * (1 - p) * phi_S, psi * trans,
                                (1 - psi) * (1 - p) * (1 - phi_S)))
  n <- rmultinom(1, n.tagged, pr)[1:6, 1]
  overlap.counts[i, ] <- n[c(2, 4, 5)]
  counts.all[i, ] <- n
  truth.all[i, ] <- truth[c(pars, "delta")]
  # Fitting: this can take a few minutes.
  fit <- safe_jags(list(n = n, N_seen = sum(n)), NULL, c(pars, "delta"), "six_cell.jags",
              n.chains = 3, n.iter = 4000, n.burnin = 1000, n.thin = 3, parallel = FALSE, verbose = FALSE)
  # If a fit errors or fails the Rhat/ESS checks, refit once with four times as many iterations
  # (same chains, same thinning, so four times the draws). Dropping hard data
  # sets without a retry would make the results look better than they are.
  if (!fit_ok(fit)) {
    retried[i] <- TRUE
    fit <- safe_jags(list(n = n, N_seen = sum(n)), NULL, c(pars, "delta"), "six_cell.jags",
                n.chains = 3, n.iter = 16000, n.burnin = 4000, n.thin = 3, parallel = FALSE, verbose = FALSE)
  }
  # Save convergence for every fit. Fits still failing after the retry are
  # left out and counted; the summary describes fits that passed diagnostics.
  rhat <- fit_rhat(fit)
  rhat.all <- c(rhat.all, rhat)
  passed[i] <- fit_ok(fit)
  # Single brackets preserve the replicate index even when diagnostics are NULL.
  diagnostics.all[i] <- list(fit_diagnostics(fit))
  if (!passed[i]) {
    issues <- diagnostic_failures(diagnostics.all[[i]])
    failed.diagnostics <- rbind(failed.diagnostics, data.frame(rep = i, issues, row.names = NULL))
    failed.fits[[as.character(i)]] <- list(rep = i, truth = truth, counts = n,
                                          N_seen = sum(n), fit = fit,
                                          diagnostics = diagnostics.all[[i]])
    cat("x", "")
    next
  }
  for (k in pars) {
    d <- fit$sims.list[[k]]
    results <- rbind(results, data.frame(rep = i, par = k, truth = truth[k], est = mean(d),
                                         lo = quantile(d, 0.025), hi = quantile(d, 0.975)))
  }
  cat(i, "")
}
cat("\n")

# 3. Summarize ----
n.kept <- sum(passed)
cat("Data sets with zero counts in downstream cells (all simulated data sets):\n")
print(colSums(overlap.counts == 0))
cat("Fits failing diagnostics after a retry:", n.reps - n.kept, "of", n.reps,
    "(results below are among the", n.kept, "fits that passed)\n")
cat("Failure rate:", round(100 * mean(!passed), 1), "%\n")
if (nrow(failed.diagnostics)) {
  cat("Failed-fit parameters and reasons (after the longer retry):\n")
  print(failed.diagnostics, row.names = FALSE)
}
# Keep the evidence even if every fit fails and the coverage summary stops.
dir.create("diagnostics", showWarnings = FALSE)
replicate.summary <- data.frame(rep = seq_len(n.reps), passed = passed, retried = retried,
                                max.rhat = rhat.all, counts.all, truth.all)
write.csv(replicate.summary, "diagnostics/09_replicate_status.csv", row.names = FALSE)
write.csv(failed.diagnostics, "diagnostics/09_failed_diagnostics.csv", row.names = FALSE)
saveRDS(list(failed.fits = failed.fits, diagnostics.all = diagnostics.all,
             replicate.summary = replicate.summary, model = readLines("six_cell.jags"),
             n.tagged = n.tagged, parameters = c(pars, "delta"),
             settings = list(n.chains = 3, initial.iter = 4000, initial.burnin = 1000,
                             retry.iter = 16000, retry.burnin = 4000, n.thin = 3),
             session = sessionInfo()), "diagnostics/09_fit_diagnostics.rds")
cat("Saved diagnostics and failed-fit data/objects under diagnostics/.\n")
if (!n.kept) stop("No fits passed diagnostics; inspect fit failures before summarizing.")
results$covered <- results$lo <= results$truth & results$truth <= results$hi
results$error <- results$est - results$truth
summary.table <- do.call(rbind, lapply(split(results, results$par), function(d)
  data.frame(par = d$par[1], mean.error = round(mean(d$error), 4),
             mean.abs.error = round(mean(abs(d$error)), 4),
             coverage.pct = round(100 * mean(d$covered), 1))))
print(summary.table[pars, ], row.names = FALSE)
cat("With", n.kept, "fits kept, coverage within about +/-", round(100 * 1.96 * sqrt(0.95 * 0.05 / n.kept)),
    "points of 95% is consistent with binomial noise.\n")

# 4. Plot: estimate vs truth for each parameter ----
results$par <- factor(results$par, levels = pars)
results$status <- ifelse(results$covered, "Interval covers truth", "Interval misses truth")
p1 <- ggplot(results, aes(truth, est, color = status)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_linerange(aes(ymin = lo, ymax = hi)) +
  geom_point(size = 1.5) +
  facet_wrap(~ par, scales = "free") +
  scale_color_manual(NULL, values = c("Interval covers truth" = "steelblue", "Interval misses truth" = "firebrick")) +
  labs(x = "Truth", y = "Estimate", title = "Estimate vs truth (dashed: 1:1 line)") +
  theme_bw(base_size = 13) + theme(legend.position = "top")
print(p1)

# Try this ----
# Predict before rerunning; reset each change before trying the next.
# 1. Set n.tagged <- 500. Which parameters become less precise? Does any
#    parameter's coverage deteriorate, and how many fits fail diagnostics?
# 2. Set n.reps <- 200. How does the Monte Carlo uncertainty around coverage
#    change, and does the ranking of parameter errors stay the same?

# Lessons ----
# - Simulate from known values, fit, and compare. If a model can't recover
#   known values, it won't recover unknown ones.
# - Check two things over many repeats: bias (does the estimate land on the
#   truth on average?) and coverage (do 95% intervals hold the truth about
#   95% of the time?).
# - With a modest number of repeats, coverage bounces a few points around 95%
#   from chance alone.
# - Low recovery makes overlap cells sparse. Examine which parameters become
#   prior-sensitive; no particular bias direction or coverage is guaranteed.
# - Report diagnostic failure rate alongside conditional bias/coverage. Inspect
#   failed.diagnostics and readRDS("diagnostics/09_fit_diagnostics.rds") rather
#   than silently discarding the hard fits or loosening diagnostic thresholds.

# Where BPA teaches this ----
# - 1.5 in the book: "Truth is known, so we know what to expect from the
#   analysis." Simulated data run through the whole book for this reason.
# - 7.9  Parameter identifiability: extrinsic problems "best studied using
#   simulation"
