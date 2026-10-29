# Day and night: why split GE before wild borrows from hatchery
# (Run 00_setup.R once first: it installs packages and checks JAGS.)
# (Model fits can take a few minutes to run; R will look busy while JAGS samples.)
#
# What this script should leave you with: a pooled rate is a blend, and the
# blend depends on the mix. Wild and hatchery can behave identically by day
# and by night and still have different pooled GE, just because they arrive
# at different hours. Borrow the parts, not the blend.
#
# Guidance differs sharply between day and night. Fish detected at Lower
# Granite (GRJ or GRS) have a passage time, so cells c1, c2, c3, c4, and c6
# each split into a day and a night version. A c5 fish was seen only at GOJ,
# so nobody knows when it passed Lower Granite; c5 stays pooled.
# 5 cells x 2 + 1 = 11 cells.
#
#   a           share of fish passing at night
#   psi_day     daytime guidance
#   psi_night   logit(psi_night) = logit(psi_day) + night_eff
#   GE          the blend (1 - a) psi_day + a psi_night
#
# Assumption: GRS detection p is the same day and night. A single group
# cannot freely identify both detection rates when c5 has unknown timing.
# Sharing guidance across rear types with different night shares can supply
# additional identifying information, but precision may still be poor.
# Equal p is a simplifying assumption here, not a universally untestable one.
#
# Why it matters: 06_hierarchical.R has wild fish borrow strength from the
# far more numerous hatchery fish. Here wild and hatchery share the same true
# day and night GE, but wild fish pass at night more often. So wild's pooled
# GE really is higher than hatchery's, only because of timing. Three ways to
# estimate wild GE for one week:
#   Wild alone:       six-cell model on wild fish only (sparse, without borrowing)
#   Borrow pooled:    wild and hatchery forced to share one pooled psi
#   Borrow split:     wild and hatchery forced to share psi_day and psi_night;
#                     each rear keeps its own night share a
# The two borrowing models are limiting cases: complete pooling, where the
# shared values are identical for wild and hatchery. 06_hierarchical.R uses
# partial pooling, where wild keeps its own values and is pulled toward the
# shared structure. The extremes make the timing problem easy to see; the
# same problem leaks into partial pooling whenever pooled GE is what's shared.
# Abundance is catch divided by GE, so the script reports the error in wild
# abundance, using the posterior median of 1/GE (02_lincoln_petersen.R).
# This remains finite when no guided wild fish are observed. In that case,
# the wild-alone model's flat GE prior gives an infinite E[1/GE]; a finite
# MCMC average would misleadingly print a finite answer. We therefore use
# the same median-based abundance summary for all three methods.
#
# Teaching version: one week at a time, no covariates. The talk adds a night
# offset to the hierarchical daytime regression (including spill and outflow).
# Here daytime GE is a single probability, so this is not that regression.
# Section 5 repeats the
# week 100 times (300 more model fits, several minutes). It's off by default:
# set run.repeats to TRUE below to run it.
# All data are simulated. Requires JAGS and the jagsUI and ggplot2 packages.

library(jagsUI)
source("helpers.R")
library(ggplot2)

# 0. Simpson's paradox: the talk's two-week toy table ----
simpson <- data.frame(week = c("A", "B"), fish = c(1000, 1000),
                      day.ge = c(0.04, 0.02), night.ge = c(0.30, 0.25),
                      night.share = c(0.10, 0.40))
simpson$pooled.ge <- with(simpson, (1 - night.share) * day.ge + night.share * night.ge)
simpson$guided <- simpson$fish * simpson$pooled.ge
print(simpson)
stopifnot(isTRUE(all.equal(simpson$pooled.ge, c(0.066, 0.112))))
# Both within-period rates fall, but pooled GE rises because the timing mix
# changes. This illustrates aggregation, not an estimated causal spill effect.

set.seed(6)

run.repeats <- FALSE   # TRUE adds section 5, the 100-week repeat (several minutes)

# 1. True values ----

psi_day <- 0.02
night_eff <- 3.0
psi_night <- plogis(qlogis(psi_day) + night_eff)
a.true <- c(wild = 0.18, hatchery = 0.09)          # wild pass at night about twice as often
n.tagged <- c(wild = 80, hatchery = 6000)          # wild is thin
trans <- 0.60; phi_S <- 0.50; delta <- 0.40; p <- 0.40
phi_B <- plogis(qlogis(phi_S) - delta)
ge.true <- (1 - a.true) * psi_day + a.true * psi_night
print(round(ge.true, 3))

cell.probs <- function(a) c(
  (1 - a) * psi_day * (1 - trans) * (1 - phi_B),          # c1 day
  (1 - a) * psi_day * (1 - trans) * phi_B,                # c2 day
  (1 - a) * (1 - psi_day) * p * (1 - phi_S),              # c3 day
  (1 - a) * (1 - psi_day) * p * phi_S,                    # c4 day
  (1 - a) * psi_day * trans,                              # c6 day
  a * psi_night * (1 - trans) * (1 - phi_B),              # c1 night
  a * psi_night * (1 - trans) * phi_B,                    # c2 night
  a * (1 - psi_night) * p * (1 - phi_S),                  # c3 night
  a * (1 - psi_night) * p * phi_S,                        # c4 night
  a * psi_night * trans,                                  # c6 night
  ((1 - a) * (1 - psi_day) + a * (1 - psi_night)) * (1 - p) * phi_S,          # c5 pooled
  ((1 - a) * (1 - psi_day) + a * (1 - psi_night)) * (1 - p) * (1 - phi_S))    # never

simulate.week <- function() {
  # Data layout: n11[rear, cell], 2 x 11; rear 1 = wild, 2 = hatchery.
  # Cell order: c1 c2 c3 c4 c6 by day (1 to 5), the same by night (6 to 10),
  # then c5 pooled (11). Day or night comes from the fish's detection time at
  # Lower Granite (GRJ or GRS); c5 fish have none, so they stay pooled.
  # From fish-level records, add a day/night column to the cell from
  # 05_six_cell.R, section 2, and tabulate in this order.
  n11 <- t(sapply(1:2, function(r) rmultinom(1, n.tagged[r], cell.probs(a.true[r]))[1:11, 1]))
  # Collapse to six pooled cells: c1 c2 c3 c4 c5 c6
  n6 <- cbind(n11[, 1] + n11[, 6], n11[, 2] + n11[, 7], n11[, 3] + n11[, 8],
              n11[, 4] + n11[, 9], n11[, 11], n11[, 5] + n11[, 10])
  list(n11 = n11, n6 = n6)
}

# 2. The three models ----

# Wild alone: the six-cell model from 05_six_cell.R, fit to wild fish only
cat(file = "wild_alone.jags", "
model {
# Priors
psi ~ dbeta(1, 1)                     # wild GE, pooled over day and night
p ~ dbeta(1, 1)                       # GRS detection
trans ~ dbeta(1, 1)                   # share of guided fish barged
phi_S ~ dbeta(1, 1)                   # spilled fish: survive to GOJ and detected there
# Reminder: JAGS's dnorm(mean, precision) uses precision = 1 / SD^2, not SD.
delta ~ dnorm(0, 0.1)                 # bypassed fish penalty, logit scale
phi_B <- ilogit(logit(phi_S) - delta)

# Cell probabilities
pi_raw[1] <- psi * (1 - trans) * (1 - phi_B)    # c1 GRJ only
pi_raw[2] <- psi * (1 - trans) * phi_B          # c2 GRJ and GOJ
pi_raw[3] <- (1 - psi) * p * (1 - phi_S)        # c3 GRS only
pi_raw[4] <- (1 - psi) * p * phi_S              # c4 GRS and GOJ
pi_raw[5] <- (1 - psi) * (1 - p) * phi_S        # c5 GOJ only
pi_raw[6] <- psi * trans                        # c6 transported

# Condition on being detected
p_seen <- sum(pi_raw[1:6])
for (cc in 1:6) { pi_obs[cc] <- pi_raw[cc] / p_seen }

# Each detected wild fish is one roll of a six-sided die
n[1:6] ~ dmulti(pi_obs[1:6], N_seen)

wild_ge <- psi
}
")

# Borrow pooled: one pooled psi forced to be the same for wild (row 1) and
# hatchery (row 2)
cat(file = "borrow_pooled.jags", "
model {
# Priors
psi ~ dbeta(1, 1)                     # one pooled GE for both rears
p ~ dbeta(1, 1)                       # GRS detection
trans ~ dbeta(1, 1)                   # share of guided fish barged
phi_S ~ dbeta(1, 1)                   # spilled fish: survive to GOJ and detected there
delta ~ dnorm(0, 0.1)                 # bypassed fish penalty, logit scale
phi_B <- ilogit(logit(phi_S) - delta)

# Cell probabilities, the same for both rears
pi_raw[1] <- psi * (1 - trans) * (1 - phi_B)    # c1 GRJ only
pi_raw[2] <- psi * (1 - trans) * phi_B          # c2 GRJ and GOJ
pi_raw[3] <- (1 - psi) * p * (1 - phi_S)        # c3 GRS only
pi_raw[4] <- (1 - psi) * p * phi_S              # c4 GRS and GOJ
pi_raw[5] <- (1 - psi) * (1 - p) * phi_S        # c5 GOJ only
pi_raw[6] <- psi * trans                        # c6 transported

# Condition on being detected
p_seen <- sum(pi_raw[1:6])
for (cc in 1:6) { pi_obs[cc] <- pi_raw[cc] / p_seen }

# Each detected fish is one roll of a six-sided die; wild and hatchery
# share the same die
n[1,1:6] ~ dmulti(pi_obs[1:6], N_seen[1])
n[2,1:6] ~ dmulti(pi_obs[1:6], N_seen[2])

wild_ge <- psi
}
")

# Borrow split: psi_day and psi_night forced to be the same for wild and
# hatchery; each rear keeps its own night share
# Model in words: one daytime GE and a night shift are shared by wild and
# hatchery. Each rear has its own night share a, estimated from its own
# detection times. Each rear's 11 cells come from those pieces, conditioned
# on being detected, and its detected fish are rolls of an 11-sided die.
# Wild GE is the shared day and night GE blended with wild's own night share.
cat(file = "borrow_split.jags", "
model {
# Priors
psi_day ~ dbeta(1, 1)                 # daytime GE, shared
night_eff ~ dnorm(0, 0.1)             # night shift on the logit scale, shared
psi_night <- ilogit(logit(psi_day) + night_eff)
p ~ dbeta(1, 1)                       # GRS detection
trans ~ dbeta(1, 1)                   # share of guided fish barged
phi_S ~ dbeta(1, 1)                   # spilled fish: survive to GOJ and detected there
delta ~ dnorm(0, 0.1)                 # bypassed fish penalty, logit scale
phi_B <- ilogit(logit(phi_S) - delta)

for (r in 1:2) {
   a[r] ~ dbeta(1, 1)                 # night share for rear r
   # day: cells 1 to 5 (c1 c2 c3 c4 c6); night: 6 to 10; c5 pooled: 11
   pi_raw[r,1] <- (1 - a[r]) * psi_day * (1 - trans) * (1 - phi_B)
   pi_raw[r,2] <- (1 - a[r]) * psi_day * (1 - trans) * phi_B
   pi_raw[r,3] <- (1 - a[r]) * (1 - psi_day) * p * (1 - phi_S)
   pi_raw[r,4] <- (1 - a[r]) * (1 - psi_day) * p * phi_S
   pi_raw[r,5] <- (1 - a[r]) * psi_day * trans
   pi_raw[r,6] <- a[r] * psi_night * (1 - trans) * (1 - phi_B)
   pi_raw[r,7] <- a[r] * psi_night * (1 - trans) * phi_B
   pi_raw[r,8] <- a[r] * (1 - psi_night) * p * (1 - phi_S)
   pi_raw[r,9] <- a[r] * (1 - psi_night) * p * phi_S
   pi_raw[r,10] <- a[r] * psi_night * trans
   pi_raw[r,11] <- ((1 - a[r]) * (1 - psi_day) + a[r] * (1 - psi_night)) * (1 - p) * phi_S
   # Condition on being detected
   p_seen[r] <- sum(pi_raw[r,1:11])
   for (cc in 1:11) { pi_obs[r,cc] <- pi_raw[r,cc] / p_seen[r] }
   # Each detected fish is one roll of an 11-sided die
   n[r,1:11] ~ dmulti(pi_obs[r,1:11], N_seen[r])
   }

# Wild GE: shared day and night GE, blended with wild's own night share
wild_ge <- (1 - a[1]) * psi_day + a[1] * psi_night
}
")

# Parameters to save. The shared parameters are saved too, so convergence is
# checked on the whole model, not just on wild GE.
params.alone  <- c("wild_ge", "psi", "p", "trans", "phi_S", "delta")
params.pooled <- c("wild_ge", "psi", "p", "trans", "phi_S", "delta")
params.split  <- c("wild_ge", "psi_day", "night_eff", "a", "p", "trans", "phi_S", "delta")

# Percent error in median expanded wild abundance: catch is held fixed.
# This isolates GE estimation, not all sources of abundance uncertainty.
abund.error <- function(fit) 100 * (median(1 / fit$sims.list$wild_ge) * ge.true["wild"] - 1)

labs3 <- c("Wild alone", "Borrow pooled", "Borrow split")
cols3 <- c("Wild alone" = "gray40", "Borrow pooled" = "darkorange", "Borrow split" = "steelblue")

# 3. One week ----

wk <- simulate.week()
cat("Zero guided wild fish this week:", sum(wk$n6[1, c(1, 2, 6)]) == 0, "\n")
cat("Detected fish, wild:", sum(wk$n11[1, ]), "  hatchery:", sum(wk$n11[2, ]), "\n")

# Data for each model
data.alone  <- list(n = wk$n6[1, ], N_seen = sum(wk$n6[1, ]))     # wild only, 6 cells
data.pooled <- list(n = wk$n6, N_seen = rowSums(wk$n6))           # wild and hatchery, 6 cells
data.split  <- list(n = wk$n11, N_seen = rowSums(wk$n11))         # wild and hatchery, 11 cells

# Fitting: this can take a few minutes.
fit.alone  <- jags(data.alone, NULL, params.alone, "wild_alone.jags",
                   n.chains = 3, n.iter = 6000, n.burnin = 2000, n.thin = 2, parallel = TRUE)
fit.pooled <- jags(data.pooled, NULL, params.pooled, "borrow_pooled.jags",
                   n.chains = 3, n.iter = 6000, n.burnin = 2000, n.thin = 2, parallel = TRUE)
fit.split  <- jags(data.split, NULL, params.split, "borrow_split.jags",
                   n.chains = 3, n.iter = 6000, n.burnin = 2000, n.thin = 2, parallel = TRUE)

print(fit_diagnostics(fit.alone))
if (!fit_ok(fit.alone)) warning("fit.alone: diagnostics need attention before interpreting results.")
print(fit_diagnostics(fit.pooled))
if (!fit_ok(fit.pooled)) warning("fit.pooled: diagnostics need attention before interpreting results.")
print(fit_diagnostics(fit.split))
if (!fit_ok(fit.split)) warning("fit.split: diagnostics need attention before interpreting results.")
cat("Max Rhat  wild alone:", round(fit_rhat(fit.alone), 3),
    "  borrow pooled:", round(fit_rhat(fit.pooled), 3),
    "  borrow split:", round(fit_rhat(fit.split), 3), "\n")

one.week <- rbind(
  wild.alone    = c(estimate = mean(fit.alone$sims.list$wild_ge),
                    lower = unname(quantile(fit.alone$sims.list$wild_ge, 0.025)),
                    upper = unname(quantile(fit.alone$sims.list$wild_ge, 0.975)),
                    abundance.error.pct = abund.error(fit.alone)),
  borrow.pooled = c(estimate = mean(fit.pooled$sims.list$wild_ge),
                    lower = unname(quantile(fit.pooled$sims.list$wild_ge, 0.025)),
                    upper = unname(quantile(fit.pooled$sims.list$wild_ge, 0.975)),
                    abundance.error.pct = abund.error(fit.pooled)),
  borrow.split  = c(estimate = mean(fit.split$sims.list$wild_ge),
                    lower = unname(quantile(fit.split$sims.list$wild_ge, 0.025)),
                    upper = unname(quantile(fit.split$sims.list$wild_ge, 0.975)),
                    abundance.error.pct = abund.error(fit.split)))
cat("Abundance errors use the posterior median of inverse GE, not its mean.\n")
cat("True wild GE:", round(ge.true["wild"], 3), "\n")
print(round(one.week, 3))

# 4. Plot the one week ----
week.df <- data.frame(method = factor(labs3, levels = rev(labs3)), as.data.frame(one.week[, 1:3]))
p1 <- ggplot(week.df, aes(y = method, color = method)) +
  geom_vline(xintercept = ge.true["wild"], linetype = "dashed") +
  geom_linerange(aes(xmin = lower, xmax = upper), linewidth = 1) +
  geom_point(aes(x = estimate), size = 3) +
  scale_color_manual(values = cols3, guide = "none") +
  expand_limits(x = 0) +
  labs(title = "One week (dashed: true wild GE)", x = "Wild guidance efficiency", y = NULL) +
  theme_bw(base_size = 13)
print(p1)

# 5. Repeat the week many times ----
# Only runs if run.repeats is TRUE (set near the top of the script). It's
# FALSE by default, so this whole section is skipped unless you turn it on.
# Errors or inadequate Rhat/ESS diagnostics trigger one longer refit.
# Weeks where any
# of the three models still fails are left out for all three, so the
# comparison stays paired, and the failure count is printed. Results describe
# fits that passed diagnostics; if many fail, the hard weeks are the ones
# missing.

if (run.repeats) {
  n.reps <- 100
  errors <- rhats <- matrix(NA, n.reps, 3,
                            dimnames = list(NULL, c("wild.alone", "borrow.pooled", "borrow.split")))
  passed <- matrix(FALSE, n.reps, 3, dimnames = dimnames(errors))
  zero.guided <- logical(n.reps)
  for (i in 1:n.reps) {
    w <- simulate.week()
    zero.guided[i] <- sum(w$n6[1, c(1, 2, 6)]) == 0
    d.alone  <- list(n = w$n6[1, ], N_seen = sum(w$n6[1, ]))
    d.pooled <- list(n = w$n6, N_seen = rowSums(w$n6))
    d.split  <- list(n = w$n11, N_seen = rowSums(w$n11))
    f.alone  <- safe_jags(d.alone, NULL, params.alone, "wild_alone.jags", n.chains = 3,
                     n.iter = 4000, n.burnin = 1000, n.thin = 3, parallel = FALSE, verbose = FALSE)
    f.pooled <- safe_jags(d.pooled, NULL, params.pooled, "borrow_pooled.jags", n.chains = 3,
                     n.iter = 4000, n.burnin = 1000, n.thin = 3, parallel = FALSE, verbose = FALSE)
    f.split  <- safe_jags(d.split, NULL, params.split, "borrow_split.jags", n.chains = 3,
                     n.iter = 4000, n.burnin = 1000, n.thin = 3, parallel = FALSE, verbose = FALSE)
    # Any fit that errors or fails diagnostics is refit with four times as many
    # iterations (same chains, same thinning, so four times the draws).
    if (!fit_ok(f.alone))
      f.alone <- safe_jags(d.alone, NULL, params.alone, "wild_alone.jags", n.chains = 3,
                      n.iter = 16000, n.burnin = 4000, n.thin = 3, parallel = FALSE, verbose = FALSE)
    if (!fit_ok(f.pooled))
      f.pooled <- safe_jags(d.pooled, NULL, params.pooled, "borrow_pooled.jags", n.chains = 3,
                       n.iter = 16000, n.burnin = 4000, n.thin = 3, parallel = FALSE, verbose = FALSE)
    if (!fit_ok(f.split))
      f.split <- safe_jags(d.split, NULL, params.split, "borrow_split.jags", n.chains = 3,
                      n.iter = 16000, n.burnin = 4000, n.thin = 3, parallel = FALSE, verbose = FALSE)
    passed[i, ] <- c(fit_ok(f.alone), fit_ok(f.pooled), fit_ok(f.split))
    if (all(passed[i, ]))
      errors[i, ] <- c(abund.error(f.alone), abund.error(f.pooled), abund.error(f.split))
    rhats[i, ] <- c(fit_rhat(f.alone),
                    fit_rhat(f.pooled),
                    fit_rhat(f.split))
    if (i %% 10 == 0) cat(i, "")
  }
  cat("\n")
  cat("Fits failing diagnostics after a retry:\n")
  print(colSums(!passed))
  cat("Zero-guided wild weeks:", sum(zero.guided), "of", n.reps,
      "(retained whenever all fits pass; medians remain defined)\n")
  bad.week <- !apply(passed, 1, all)
  cat("Weeks left out for all three models:", sum(bad.week), "of", n.reps,
      "(results below are among weeks where every fit passed)\n")
  errors[bad.week, ] <- NA
  if (all(bad.week)) stop("No paired weeks passed diagnostics; inspect fit failures before summarizing.")

  summ <- t(apply(errors, 2, function(e) round(c(n = sum(!is.na(e)), median = median(e, na.rm = TRUE),
                                                 quantile(e, c(0.05, 0.95), na.rm = TRUE)), 1)))
  colnames(summ) <- c("weeks used", "median % error", "5%", "95%")
  print(summ)

  # Error in wild abundance over the repeats, as violins.
  # Each violin covers the middle 95% of its weeks (2.5% to 97.5%), so rare
  # extreme weeks don't stretch the axis. Thin wild data occasionally give a
  # near-zero GE draw, and 1 / GE blows up; those weeks are in the printed table.
  # The axis stops at 150%; anything past it is clipped and noted in text.
  x.cap <- 150
  err.df <- do.call(rbind, lapply(1:3, function(k) {
    e <- errors[!is.na(errors[, k]), k]
    keep <- e >= quantile(e, 0.025) & e <= quantile(e, 0.975)
    data.frame(method = labs3[k], error = e[keep])
  }))
  err.df$method <- factor(err.df$method, levels = rev(labs3))
  med.df <- data.frame(method = factor(labs3, levels = rev(labs3)), median = apply(errors, 2, median, na.rm = TRUE))
  tail.df <- aggregate(error ~ method, err.df, max)
  tail.df <- tail.df[tail.df$error > x.cap, ]

  p2 <- ggplot(err.df, aes(error, method, fill = method, color = method)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_violin(alpha = 0.35, scale = "width", trim = TRUE) +
    geom_point(data = med.df, aes(median, method), shape = 124, size = 10, inherit.aes = FALSE,
               color = cols3[as.character(med.df$method)]) +
    geom_text(data = tail.df, aes(x = x.cap, y = method, label = paste0("tail to ", round(error), "%")),
              inherit.aes = FALSE, hjust = 1, vjust = -2.2, size = 3.2, color = "gray30") +
    scale_fill_manual(values = cols3, guide = "none") +
    scale_color_manual(values = cols3, guide = "none") +
    coord_cartesian(xlim = c(min(err.df$error, 0), x.cap)) +
    labs(title = paste(n.reps, "simulated weeks (bar: median; dashed: no error)"),
         x = "Error in median expanded wild abundance (%), middle 95% of weeks", y = NULL) +
    theme_bw(base_size = 13)
  print(p2)
}

# Try this ----
# Predict before rerunning; reset each change before trying the next.
# 1. Set a.true <- c(wild = 0.09, hatchery = 0.09). Does borrowing pooled
#    GE still cause the same systematic mismatch?
# 2. Set night_eff <- 0 in the true values. Does a difference in night share
#    still create a difference in true pooled GE?
# 3. Set n.tagged <- c(wild = 20, hatchery = 6000). What happens to wild-only
#    uncertainty and zero-guided weeks? Check diagnostics before comparisons.

# Lessons ----
# - Guidance can differ sharply by time of day. Pooled GE mixes guidance with
#   timing: it moves whenever the night share moves.
# - Wild borrowing hatchery's pooled GE carries hatchery's timing into the
#   wild estimate and biases wild abundance. Borrowing day and night GE,
#   blended with wild's own night share, addresses that timing mismatch.
# - The borrowing models here are complete pooling, the extreme case. The
#   hierarchical model in 06 pools partially, but the same timing problem
#   leaks in whenever pooled GE is what gets shared.
# - Splitting addresses timing bias; it does not guarantee an unbiased
#   abundance estimate. Wild still estimates its night share from few fish.
# - Wild alone can have a long right tail. With no guided fish, its
#   posterior mean of inverse GE is infinite under this prior. The median
#   remains defined; it is a different summary, not a cure for sparse data.
# - Equal GRS detection day and night is an assumption. Relaxing it needs
#   additional information or sharing assumptions; assess identifiability
#   and precision for the particular model.

# Where BPA teaches this ----
# BPA has no day/night guidance example. The closest pieces:
# - 9.2  Multistate model structure: splitting a state adds cells
# - 7.6  Models with time and group effects: what groups share and what they
#        keep separate
