# Hierarchical guidance: borrowing strength across weeks and rear types
# (Run 00_setup.R once first: it installs packages and checks JAGS.)
# (Model fits can take a few minutes to run; R will look busy while JAGS samples.)
#
# What this script should leave you with: a week with 10 wild fish can't
# stand on its own, but it doesn't have to. If wild and hatchery respond to
# the river in similar ways, the hatchery fish tell you most of what a thin
# wild week would have told you, and the wild fish only have to say how they
# differ. The price is that an unusual week gets pulled toward normal.
#
# Wild fish are the target, but most PIT-tagged fish are hatchery. Fitting
# each wild week alone leaves thin weeks with wide, jumpy estimates. A
# hierarchical model shares what wild and hatchery fish have in common and
# keeps what differs:
#
#   logit(psi[r,s]) ~ Normal(mu_strat[g(s)] + rear_psi[r] + beta * spill[s], sigma_psi)
#   mu_strat[g]     ~ Normal(alpha, sigma_strat)
#
#   r = rear type (1 wild, 2 hatchery), s = week, g(s) = the week's stratum
#   rear_psi: fixed offsets, centered to sum to zero (two levels are too few
#             for a random effect; BPA 4.1.4)
#   beta <= 0: spill can only lower guidance
#
# Weeks shrink toward their stratum; strata shrink toward the season average.
# Hatchery fish pin down the strata, the spill slope, and how much weeks
# wiggle, so a thin wild week gets a sensible starting point instead of a
# flat prior. A week with many fish is set by its own fish; a thin week gets
# pulled in. Each rear-week's data come from the six-cell likelihood in
# 05_six_cell.R.
#
# The script fits two models and compares wild weekly GE:
#   Model 1: every rear-week on its own (flat prior on each psi)
#   Model 2: hierarchical
# Both models share detection, transport, and survival across all rear-weeks;
# only the guidance part differs. Teaching version: those shared pieces are
# held constant over the season.
#
# Section 6 repeats the season 20 times so the gain isn't one lucky draw.
# Section 6 is off by default: set run.repeats to TRUE below to run it.
# All data are simulated.
# Requires JAGS and the jagsUI and ggplot2 packages.

library(jagsUI)
library(ggplot2)

set.seed(6)

run.repeats <- FALSE   # TRUE adds section 6, the 20-season repeat (a few minutes)

# 1. Simulate a season ----

R <- 2                                   # 1 = wild, 2 = hatchery
S <- 16                                  # weeks
parent <- rep(1:4, each = 4)             # stratum of each week
n_strat <- 4

alpha <- -2.0
sigma_strat <- 0.5
rear_psi <- c(-0.2, 0.2)                 # centered: sums to zero
beta <- -0.8
sigma_psi <- 0.3
trans <- 0.60; phi_S <- 0.50; delta <- 0.40; p <- 0.40
phi_B <- plogis(qlogis(phi_S) - delta)

simulate.season <- function() {
  mu_strat <- rnorm(n_strat, alpha, sigma_strat)
  spill <- runif(S, 0.2, 0.6)
  spill_std <- (spill - mean(spill)) / sd(spill)
  logit_psi <- matrix(NA, R, S)
  for (r in 1:R) for (s in 1:S) {
    logit_psi[r, s] <- rnorm(1, mu_strat[parent[s]] + rear_psi[r] + beta * spill_std[s], sigma_psi)
  }
  psi <- plogis(logit_psi)
  n.tagged <- rbind(round(exp(runif(S, log(8), log(150)))),   # wild: thin, some very thin
                    round(runif(S, 400, 1500)))              # hatchery: plenty
  # Data layout: n[rear, week, cell], dimensions R x S x 6.
  #   rear: 1 = wild, 2 = hatchery; week: 1 to S; cell: c1 c2 c3 c4 c5 c6
  # From real fish-level records (one row per tag, with rear, week, and the
  # cell from 05_six_cell.R, section 2), the same array is:
  #   table(factor(rear, 1:R), factor(week, 1:S),
  #         factor(cell, c("c1", "c2", "c3", "c4", "c5", "c6")))
  # N_seen[rear, week] is the sum over cells.
  n <- array(0, c(R, S, 6))
  for (r in 1:R) for (s in 1:S) {
    pr <- c(psi[r, s] * (1 - trans) * (1 - phi_B),
            psi[r, s] * (1 - trans) * phi_B,
            (1 - psi[r, s]) * p * (1 - phi_S),
            (1 - psi[r, s]) * p * phi_S,
            (1 - psi[r, s]) * (1 - p) * phi_S,
            psi[r, s] * trans,
            (1 - psi[r, s]) * (1 - p) * (1 - phi_S))
    n[r, s, ] <- rmultinom(1, n.tagged[r, s], pr)[1:6, 1]
  }
  list(n = n, N_seen = apply(n, c(1, 2), sum), psi = psi, spill_std = spill_std, n.tagged = n.tagged)
}

# 2. Model 1: every rear-week on its own ----
# Model in words: each rear-week gets its own GE with a flat prior; nothing
# about guidance is shared. Detection, transport, and GOJ are shared.

cat(file = "ge_separate.jags", "
model {
# Guidance: every rear-week on its own
for (r in 1:R) {
   for (s in 1:S) {
      psi[r,s] ~ dbeta(1, 1)
      }
   }

# Observation model: the same six cells as 05_six_cell.R, for every
# rear-week. Detection, transport, and GOJ (survival and detection) are
# shared by all rear-weeks.
trans ~ dbeta(1, 1)
phi_S ~ dbeta(1, 1)
# Reminder: JAGS's dnorm(mean, precision) uses precision = 1 / SD^2, not SD.
delta ~ dnorm(0, 0.1)
phi_B <- ilogit(logit(phi_S) - delta)
p ~ dbeta(1, 1)
for (r in 1:R) {
   for (s in 1:S) {
      pi_raw[r,s,1] <- psi[r,s] * (1 - trans) * (1 - phi_B)    # c1 GRJ only
      pi_raw[r,s,2] <- psi[r,s] * (1 - trans) * phi_B          # c2 GRJ and GOJ
      pi_raw[r,s,3] <- (1 - psi[r,s]) * p * (1 - phi_S)        # c3 GRS only
      pi_raw[r,s,4] <- (1 - psi[r,s]) * p * phi_S              # c4 GRS and GOJ
      pi_raw[r,s,5] <- (1 - psi[r,s]) * (1 - p) * phi_S        # c5 GOJ only
      pi_raw[r,s,6] <- psi[r,s] * trans                        # c6 transported
      # condition on being detected
      p_seen[r,s] <- sum(pi_raw[r,s,1:6])
      for (cc in 1:6) { pi_obs[r,s,cc] <- pi_raw[r,s,cc] / p_seen[r,s] }
      # each detected fish is one roll of a six-sided die
      n[r,s,1:6] ~ dmulti(pi_obs[r,s,1:6], N_seen[r,s])
      }
   }
}
")

# 3. Model 2: hierarchical ----
# Model in words: each rear-week's GE, on the logit scale, starts from its
# stratum's average, shifts by a fixed wild or hatchery offset and by the
# spill effect, then gets its own random wiggle. Stratum averages are
# themselves drawn around a season average. The detection part is the same
# six cells as Model 1. How far a week strays from its stratum is set by
# sigma_psi, which wild and hatchery share.

cat(file = "ge_hierarchical.jags", "
model {
# Guidance: hierarchical
alpha ~ dnorm(0, 0.1)
sigma_strat ~ dunif(0, 3)
sigma_psi ~ dunif(0, 3)
tau_strat <- pow(sigma_strat, -2)
tau_psi <- pow(sigma_psi, -2)
beta ~ dnorm(0, 0.1) T(, 0)                  # spill can only lower GE

# rear offsets: fixed, centered to sum to zero
for (r in 1:R) { rear_raw[r] ~ dnorm(0, 0.1) }
for (r in 1:R) { rear_psi[r] <- rear_raw[r] - mean(rear_raw[]) }

# strata shrink toward the season average
for (g in 1:n_strat) { mu_strat[g] ~ dnorm(alpha, tau_strat) }

# weeks shrink toward their stratum
for (r in 1:R) {
   for (s in 1:S) {
      logit_psi[r,s] ~ dnorm(mu_strat[parent[s]] + rear_psi[r] + beta * spill_std[s], tau_psi)
      psi[r,s] <- ilogit(logit_psi[r,s])     # convert logit scale to probability scale
      }
   }

# Observation model: the same six cells as 05_six_cell.R, for every
# rear-week. Detection, transport, and GOJ (survival and detection) are
# shared by all rear-weeks.
trans ~ dbeta(1, 1)
phi_S ~ dbeta(1, 1)
delta ~ dnorm(0, 0.1)
phi_B <- ilogit(logit(phi_S) - delta)
p ~ dbeta(1, 1)
for (r in 1:R) {
   for (s in 1:S) {
      pi_raw[r,s,1] <- psi[r,s] * (1 - trans) * (1 - phi_B)    # c1 GRJ only
      pi_raw[r,s,2] <- psi[r,s] * (1 - trans) * phi_B          # c2 GRJ and GOJ
      pi_raw[r,s,3] <- (1 - psi[r,s]) * p * (1 - phi_S)        # c3 GRS only
      pi_raw[r,s,4] <- (1 - psi[r,s]) * p * phi_S              # c4 GRS and GOJ
      pi_raw[r,s,5] <- (1 - psi[r,s]) * (1 - p) * phi_S        # c5 GOJ only
      pi_raw[r,s,6] <- psi[r,s] * trans                        # c6 transported
      # condition on being detected
      p_seen[r,s] <- sum(pi_raw[r,s,1:6])
      for (cc in 1:6) { pi_obs[r,s,cc] <- pi_raw[r,s,cc] / p_seen[r,s] }
      # each detected fish is one roll of a six-sided die
      n[r,s,1:6] ~ dmulti(pi_obs[r,s,1:6], N_seen[r,s])
      }
   }
}
")

# 4. One season ----

season <- simulate.season()
print(season$N_seen)

# Data for each model
jags.data1 <- list(n = season$n, N_seen = season$N_seen, R = R, S = S)
jags.data2 <- list(n = season$n, N_seen = season$N_seen, R = R, S = S,
                   parent = parent, n_strat = n_strat, spill_std = season$spill_std)

# Initial values
inits1 <- function() list(psi = matrix(runif(R * S, 0.05, 0.3), R, S), p = runif(1, 0.2, 0.8),
                          trans = runif(1, 0.2, 0.8), phi_S = runif(1, 0.2, 0.8))
inits2 <- function() list(alpha = rnorm(1, -2, 0.5), beta = -runif(1, 0, 1), p = runif(1, 0.2, 0.8),
                          trans = runif(1, 0.2, 0.8), phi_S = runif(1, 0.2, 0.8),
                          sigma_strat = runif(1, 0.1, 1), sigma_psi = runif(1, 0.1, 1))

# Parameters to save. The shared and variance parameters are saved too, so
# convergence is checked on the whole model, not just on psi.
# Save rear_psi, not rear_raw. Centering cancels the common mean of rear_raw,
# so the data say nothing about it: its posterior is just its broad prior.
# The sampler can move slowly along that prior-only direction, so rear_raw's
# Rhat can look poor even when every quantity the model uses has converged.
params1 <- c("psi", "p", "trans", "phi_S", "delta")
params2 <- c("psi", "beta", "alpha", "mu_strat", "rear_psi", "sigma_psi", "sigma_strat",
             "p", "trans", "phi_S", "delta")

# Fitting: this can take a few minutes.
fit1 <- jags(jags.data1, inits1, params1, "ge_separate.jags",
             n.chains = 3, n.iter = 10000, n.burnin = 5000, n.thin = 5, parallel = TRUE)
fit2 <- jags(jags.data2, inits2, params2, "ge_hierarchical.jags",
             n.chains = 3, n.iter = 10000, n.burnin = 5000, n.thin = 5, parallel = TRUE)

cat("Max Rhat, separate:    ", round(max(unlist(fit1$Rhat), na.rm = TRUE), 3), "\n")
cat("Max Rhat, hierarchical:", round(max(unlist(fit2$Rhat), na.rm = TRUE), 3), "\n")

# Wild weekly GE: posterior mean and 95% CrI for each week
wild1 <- fit1$sims.list$psi[, 1, ]     # draws x weeks, wild only
wild2 <- fit2$sims.list$psi[, 1, ]
w1 <- data.frame(mean = colMeans(wild1), lo = apply(wild1, 2, quantile, 0.025), hi = apply(wild1, 2, quantile, 0.975))
w2 <- data.frame(mean = colMeans(wild2), lo = apply(wild2, 2, quantile, 0.025), hi = apply(wild2, 2, quantile, 0.975))
truth <- season$psi[1, ]

cat("\nWild weekly GE, mean absolute error:  separate", round(mean(abs(w1$mean - truth)), 3),
    "  hierarchical", round(mean(abs(w2$mean - truth)), 3), "\n")
cat("Wild weekly GE, mean 95% CrI width:   separate", round(mean(w1$hi - w1$lo), 3),
    "  hierarchical", round(mean(w2$hi - w2$lo), 3), "\n")
cat("Spill slope beta: truth", beta, " estimate", round(fit2$mean$beta, 2), "\n")

# 5. Plot the one season ----
wk.df <- rbind(data.frame(week = 1:S, model = "Each week on its own", w1),
               data.frame(week = 1:S, model = "Hierarchical", w2))
dodge <- position_dodge(width = 0.5)
p1 <- ggplot(wk.df, aes(factor(week), mean, color = model)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 1, position = dodge) +
  geom_point(size = 2.5, position = dodge) +
  geom_point(data = data.frame(week = 1:S, truth = truth), aes(factor(week), truth),
             inherit.aes = FALSE, shape = 4, size = 3, stroke = 1.2) +
  annotate("text", x = 1:S, y = max(wk.df$hi, truth) * 1.06, label = season$n.tagged[1, ], size = 3, color = "gray40") +
  scale_color_manual(NULL, values = c("Each week on its own" = "darkorange", "Hierarchical" = "steelblue")) +
  labs(title = "One season: wild weekly GE (x = truth; numbers = wild fish tagged)",
       x = "Week", y = "Wild guidance efficiency") +
  theme_bw(base_size = 13) + theme(legend.position = "top")
print(p1)

# 6. Repeat the season ----
# Only runs if run.repeats is TRUE (set near the top of the script). It's
# FALSE by default, so this whole section is skipped unless you turn it on.
# Pool every wild week from every season and ask how each model does as a
# function of how many wild fish were tagged that week. A fit with any Rhat
# above 1.1 is refit with four times as many iterations. Seasons where either model still
# fails are left out of BOTH models, so the comparison stays paired, and the
# failure rate is printed. Results describe fits that passed diagnostics; if
# many fail, the hard data sets are the ones missing, and the summary looks
# better than the truth.

if (run.repeats) {
  n.reps <- 20
  rhats <- matrix(NA, 0, 2, dimnames = list(NULL, c("separate", "hierarchical")))
  all.weeks <- data.frame()
  for (i in 1:n.reps) {
    ssn <- simulate.season()
    d1 <- list(n = ssn$n, N_seen = ssn$N_seen, R = R, S = S)
    d2 <- list(n = ssn$n, N_seen = ssn$N_seen, R = R, S = S,
               parent = parent, n_strat = n_strat, spill_std = ssn$spill_std)
    f1 <- jags(d1, inits1, params1, "ge_separate.jags", n.chains = 3,
               n.iter = 6000, n.burnin = 3000, n.thin = 3, parallel = FALSE, verbose = FALSE)
    f2 <- jags(d2, inits2, params2, "ge_hierarchical.jags", n.chains = 3,
               n.iter = 6000, n.burnin = 3000, n.thin = 3, parallel = FALSE, verbose = FALSE)
    # If either fit has any Rhat above 1.1, refit it once with four times as
    # many iterations (same chains, same thinning, so four times the draws).
    if (max(unlist(f1$Rhat), na.rm = TRUE) > 1.1)
      f1 <- jags(d1, inits1, params1, "ge_separate.jags", n.chains = 3,
                 n.iter = 24000, n.burnin = 12000, n.thin = 3, parallel = FALSE, verbose = FALSE)
    if (max(unlist(f2$Rhat), na.rm = TRUE) > 1.1)
      f2 <- jags(d2, inits2, params2, "ge_hierarchical.jags", n.chains = 3,
                 n.iter = 24000, n.burnin = 12000, n.thin = 3, parallel = FALSE, verbose = FALSE)
    r1 <- max(unlist(f1$Rhat), na.rm = TRUE); r2 <- max(unlist(f2$Rhat), na.rm = TRUE)
    rhats <- rbind(rhats, c(r1, r2))
    if (max(r1, r2) > 1.1) { cat("x", ""); next }
    for (m in 1:2) {
      d <- if (m == 1) f1$sims.list$psi[, 1, ] else f2$sims.list$psi[, 1, ]
      w <- data.frame(mean = colMeans(d), lo = apply(d, 2, quantile, 0.025), hi = apply(d, 2, quantile, 0.975))
      all.weeks <- rbind(all.weeks, data.frame(
        season = i, model = c("Each week on its own", "Hierarchical")[m],
        n.tagged = ssn$n.tagged[1, ], truth = ssn$psi[1, ], w))
    }
    cat(i, "")
  }
  cat("\n")
  cat("Seasons failing diagnostics after a retry:", sum(apply(rhats, 1, max) > 1.1), "of", n.reps,
      "(left out of both models; results below are among fits passing diagnostics)\n")

  all.weeks$abs.error <- abs(all.weeks$mean - all.weeks$truth)
  all.weeks$width <- all.weeks$hi - all.weeks$lo
  all.weeks$covered <- all.weeks$lo <= all.weeks$truth & all.weeks$truth <= all.weeks$hi
  all.weeks$size.class <- cut(all.weeks$n.tagged, c(0, 20, 50, Inf),
                              labels = c("8 to 20 tagged", "21 to 50", "51 to 150"))

  tab <- aggregate(cbind(abs.error, width, covered) ~ model + size.class, all.weeks, mean)
  tab$covered <- round(100 * tab$covered, 1)
  tab$abs.error <- round(tab$abs.error, 3); tab$width <- round(tab$width, 3)
  names(tab) <- c("model", "wild fish tagged", "mean abs error", "mean CrI width", "coverage %")
  print(tab[order(tab[[2]], tab$model), ], row.names = FALSE)

  p2 <- ggplot(all.weeks, aes(n.tagged, abs.error, color = model)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 1.2, formula = y ~ x) +
    scale_x_log10() +
    scale_color_manual(NULL, values = c("Each week on its own" = "darkorange", "Hierarchical" = "steelblue")) +
    labs(title = paste(n.reps, "simulated seasons: error in wild weekly GE"),
         x = "Wild fish tagged that week (log scale)", y = "Absolute error in GE") +
    theme_bw(base_size = 13) + theme(legend.position = "top")
  print(p2)
}

# Lessons ----
# - Thin groups borrow strength from the rest: a thin wild week gets pulled
#   toward its stratum and the spill effect, which hatchery fish pin down.
# - The gain is biggest where data are thinnest. Well-sampled weeks are set
#   by their own fish, and the two models agree there.
# - Shrinkage has a cost: a week that really is unusual gets pulled toward
#   its stratum too. Check coverage to see whether the trade is honest.
# - Two rear types are too few for a random effect, so rear offsets are fixed
#   and centered (BPA 4.1.4).

# Where BPA teaches this ----
# - 4.1  Introduction to random effects, including when a factor has too few
#        levels to treat as random (4.1.4 in the book)
# - 4.3  Mixed models with random effects for variability among groups
# - 7.4.2  Random time effects (weeks here)
# - 7.6.2  Fixed group and random time effects (rear type fixed, weeks random)
