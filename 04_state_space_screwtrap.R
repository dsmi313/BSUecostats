# State-space model for daily screw trap abundance
# (Run 00_setup.R once first: it installs packages and checks JAGS.)
# (Model fits can take a few minutes to run; R will look busy while JAGS samples.)
#
# Builds on Kery and Schaub (2012) Bayesian Population Analysis, Chapter 5.
# The process model is the log-scale model from BPA 5.4 (house martins).
# BPA 5.3 and 5.5 make the key point: counts alone estimate an index (N * p),
# not abundance. To get N you need separate information on detection (p).
# At a screw trap, that information is the efficiency trials: marked fish
# released upstream and recaptured. This script shows both halves:
#   Model 1: counts only (adapted from BPA Chapter 5). It tracks N * p.
#   Model 2: counts plus efficiency trials. It recovers N.
#
# This is a teaching example, not the model IDFG uses to estimate screw trap
# abundance. The concepts carry over: a process model for true abundance, an
# observation model for what the trap sees, and separate data on detection.
# IDFG's production model (J. McCormick, unpublished 2023; Smith et al. 2024,
# Appendix I) has the same structure, with two differences:
# - It writes efficiency as logit(p[t]) <- logit(p[t-1] + g[t]). The two
#   logits cancel, so p takes random-walk steps on the probability scale.
#   Model 2 here steps on the logit scale, which keeps p between 0 and 1.
# - It credits each recapture to the day the fish was marked:
#   R[t] ~ dbin(p[t], M[t]), with M[t] the fish marked on day t. Here,
#   trials are separate releases every other day.
#
# All data are simulated, so we know the truth and can check the answer.
# Requires JAGS (https://sourceforge.net/projects/mcmc-jags/) and the jagsUI package.

library(jagsUI)
source("helpers.R")
library(ggplot2)

set.seed(2026)

# 1. Simulate a season of daily passage ----

n.days <- 80

# State process (BPA 5.4): log(N[t+1]) = log(N[t]) + r[t], r[t] ~ Normal(mean.r, sigma.proc)
logN1 <- log(150)        # Fish passing on day 1
mean.r <- 0.01           # Average daily change on the log scale
sigma.proc <- 0.25       # Day-to-day process variation

logN <- numeric(n.days)
logN[1] <- logN1
r <- rnorm(n.days - 1, mean.r, sigma.proc)
for (t in 1:(n.days - 1)) {
  logN[t + 1] <- logN[t] + r[t]
}
N <- round(exp(logN))    # True daily passage (integers)

# Trap efficiency p wanders from day to day on the logit scale.
# Same random walk idea as the state process, one level down.
sigma.p <- 0.15
lp <- numeric(n.days)
lp[1] <- qlogis(0.20)
for (t in 1:(n.days - 1)) {
  lp[t + 1] <- lp[t] + rnorm(1, 0, sigma.p)
}
p <- plogis(lp)

# Observation process, part 1: daily catch. Each passing fish is caught with probability p.
C <- rbinom(n.days, N, p)

# The trap was pulled for debris on days 41 to 45. Those days become NA.
# As in BPA 5.4, JAGS treats NA data as unknowns and estimates N for them anyway.
trap.out <- 41:45
C[trap.out] <- NA

# Observation process, part 2: efficiency trials every other day.
# M fish are marked, released upstream, and R of them are recaptured.
trial.days <- seq(2, n.days, by = 2)
trial.days <- setdiff(trial.days, trap.out)
M <- rep(30, length(trial.days))
R <- rbinom(length(trial.days), M, p[trial.days])

# 2. Model 1: counts only (adapted from BPA Chapter 5) ----
# Poisson counts around the latent abundance, with no detection term.
# This is exactly the situation BPA 5.3 warns about: with p < 1,
# the estimated "N" is really N * p.

cat(file = "ssm_counts_only.jags", "
model {
# Priors
logN.est[1] ~ dunif(0, 12)           # Initial log abundance
# Reminder: JAGS's dnorm(mean, precision) uses precision = 1 / SD^2, not SD.
mean.r ~ dnorm(0, 0.01)              # Mean daily growth rate on the log scale
sigma.proc ~ dunif(0, 2)             # SD of the state process
tau.proc <- pow(sigma.proc, -2)

# State process
for (t in 1:(T - 1)) {
   r[t] ~ dnorm(mean.r, tau.proc)
   logN.est[t + 1] <- logN.est[t] + r[t]
   }

# Observation process
for (t in 1:T) {
   N.est[t] <- exp(logN.est[t])
   C[t] ~ dpois(N.est[t])
   }

N.total <- sum(N.est[])
}
")

jags.data1 <- list(C = C, T = n.days)
inits1 <- function() list(logN.est = c(log(mean(C, na.rm = TRUE)), rep(NA, n.days - 1)),
                          mean.r = 0, sigma.proc = runif(1, 0.1, 0.5))
params1 <- c("N.est", "N.total", "mean.r", "sigma.proc")

ni <- 60000; nb <- 30000; nt <- 30; nc <- 3

# Fitting: this can take a few minutes.
fit1 <- jags(jags.data1, inits1, params1, "ssm_counts_only.jags",
             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)

# 3. Model 2: counts plus efficiency trials ----
# The trials give the model direct information on p, so N and p separate.
# Catch is Poisson(p * N), an approximation to binomial thinning
# when N is a continuous latent state.
#
# Model in words: each day's true passage is yesterday's times a random daily
# change (on the log scale). Each day's trap efficiency is yesterday's plus a
# small random step (on the logit scale). The catch is the passage times the
# efficiency, plus Poisson noise. The trials are binomial: of M marked fish,
# each is recaught with that day's efficiency.

cat(file = "ssm_with_trials.jags", "
model {
# Priors
logN.est[1] ~ dunif(0, 12)           # Initial log abundance
mean.r ~ dnorm(0, 0.01)              # Mean daily growth rate on the log scale
sigma.proc ~ dunif(0, 2)             # SD of the state process
tau.proc <- pow(sigma.proc, -2)

lp[1] ~ dlogis(0, 1)                 # Same as p[1] ~ Uniform(0, 1)
sigma.p ~ dunif(0, 1)                # SD of the daily change in logit p
tau.p <- pow(sigma.p, -2)

# State process: abundance
for (t in 1:(T - 1)) {
   r[t] ~ dnorm(mean.r, tau.proc)
   logN.est[t + 1] <- logN.est[t] + r[t]
   }

# Trap efficiency: random walk on the logit scale
for (t in 1:(T - 1)) {
   lp[t + 1] ~ dnorm(lp[t], tau.p)
   }

# Observation process, part 1: daily catch
for (t in 1:T) {
   p[t] <- ilogit(lp[t])
   N.est[t] <- exp(logN.est[t])
   C[t] ~ dpois(p[t] * N.est[t])
   }

# Observation process, part 2: efficiency trials
for (k in 1:K) {
   R[k] ~ dbin(p[trial.day[k]], M[k])
   }

N.total <- sum(N.est[])
}
")

jags.data2 <- list(C = C, T = n.days, R = R, M = M, K = length(trial.days), trial.day = trial.days)
inits2 <- function() list(logN.est = c(log(mean(C, na.rm = TRUE) / 0.2), rep(NA, n.days - 1)),
                          mean.r = 0, sigma.proc = runif(1, 0.1, 0.5),
                          sigma.p = runif(1, 0.05, 0.3))
params2 <- c("N.est", "N.total", "p", "mean.r", "sigma.proc", "sigma.p")

fit2 <- jags(jags.data2, inits2, params2, "ssm_with_trials.jags",
             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)

# 4. Check convergence ----
# Rhat near 1 (below 1.01) for every monitored node.
print(fit_diagnostics(fit1))
if (!fit_ok(fit1)) warning("fit1: diagnostics need attention before interpreting results.")
print(fit_diagnostics(fit2))
if (!fit_ok(fit2)) warning("fit2: diagnostics need attention before interpreting results.")
cat("Model 1 max Rhat:", round(fit_rhat(fit1), 3), "\n")
cat("Model 2 max Rhat:", round(fit_rhat(fit2), 3), "\n")

# 5. Compare to the truth ----
ci <- function(x) round(quantile(x, c(0.025, 0.975)))
cat("\nTrue season total:           ", sum(N), "\n")
cat("True total of N * p:         ", round(sum(N * p)), "\n")
cat("Model 1 (counts only):       ", round(fit1$mean$N.total),
    " 95% CrI", ci(fit1$sims.list$N.total), "\n")
cat("Model 2 (counts plus trials):", round(fit2$mean$N.total),
    " 95% CrI", ci(fit2$sims.list$N.total), "\n")

# 6. Plot ----
# Each plot prints on its own.
band.df <- function(fit, label) {
  data.frame(day = 1:n.days, model = label, mean = fit$mean$N.est,
             lo = apply(fit$sims.list$N.est, 2, quantile, 0.025),
             hi = apply(fit$sims.list$N.est, 2, quantile, 0.975))
}
bands <- rbind(band.df(fit1, "Counts only"), band.df(fit2, "Counts plus trials"))
truth.df <- rbind(data.frame(day = 1:n.days, value = N, series = "True N"),
                  data.frame(day = 1:n.days, value = N * p, series = "True N * p"))

p1 <- ggplot() +
  annotate("rect", xmin = min(trap.out) - 0.5, xmax = max(trap.out) + 0.5, ymin = -Inf, ymax = Inf,
           fill = "gray90") +
  geom_ribbon(data = bands, aes(day, ymin = lo, ymax = hi, fill = model), alpha = 0.25) +
  geom_line(data = bands, aes(day, mean, color = model), linewidth = 1) +
  geom_line(data = truth.df, aes(day, value, linetype = series)) +
  geom_point(data = data.frame(day = 1:n.days, C = C), aes(day, C), size = 1, color = "gray40", na.rm = TRUE) +
  scale_color_manual(NULL, values = c("Counts only" = "darkorange", "Counts plus trials" = "steelblue")) +
  scale_fill_manual(NULL, values = c("Counts only" = "darkorange", "Counts plus trials" = "steelblue")) +
  scale_linetype_manual(NULL, values = c("True N" = "solid", "True N * p" = "dashed")) +
  labs(title = "Daily abundance (gray band: trap out; points: catch)", x = "Day of season", y = "Fish passing") +
  theme_bw(base_size = 13) + theme(legend.position = "top")
print(p1)

p.df <- data.frame(day = 1:n.days, truth = p, mean = fit2$mean$p,
                   lo = apply(fit2$sims.list$p, 2, quantile, 0.025),
                   hi = apply(fit2$sims.list$p, 2, quantile, 0.975))
p2 <- ggplot(p.df, aes(day)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = mean, color = "Estimated p"), linewidth = 1) +
  geom_line(aes(y = truth, color = "True p"), linewidth = 1) +
  geom_point(data = data.frame(day = trial.days, rate = R / M), aes(day, rate, color = "Trial R / M")) +
  scale_color_manual(NULL, values = c("Estimated p" = "steelblue", "True p" = "black", "Trial R / M" = "gray40")) +
  guides(color = guide_legend(override.aes = list(linetype = c(1, 0, 1), shape = c(NA, 16, NA)))) +
  labs(title = "Trap efficiency", x = "Day of season", y = "Capture probability") +
  theme_bw(base_size = 13) + theme(legend.position = "top")
print(p2)

# Try this ----
# Predict before rerunning; reset each change before trying the next.
# 1. Change trap.out <- 41:45 to trap.out <- 31:55. How do the passage
#    intervals inside the outage compare with the observed days?
# 2. Change M <- rep(30, length(trial.days)) to use 10 marked fish per trial.
#    Which quantities become less precise, and do both fits pass diagnostics?

# Lessons ----
# - A state-space model separates the process (how many fish really pass each
#   day) from the observation (how many the trap catches).
# - Counts alone estimate N * p, an index, not N (BPA 5.3 and 5.5). If p
#   drifts over the season, the index drifts with it and looks like a change
#   in abundance.
# - Adding separate data on detection (efficiency trials) is what turns the
#   index into abundance.
# - Days with the trap out can be left as NA. The model fills them in from
#   the days around them, with uncertainty under the assumed model.
#   The binomial simulation and Poisson observation model differ: Poisson
#   variance is N*p rather than N*p*(1-p). Check coverage before adaptation.

# Where BPA teaches this ----
# - 5.2  A simple state-space model
# - 5.3  Systematic bias in the observation process: counts track N * p
# - 5.4  House martin example: the log-scale model used here, and NA data
#        estimated as unknowns
# - 5.5  Summary and outlook: counts alone give an index; getting N needs
#        separate information on detection

# The six-cell multistate model at Lower Granite Dam
# (Run 00_setup.R once first: it installs packages and checks JAGS.)
# (Model fits can take a few minutes to run; R will look busy while JAGS samples.)
#
# What this script should leave you with: you never see a fish's route, only
# where it was detected, and each detection site exists to identify one
# piece of the puzzle. Take a site away and the model still runs, but it
# can't tell two answers apart. Feed it the wrong fish and it confidently
# gives the wrong answer.
#
# Each PIT-tagged fish reaches Lower Granite and takes one route:
#   guided (psi) into the juvenile bypass, read at GRJ, or
#   not guided (1 - psi), passing the spillway or turbines.
# Guided fish are barged (trans) or returned to the river. Nonguided fish are
# read at the spillway array GRS with probability p. Some are read again
# downstream at Little Goose (GOJ).
#   phi_B, phi_S: the chance a bypassed (B) or spilled (S) fish survives to GOJ
#   AND is detected there. Survival and GOJ detection can't be separated with
#   these data, so phi is their product, not survival alone.
#
# The route is the hidden state; detections are what we see (BPA Chapter 9).
# Six detectable cells plus one never-detected leaf:
#   c1 GRJ only        psi (1 - trans)(1 - phi_B)
#   c2 GRJ and GOJ     psi (1 - trans) phi_B
#   c3 GRS only        (1 - psi) p (1 - phi_S)
#   c4 GRS and GOJ     (1 - psi) p phi_S
#   c5 GOJ only        (1 - psi)(1 - p) phi_S
#   c6 transported     psi trans
#   never detected     (1 - psi)(1 - p)(1 - phi_S)
#
# Never-detected fish can't be counted, so the likelihood conditions on
# being detected: divide each cell by the chance of being seen at all.
#
# Sections 7 and 8 show two ways this goes wrong:
#   7. Fish tagged after the decision: marked fish that aren't a fair sample
#   8. No downstream site: without GOJ, p and psi can't be separated
#
# A teaching version, simplified from smoltEASE: one group, no covariates.
# All data are simulated. Requires JAGS and the jagsUI package.

library(jagsUI)
source("helpers.R")
library(ggplot2)

# 0. The row method: read efficiency within spillbay 1 ----
# Toy row-read probabilities from "Step 2: The Row Method". A fish in the
# monitored route crosses all three rows. Independent, equally readable
# fish give d = 1 - product(1 - row probabilities).
row.p <- c(0.64, 0.61, 0.53)
row.d <- 1 - prod(1 - row.p)
cat("Talk example, array read probability within bay 1:", row.d, "\n")
# Observe only seven histories: 111, 110, 101, 011, 100, 010, 001.
# The unseen 000 history is absent from the field data.
set.seed(12)
row.reads <- sapply(row.p, function(prob) rbinom(10000, 1, prob))
row.seen <- row.reads[rowSums(row.reads) > 0, , drop = FALSE]
print(table(apply(row.seen, 1, paste0, collapse = "")))
# Row 2's reads supply the recapture sample for estimating row 1's p.
# Other overlaps estimate rows 2 and 3. Conditional independence makes
# these fractions valid despite the missing 000 fish.
row.p.hat <- c(mean(row.seen[row.seen[, 2] == 1, 1]),
               mean(row.seen[row.seen[, 1] == 1, 2]),
               mean(row.seen[row.seen[, 2] == 1, 3]))
row.d.hat <- 1 - prod(1 - row.p.hat)
cat("Simulated row estimates:", row.p.hat, "  array d:", row.d.hat, "\n")
# The rows estimate d, not the share q of nonguided fish using bay 1.
# Whole-route detection p = q*d. Shared unreadability breaks independence
# and can make the row method overstate d. GOJ supplies additional information.

set.seed(4)

# 1. Simulate one record per tagged fish ----
# Real data arrive one row per PIT tag, with a flag for each detection site.
# Simulate data in that same shape so section 2 can show the formatting step.

n.tagged <- 6000     # Tagged fish reaching Lower Granite (includes never-detected)
psi <- 0.10          # Guidance efficiency, the target
trans <- 0.60        # Share of guided fish that are barged
phi_S <- 0.50        # Survive to GOJ and be detected there, after spill or turbine passage
delta <- 0.40        # Bypassed fish are less likely to be seen at GOJ, on the logit scale
phi_B <- plogis(qlogis(phi_S) - delta)
p <- 0.40            # GRS detection of nonguided fish

is.guided <- rbinom(n.tagged, 1, psi)
barged <- ifelse(is.guided == 1, rbinom(n.tagged, 1, trans), 0)
fish <- data.frame(
  tag = sprintf("3D9.%06d", 1:n.tagged),
  GRJ = is.guided,                                                   # read in the bypass
  barged = barged,                                                # loaded on a barge
  GRS = ifelse(is.guided == 0, rbinom(n.tagged, 1, p), 0),           # read at the spillway
  GOJ = ifelse(is.guided == 1,
               ifelse(barged == 1, 0, rbinom(n.tagged, 1, phi_B)),
               rbinom(n.tagged, 1, phi_S)))                       # read at Little Goose
# Real data never contain fish detected nowhere; drop them as PTAGIS would
fish <- fish[fish$GRJ + fish$GRS + fish$GOJ > 0, ]
head(fish)

# 2. Format: detection histories to cell counts ----
# This is the step where most real analyses trip up. Each fish goes to
# exactly one cell, and the model wants a vector of 6 counts, n[cell], in
# this order: c1 c2 c3 c4 c5 c6. With groups and weeks it becomes an array
# n[group, week, cell] (see 06_hierarchical.R); build it with table() on the
# same cell column plus group and week columns.

cell.of <- function(GRJ, GRS, GOJ, barged) {
  flags <- list(GRJ, GRS, GOJ, barged)
  if (length(unique(lengths(flags))) != 1L ||
      anyNA(unlist(flags)) || !all(unlist(flags) %in% c(0, 1)))
    stop("Detection and transport flags must be equally long, nonmissing 0/1 vectors.")
  if (any(GRJ == 1 & GRS == 1) || any(barged == 1 & (GRJ != 1 | GRS == 1 | GOJ == 1)))
    stop("History conflicts with the assumed routes. Resolve passage episodes before tabulating.")
  ifelse(barged == 1,             "c6",     # transported
  ifelse(GRJ == 1 & GOJ == 1,     "c2",     # GRJ and GOJ
  ifelse(GRJ == 1,                "c1",     # GRJ only
  ifelse(GRS == 1 & GOJ == 1,     "c4",     # GRS and GOJ
  ifelse(GRS == 1,                "c3",     # GRS only
  ifelse(GOJ == 1,                "c5",     # GOJ only
                                  NA))))))
}
stopifnot(!anyNA(fish$tag), !anyDuplicated(fish$tag))
fish$cell <- cell.of(fish$GRJ, fish$GRS, fish$GOJ, fish$barged)
stopifnot(!any(is.na(fish$cell)))            # every detected fish lands in a cell

n <- table(factor(fish$cell, levels = c("c1", "c2", "c3", "c4", "c5", "c6")))
n <- setNames(as.vector(n), names(n))        # named vector of 6 counts
N_seen <- sum(n)
print(n)

# True cell probabilities, for reference
cell.prob <- c(c1 = psi * (1 - trans) * (1 - phi_B),
               c2 = psi * (1 - trans) * phi_B,
               c3 = (1 - psi) * p * (1 - phi_S),
               c4 = (1 - psi) * p * phi_S,
               c5 = (1 - psi) * (1 - p) * phi_S,
               c6 = psi * trans,
               never = (1 - psi) * (1 - p) * (1 - phi_S))
stopifnot(abs(sum(cell.prob) - 1) < 1e-12)

# 3. Back-of-envelope estimates ----
# c4 and c5 both skipped the bypass and both were detected at GOJ. They differ
# only in whether GRS read them, so phi_S (survival and GOJ detection) cancels
# (Lincoln-Petersen logic). It assumes GRS-read and GRS-missed fish reach and
# are detected at GOJ at the same rate.
p.hat <- n["c4"] / (n["c4"] + n["c5"])
guided <- n["c1"] + n["c2"] + n["c6"]
nonguided.read <- n["c3"] + n["c4"]
ge.hat <- guided / (guided + nonguided.read / p.hat)
cat("Closed form: p-hat =", round(p.hat, 3), " GE =", round(ge.hat, 3), "\n")

# 3b. The talk's MY2025 wild counts: GE before the multinomial ----
# These are the season counts shown in the talk, not the simulated fish above.
# Keep slide.* names so the JAGS example below still fits its simulated data.
slide.counts <- c(c1 = 320, c2 = 35, c3 = 3662, c4 = 123, c5 = 195, c6 = 156)
slide.guided <- sum(slide.counts[c("c1", "c2", "c6")])   # 511 GRJ / transported
slide.grs <- sum(slide.counts[c("c3", "c4")])             # 3,785 GRS
# Step 1: count only fish seen at GRJ or GRS.
slide.ge.naive <- slide.guided / (slide.guided + slide.grs)
# Step 2: correct read efficiency in spillbay 1 using the three-row example.
# This still misses nonguided fish using other routes.
slide.row.d <- 1 - prod(1 - c(0.64, 0.61, 0.53))
slide.ge.rows <- slide.guided / (slide.guided + slide.grs / slide.row.d)
# Step 3: GOJ-only fish reveal nonguided passage missed by GRS altogether.
slide.p <- unname(slide.counts["c4"] / sum(slide.counts[c("c4", "c5")]))
slide.ge.goj <- slide.guided / (slide.guided + slide.grs / slide.p)
print(data.frame(step = c("GRJ and GRS only", "Correct bay-1 row reads", "Add GOJ recovery"),
                 GE = c(slide.ge.naive, slide.ge.rows, slide.ge.goj)))
cat("GRS detection from GOJ overlaps:", slide.p, "\n")
# Step 4: distinguish transport from return to river and route recovery.
cat("Transport fraction:", slide.counts["c6"] / slide.guided,
    "  bypass recovery:", slide.counts["c2"] / sum(slide.counts[c("c1", "c2")]),
    "  nonguided recovery:", slide.counts["c4"] / slide.grs, "\n")
stopifnot(slide.guided == 511, slide.grs == 3785,
          isTRUE(all.equal(round(c(slide.ge.naive, slide.ge.rows, slide.ge.goj), 2),
                           c(0.12, 0.11, 0.05))))
# The overlap calculation assumes GRS-read and GRS-missed nonguided fish
# have the same downstream recovery probability. These are season-level
# closed-form calculations, not the hierarchical day/night model estimates.

# 4. JAGS model ----
# Model in words: five probabilities (psi, p, trans, phi_S, and the bypass
# penalty delta) multiply along the branches of the detection tree to give
# each cell's chance. Never-detected fish can't be counted, so each cell's
# chance is divided by the chance of being detected at all. Then each
# detected fish is one roll of a six-sided die with those chances.

cat(file = "six_cell.jags", "
model {
# Priors: flat on the probability scale
psi ~ dbeta(1, 1)
p ~ dbeta(1, 1)
trans ~ dbeta(1, 1)
phi_S ~ dbeta(1, 1)
# Reminder: JAGS's dnorm(mean, precision) uses precision = 1 / SD^2, not SD.
delta ~ dnorm(0, 0.1)
phi_B <- ilogit(logit(phi_S) - delta)

# Cell probabilities before conditioning
pi_raw[1] <- psi * (1 - trans) * (1 - phi_B)    # c1 GRJ only
pi_raw[2] <- psi * (1 - trans) * phi_B          # c2 GRJ and GOJ
pi_raw[3] <- (1 - psi) * p * (1 - phi_S)        # c3 GRS only
pi_raw[4] <- (1 - psi) * p * phi_S              # c4 GRS and GOJ
pi_raw[5] <- (1 - psi) * (1 - p) * phi_S        # c5 GOJ only
pi_raw[6] <- psi * trans                        # c6 transported

# Condition on being detected: only fish we can see
p_seen <- sum(pi_raw[1:6])
for (cc in 1:6) {
    pi_obs[cc] <- pi_raw[cc] / p_seen
    }

# Each detected fish is one roll of a six-sided die
n[1:6] ~ dmulti(pi_obs[1:6], N_seen)

# Expansion holding N_seen fixed: detection-probability uncertainty only.
# This is not a full posterior for the number of tagged arrivals. That would
# require an abundance model and a likelihood for the total detected count.
N_total <- N_seen / p_seen
}
")

jags.data <- list(n = as.vector(n), N_seen = N_seen)
inits <- function() list(psi = runif(1, 0.05, 0.5), p = runif(1, 0.2, 0.8),
                         trans = runif(1, 0.2, 0.8), phi_S = runif(1, 0.2, 0.8), delta = rnorm(1, 0, 0.5))
params <- c("psi", "p", "trans", "phi_S", "phi_B", "delta", "N_total")

# Fitting: this can take a few minutes.
fit <- jags(jags.data, inits, params, "six_cell.jags",
            n.chains = 3, n.iter = 6000, n.burnin = 2000, n.thin = 2, parallel = TRUE)

# 5. Compare to the truth ----
print(fit_diagnostics(fit))
if (!fit_ok(fit)) warning("fit: diagnostics need attention before interpreting results.")
cat("Max Rhat:", round(fit_rhat(fit), 3), "\n")
truth <- c(psi = psi, p = p, trans = trans, phi_S = phi_S, phi_B = phi_B)
est <- sapply(names(truth), function(k) fit$mean[[k]])
lo <- sapply(names(truth), function(k) quantile(fit$sims.list[[k]], 0.025))
hi <- sapply(names(truth), function(k) quantile(fit$sims.list[[k]], 0.975))
print(round(cbind(truth, estimate = est, lower = lo, upper = hi), 3))
cat("Expanded tagged passage (detection-probability uncertainty only): truth", n.tagged, " estimate", round(fit$mean$N_total),
    " 95% expansion interval", round(quantile(fit$sims.list$N_total, c(0.025, 0.975))), "\n")

# 6. Plot ----
rec.df <- data.frame(par = factor(names(truth), levels = rev(names(truth))),
                     truth = truth, est = est, lo = lo, hi = hi)
p1 <- ggplot(rec.df, aes(y = par)) +
  geom_linerange(aes(xmin = lo, xmax = hi), color = "steelblue", linewidth = 1) +
  geom_point(aes(x = est, shape = "Posterior mean, 95% CrI"), color = "steelblue", size = 2.5) +
  geom_point(aes(x = truth, shape = "Truth"), size = 3, stroke = 1.2) +
  scale_shape_manual(NULL, values = c("Posterior mean, 95% CrI" = 16, "Truth" = 4)) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "Probability", y = NULL) +
  theme_bw(base_size = 13) + theme(legend.position = "top")
print(p1)

# 7. Gotcha: fish tagged after the decision ----
# Some fish are PIT tagged inside the Lower Granite facility, after they were
# already guided into the bypass. They are in the data only because they were
# guided, so they aren't a fair sample of fish arriving at the dam. Most are
# barged (c6); the rest are returned to the river (c1, c2).
# Mixing them in makes guidance look far better than it is.

n.facility <- 1500
facility.fish <- rmultinom(1, n.facility, c(c1 = (1 - trans) * (1 - phi_B),
                                            c2 = (1 - trans) * phi_B,
                                            c6 = trans))[, 1]
n.mixed <- n
n.mixed[c("c1", "c2", "c6")] <- n.mixed[c("c1", "c2", "c6")] + facility.fish

fit.mixed <- jags(list(n = as.vector(n.mixed), N_seen = sum(n.mixed)), inits, params, "six_cell.jags",
                  n.chains = 3, n.iter = 6000, n.burnin = 2000, n.thin = 2, parallel = TRUE)

ci <- function(x) round(quantile(x, c(0.025, 0.975)), 3)
cat("\nGE truth:", psi, "\n")
cat("  Fish tagged upstream only:      ", round(fit$mean$psi, 3), " 95% CrI", ci(fit$sims.list$psi), "\n")
cat("  Plus fish tagged in the facility:", round(fit.mixed$mean$psi, 3), " 95% CrI", ci(fit.mixed$sims.list$psi), "\n")
cat("Fix: drop fish whose tagging site guarantees a route before fitting.\n")

# 8. Identifiability: no downstream site ----
# Remove GOJ. Every detected fish is now either guided (GRJ, including barged)
# or nonguided and read at GRS. GOJ-only fish (c5) become never-detected.
# One observed ratio is left:
#   guided / (guided + GRS) = psi / (psi + (1 - psi) p)
# One equation, two unknowns: low psi with poor GRS detection fits the same
# counts as higher psi with good GRS detection. The sampler still runs and
# Rhat still looks fine. The warning sign: p's posterior spreads over most of
# 0 to 1, and its shape comes from the priors, not the data (here the flat
# prior on psi tilts it toward high p). Compare posterior to prior (BPA 7.9).

cat(file = "no_goj.jags", "
model {
psi ~ dbeta(1, 1)
p ~ dbeta(1, 1)
p_seen <- psi + (1 - psi) * p
guided ~ dbin(psi / p_seen, N_seen)
}
")

guided <- sum(n[c("c1", "c2", "c6")])
grs <- sum(n[c("c3", "c4")])
fit.nogoj <- jags(list(guided = guided, N_seen = guided + grs), NULL, c("psi", "p"), "no_goj.jags",
                  n.chains = 3, n.iter = 20000, n.burnin = 5000, n.thin = 5, parallel = TRUE)
cat("\nMax Rhat without GOJ:", round(fit_rhat(fit.nogoj), 3), "\n")
cat("p:   truth", p, "  with GOJ", round(fit$mean$p, 3), ci(fit$sims.list$p),
    "  without GOJ", round(fit.nogoj$mean$p, 3), ci(fit.nogoj$sims.list$p), "\n")
cat("psi: truth", psi, "  with GOJ", round(fit$mean$psi, 3), ci(fit$sims.list$psi),
    "  without GOJ", round(fit.nogoj$mean$psi, 3), ci(fit.nogoj$sims.list$psi), "\n")

# Posterior of p without GOJ vs its flat prior
p.true <- p
p2 <- ggplot(data.frame(p = fit.nogoj$sims.list$p), aes(p)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "gray80", color = "white",
                 boundary = 0) +
  geom_hline(aes(yintercept = 1, linetype = "Flat prior")) +
  geom_vline(aes(xintercept = p.true, linetype = "Truth")) +
  scale_linetype_manual(NULL, values = c("Flat prior" = "dashed", "Truth" = "solid")) +
  labs(title = "Without GOJ, the posterior for p", x = "GRS detection p", y = "Density") +
  theme_bw(base_size = 13) + theme(legend.position = "top")
print(p2)

# Joint draws trace a ridge: every pair on it fits the data equally well
r.obs <- guided / (guided + grs)
ridge <- data.frame(p = seq(0.01, 1, length.out = 200))
ridge$psi <- r.obs * ridge$p / (1 - r.obs + r.obs * ridge$p)    # solves psi / (psi + (1 - psi) p) = r.obs
p3 <- ggplot(data.frame(p = fit.nogoj$sims.list$p, psi = fit.nogoj$sims.list$psi), aes(p, psi)) +
  geom_point(alpha = 0.15, size = 0.8, color = "steelblue") +
  geom_line(data = ridge, aes(linetype = "Same observed ratio"), linewidth = 1) +
  geom_point(data = data.frame(p = p, psi = psi), aes(shape = "Truth"), size = 4, stroke = 1.5) +
  scale_linetype_manual(NULL, values = c("Same observed ratio" = "dashed")) +
  scale_shape_manual(NULL, values = c("Truth" = 4)) +
  labs(title = "Without GOJ: psi and p trade off", x = "GRS detection p", y = "Guidance efficiency psi") +
  theme_bw(base_size = 13) + theme(legend.position = "top")
print(p3)

# 9. Incomplete transport flags: why the transport cell is there ----
# Suppose a barged fish retains its GRJ record but loses its transport flag.
# It moves from c6 to c1, so it is still guided. This differs from losing the
# fish entirely or assigning it to a different passage period.
transport.summary <- function(counts) {
  guided <- sum(counts[c("c1", "c2", "c6")])
  phat <- counts["c4"] / (counts["c4"] + counts["c5"])
  c(GE = unname(guided / (guided + sum(counts[c("c3", "c4")]) / phat)),
    transport = unname(counts["c6"] / guided),
    bypass.recovery = unname(counts["c2"] / sum(counts[c("c1", "c2")])))
}
n.flags.missing <- n
n.flags.missing["c1"] <- n["c1"] + n["c6"]
n.flags.missing["c6"] <- 0
flag.check <- rbind(complete = transport.summary(n),
                    missing = transport.summary(n.flags.missing))
print(flag.check)
stopifnot(isTRUE(all.equal(flag.check[1, "GE"], flag.check[2, "GE"])))
# That's the point: GE doesn't change, bypass recovery does. Unflagged barged
# fish look like bypassed fish that never reached GOJ, so bypass recovery
# drops. The transport cell keeps bypass recovery close to what the data
# actually do. Audit transport records and passage versus loading dates.

# Try this ----
# Predict before rerunning; reset each change before trying the next.
# 1. Set n.facility <- 5000. How far does estimated GE move from truth?
#    Does the credible interval alert you to the tagging-site problem?
# 2. In section 8 ONLY, change psi ~ dbeta(1, 1) to psi ~ dbeta(1, 9).
#    Does p's posterior move even though its own prior and the data are unchanged?
# 3. Set n.tagged <- 600. Which probability has the widest interval relative
#    to its true value, and which detection histories became sparse?

# Lessons ----
# - The route a fish takes at the dam is hidden; the detection history is what
#   we see. That makes this a multistate model.
# - Never-detected fish can't be counted, so the likelihood conditions on
#   being detected: divide each cell by the chance of being seen at all.
# - p is identified by c4 and c5: fish that skipped the bypass and were
#   detected at GOJ, read or missed at GRS. phi_S cancels.
# - Marked fish have to be a fair sample. Fish tagged after the guidance
#   decision make GE look much higher than it is, with a tight interval that
#   misses the truth.
# - Without a downstream site, only the combination of GE and p is
#   identified. The sampler still runs and Rhat still looks fine, so check
#   posteriors against priors and simulate with known truth.

# Where BPA teaches this ----
# - 9.2  Multistate model: estimation of movement between two sites (model
#        description, simulated data, analysis). Route at the dam is the state.
# - 7.10  The multinomial likelihood (m-array format): counts of detection
#         histories as one multinomial draw, as in dmulti here
# - 7.9  Parameter identifiability. Its example is the same problem as no
#        GOJ: in the time-dependent CJS model, the last survival and last
#        recapture can't be separated, only their product. It shows the
#        check used here, comparing each posterior to its flat prior
#        (Gimenez et al. 2009b), and warns the problem is greater in
#        multistate models.

# Glossary

Short definitions of the fish and statistics terms in these scripts. If a script uses a word you don't know, look here first.

## The setting in one paragraph

Young salmon and steelhead (smolts) swim down the Snake River to the ocean in spring. At Lower Granite Dam, each fish passes one of three ways: it's guided into a juvenile bypass channel, or it goes over the spillway, or it passes through the turbines. A little of each route is read at PIT-tag antennas; the models combine missing and observed detections to estimate how many fish and what share were guided.

## Fish and fisheries terms

**Barge (transport).** Some fish guided into the bypass are loaded on barges and carried past the downstream dams. A barged fish can't be read at GOJ.

**Bypass.** A channel at the dam that guides fish away from the turbines. Fish in it are read at GRJ.

**Detection history.** The set of sites that read a tagged fish, for example "GRS and GOJ". Each history is one cell in the models.

**Efficiency trial.** Marking a known number of fish, releasing them upstream of a trap, and counting how many the trap catches again. The share recaught estimates trap efficiency.

**GE (guidance efficiency, psi).** The chance a fish arriving at the dam is guided into the bypass. See also Bypass.

**GOJ.** PIT detection site in the juvenile bypass at Little Goose Dam, downstream of Lower Granite.

**GRJ.** PIT detection site in the juvenile bypass at Lower Granite Dam.

**GRS.** PIT detection array in spillbay 1 at Lower Granite. Its rows estimate read probability within that route, not the fraction of all nonguided fish using it.

**Hatchery and wild (rear type).** Fish raised in a hatchery vs fish born in the river. Wild fish are usually the target, but most tagged fish are hatchery.

**Night share (a).** The fraction of fish that pass the dam at night.

**PIT tag.** A small passive integrated transponder tag with a unique code, read by antennas without catching the fish.

**PTAGIS.** The public database of PIT tag detections in the Columbia River basin: https://www.ptagis.org/

**RSW (removable spillway weir).** A structure providing surface passage through a spillbay. At Lower Granite, it directs surface-oriented passage toward spillbay 1.

**Sample rate.** Here, the fraction of time the trap sampling gates are open. Expanding by this rate assumes sampled periods represent passage during unsampled periods.

**SCRAPI / SCRAPI2.** The smolt abundance and composition analysis workflows discussed in the talk. The teaching comparison holds GE fixed for SCRAPI-style bootstraps and carries GE posterior draws for SCRAPI2-style expansion. The examples omit several production steps.

**Screw trap.** A floating, rotating cone in a river that catches a fraction of the fish passing it.

**Smolt.** A young salmon or steelhead migrating to the ocean.

**smoltEASE.** The modeling code used for the talk's PIT-based guidance estimates, including hierarchical guidance, detection, recovery, and transport components. The scripts here are simplified teaching examples.

**Spill.** Water released over the spillway instead of through the turbines. More spill tends to send more fish over the spillway and fewer into the bypass.

**Stratum (strata).** A block of consecutive weeks grouped together.

**Trap efficiency (p).** The chance a passing fish is caught by the trap. See also Efficiency trial. (In the dam scripts, p is the chance GRS reads a fish given that it wasn't guided into the bypass, whether it then went over the spillway or through the turbines. It is not the chance of being read given spillway passage.)

## Statistics terms

**Bernoulli trial.** One yes-or-no outcome with probability p, like one fish guided or not.

**Bias.** How far an estimate lands from the truth, on average over many repeats.

**Binomial.** The number of yeses in a fixed number of independent Bernoulli trials, like recaptures out of M marked fish.

**Bootstrap.** Estimating uncertainty by recomputing an answer many times on randomly redrawn data.

**Burn-in.** The first MCMC draws, thrown away while the chains move away from their starting points.

**Chain.** One run of the MCMC sampler. Running several chains from different starts and checking they agree is how you judge convergence. See also Rhat and Trace plot.

**Complete pooling vs partial pooling.** Complete pooling forces groups to share one value. Partial pooling (a hierarchical model) lets each group have its own value but pulls it toward the group average, sharing strength without pretending all groups are identical.

**Conditioning on detection.** Fish never detected anywhere aren't in the data, so the model divides each cell's probability by the chance of being detected at all.

**Convergence.** When MCMC chains have settled into the posterior. Checked with Rhat and trace plots (see both below).

**Convex.** Curving upward, like 1 / p. Averaging a convex function gives more than applying it to the average (Jensen's inequality).

**Coverage.** Over many simulated data sets, the share of 95% intervals that contain the truth. It measures calibration under the simulated conditions. Nominal 95% Bayesian intervals do not guarantee 95% coverage at every parameter value.

**Credible interval.** A Bayesian interval: given the model and the data, the parameter is inside it with the stated probability. See also Posterior and Coverage.

**Hierarchical model.** A model where group-level values (each week's GE, say) are themselves drawn from a distribution with its own parameters.

**Identifiability.** Whether the data can pin down a parameter. If two parameters can trade off and fit the data equally well, they aren't separately identified.

**Jensen's inequality.** For a curved function, the average of the function isn't the function of the average. Here: the average of 1 / p is larger than 1 / (average p). See also Convex.

**Latent state.** A true quantity that isn't observed directly, like daily passage or a fish's route.

**Likelihood.** How probable the observed data are for each possible parameter value.

**Lincoln-Petersen.** A capture-recapture estimator using the overlap between marked and subsequently sampled animals: the basic form is M*C/R. Corrections depend on the sampling design. The talk's Bailey modification is M*(C+1)/(R+1).

**Logit.** log(p / (1 − p)). It maps probabilities (0 to 1) onto the whole number line, so effects can be added without pushing a probability past 0 or 1.

**MCMC (Markov chain Monte Carlo).** The simulation method JAGS uses to draw samples from the posterior when there is no exact formula.

**Multinomial.** The joint counts across mutually exclusive categories in a fixed number of independent trials with common category probabilities. Each fish has one categorical outcome; the counts are one multinomial draw.

**n.eff.** How many independent draws the correlated MCMC draws are worth.

**Posterior.** What you believe about a parameter after seeing the data: prior times likelihood, rescaled.

**Precision.** 1 / variance. JAGS's dnorm uses precision, not standard deviation: dnorm(0, 0.1) has SD about 3.2.

**Prior.** What you believe about a parameter before seeing the data.

**Rhat.** Compares variation between chains to variation within chains. Near 1 supports agreement. These lessons also check rank-normalized split Rhat <= 1.01, bulk/tail effective sample sizes >= 400, and finite Monte Carlo error for key outputs.

**Sampling distribution.** All the estimates you could have gotten from repeating the same study under the same truth.

**Shrinkage.** In a hierarchical model, a group's estimate pulled toward the average of all groups. See also Hierarchical model and Complete pooling vs partial pooling.

**Simpson's paradox.** An aggregate comparison reverses the comparisons within its component groups because their proportions differ. In the talk, day and night GE both fall while pooled GE rises as night passage becomes more common.

**State-space model.** A model with two parts: a process model for the true state over time, and an observation model for what you actually measure.

**Thinning.** Keeping only every k-th MCMC draw to save memory.

**Trace plot.** Each chain's draws plotted in order. Converged chains look like grassy, overlapping fuzz.

# How Many Fish Are There? Mark-Recapture with JAGS

**First, download and install JAGS:** https://sourceforge.net/projects/mcmc-jags/ (JAGS is a separate program from R; install it before the R packages, then restart R.)

Companion code for the BSU ecological statistics group talk on IDFG mark-recapture, from screw traps to Lower Granite Dam.

**[View the presentation slides](https://docs.google.com/presentation/d/1t2OOFYyH-QrZKkhVpIfSck41ejtzsgCnnv2lWwjo-_E/edit?slide=id.p1#slide=id.p1)**

Most scripts simulate data from known values, fit a model, and check the answer against the truth. Start at 01 and work up: each step adds one idea.

New to Bayesian statistics? Start with [01_first_bayesian_model.R](01_first_bayesian_model.R), and keep [GLOSSARY.md](GLOSSARY.md) open: it explains the fish terms (GE, GRJ, PIT tag, spill) and the statistics terms (prior, likelihood, posterior, Rhat).

**These are teaching examples, not the models IDFG uses for its estimates.** They are simplified so each idea stands on its own. The concepts carry over to the production models. Written with help from several fisheries and statistical collaborators.

## Requirements

- R (4.0 or later)
- JAGS 4.x, installed separately from R: https://sourceforge.net/projects/mcmc-jags/
- The jagsUI, ggplot2 (>= 3.4.0), and posterior packages: `install.packages(c("jagsUI", "ggplot2", "posterior"))`

On Windows, install JAGS before jagsUI. If R can't find JAGS, restart R after installing it.

Run [00_setup.R](00_setup.R) first. It installs any missing packages, checks that R can find JAGS, fits a tiny test model, and prints your R, JAGS, and package versions. Save that output with your notes when you compare runs or report a problem.

## Quick start

1. Download this repository using **Code → Download ZIP**, then extract it.
2. Open R or RStudio and set the working directory to the extracted folder.
3. Run the setup check:

   ```r
   source("00_setup.R")
   ```

4. Begin with [01_first_bayesian_model.R](01_first_bayesian_model.R). Run one numbered lesson at a time, reading the comments and inspecting the plots.
5. Keep the [glossary](GLOSSARY.md) handy; see [Using This Code](USING_THIS_CODE.md) to adapt an example.

The numbered R scripts and `helpers.R` stay together so they can be run directly from this folder. The helper checks rank-normalized split Rhat (<= 1.01), bulk/tail effective sample sizes (>= 400), and finite Monte Carlo error for summaries like inverse GE.

Simulation fits that error or fail diagnostics are retried once and counted if they still fail. Comparisons remain paired. Lesson 07 compares posterior-median abundance expansions because the wild-by-hatchery pooled GE is not the quantity to borrow.

## Levels

- **Start here, no Bayesian background needed:** `01`, `02`, `03`. The basics: prior, likelihood, posterior, MCMC, and the logit.
- **Next steps:** `04`, `05`. State-space and multistate models, conditioning on detection, and identifiability.
- **Graduate level:** `06`, `07`, `08`, `09`. Hierarchical models, borrowing strength, propagating uncertainty, and simulation studies.

## Scripts

| Script | Idea | Talk section |
|---|---|---|
| [00_setup.R](00_setup.R) | Installs packages, checks JAGS, runs a test model, and records versions. | |
| [01_first_bayesian_model.R](01_first_bayesian_model.R) | Prior, likelihood, and posterior for one efficiency trial, solved exactly by hand, then the same model in JAGS. How to read JAGS output, check convergence, and interpret a 95% interval. | Why Bayesian models help |
| [02_lincoln_petersen.R](02_lincoln_petersen.R) | The slide's Bailey example and summing strata, then a separate efficiency-trial simulation with its own +1 correction and Jensen's inequality. Plots show why the naive ratio can be biased. | We control the marks |
| [03_bernoulli_logit.R](03_bernoulli_logit.R) | A logit GLM in JAGS: each fish is a Bernoulli trial, and spill shifts guidance on the logit scale. | Why model GE on the logit scale |
| [04_state_space_screwtrap.R](04_state_space_screwtrap.R) | A state-space model for daily screw trap passage. Counts alone estimate N × p; adding efficiency trials recovers N. Days with the trap out can be left as missing. | State-space model for daily passage |
| [05_six_cell.R](05_six_cell.R) | The multistate model at Lower Granite: six detection cells, conditioning on being detected, and a multinomial likelihood. Then two ways it goes wrong: fish tagged after the decision and no downstream site. | Six-cell route model |
| [06_hierarchical.R](06_hierarchical.R) | Borrowing strength across weeks and rear types. Compares each wild week on its own to the hierarchical fit. | Borrow strength; the hierarchical GE regression |
| [07_day_night.R](07_day_night.R) | Day and night: five cells split by passage time, c5 stays pooled, 11 cells in all. Shows why wild should borrow hatchery's day and night GE, not its pooled GE, and how Simpson's paradox appears. | Day and night GE |
| [08_bootstrap.R](08_bootstrap.R) | The slide's one-day fixed-GE bootstrap, then weekly fixed GE versus GE draws. Coverage of each over 300 simulated weeks. Plain R, no JAGS. | Where bootstrap intervals come from |
| [09_recover_truth.R](09_recover_truth.R) | A small simulation study: bias and 95% interval coverage over 50 simulated data sets. | Does it recover GE? |

Every lesson script runs top to bottom and prints its results and a plot. Scripts that use JAGS write their models to `.jags` files in the working directory. Run times are a few seconds to a few minutes depending on your machine.

## Real data to try

The scripts use simulated data, but the same kinds of data are public:

- **PIT tag detections:** PTAGIS, https://www.ptagis.org/
- **River flow and spill at the dams:** Fish Passage Center historical flow and spill query, https://www.fpc.org/smolt/smolt_queries/Q_river_historicflowspill_dataquery.php

[USING_THIS_CODE.md](USING_THIS_CODE.md) explains how to turn fish-level detection records into the cell counts the models need.

## What to read alongside

Kéry, M., and M. Schaub. 2012. *Bayesian Population Analysis using WinBUGS: A Hierarchical Perspective.* Academic Press. The book's code was written for WinBUGS; the authors also provide a translation to JAGS and a good discussion of identifiability and complex observation models.

Each script ends with a `Where BPA teaches this` section pointing to the matching book sections.

- Chapter 1: the binomial observation process, and why simulate (scripts 02, 09)
- Chapter 2: Bayes rule, priors, posteriors, MCMC, and convergence (script 01)
- Chapter 3: GLMs, link functions, and the logit (scripts 03, 06)
- Chapter 4: random effects (script 06)
- Chapter 5: state-space models for counts (script 04)
- Chapter 6: closed-population size from capture-recapture (script 02)
- Chapter 7: time and group effects, identifiability, and the multinomial likelihood (scripts 05, 06, 09)
- Chapter 9: multistate models (scripts 05, 07)

See [USING_THIS_CODE.md](USING_THIS_CODE.md) for how to adapt these models to your own data, including how to format raw detection histories into the arrays the models need, and a list of gotchas to check before you interpret a fit.

**Before using any adapted model on real data, run it through [09_recover_truth.R](09_recover_truth.R).** It's a template, not just an exercise: if a model can't recover known values from simulated data, it won't recover unknown ones from real data.

## Regression checks

After setup, run `Rscript tests/audit_checks.R` from the repository root. These focused checks cover diagnostic failures, detection-history validation, and the finite median used for inverse-GE expansions.

## Matching the talk examples

Use slide titles because slide numbers change as the deck develops.

| Talk topic | Companion example and scope |
|---|---|
| We Control the Marks | `02`, section 0: Bailey's 9,533 example and sum across strata. The following independent-trial simulation uses a different correction. |
| Where Bootstrap Intervals Come From | `08`, section 0: 130 caught, sample rate 0.25, GE 0.13, 2,000 replicates, 90% interval. The slide's 3,446–4,554 endpoints are from a particular simulation. |
| Step 2: The Row Method | `05`, section 0: three rows with read probabilities 0.64, 0.61, 0.53 and a combined 0.934012. Row overlaps illustrate estimation under independent, equally readable fish. |
| Most of the PIT Fish Are Hatchery | `06`, section 1b: an illustrative adjacent-week pooling rule requiring both c4+c5 > 10 and c3+c4 > 10, with a trailing group merged backward. The fitted simulation retains its fixed teaching strata. |
| Simpson's Paradox | `07`, section 0: the toy reversal from 0.066 to 0.112. |
| Lessons for Your Own Models: transport flags | `05`, section 9: missing flags leave GE unchanged but drag bypass recovery down, because barged fish look like bypassed fish that never reached GOJ. |

The talk's regression includes outflow; `06` uses spill only. The talk's day/night model extends that regression; `07` uses one week without covariates. The talk carries genetic stock-assignment uncertainty; `08` leaves it out along with composition resampling.

Each lesson (`01`–`09`) ends with two or three **Try this** exercises just before its lessons summary. Predict the result, change one setting, rerun, and inspect diagnostics before interpreting it. Reset the change before trying the next exercise. Answers are intentionally omitted.

Lesson `05`, section 3b, walks through the talk's MY2025 wild counts with plain division before introducing the multinomial: naive GE, a row-read correction, then the GOJ overlap correction. These season calculations are kept separate from the simulated JAGS data and from the production hierarchical estimates.

# Using This Code

The scripts are built around fish passing a dam, but the models are general. If you mark animals and detect them again, count animals with imperfect detection, or estimate a probability that changes with covariates or time, the same steps apply: define latent states, write the observation model, and fit it with MCMC.

## What each model is good for

**01 First Bayesian model.** Any probability estimated from successes out of trials: survival, detection, occupancy, nest success. It's the template for reading any JAGS output: check Rhat and trace plots before interpreting the mean.

**02 Trap efficiency.** Any study where you release a known number of marked animals and count how many come back: trap efficiency trials, mark-resight, survey recapture. Use it to see how sampling probabilities and ratios interact, and how a simple +1 correction changes the answer.

**03 Logit GLM.** Any yes/no outcome with covariates: survival to a site, presence at a plot, nest success, a fish taking one route or another. The pattern (Bernoulli trial, logit link, linear predictor) is the backbone of many modern ecological models.

**04 State-space model.** Repeated counts over time where detection is below 1: daily trap catches, weekly redd counts, annual surveys. Two lessons carry over:

- Counts alone estimate an index (N × p), not N. To get N you need separate data on detection, such as efficiency trials, double observers, or marked animals.
- Missing days or years can be left as `NA`. The model estimates them from the days around them.

**05 Six-cell multistate model.** Any system where animals take one of several hidden routes or states, and you see them at some detection sites and not others: dam passage, river junctions, breeding stages, or habitat use. The detection history is the data; the state is the process.

**06 Hierarchical model.** Many small groups that are alike but not identical: weeks, sites, streams, years. Partial pooling borrows strength from the rest, and well-sampled groups mostly speak for themselves.

**07 Day and night.** What to do when one state splits further (day and night here), and why it matters when groups borrow from each other. If two groups share the same underlying rates but differ in their mix of time periods, a pooled estimate can be misleading.

**08 Bootstrap.** Any estimate built by expanding a count by an estimated rate (catch / efficiency, counts / detection). If the rate is uncertain, the interval has to carry that uncertainty too, or it will be too narrow.

**09 Recover the truth.** Not an optional exercise: a template for every model you adapt. Before trusting a modified model on real data, simulate from it with known values and confirm it gets the right answer.

## Adapting a script to your data

1. **Run the script unchanged first.** Make sure JAGS works and you get the printed results and plot.
2. **Find the data the model needs.** Each JAGS call takes a `jags.data` list. That list is the only thing your real data have to match. For example:
   - `04`: a vector of daily counts `C` (with `NA` for missed days), plus trial releases `M` and recaptures `R` with the day each trial ran.
   - `05`: a vector of six cell counts `n` and their total `N_seen`.
   - `06`: an array `n[group, time, cell]`, plus covariates and a stratum index.
   - `07`: a matrix `n11[group, cell]` with 11 cells in a fixed order (day cells, night cells, then c5 pooled).
3. **Format your data into that shape.** This is where most analyses go wrong, not the math. `05_six_cell.R` section 2 starts from one row per tagged fish, assigns each fish to exactly one cell, and tabulates the counts.
4. **Replace the simulation section with your data.** Keep the simulation code: you'll need it for step 6.
5. **Change the model to match your system.** Redraw the detection tree for your sites, write one `pi_raw` line per detectable history, and keep the conditioning lines (`p_seen` and `pi_obs`) as the last step before the likelihood.
6. **Simulate from your new model and refit it** using `09` as the template, before using real data. If it can't recover known values, it won't recover unknown ones. Do this every time you change the structure.

## Where to get real data

- **PIT tag detections:** PTAGIS, https://www.ptagis.org/. Detection records come one row per tag and site; step 3 above turns them into cell counts.
- **River flow and spill:** Fish Passage Center historical flow and spill query, https://www.fpc.org/smolt/smolt_queries/Q_river_historicflowspill_dataquery.php. Covariates like the spill effect are available there.

## Gotchas that bite real studies

- **Marked animals that aren't a fair sample.** Animals tagged after the outcome you're estimating (fish tagged inside the bypass, after they were already guided) make the rate look far better than it is, with a tight interval that misses the truth. Check where and when every tag was applied, and drop tags whose release site guarantees a route (`05`, section 7).
- **Incomplete transport flags.** Check transport records and distinguish passage time from loading time. Missing flags leave GE unchanged but drag bypass recovery down, because barged fish look like bypassed fish that never reached GOJ (`05`, section 9).

- **A parameter the data can't separate.** Without the right detection site, two parameters trade off along a ridge. The sampler still runs and Rhat still looks fine (`05`, section 8). Ask what part of the data identifies each parameter before you fit.
- **Expanding by an uncertain rate.** 1 / rate is curved, so expanding by the average rate and averaging the expansions give different answers, and neither is guaranteed to be unbiased (`02` and `08`). Carry the rate's uncertainty into the interval, and check the point estimate against simulated truth.
- **An inverse with no finite mean.** With no guided fish, the wild-only flat-GE prior in `07` yields an infinite posterior mean of `1 / GE`. A finite MCMC average hides this. Lesson `07` uses posterior draws and a finite median before inversion.
- **A conditional expansion labeled as abundance.** `05` holds `N_seen` fixed and propagates detection-probability uncertainty through `N_seen / p_seen`. Its expansion interval is not a full abundance posterior.
- **An interval that ignores an input.** A bootstrap that holds an estimated rate fixed gives an interval that's far too narrow (`08`).
- **Borrowing a mixed quantity.** If groups differ in timing or habitat use, a pooled rate mixes the thing you want with the mix itself. Borrow the components, not the blend (`07`).
- **Monitoring the wrong parameter.** When offsets are centered to sum to zero, their common mean cancels and the data say nothing about it; its posterior is just its prior. Rhat on the raw offset is not enough; inspect the parameter itself.
- **Dropping fits that fail diagnostics.** In a simulation study, the fits that fail are often the hard data sets. Retry them with longer chains, keep comparisons paired, and report the failure rate.

## Checks to run every time

- **Convergence.** Use `fit_diagnostics()` from `helpers.R`: rank-normalized split Rhat <= 1.01 and bulk/tail effective sample sizes >= 400 are initial screens. Missing diagnostics fail the check.
- **Prior vs posterior.** If a parameter's posterior looks just like its prior, the data aren't informing it. In a Bayesian fit an unidentified parameter still gets a posterior (BPA 7.9).
- **Coverage.** Measure how often intervals contain the truth under the conditions you simulate. Nominal 95% Bayesian intervals do not guarantee 95% coverage at every fixed parameter value. Report the simulation design, not just the nominal interval.
- **Use the draws, not the summary.** When an answer is a nonlinear function of parameters (1 / p, a product, a ratio), compute it for every posterior draw and then summarize. That carries the uncertainty through the transform.

## Common extensions

- Covariates on any probability: add terms to its logit linear predictor, as in `03` and `06`.
- Random effects for time or site: replace a fixed parameter with draws from a shared distribution, as in `06`.
- More detection sites: add branches to the tree and cells to the multinomial.
- Combining data sets: add another likelihood line to the same model. The parameters they share are estimated from both.

## Questions

David Smith, Idaho Department of Fish and Game, david.smith@idfg.idaho.gov

# Shared checks for the JAGS lessons. Keep this file beside the scripts.
# Modern rank-normalized split Rhat and bulk/tail ESS supplement the
# traditional diagnostics printed by jagsUI. Thresholds are initial screens;
# inspect traces and Monte Carlo error for the actual summaries of interest.

fit_diagnostics <- function(fit) {
  if (is.null(fit) || is.null(fit$samples)) return(NULL)
  tryCatch({
    chains <- lapply(fit$samples, as.matrix)
    if (length(chains) < 2L || any(vapply(chains, nrow, integer(1)) < 4L))
      return(NULL)
    if (!all(vapply(chains, function(x) identical(dim(x), dim(chains[[1]])), logical(1))))
      return(NULL)
    draws <- array(NA_real_, c(nrow(chains[[1]]), length(chains), ncol(chains[[1]])),
                   dimnames = list(NULL, NULL, colnames(chains[[1]])))
    for (i in seq_along(chains)) draws[, i, ] <- chains[[i]]
    if (any(!is.finite(draws))) return(NULL)
    # Inverse GE is reported in lesson 07: check its sampling precision too.
    if ("wild_ge" %in% dimnames(draws)[[3]]) {
      inverse <- 1 / draws[, , "wild_ge"]
      if (any(!is.finite(inverse))) return(NULL)
      extended <- array(NA_real_, dim(draws) + c(0L, 0L, 1L),
                        dimnames = list(NULL, NULL, c(dimnames(draws)[[3]], "inverse_ge")))
      extended[, , seq_len(dim(draws)[3])] <- draws
      extended[, , dim(extended)[3]] <- inverse
      draws <- extended
    }
    # Pass the functions themselves: measure names given as strings are looked
    # up by name, so a variable called rhat in the workspace would break them.
    posterior::summarise_draws(posterior::as_draws_array(draws),
                              rhat = posterior::rhat,
                              ess_bulk = posterior::ess_bulk,
                              ess_tail = posterior::ess_tail,
                              mcse_median = posterior::mcse_median)
  }, error = function(e) NULL)
}

fit_ok <- function(fit) {
  d <- fit_diagnostics(fit)
  !is.null(d) && nrow(d) > 0L &&
    all(is.finite(as.matrix(d[, c("rhat", "ess_bulk", "ess_tail", "mcse_median")]))) &&
    all(d$rhat <= 1.01 & d$ess_bulk >= 400 & d$ess_tail >= 400)
}

fit_rhat <- function(fit) {
  d <- fit_diagnostics(fit)
  if (is.null(d) || !nrow(d) || any(!is.finite(d$rhat))) return(Inf)
  max(d$rhat)
}

# Simulation loops retry NULL just as they retry poor mixing. Keep and count
# the failure instead of aborting a study or treating absent diagnostics as OK.
safe_jags <- function(...) {
  tryCatch(jagsUI::jags(...), error = function(e) {
    message("Fit failed: ", conditionMessage(e))
    NULL
  })
}

# state-space model file ends here
