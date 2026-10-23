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
#   Model 1: counts only (BPA Chapter 5 as written). It tracks N * p.
#   Model 2: counts plus efficiency trials. It recovers N.
#
# This is a teaching example, not the model IDFG uses to estimate screw trap
# abundance. The concepts carry over: a process model for true abundance, an
# observation model for what the trap sees, and separate data on detection.
#
# All data are simulated, so we know the truth and can check the answer.
# Requires JAGS (https://sourceforge.net/projects/mcmc-jags/) and the jagsUI package.

library(jagsUI)
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

# 2. Model 1: counts only (BPA Chapter 5 as written) ----
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
# Rhat near 1 (below 1.1) for every monitored node.
cat("Model 1 max Rhat:", round(max(unlist(fit1$Rhat), na.rm = TRUE), 3), "\n")
cat("Model 2 max Rhat:", round(max(unlist(fit2$Rhat), na.rm = TRUE), 3), "\n")

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

# Lessons ----
# - A state-space model separates the process (how many fish really pass each
#   day) from the observation (how many the trap catches).
# - Counts alone estimate N * p, an index, not N (BPA 5.3 and 5.5). If p
#   drifts over the season, the index drifts with it and looks like a change
#   in abundance.
# - Adding separate data on detection (efficiency trials) is what turns the
#   index into abundance.
# - Days with the trap out can be left as NA. The model fills them in from
#   the days around them, with honest uncertainty.

# Where BPA teaches this ----
# - 5.2  A simple state-space model
# - 5.3  Systematic bias in the observation process: counts track N * p
# - 5.4  House martin example: the log-scale model used here, and NA data
#        estimated as unknowns
# - 5.5  Summary and outlook: counts alone give an index; getting N needs
#        separate information on detection
