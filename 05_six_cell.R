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
library(ggplot2)

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
  ifelse(barged == 1,             "c6",     # transported
  ifelse(GRJ == 1 & GOJ == 1,     "c2",     # GRJ and GOJ
  ifelse(GRJ == 1,                "c1",     # GRJ only
  ifelse(GRS == 1 & GOJ == 1,     "c4",     # GRS and GOJ
  ifelse(GRS == 1,                "c3",     # GRS only
  ifelse(GOJ == 1,                "c5",     # GOJ only
                                  NA))))))
}
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

# Derived: how many tagged fish reached the dam, detected or not
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
cat("Max Rhat:", round(max(unlist(fit$Rhat), na.rm = TRUE), 3), "\n")
truth <- c(psi = psi, p = p, trans = trans, phi_S = phi_S, phi_B = phi_B)
est <- sapply(names(truth), function(k) fit$mean[[k]])
lo <- sapply(names(truth), function(k) quantile(fit$sims.list[[k]], 0.025))
hi <- sapply(names(truth), function(k) quantile(fit$sims.list[[k]], 0.975))
print(round(cbind(truth, estimate = est, lower = lo, upper = hi), 3))
cat("Tagged fish reaching the dam: truth", n.tagged, " estimate", round(fit$mean$N_total),
    " 95% CrI", round(quantile(fit$sims.list$N_total, c(0.025, 0.975))), "\n")

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
cat("\nMax Rhat without GOJ:", round(max(unlist(fit.nogoj$Rhat), na.rm = TRUE), 3), "\n")
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
