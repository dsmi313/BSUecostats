# Trap efficiency at a screw trap
# (Run 00_setup.R once first: it installs packages and checks JAGS.)
#
# One stratum at a rotary screw trap. Some passing fish are caught. Separately,
# M fish are marked and released upstream (an efficiency trial), and R of them
# are caught again. R / M estimates trap efficiency, and the unmarked catch C
# divided by efficiency estimates how many fish passed.
#
# Two estimators of passage N:
#   Simple expansion:   C * M / R          (catch divided by R / M; it has the
#                                           Lincoln-Petersen form)
#   +1 adjusted:        C * (M + 1) / (R + 1)
# The +1 version is the same idea as Chapman's correction, written for this
# setup: here C is the unmarked catch and the marked fish come from a
# separate trial, which is not the two-sample design Chapman's formula
# assumes. For R ~ Binomial(M, p), E[(M + 1) / (R + 1)] = (1 - (1 - p)^(M + 1)) / p,
# so the adjusted estimator averages N * (1 - (1 - p)^(M + 1)): nearly
# unbiased unless M is tiny or p is very small.
#
# Plain R, no JAGS. All data are simulated, so we know the truth.
# This is a teaching example, not the estimator IDFG uses.

library(ggplot2)

# 0. The talk's exact toy example: Bailey's modification ----
# Keep the slide's estimator and numbers separate from the independent-trial
# simulation below. Swapping M and C inside the +1 terms changes the estimator.
bailey <- function(M, C, R) M * (C + 1) / (R + 1)
slide.M <- 200; slide.C <- 1000; slide.R <- 20
cat("Talk example, Bailey:", round(bailey(slide.M, slide.C, slide.R)),
    "  plain Petersen:", slide.M * slide.C / slide.R, "\n")
stopifnot(round(bailey(slide.M, slide.C, slide.R)) == 9533)
# With multiple disjoint strata, apply the formula to each stratum and sum.
# This second row is an additional hypothetical stratum, not a slide value.
strata.example <- data.frame(stratum = 1:2, M = c(200, 100),
                            C = c(1000, 400), R = c(20, 10))
strata.example$N.bailey <- with(strata.example, bailey(M, C, R))
print(strata.example)
cat("Sum of stratum estimates:", sum(strata.example$N.bailey), "\n")
# The remaining lesson uses a different design: C is an independent unmarked
# binomial catch and R comes from a separate release trial. Its adjustment is
# C*(M+1)/(R+1), which gives 9,571 for the slide inputs, NOT Bailey's 9,533.
# The expectation identity below applies to that simulation's adjustment.

set.seed(1)

# 1. One stratum ----

N <- 4000      # True number of fish passing the trap
p <- 0.15      # True trap efficiency (chance a passing fish is caught)
M <- 50        # Marked fish released upstream (few, so efficiency is uncertain)

C <- rbinom(1, N, p)   # Unmarked catch
R <- rbinom(1, M, p)   # Marked fish recaptured

p.hat <- R / M
N.simple <- C * M / R
N.adj <- C * (M + 1) / (R + 1)

cat("Catch C =", C, " Recaptures R =", R, "of M =", M, "\n")
cat("Efficiency estimate R / M =", round(p.hat, 3), " (truth", p, ")\n")
cat("Simple expansion N =", round(N.simple), "\n")
cat("+1 adjusted N      =", round(N.adj), "\n")
cat("Truth              =", N, "\n")

# 2. Repeat the same stratum many times ----
# Rerun this one stratum 10,000 times with the same truth (N, p, M). Only
# chance changes: which fish get caught and which marked fish come back.
# This is the sampling distribution: everything you could have gotten. In
# real life you get one of these. Not 10,000 different days or strata.

n.sims <- 10000
C.sim <- rbinom(n.sims, N, p)
R.sim <- rbinom(n.sims, M, p)
simple.sim <- C.sim * M / R.sim
adj.sim <- C.sim * (M + 1) / (R.sim + 1)

cat("\nAcross", n.sims, "simulated repeats of the same stratum:\n")
# With few marked fish, some repeats recapture none. The simple expansion
# then divides by zero; the +1 terms keep the adjusted version finite.
cat("Repeats with zero recaptures:", sum(R.sim == 0), "(simple expansion is infinite there)\n")
ok <- R.sim > 0
cat("Mean simple expansion (R > 0 only):", round(mean(simple.sim[ok])),
    " bias", round(100 * (mean(simple.sim[ok]) / N - 1), 1), "%\n")
cat("Mean +1 adjusted:                  ", round(mean(adj.sim)),
    " bias", round(100 * (mean(adj.sim) / N - 1), 1), "%\n")
cat("Expected from the formula:         ", round(N * (1 - (1 - p)^(M + 1))), "\n")
cat("Middle 95% of +1 adjusted:", round(quantile(adj.sim, c(0.025, 0.975))), "\n")

# 3. Jensen's inequality ----
# Efficiency is uncertain. Expanding by 1 / (average efficiency) is not the
# same as averaging the expansions. Draw efficiency from its posterior
# under a flat prior, Beta(R + 1, M - R + 1), and compare.

p.draws <- rbeta(10000, R + 1, M - R + 1)
cat("\nC / mean(p)  =", round(C / mean(p.draws)), "\n")
cat("mean(C / p)  =", round(mean(C / p.draws)), "\n")
cat("1 / p is convex, so the second is always larger. The gap grows as p gets small.\n")

# 4. Plot ----
# Each plot prints on its own.

p1 <- ggplot(data.frame(est = adj.sim), aes(est)) +
  geom_histogram(bins = 60, fill = "gray80", color = "white") +
  geom_vline(aes(xintercept = N, linetype = "Truth")) +
  geom_vline(aes(xintercept = mean(adj.sim), linetype = "Mean estimate"), color = "steelblue") +
  scale_linetype_manual(NULL, values = c("Truth" = "solid", "Mean estimate" = "dashed")) +
  labs(title = "+1 adjusted estimates, 10,000 simulated repeats of the same stratum", x = "Estimated N", y = "Count") +
  theme_bw(base_size = 13) + theme(legend.position = "top")
print(p1)

p2 <- ggplot(data.frame(x = C / p.draws), aes(x)) +
  geom_histogram(bins = 60, fill = "gray80", color = "white") +
  geom_vline(aes(xintercept = C / mean(p.draws), color = "C / mean(p)"), linewidth = 1) +
  geom_vline(aes(xintercept = mean(C / p.draws), color = "mean(C / p)"), linewidth = 1) +
  scale_color_manual(NULL, values = c("C / mean(p)" = "darkorange", "mean(C / p)" = "steelblue")) +
  coord_cartesian(xlim = c(0, quantile(C / p.draws, 0.995))) +   # trim the rare extreme draws from the axis
  labs(title = "Expansion C / p over efficiency draws", x = "C / p", y = "Count") +
  theme_bw(base_size = 13) + theme(legend.position = "top")
print(p2)

# Try this ----
# Predict before rerunning; reset each change before trying the next.
# 1. Set M <- 15 in the simulation. How many repeats have zero recaptures?
#    What happens to the simple estimator and to the adjusted estimates?
# 2. Set p <- 0.03. Does the +1 correction still give a mean close to truth?
#    How does the printed expectation help interpret the simulation?

# Lessons ----
# - Trap efficiency comes from marked fish: R / M. Abundance is catch divided
#   by efficiency, so every error in efficiency carries straight into N.
# - The simple expansion runs high when few marked fish are recaptured, and
#   is infinite when none are. Adding 1 to M and R corrects nearly all of that.
# - Match the estimator to the design: Chapman's formula is for two samples
#   from the same population; an efficiency trial is a different setup.
# - Efficiency is uncertain, and 1 / efficiency is curved. Averaging the
#   expansions, mean(C / p), is always larger than expanding by the average,
#   C / mean(p). That gap is Jensen's inequality, and it grows as p gets small.

# Where BPA teaches this ----
# Kery and Schaub (2012), section numbers as in the book and its JAGS code:
# - 1.3  The binomial distribution as a canonical description of the
#        observation process (catch and recaptures here are binomial draws)
# - 6.2.1  Closed-population size from capture-recapture, the simplest case
#          (model M0): estimating N from marked animals in a Bayesian model.
