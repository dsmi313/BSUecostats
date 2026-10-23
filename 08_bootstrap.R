# From GE to abundance: whose uncertainty does the bootstrap carry?
# (Run 00_setup.R once first: it installs packages and checks JAGS.)
#
# Lower Granite abundance expands the trap catch by the chance a passing fish
# is caught: the sampling rate times guidance efficiency (GE).
#   daily passage = daily catch / (sampling rate x GE)
#
# Two bootstraps for a week's total:
#   Fixed GE (SCRAPI-style): each replicate redraws each day's catch from a
#     binomial, rbinom(1, dailypass, Ptrue), with GE held at one value. The
#     interval carries trap-sampling noise but nothing about GE, even though
#     GE is the shakiest input.
#   GE draws (SCRAPI2-style): each replicate also gets its own GE draw from
#     the GE posterior, and the point estimate is the mean of the replicates.
#     Because 1 / GE is convex, the exact average of the expansion over the GE
#     posterior is above catch / (rate x mean GE) (Jensen's inequality,
#     02_lincoln_petersen.R). A finite bootstrap adds catch-redraw noise, so a
#     single week's replicate mean isn't guaranteed to land above it.
#
# This is a hybrid: catch noise comes from the bootstrap, GE uncertainty from
# posterior draws. It is not a full Bayesian posterior for abundance.
#
# Teaching version: plain R, no JAGS. The GE posterior here is a Beta from a
# simple PIT sample, standing in for the posterior draws from 05 to 07.
# Real SCRAPI also resamples measured fish for the composition splits; that
# part is left out. All data are simulated.

library(ggplot2)

set.seed(8)

# 1. True values ----

N.day <- c(3000, 4500, 6000, 8000, 7000, 5000, 3500)   # true daily passage
N.week <- sum(N.day)
rate <- 0.25          # sampling rate
GE <- 0.05            # true guidance efficiency
n.pit <- 600          # PIT-tagged fish informing GE this week
B <- 2000             # bootstrap replicates

# 2. One week, both bootstraps ----

one.week <- function() {
  C <- rbinom(length(N.day), N.day, rate * GE)            # daily trap catch
  g <- rbinom(1, n.pit, GE)                               # PIT fish guided
  ge.draws <- rbeta(B, g + 1, n.pit - g + 1)              # GE posterior draws
  ge.mean <- mean(ge.draws)

  # Fixed GE: one GE for every replicate
  pass.fixed <- C / (rate * ge.mean)
  plug.in <- sum(pass.fixed)
  boot.fixed <- replicate(B, sum(rbinom(length(C), round(pass.fixed), rate * ge.mean) / (rate * ge.mean)))

  # GE draws: each replicate carries its own GE
  boot.draws <- sapply(ge.draws, function(ge) {
    pass <- C / (rate * ge)
    sum(rbinom(length(C), round(pass), rate * ge) / (rate * ge))
  })

  list(plug.in = plug.in, boot.fixed = boot.fixed, boot.draws = boot.draws,
       est.draws = mean(boot.draws), med.draws = median(boot.draws))
}

wk <- one.week()
ci <- function(x) round(quantile(x, c(0.025, 0.975)))
cat("True weekly passage:", N.week, "\n")
cat("Fixed GE:  plug-in estimate", round(wk$plug.in), "  95% interval", ci(wk$boot.fixed), "\n")
cat("GE draws:  mean of replicates", round(wk$est.draws), "  95% interval", ci(wk$boot.draws), "\n")
cat("Plug-in sits at the", round(100 * mean(wk$boot.draws < wk$plug.in)),
    "th percentile of the GE-draw replicates\n")

# 3. Repeat over many weeks ----

n.weeks <- 300
res <- t(replicate(n.weeks, {
  w <- one.week()
  c(plug.in = w$plug.in, est.draws = w$est.draws, med.draws = w$med.draws,
    lo.fixed = unname(quantile(w$boot.fixed, 0.025)), hi.fixed = unname(quantile(w$boot.fixed, 0.975)),
    lo.draws = unname(quantile(w$boot.draws, 0.025)), hi.draws = unname(quantile(w$boot.draws, 0.975)))
}))
res <- as.data.frame(res)

cover.fixed <- mean(res$lo.fixed <= N.week & N.week <= res$hi.fixed)
cover.draws <- mean(res$lo.draws <= N.week & N.week <= res$hi.draws)
cat("\nOver", n.weeks, "simulated weeks:\n")
cat("95% interval coverage   fixed GE:", round(100 * cover.fixed, 1), "%   GE draws:", round(100 * cover.draws, 1), "%\n")
cat("Median interval width   fixed GE:", round(median(res$hi.fixed - res$lo.fixed)),
    "   GE draws:", round(median(res$hi.draws - res$lo.draws)), "\n")
cat("Mean % error vs truth   plug-in:", round(100 * (mean(res$plug.in) / N.week - 1), 1),
    "%   mean of GE-draw replicates:", round(100 * (mean(res$est.draws) / N.week - 1), 1),
    "%   median of GE-draw replicates:", round(100 * (mean(res$med.draws) / N.week - 1), 1), "%\n")
cat("Weeks where the plug-in falls below the mean of the GE-draw replicates:",
    round(100 * mean(res$plug.in < res$est.draws)), "%\n")

# 4. Plot ----
# Each plot prints on its own.

boot.df <- rbind(data.frame(method = "Fixed GE", total = wk$boot.fixed),
                 data.frame(method = "GE draws", total = wk$boot.draws))
lines.df <- data.frame(what = c("Truth", "Plug-in (catch / mean GE)", "Mean of GE-draw replicates"),
                       x = c(N.week, wk$plug.in, wk$est.draws))
p1 <- ggplot(boot.df, aes(total, fill = method)) +
  geom_density(alpha = 0.4, color = NA) +
  geom_vline(data = lines.df, aes(xintercept = x, linetype = what)) +
  scale_fill_manual(NULL, values = c("Fixed GE" = "darkorange", "GE draws" = "steelblue")) +
  scale_linetype_manual(NULL, values = c("Truth" = "solid", "Plug-in (catch / mean GE)" = "dotted",
                                         "Mean of GE-draw replicates" = "dashed")) +
  labs(title = "One week: bootstrap replicates of weekly passage", x = "Weekly passage", y = "Density") +
  theme_bw(base_size = 13) + theme(legend.position = "top", legend.box = "vertical")
print(p1)

show <- 1:40
int.df <- rbind(data.frame(week = show, method = "Fixed GE", lo = res$lo.fixed[show], hi = res$hi.fixed[show]),
                data.frame(week = show, method = "GE draws", lo = res$lo.draws[show], hi = res$hi.draws[show]))
int.df$covers <- ifelse(int.df$lo <= N.week & N.week <= int.df$hi, "Covers truth", "Misses truth")
p2 <- ggplot(int.df, aes(y = week, xmin = 100 * (lo / N.week - 1), xmax = 100 * (hi / N.week - 1),
                         color = covers)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_linerange(linewidth = 0.8) +
  facet_wrap(~ method) +
  scale_color_manual(NULL, values = c("Covers truth" = "steelblue", "Misses truth" = "firebrick")) +
  labs(title = "95% intervals for 40 simulated weeks", x = "Interval relative to truth (%)", y = "Week") +
  theme_bw(base_size = 13) + theme(legend.position = "top")
print(p2)

# Lessons ----
# - An interval has to carry the uncertainty of every estimated input. Holding
#   GE fixed gives an interval that is far too narrow and misses the truth
#   most of the time.
# - Giving every bootstrap replicate its own GE draw improved coverage to
#   about 95% in this simulation.
# - Propagating uncertainty and removing bias are different jobs. Carrying GE
#   draws widened the intervals to an honest width, but the inverse also
#   raised the mean of the replicates: in this simulation it ran a few percent
#   high, while the plug-in and the median landed closer to the truth.
# - The mean of the replicates sat above catch / (rate x mean GE) in every
#   simulated week here. Jensen's inequality guarantees that ordering for the
#   exact average over GE; the catch redraws add noise, so it isn't
#   guaranteed week by week.
# - Which center to report is a choice. Simulate your own setup to see the
#   trade.

# Where BPA teaches this ----
# BPA doesn't use the bootstrap. The Bayesian habit it teaches throughout is
# the same idea as giving each replicate its own GE draw: compute derived
# quantities from every posterior draw, then summarize. Example:
# - 5.4  The probability that the house martin population in 2015 is below
#        2009, computed from the posterior draws of N
