# Your first Bayesian model: how efficient is the trap?
# (Run 00_setup.R once first: it installs packages and checks JAGS.)
# (The model fit in section 3 takes a few seconds to run.)
#
# What this script should leave you with: a Bayesian analysis is three
# pieces. A prior says what you believed before the data. A likelihood says
# how likely the data are for each possible answer. The posterior combines
# the two into what you believe after the data. For this tiny model we can
# do it exactly by hand, then let JAGS do it by simulation and check that
# the two agree. Every later script is the same idea with a bigger model.
#
# The question: a screw trap (a floating, rotating cone that catches
# downstream-migrating juvenile fish) catches some fraction p of the fish
# passing it. To estimate p, biologists mark M fish, release them upstream,
# and count how many, R, the trap catches again. That is an efficiency trial.
#
# Teaching example; the numbers are made up.

library(jagsUI)
library(ggplot2)

set.seed(1)

# 1. The data ----

M <- 50    # marked fish released upstream
R <- 8     # of those, caught again by the trap

cat("Recaptured", R, "of", M, "marked fish. Simple estimate R / M =", R / M, "\n")

# 2. Bayes by hand ----
#
# Prior: before the trial, every efficiency from 0 to 1 is equally plausible.
#   That is a Beta(1, 1) distribution, which is flat.
# Likelihood: each marked fish is caught or not, with probability p, so the
#   number recaptured is Binomial(M, p). For a given p, dbinom(R, M, p) is
#   how likely it was to see exactly R recaptures.
# Posterior: proportional to prior x likelihood. Here the math works out
#   exactly: a Beta(1, 1) prior with binomial data gives a
#   Beta(1 + R, 1 + M - R) posterior.
#
# The grid below does the same thing numerically: try 1,000 values of p,
# multiply prior by likelihood for each, and rescale so the total is 1.

p.grid <- seq(0.0005, 0.9995, length.out = 1000)
prior <- dbeta(p.grid, 1, 1)
likelihood <- dbinom(R, M, p.grid)
posterior.grid <- prior * likelihood
posterior.grid <- posterior.grid / sum(posterior.grid * diff(p.grid)[1])   # rescale to a density
posterior.exact <- dbeta(p.grid, 1 + R, 1 + M - R)

cat("Largest gap between grid and exact posterior:",
    signif(max(abs(posterior.grid - posterior.exact)), 2), "(essentially zero)\n")

post.mean <- (1 + R) / (2 + M)
post.ci <- qbeta(c(0.025, 0.975), 1 + R, 1 + M - R)
cat("Posterior mean of p:", round(post.mean, 3), "\n")
cat("95% credible interval:", round(post.ci, 3), "\n")
# A 95% credible interval means: given the model and these data, there is
# a 95% probability that p lies in this range. (A frequentist confidence
# interval does not mean that; see BPA 2.4.)

bayes.df <- rbind(
  data.frame(p = p.grid, density = prior, piece = "Prior (flat)"),
  data.frame(p = p.grid, density = likelihood / sum(likelihood * diff(p.grid)[1]),
             piece = "Likelihood (rescaled)"),
  data.frame(p = p.grid, density = posterior.exact, piece = "Posterior"))
# Draw the likelihood last and thin, so it shows on top of the posterior
bayes.df$piece <- factor(bayes.df$piece, levels = c("Prior (flat)", "Posterior", "Likelihood (rescaled)"))
p1 <- ggplot(bayes.df, aes(p, density, color = piece, linetype = piece, linewidth = piece)) +
  geom_line() +
  scale_linewidth_manual(NULL, values = c("Prior (flat)" = 1, "Posterior" = 2.5, "Likelihood (rescaled)" = 0.8)) +
  scale_color_manual(NULL, values = c("Prior (flat)" = "gray50", "Likelihood (rescaled)" = "darkorange",
                                      "Posterior" = "steelblue")) +
  scale_linetype_manual(NULL, values = c("Prior (flat)" = "dashed", "Likelihood (rescaled)" = "dotted",
                                         "Posterior" = "solid")) +
  coord_cartesian(xlim = c(0, 0.5)) +
  labs(title = "Prior x likelihood = posterior (with a flat prior, the two lines overlap)",
       x = "Trap efficiency p", y = "Density") +
  theme_bw(base_size = 13) + theme(legend.position = "top")
print(p1)
# With a flat prior, the posterior has the same shape as the likelihood:
# the data do all the talking.

# 3. The same model in JAGS ----
#
# Most models have no exact answer, so we let JAGS draw samples from the
# posterior instead (Markov chain Monte Carlo, MCMC). If JAGS works, its
# draws should match the exact Beta posterior above.
#
# Model in words: p has a flat prior; R is binomial with M tries and
# success probability p. That's the whole model.
#
# Two JAGS habits to know:
# - dbin(p, M) takes the probability FIRST and the number of tries second.
#   R's rbinom() and dbinom() put the number of tries first. Easy to mix up.
# - In later scripts, dnorm(mean, precision) takes a PRECISION, not a
#   standard deviation: precision = 1 / SD^2. dnorm(0, 0.1) has SD about 3.2.

cat(file = "first_model.jags", "
model {
p ~ dbeta(1, 1)      # prior: flat from 0 to 1
R ~ dbin(p, M)       # likelihood: R recaptures out of M, each with probability p
}
")

jags.data <- list(R = R, M = M)    # names must match the names in the model
params <- "p"                        # what to save

# MCMC settings:
# - n.chains: run 3 independent chains from different starting points. If
#   they all end up in the same place, that's evidence the sampler found the
#   posterior.
# - n.burnin: throw away the first draws, while chains wander in from their
#   starting points.
# - n.iter: total draws per chain, including burn-in.
# - n.thin: keep every k-th draw. 1 keeps them all; larger values save memory
#   when draws are strongly correlated.
fit <- jags(jags.data, inits = NULL, params, "first_model.jags",
            n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)

# 4. Reading the JAGS output ----
print(fit)
# What the columns mean:
# - mean, sd: the posterior mean and standard deviation of the draws
# - 2.5%, 97.5%: the 95% credible interval from the draws
# - Rhat: compares the spread between chains to the spread within chains.
#   Near 1 means the chains agree. Above about 1.1 means run longer or fix
#   the model (BPA 2.5).
# - n.eff: how many independent draws the correlated draws are worth
# - deviance: a measure of fit that jagsUI adds; ignore it for now
# jagsUI also prints overlap0, f, and DIC with its own short explanations.
# They're useful later for regression slopes and model comparison; you can
# skip them here.

traceplot(fit, parameters = "p")
# A trace plot shows each chain's draws in order. Good chains look like
# "grassy" fuzz around a flat level, with the chains overlapping (BPA 3.3.1).
# A trend, or chains in different places, means not converged.

draws <- fit$sims.list$p     # every saved draw of p, all chains together
# These should match the exact answer from section 2, up to small MCMC noise.
# JAGS isn't a different method: it's a numerical engine for the same
# prior x likelihood calculation, for models too big to solve by hand.
cat("\nJAGS posterior mean:", round(mean(draws), 3), "  exact:", round(post.mean, 3), "\n")
cat("JAGS 95% interval:  ", round(quantile(draws, c(0.025, 0.975)), 3),
    "  exact:", round(post.ci, 3), "\n")

p2 <- ggplot(data.frame(p = draws), aes(p)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "gray80", color = "white") +
  geom_line(data = data.frame(p = p.grid, density = posterior.exact), aes(p, density),
            color = "steelblue", linewidth = 1) +
  coord_cartesian(xlim = c(0, 0.5)) +
  labs(title = "JAGS draws (bars) vs the exact posterior (line)", x = "Trap efficiency p", y = "Density") +
  theme_bw(base_size = 13)
print(p2)

# 5. When the prior matters ----
# Suppose earlier trials at this trap suggested efficiency near 0.15. A
# Beta(3, 17) prior says that (mean 0.15), with the weight of about 20 fish.
# With 50 marked fish, the prior pulls the answer a little. With 500, the
# data swamp it. Exact Beta answers, no JAGS needed.

prior.table <- function(M, R) {
  data.frame(marked = M, recaptured = R,
             flat.prior.mean = round((1 + R) / (2 + M), 3),
             informative.prior.mean = round((3 + R) / (20 + M), 3),
             flat.ci.width = round(diff(qbeta(c(0.025, 0.975), 1 + R, 1 + M - R)), 3),
             informative.ci.width = round(diff(qbeta(c(0.025, 0.975), 3 + R, 17 + M - R)), 3))
}
print(rbind(prior.table(50, 8), prior.table(500, 80)), row.names = FALSE)

# 6. From efficiency to abundance ----
# The point of the trial: the trap caught C unmarked fish, so passage is about
# C / p. Every draw of p gives a draw of C / p: expanded passage carrying the
# uncertainty in efficiency. Catch is held fixed here, so this is not a full
# posterior for abundance: it leaves out the chance in how many fish the trap
# happened to catch. 02_lincoln_petersen.R picks up here.

C <- 600
N.draws <- C / draws     # expanded passage, carrying efficiency uncertainty only
cat("\nUnmarked catch:", C, "\n")
cat("Expanded passage (efficiency uncertainty only), median:", round(median(N.draws)), "  95% interval:",
    round(quantile(N.draws, c(0.025, 0.975))), "\n")

# Lessons ----
# - Prior x likelihood = posterior. A flat prior lets the data speak; an
#   informative prior adds information and matters most when data are thin.
# - JAGS doesn't solve the equations; it draws samples from the posterior.
#   Summaries of the draws (mean, quantiles) are the estimates.
# - Check convergence every time: Rhat near 1 and grassy, overlapping trace
#   plots. More chains that agree is the evidence.
# - A 95% credible interval is a probability statement about the parameter,
#   given the model and data.
# - Anything computed from a parameter (like C / p) gets its own distribution:
#   compute it for every draw. That carries the parameter's uncertainty, and
#   only that; anything held fixed (like the catch here) adds none.
# - JAGS gotchas: dbin(p, n) puts the probability first; dnorm uses
#   precision, not SD.

# Where BPA teaches this ----
# - Chapter 2: Bayes rule, priors, and the posterior; 2.4 compares Bayesian
#   credible intervals with frequentist confidence intervals; 2.5 covers
#   MCMC, burn-in, Rhat, and thinning
# - 3.3.1: judging convergence from trace plots and Rhat on a first model
