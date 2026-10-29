# A logit GLM in JAGS: does spill change guidance?
# (Run 00_setup.R once first: it installs packages and checks JAGS.)
# (Model fits can take a few minutes to run; R will look busy while JAGS samples.)
#
# Each PIT-tagged fish arriving at the dam is either guided into the bypass
# (1) or not (0). That is a Bernoulli trial (BPA 3.5). A link function lets
# a linear predictor model a non-normal response (BPA 3.2.3); the logit link
# keeps the guidance probability between 0 and 1 while the linear predictor
# adds covariate effects (BPA 3.5).
#
#   guided[i] ~ Bernoulli(psi[i])
#   logit(psi[i]) = b0 + b1 * spill[i]
#
# This version pretends we see every fish's route. Real data don't allow
# that, which is why the later scripts need the multistate cells.
# All data are simulated. Requires JAGS and the jagsUI package.

library(jagsUI)
library(ggplot2)

set.seed(2)

# 1. Simulate fish ----

n.fish <- 400
spill <- runif(n.fish, 0.2, 0.6)                 # Spill share on the day each fish arrived
spill.std <- (spill - mean(spill)) / sd(spill)   # Standardize the covariate

b0 <- -1.0      # Logit guidance at average spill
b1 <- -0.8      # More spill, less guidance
psi <- plogis(b0 + b1 * spill.std)
guided <- rbinom(n.fish, 1, psi)

# 2. JAGS model ----

cat(file = "glm_logit.jags", "
model {
# Priors
# Reminder: JAGS's dnorm(mean, precision) uses precision = 1 / SD^2, not SD.
b0 ~ dnorm(0, 0.1)
b1 ~ dnorm(0, 0.1)

# Likelihood
for (i in 1:n) {
   guided[i] ~ dbern(psi[i])
   logit(psi[i]) <- b0 + b1 * spill[i]
   }
}
")

jags.data <- list(guided = guided, spill = spill.std, n = n.fish)
inits <- function() list(b0 = rnorm(1), b1 = rnorm(1))
params <- c("b0", "b1")

# Fitting: this can take a few minutes.
fit <- jags(jags.data, inits, params, "glm_logit.jags",
            n.chains = 3, n.iter = 6000, n.burnin = 2000, n.thin = 2, parallel = TRUE)

# 3. Results ----
# The full JAGS summary: posterior mean, SD, 95% interval, Rhat, and n.eff
# for each saved parameter (01_first_bayesian_model.R explains the columns).
print(fit)
# Trace plots: each chain's draws in order. Look for grassy, overlapping fuzz.
traceplot(fit, parameters = c("b0", "b1"))

cat("Max Rhat:", round(max(unlist(fit$Rhat), na.rm = TRUE), 3), "\n")
cat("b0: truth", b0, " estimate", round(fit$mean$b0, 2),
    " 95% CrI", round(quantile(fit$sims.list$b0, c(0.025, 0.975)), 2), "\n")
cat("b1: truth", b1, " estimate", round(fit$mean$b1, 2),
    " 95% CrI", round(quantile(fit$sims.list$b1, c(0.025, 0.975)), 2), "\n")

# Adding 1 on the logit scale multiplies the odds by e (about 2.7).
cat("A 1 SD increase in spill multiplies the odds of guidance by",
    round(median(exp(fit$sims.list$b1)), 2), "(posterior median)\n")

# 4. Plot the S-curve ----
# The model is linear on the logit scale and curved on the probability scale.
x <- seq(0.2, 0.6, length.out = 100)
x.std <- (x - mean(spill)) / sd(spill)
pred <- sapply(x.std, function(z) plogis(fit$sims.list$b0 + fit$sims.list$b1 * z))
curve.df <- data.frame(spill = x, mean = colMeans(pred),
                       lo = apply(pred, 2, quantile, 0.025), hi = apply(pred, 2, quantile, 0.975),
                       truth = plogis(b0 + b1 * x.std))

p1 <- ggplot(curve.df, aes(spill)) +
  geom_jitter(data = data.frame(spill = spill, guided = guided), aes(y = guided),
              width = 0, height = 0.03, size = 0.8, color = "gray60") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = mean, linetype = "Estimated, 95% CrI"), color = "steelblue", linewidth = 1) +
  geom_line(aes(y = truth, linetype = "True guidance"), linewidth = 1) +
  scale_linetype_manual(NULL, values = c("Estimated, 95% CrI" = "solid", "True guidance" = "dashed")) +
  labs(x = "Spill share", y = "Guided (1) or not (0)") +
  theme_bw(base_size = 13) + theme(legend.position = "top")
print(p1)

# Try this ----
# Predict before rerunning; reset each change before trying the next.
# 1. Set n.fish <- 50. Does the interval for b1 still exclude zero?
# 2. Set b1 <- 0. What should happen to the fitted curve and its uncertainty?

# Lessons ----
# - Each fish is a Bernoulli trial: guided or not. That is the noise part of
#   the GLM.
# - The logit link lets covariates add on an unbounded scale while the
#   probability stays between 0 and 1.
# - The model is linear on the logit scale and curved on the probability
#   scale. Adding 1 on the logit scale multiplies the odds by about 2.7.
# - Simulating from known values and fitting shows whether the model can get
#   the answer back. Every script in this repo uses that pattern.

# Where BPA teaches this ----
# - 3.2  Statistical models: Response = Signal + Noise (3.2.1 noise, 3.2.2 signal)
# - 3.5  Binomial GLM for bounded counts or proportions; 3.5.1 simulates
#        data and fits it, the same pattern as this script
