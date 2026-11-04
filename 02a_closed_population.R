# Closed-population size from capture histories: M0, Mt, Mb, and Mh
# (Run 00_setup.R once first: it installs packages and checks JAGS.)
# (Model fits can take a few minutes to run; R will look busy while JAGS samples.)
#
# Straight from Kéry and Schaub (2012), chapter 6. The simulation functions
# and JAGS models are the book's, from the authors' JAGS translation
# (M. Schaub; jagsUI version by M. Kéry, 2025). Changes: functions are renamed
# so sections don't overwrite each other, seeds are set, and the repo's
# diagnostic checks and comments are added.
#
# Lesson 02 used two occasions: mark M fish, then count R recaptures in a
# second sample. Here every animal caught gets a capture history over T
# occasions, e.g. 0 1 1 means missed on occasion 1 and caught on 2 and 3.
# Animals never caught (0 0 0) are missing from the data. How many of them
# there are is the whole question.
#
# Data augmentation (BPA 6.2.1) turns that into an occupancy-like problem: add
# nz all-zero rows, give every row an indicator z (1 = a real member of the
# population), and estimate N = sum(z). The models differ only in what drives
# capture probability p:
#   M0  one p for everyone on every occasion
#   Mt  p changes by occasion. Lincoln-Petersen is Mt with two occasions.
#   Mb  p changes after capture (trap-happy or trap-shy fish)
#   Mh  p differs among individuals (bigger fish are easier to catch, say)
# All data are simulated. Requires JAGS and the jagsUI package.

library(jagsUI)
source("helpers.R")
library(ggplot2)

set.seed(6)

# 1. Model M0 with data augmentation (BPA 6.2.1) ----

sim.M0 <- function(N = 100, p = 0.5, T = 3){
   yfull <- yobs <- array(NA, dim = c(N, T))
   for (j in 1:T){
      yfull[,j] <- rbinom(n = N, size = 1, prob = p)
      }
   ever.detected <- apply(yfull, 1, max)
   C <- sum(ever.detected)
   yobs <- yfull[ever.detected == 1,]
   cat(C, "out of", N, "animals present were detected.\n")
   return(list(N = N, p = p, C = C, T = T, yfull = yfull, yobs = yobs))
   }

data <- sim.M0()

# Augment the data set by 150 potential individuals that were never caught.
nz <- 150
yaug <- rbind(data$yobs, array(0, dim = c(nz, data$T)))

cat(file = "M0.jags", "
model {

# Priors
omega ~ dunif(0, 1)
p ~ dunif(0, 1)

# Likelihood
for (i in 1:M){
   z[i] ~ dbern(omega)          # Inclusion indicators
   for (j in 1:T){
      yaug[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i] * p    # Can only be detected if z=1
      } #j
   } #i

# Derived quantities
N <- sum(z[])
}
")

ni <- 2500; nt <- 2; nb <- 500; nc <- 3

# The same M0 fit is reused below with other augmented data sets.
fit.M0 <- function(yaug) {
  inits <- function() list(z = rep(1, nrow(yaug)), p = runif(1, 0, 1))
  jags(list(yaug = yaug, M = nrow(yaug), T = ncol(yaug)), inits, c("N", "p", "omega"), "M0.jags",
       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
}
out150 <- fit.M0(yaug)
print(out150, dig = 3)
if (!fit_ok(out150)) warning("M0 fit: diagnostics need attention before interpreting results.")

# 2. How much to augment (BPA 6.2.1) ----
# M = C + nz is a ceiling: N can't exceed it. With nz = 5 the posterior piles
# up against the ceiling. The diagnostic screen flags this fit only because
# N's tail ESS can't be computed when draws pile up at one value. Don't count
# on that: if N's posterior touches M, augment more and refit.

yaug5 <- rbind(data$yobs, array(0, dim = c(5, data$T)))
out5 <- fit.M0(yaug5)
if (!fit_ok(out5)) warning("M0 fit, nz = 5: fails the diagnostic screen, as expected when N piles up at the ceiling.")
for (fit in list(list(out5, nrow(yaug5)), list(out150, nrow(yaug)))) {
  cat("Augmented rows M =", fit[[2]], ": posterior median N", median(fit[[1]]$sims.list$N),
      " 95% CrI", quantile(fit[[1]]$sims.list$N, c(0.025, 0.975)),
      " share of draws at the ceiling", round(mean(fit[[1]]$sims.list$N == fit[[2]]), 3), "\n")
}
cat("Truth N =", data$N, "; animals caught C =", data$C, "\n")

aug.df <- data.frame(N = c(out5$sims.list$N, out150$sims.list$N),
                     fit = rep(c(paste0("Augmented by 5 (M = ", nrow(yaug5), ")"),
                                 paste0("Augmented by 150 (M = ", nrow(yaug), ")")),
                               c(length(out5$sims.list$N), length(out150$sims.list$N))))
p1 <- ggplot(aug.df, aes(N)) +
  geom_histogram(binwidth = 1, fill = "gray60") +
  geom_vline(xintercept = data$N, color = "red", linewidth = 1) +
  geom_vline(xintercept = data$C, linetype = "dashed") +
  facet_wrap(~ fit, ncol = 1) +
  labs(x = "Population size N (red: truth; dashed: animals caught)", y = "Posterior draws") +
  theme_bw(base_size = 13)
print(p1)

# 3. Model Mt: time effects (BPA 6.2.2) ----

sim.Mt <- function(N = 100, mean.p = 0.5, T = 3, time.eff = runif(T, -2, 2)){
   yfull <- yobs <- array(NA, dim = c(N, T) )
   p.vec <- array(NA, dim = T)
   for (j in 1:T){
      p <- plogis(log(mean.p / (1-mean.p)) + time.eff[j])
      yfull[,j] <- rbinom(n = N, size = 1, prob = p)
      p.vec[j] <- p
      }
   ever.detected <- apply(yfull, 1, max)
   C <- sum(ever.detected)
   yobs <- yfull[ever.detected == 1,]
   cat(C, "out of", N, "animals present were detected.\n")
   return(list(N = N, p.vec = p.vec, C = C, T = T, yfull = yfull, yobs = yobs))
   }

data.t <- sim.Mt(T = 3)
nz <- 150
yaug.t <- rbind(data.t$yobs, array(0, dim = c(nz, data.t$T)))

cat(file = "Mt.jags", "
model {
# Priors
omega ~ dunif(0, 1)
for (i in 1:T){
   p[i] ~ dunif(0, 1)
   }

# Likelihood
for (i in 1:M){
   z[i] ~ dbern(omega)
   for (j in 1:T){
      yaug[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i] * p[j]
      } #j
   } #i

# Derived quantities
N <- sum(z[])
} # end model
")

out.t <- jags(list(yaug = yaug.t, M = nrow(yaug.t), T = ncol(yaug.t)),
              function() list(z = rep(1, nrow(yaug.t)), p = runif(data.t$T, 0, 1)),
              c("N", "p", "omega"), "Mt.jags",
              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
print(out.t, dig = 3)
if (!fit_ok(out.t)) warning("Mt fit: diagnostics need attention before interpreting results.")

# Lincoln-Petersen from the first two occasions only, for comparison.
n1 <- sum(data.t$yobs[, 1]); n2 <- sum(data.t$yobs[, 2])
m2 <- sum(data.t$yobs[, 1] * data.t$yobs[, 2])
cat("Mt: truth N", data.t$N, " posterior median", median(out.t$sims.list$N),
    " 95% CrI", quantile(out.t$sims.list$N, c(0.025, 0.975)), "\n")
cat("Lincoln-Petersen from occasions 1 and 2 (n1 * n2 / m2):", round(n1 * n2 / m2), "\n")
cat("True p by occasion:", round(data.t$p.vec, 2),
    " estimated:", round(out.t$mean$p, 2), "\n")

# 4. Model Mb: behavioral (trap) response (BPA 6.2.3) ----
# Capture probability is c if the animal was caught on the previous occasion
# and p if not. This is an immediate trap response: the change lasts one
# occasion. The book also describes a permanent response, where p changes for
# good after first capture. Fish handled at a trap can become trap-shy; baited
# traps make animals trap-happy.

sim.Mb <- function(N = 200, T = 5, p = 0.3, c = 0.4){
   yfull <- yobs <- array(NA, dim = c(N, T) )
   p.eff <- array(NA, dim = N)

   # First capture occasion
   yfull[,1] <- rbinom(n = N, size = 1, prob = p)

   # Later capture occasions
   for (j in 2:T){
      p.eff <- (1 - yfull[,(j-1)]) * p + yfull[,(j-1)] * c
      yfull[,j] <- rbinom(n = N, size = 1, prob = p.eff)
      }

   ever.detected <- apply(yfull, 1, max)
   C <- sum(ever.detected)
   yobs <- yfull[ever.detected == 1,]
   cat(C, "out of", N, "animals present were detected.\n")
   return(list(N = N, p = p, c = c, C = C, T = T, yfull = yfull, yobs = yobs))
   }

data.b <- sim.Mb(N = 200)
nz <- 150
yaug.b <- rbind(data.b$yobs, array(0, dim = c(nz, data.b$T)))

cat(file = "Mb.jags", "
model {
# Priors
omega ~ dunif(0, 1)
p ~ dunif(0, 1)     # Cap prob when not caught at t-1
c ~ dunif(0, 1)     # Cap prob when caught at t-1

# Likelihood
for (i in 1:M){
   z[i] ~ dbern(omega)

   # First occasion
   yaug[i,1] ~ dbern(p.eff[i,1])
   p.eff[i,1] <- z[i] * p

   # All subsequent occasions
   for (j in 2:T){
      yaug[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i] * ( (1-yaug[i,(j-1)]) * p + yaug[i,(j-1)] * c )
      } #j
   } #i

# Derived quantities
N <- sum(z[])
trap.response <- c - p
} # end model
")

out.b <- jags(list(yaug = yaug.b, M = nrow(yaug.b), T = ncol(yaug.b)),
              function() list(z = rep(1, nrow(yaug.b)), p = runif(1, 0, 1)),
              c("N", "p", "c", "trap.response", "omega"), "Mb.jags",
              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
print(out.b, dig = 3)
if (!fit_ok(out.b)) warning("Mb fit: diagnostics need attention before interpreting results.")
cat("Mb: truth N", data.b$N, " posterior median", median(out.b$sims.list$N),
    " 95% CrI", quantile(out.b$sims.list$N, c(0.025, 0.975)), "\n")
cat("Trap response c - p: truth", data.b$c - data.b$p, " estimate",
    round(out.b$mean$trap.response, 2), " 95% CrI",
    round(quantile(out.b$sims.list$trap.response, c(0.025, 0.975)), 2), "\n")

# 5. Model Mh: individual heterogeneity (BPA 6.2.4) ----
# Each animal gets its own p, drawn around mean.p on the logit scale. Animals
# with low p are the ones most likely to be missed entirely, so a model that
# assumes one p (M0) underestimates N. Lincoln-Petersen carries the same
# equal-catchability assumption.

sim.Mh <- function(N = 100, mean.p = 0.4, T = 5, sd = 1){
   yfull <- yobs <- array(NA, dim = c(N, T) )
   mean.lp <- log(mean.p / (1-mean.p))
   p.vec <- plogis(mean.lp+ rnorm(N, 0, sd))

   for (i in 1:N){
      yfull[i,] <- rbinom(n = T, size = 1, prob = p.vec[i])
      }

   ever.detected <- apply(yfull, 1, max)
   C <- sum(ever.detected)
   yobs <- yfull[ever.detected == 1,]
   cat(C, "out of", N, "animals present were detected.\n")
   return(list(N = N, p.vec = p.vec, mean.lp = mean.lp, C = C, T = T, yfull = yfull, yobs = yobs))
   }

data.h <- sim.Mh(sd = 1)

# Mh needs only each animal's number of captures out of T.
y <- sort(apply(data.h$yobs, 1, sum), decreasing = TRUE)
nz <- 300
yaug.h <- c(y, rep(0, nz))

cat(file = "Mh.jags", "
model {

# Priors
omega ~ dunif(0, 1)
mean.lp <- logit(mean.p)
mean.p ~ dunif(0, 1)
tau <- 1 / (sd * sd)
sd ~ dunif(0, 5)

# Likelihood
for (i in 1:M){
   z[i] ~ dbern(omega)
   logit(p[i]) <- eps[i]
   eps[i] ~ dnorm(mean.lp, tau)T(-16, 16)   # See web appendix A in Royle (2009)
   p.eff[i] <- z[i] * p[i]
   y[i] ~ dbin(p.eff[i], T)
   }

# Derived quantities
N <- sum(z[])
}
")

# The book runs 25,000 iterations. Mh mixes slowly, so this runs five times
# longer to pass the repo's diagnostic screen. Expect a couple of minutes.
ni.h <- 125000; nt.h <- 10; nb.h <- 25000
out.h <- jags(list(y = yaug.h, M = length(yaug.h), T = ncol(data.h$yobs)),
              function() list(z = rep(1, length(yaug.h)), sd = runif(1, 0.1, 0.9)),
              c("N", "mean.p", "sd", "omega"), "Mh.jags",
              n.chains = nc, n.thin = nt.h, n.iter = ni.h, n.burnin = nb.h, parallel = TRUE)
print(out.h, dig = 3)
if (!fit_ok(out.h)) warning("Mh fit: diagnostics need attention before interpreting results.")

# Now fit M0 to the same heterogeneous data.
yaug.h0 <- rbind(data.h$yobs, array(0, dim = c(nz, data.h$T)))
out.h0 <- fit.M0(yaug.h0)
if (!fit_ok(out.h0)) warning("M0 fit to Mh data: diagnostics need attention before interpreting results.")

cat("\nHeterogeneous catchability, truth N =", data.h$N, "; caught C =", data.h$C, "\n")
cat("  M0 (one p): posterior median", median(out.h0$sims.list$N),
    " 95% CrI", quantile(out.h0$sims.list$N, c(0.025, 0.975)), "\n")
cat("  Mh (own p): posterior median", median(out.h$sims.list$N),
    " 95% CrI", quantile(out.h$sims.list$N, c(0.025, 0.975)),
    " share of draws at the ceiling", round(mean(out.h$sims.list$N == length(yaug.h)), 3), "\n")

h.df <- data.frame(N = c(out.h0$sims.list$N, out.h$sims.list$N),
                   model = rep(c("M0: one p for every animal", "Mh: each animal has its own p"),
                               c(length(out.h0$sims.list$N), length(out.h$sims.list$N))))
p2 <- ggplot(h.df, aes(N)) +
  geom_histogram(binwidth = 2, fill = "gray60") +
  geom_vline(xintercept = data.h$N, color = "red", linewidth = 1) +
  geom_vline(xintercept = data.h$C, linetype = "dashed") +
  facet_wrap(~ model, ncol = 1) +
  labs(x = "Population size N (red: truth; dashed: animals caught)", y = "Posterior draws") +
  theme_bw(base_size = 13)
print(p2)

# Try this ----
# Predict before rerunning; reset each change before trying the next.
# 1. In section 5, change sim.Mh(sd = 1) to sim.Mh(sd = 2). How far low does M0
#    land now, and what happens to the width of the Mh interval?
# 2. In section 3, change sim.Mt(T = 3) to sim.Mt(T = 2). How close is the
#    posterior median of N to the Lincoln-Petersen estimate?
# 3. In section 4, change sim.Mb(N = 200) to sim.Mb(N = 200, c = 0.3), so there
#    is no trap response. Does the interval for c - p include zero?

# Lessons ----
# - Animals never caught are invisible. Data augmentation adds all-zero rows
#   and lets the model decide how many of them are real.
# - Augment enough: if N's posterior reaches the ceiling M, the answer is
#   wrong even when diagnostics pass.
# - The four models differ only in what changes capture probability: nothing,
#   occasion, previous capture, or the individual.
# - Heterogeneous catchability biases N low when the model ignores it. The
#   hardest-to-catch animals are the ones missing from the data.
# - Lincoln-Petersen (lesson 02) is model Mt with two occasions. Its
#   equal-catchability assumption is the one Mh relaxes.

# Where BPA teaches this ----
# - 6.1  Introduction
# - 6.2.1  Data augmentation and model M0, including how much to augment
# - 6.2.2  Time effects, model Mt
# - 6.2.3  Behavioral or memory effects, model Mb
# - 6.2.4  Individual (random) effects, model Mh
# - 6.2.5  Combined effects, model Mth; 6.3 adds a behavioral term (Mtbh) for
#          real data on species richness
# - 6.4  Individual covariates, model Mt+X: p depends on a measured trait
#        such as body size
