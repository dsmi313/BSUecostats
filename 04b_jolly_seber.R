# Survival, recruitment, and population size: the Jolly-Seber (JS) model
# (Run 00_setup.R once first: it installs packages and checks JAGS.)
# (Run time: about 3 to 8 minutes. Times vary up to 2x between computers.)
#
# Straight from Kéry and Schaub (2012), chapter 10. The simulation function
# and JAGS model are the book's, from the authors' JAGS translation
# (M. Schaub; jagsUI version by M. Kéry, 2025). Changes: a seed is set, fewer
# augmented rows are used so it runs faster, and the repo's diagnostic checks
# and comments are added.
#
# Lesson 02a assumed a closed population: no births, deaths, or movement.
# Lesson 04a let animals die but conditioned on first capture, so it could
# not estimate N. Jolly-Seber does both: animals enter the population
# (births or immigration), survive or not, and are caught or not. It does not
# condition on first capture, so the first capture also carries information
# about how many animals were never caught.
#
# The book fits JS three ways (BPA 10.3). This lesson uses the first: a
# restricted dynamic occupancy model on augmented data, as in 02a. Each row
# is a potential animal; z[i,t] = 1 if it is alive and in the population at
# occasion t. gamma[t] is the chance a not-yet-entered row enters at t.
# All data are simulated. Requires JAGS and the jagsUI package.

library(jagsUI)
source("helpers.R")
library(ggplot2)

set.seed(1)

# 1. Simulate an open population (BPA 10.4) ----

n.occasions <- 7                         # Number of capture occasions
N <- 400                                 # Superpopulation size
phi <- rep(0.7, n.occasions-1)           # Survival probabilities
b <- c(0.34, rep(0.11, n.occasions-1))   # Entry probabilities
p <- rep(0.5, n.occasions)               # Capture probabilities

PHI <- matrix(rep(phi, N), ncol = n.occasions-1, nrow = N, byrow = T)
P <- matrix(rep(p, N), ncol = n.occasions, nrow = N, byrow = T)

# Function to simulate capture-recapture data under the JS model
simul.js <- function(PHI, P, b, N){
   B <- rmultinom(1, N, b) # Generate no. of entering ind. per occasion
   n.occasions <- dim(PHI)[2] + 1
   CH.sur <- CH.p <- matrix(0, ncol = n.occasions, nrow = N)
   # Define a vector with the occasion of entering the population
   ent.occ <- numeric()
   for (t in 1:n.occasions){
      ent.occ <- c(ent.occ, rep(t, B[t]))
      }
   # Simulating survival
   for (i in 1:N){
      CH.sur[i, ent.occ[i]] <- 1   # Write 1 when ind. enters the pop.
      if (ent.occ[i] == n.occasions) next
      for (t in (ent.occ[i]+1):n.occasions){
         # Bernoulli trial: has individual survived occasion?
         sur <- rbinom(1, 1, PHI[i,t-1])
         ifelse (sur==1, CH.sur[i,t] <- 1, break)
         } #t
      } #i
   # Simulating capture
   for (i in 1:N){
      CH.p[i,] <- rbinom(n.occasions, 1, P[i,])
      } #i
   # Full capture-recapture matrix
   CH <- CH.sur * CH.p

   # Remove individuals never captured
   cap.sum <- rowSums(CH)
   never <- which(cap.sum == 0)
   CH <- CH[-never,]
   Nt <- colSums(CH.sur)    # Actual population size
   return(list(CH=CH, B=B, N=Nt))
   }

sim <- simul.js(PHI, P, b, N)
CH <- sim$CH
cat(nrow(CH), "of", N, "animals in the superpopulation were caught at least once.\n")

# 2. Fit the JS model as a restricted occupancy model (BPA 10.3.1, 10.4.1) ----

cat(file = "js-rest.occ.jags", "
model {
# Priors and constraints
for (i in 1:M){
   for (t in 1:(n.occasions-1)){
      phi[i,t] <- mean.phi
      } #t
   for (t in 1:n.occasions){
      p[i,t] <- mean.p
      } #t
   } #i
mean.phi ~ dunif(0, 1)
mean.p ~ dunif(0, 1)

for (t in 1:n.occasions){
   gamma[t] ~ dunif(0, 1)
   } #t

# Likelihood
for (i in 1:M){
   # First occasion
   # State process
   z[i,1] ~ dbern(gamma[1])
   # Observation process
   mu1[i] <- z[i,1] * p[i,1]
   y[i,1] ~ dbern(mu1[i])
   # Subsequent occasions
   for (t in 2:n.occasions){
      # State process
      q[i,t-1] <- 1-z[i,t-1]      # Availability for recruitment
      mu2[i,t] <- phi[i,t-1] * z[i,t-1] + gamma[t] * prod(q[i,1:(t-1)])
      z[i,t] ~ dbern(mu2[i,t])
      # Observation process
      mu3[i,t] <- z[i,t] * p[i,t]
      y[i,t] ~ dbern(mu3[i,t])
      } #t
   } #i

# Calculate derived population parameters
for (t in 1:n.occasions){
   qgamma[t] <- 1-gamma[t]
   }
cprob[1] <- gamma[1]
for (t in 2:n.occasions){
   cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
   } #t
psi <- sum(cprob[])            # Inclusion probability
for (t in 1:n.occasions){
   b[t] <- cprob[t] / psi      # Entry probability
   } #t
for (i in 1:M){
   recruit[i,1] <- z[i,1]
   for (t in 2:n.occasions){
      recruit[i,t] <- (1-z[i,t-1]) * z[i,t]
      } #t
   } #i
for (t in 1:n.occasions){
   N[t] <- sum(z[1:M,t])        # Actual population size
   B[t] <- sum(recruit[1:M,t])  # Number of entries
   } #t
for (i in 1:M){
   Nind[i] <- sum(z[i,1:n.occasions])
   Nalive[i] <- 1-equals(Nind[i], 0)
   } #i
Nsuper <- sum(Nalive[])         # Superpopulation size
}
")

# Augment the capture histories. The book adds 500 rows; 300 is enough here
# and runs faster. Check that Nsuper's posterior stays below the ceiling M.
nz <- 300
CH.aug <- rbind(CH, matrix(0, ncol = dim(CH)[2], nrow = nz))

jags.data <- list(y = CH.aug, n.occasions = dim(CH.aug)[2], M = dim(CH.aug)[1])

# Good initial values for z are needed in JAGS. The simplest that works is a
# 1 everywhere.
z.init <- CH.aug
z.init[z.init==0] <- 1
inits <- function(){list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), z = z.init)}

parameters <- c("mean.p", "mean.phi", "b", "Nsuper", "N", "B")
# The book runs 5,000 iterations. Three times longer passes the diagnostic
# screen. Expect several minutes.
ni <- 15000; nt <- 10; nb <- 5000; nc <- 3

js.occ <- jags(jags.data, inits, parameters, "js-rest.occ.jags",
               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
print(js.occ, digits = 3)
if (!fit_ok(js.occ)) warning("JS fit: diagnostics need attention before interpreting results.")

# 3. Results ----

q95 <- function(x) round(quantile(x, c(0.025, 0.975)))
cat("Superpopulation: truth", N, " posterior median", median(js.occ$sims.list$Nsuper),
    " 95% CrI", q95(js.occ$sims.list$Nsuper),
    " share of draws at the ceiling", round(mean(js.occ$sims.list$Nsuper == nrow(CH.aug)), 3), "\n")
cat("Survival: truth", phi[1], " estimate", round(js.occ$mean$mean.phi, 3), " 95% CrI",
    round(quantile(js.occ$sims.list$mean.phi, c(0.025, 0.975)), 3), "\n")
cat("Capture: truth", p[1], " estimate", round(js.occ$mean$mean.p, 3), " 95% CrI",
    round(quantile(js.occ$sims.list$mean.p, c(0.025, 0.975)), 3), "\n")
print(data.frame(occasion = 1:n.occasions,
                 N.true = sim$N, N.est = round(js.occ$mean$N),
                 N.lo = js.occ$q2.5$N, N.hi = js.occ$q97.5$N,
                 B.true = as.vector(sim$B), B.est = round(js.occ$mean$B)))

n.df <- data.frame(occasion = 1:n.occasions, truth = sim$N, mean = js.occ$mean$N,
                   lo = js.occ$q2.5$N, hi = js.occ$q97.5$N, caught = colSums(CH))
p1 <- ggplot(n.df, aes(occasion)) +
  geom_pointrange(aes(y = mean, ymin = lo, ymax = hi, color = "Estimate, 95% CrI")) +
  geom_point(aes(y = truth, color = "Truth"), shape = 4, size = 3, stroke = 1.2) +
  geom_point(aes(y = caught, color = "Animals caught"), shape = 1, size = 2.5) +
  scale_color_manual(NULL, values = c("Estimate, 95% CrI" = "steelblue", "Truth" = "red",
                                      "Animals caught" = "gray40")) +
  labs(x = "Occasion", y = "Population size N[t]") +
  theme_bw(base_size = 13) + theme(legend.position = "top")
print(p1)

# Try this ----
# Predict before rerunning; reset each change before trying the next.
# 1. Set p <- rep(0.2, n.occasions). How many animals are caught, and how much
#    wider do the intervals for N[t] and Nsuper get?
# 2. Set nz <- 30. Where does Nsuper's posterior end up, and what does the
#    ceiling share tell you?

# Lessons ----
# - Jolly-Seber adds entries to CJS: animals are born or arrive, survive, and
#   are caught or not. It estimates population size at each occasion, the
#   number of entries, and the superpopulation that was ever present.
# - It doesn't condition on first capture, so first captures inform N too.
# - Data augmentation works the same way as in 02a: add all-zero rows and
#   check that the posterior stays below the ceiling.
# - N[t] and B[t] are derived from the latent states z, so their uncertainty
#   comes along for free.

# Where BPA teaches this ----
# - 10.1  Introduction
# - 10.2  The JS model as a state-space model
# - 10.3  Fitting JS with data augmentation: as a restricted dynamic occupancy
#         model (10.3.1, used here), as a multistate model (10.3.2), and under
#         the superpopulation parameterization (10.3.3)
# - 10.4  Constant survival and time-dependent entry, the example here,
#         fitted all three ways and compared (10.4.4)
# - 10.5  Individual capture heterogeneity
# - 10.7  Real data: survival, recruitment, and population size of Leisler's bats
