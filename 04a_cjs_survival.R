# Survival from recaptures: the Cormack-Jolly-Seber (CJS) model
# (Run 00_setup.R once first: it installs packages and checks JAGS.)
# (Run time: about 30 to 60 seconds. Times vary up to 2x between computers.)
#
# Straight from Kéry and Schaub (2012), chapter 7. The simulation functions
# and JAGS models are the book's, from the authors' JAGS translation
# (M. Schaub; jagsUI version by M. Kéry, 2025). Changes: seeds are set, one
# model is written as a constant-parameter version of another, and the repo's
# diagnostic checks and comments are added.
#
# Lessons 02 and 02a estimate how many animals there are. CJS asks a different
# question: of the animals marked, how many survive from one occasion to the
# next? It conditions on each animal's first capture, so it never estimates
# N. It estimates apparent survival phi (survived and stayed in the study
# area) and recapture probability p.
#
# It is a state-space model, like lesson 04, but for each individual:
#   state:        z[i,t] ~ Bernoulli(phi * z[i,t-1])    alive (1) or not (0)
#   observation:  y[i,t] ~ Bernoulli(p * z[i,t])        seen only if alive
#
# Reach survival of PIT-tagged juvenile salmon through the Snake River dams is
# commonly estimated this way, with dams as the recapture occasions. Lesson 05
# extends the idea: a fish's route at Lower Granite is the hidden state.
# All data are simulated. Requires JAGS and the jagsUI package.

library(jagsUI)
source("helpers.R")
library(ggplot2)

set.seed(1)

# 1. Constant survival and recapture (BPA 7.3) ----

n.occasions <- 6                   # Number of capture occasions
marked <- rep(50, n.occasions-1)   # Annual number of newly marked individuals
phi <- rep(0.65, n.occasions-1)
p <- rep(0.4, n.occasions-1)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked))
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Define function to simulate a capture-history (CH) matrix
simul.cjs <- function(PHI, P, marked){
   n.occasions <- dim(PHI)[2] + 1
   CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))
   # Define a vector with the occasion of marking
   mark.occ <- rep(1:length(marked), marked[1:length(marked)])
   # Fill the CH matrix
   for (i in 1:sum(marked)){
      CH[i, mark.occ[i]] <- 1       # Write an 1 at the release occasion
      if (mark.occ[i]==n.occasions) next
      for (t in (mark.occ[i]+1):n.occasions){
         # Bernoulli trial: does individual survive occasion?
         sur <- rbinom(1, 1, PHI[i,t-1])
         if (sur==0) break          # If dead, move to next individual
         # Bernoulli trial: is individual recaptured?
         rp <- rbinom(1, 1, P[i,t-1])
         if (rp==1) CH[i,t] <- 1
         } #t
      } #i
   return(CH)
   }

CH <- simul.cjs(PHI, P, marked)
print(head(CH))

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

cat(file = "cjs-c-c.jags", "
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- mean.phi
      p[i,t] <- mean.p
      } #t
   } #i

mean.phi ~ dunif(0, 1)         # Prior for mean survival
mean.p ~ dunif(0, 1)           # Prior for mean recapture

# Likelihood
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}
")

# An animal seen at occasions 2 and 5 was alive at 3 and 4 too. Passing those
# known states as data speeds up the fit (BPA 7.3.1).
known.state.cjs <- function(ch){
   state <- ch
   for (i in 1:dim(ch)[1]){
      n1 <- min(which(ch[i,]==1))
      n2 <- max(which(ch[i,]==1))
      state[i,n1:n2] <- 1
      state[i,n1] <- NA
      }
   state[state==0] <- NA
   return(state)
   }

# Initial values for the unknown z: those not supplied as data
cjs.init.z <- function(ch,f){
   for (i in 1:dim(ch)[1]){
      if (sum(ch[i,])==1) next
      n2 <- max(which(ch[i,]==1))
      ch[i,f[i]:n2] <- NA
      }
   for (i in 1:dim(ch)[1]){
   ch[i,1:f[i]] <- NA
   }
   return(ch)
   }

jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH))
inits <- function(){list(z = cjs.init.z(CH, f), mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1))}
ni <- 10000; nt <- 6; nb <- 5000; nc <- 3

cjs.c.c <- jags(jags.data, inits, c("mean.phi", "mean.p"), "cjs-c-c.jags",
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
print(cjs.c.c, digits = 3)
if (!fit_ok(cjs.c.c)) warning("CJS constant fit: diagnostics need attention before interpreting results.")
cat("phi: truth", phi[1], " estimate", round(cjs.c.c$mean$mean.phi, 3), " 95% CrI",
    round(quantile(cjs.c.c$sims.list$mean.phi, c(0.025, 0.975)), 3), "\n")
cat("p:   truth", p[1], " estimate", round(cjs.c.c$mean$mean.p, 3), " 95% CrI",
    round(quantile(cjs.c.c$sims.list$mean.p, c(0.025, 0.975)), 3), "\n")

# 2. Time-varying survival, the m-array, and identifiability (BPA 7.9, 7.10) ----
# Survival and recapture now change every year. Instead of one row per
# animal, the m-array counts how many animals released at each occasion were
# next recaptured at each later occasion, plus how many were never seen
# again. Each row is one multinomial draw (dmulti), as in lesson 05.

n.occasions <- 12                  # Number of capture occasions
marked <- rep(30, n.occasions-1)   # Annual number of newly marked individuals
phi <- c(0.6, 0.5, 0.55, 0.6, 0.5, 0.4, 0.6, 0.5, 0.55, 0.6, 0.7)
p <- c(0.4, 0.65, 0.4, 0.45, 0.55, 0.68, 0.66, 0.28, 0.55, 0.45, 0.35)

PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
CH <- simul.cjs(PHI, P, marked)

# Function to create a m-array based on capture-histories (CH)
marray <- function(CH){
   nind <- dim(CH)[1]
   n.occasions <- dim(CH)[2]
   m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
   # Calculate the number of released individuals at each time period
   for (t in 1:n.occasions){
      m.array[t,1] <- sum(CH[,t])
      }
   for (i in 1:nind){
      pos <- which(CH[i,]!=0)
      g <- length(pos)
      for (z in 1:(g-1)){
         m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
         } #z
      } #i
   # Calculate the number of individuals that is never recaptured
   for (t in 1:n.occasions){
      m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
      }
   out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
   return(out)
   }

marr <- marray(CH)
print(marr)   # Rows: release occasion. Columns: next recapture, then never seen.

cat(file = "cjs-mnl.jags", "
model {
# Priors and constraints
for (t in 1:(n.occasions-1)){
   phi[t] ~ dunif(0, 1)         # Priors for survival
   p[t] ~ dunif(0, 1)           # Priors for recapture
   }
# Define the multinomial likelihood
for (t in 1:(n.occasions-1)){
   marr[t,1:n.occasions] ~ dmulti(pr[t, ], rel[t])
   }
# Define the cell probabilities of the m-array
# Main diagonal
for (t in 1:(n.occasions-1)){
   q[t] <- 1-p[t]                # Probability of non-recapture
   pr[t,t] <- phi[t]*p[t]
   # Above main diagonal
   for (j in (t+1):(n.occasions-1)){
      pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t,j] <- 0
      } #j
   } #t
# Last column: probability of non-recapture
for (t in 1:(n.occasions-1)){
   pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
   } #t

# Assess model fit using Freeman-Tukey statistic
# Compute fit statistics for observed data
for (t in 1:(n.occasions-1)){
   for (j in 1:n.occasions){
      expmarr[t,j] <- rel[t]*pr[t,j]
      E.org[t,j] <- pow((pow(marr[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      } #j
   } #t
# Generate replicate data and compute fit stats from them
for (t in 1:(n.occasions-1)){
   marr.new[t,1:n.occasions] ~ dmulti(pr[t, ], rel[t])
   for (j in 1:n.occasions){
      E.new[t,j] <- pow((pow(marr.new[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      } #j
   } #t
fit <- sum(E.org[,])
fit.new <- sum(E.new[,])
}
")

jags.data <- list(marr = marr, n.occasions = dim(marr)[2], rel = rowSums(marr))
# JAGS draws starting values for p from its prior.
inits <- function(){list(phi = runif(dim(marr)[2]-1, 0, 1))}
ni <- 10000; nt <- 3; nb <- 5000; nc <- 3

cjs <- jags(jags.data, inits, c("phi", "p", "fit", "fit.new"), "cjs-mnl.jags",
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
print(cjs, digits = 3)

if (!fit_ok(cjs)) warning("CJS m-array fit: diagnostics need attention before interpreting results.")

# The last survival and the last recapture probability aren't separately
# identified: an animal missing at occasion 12 either died or survived and
# was missed, and no later occasion can tell which. Only their product is
# identified. Diagnostics pass anyway. The symptom is a posterior nearly as
# wide as the flat prior (BPA 7.9).

phi.last <- cjs$sims.list$phi[, n.occasions - 1]
p.last <- cjs$sims.list$p[, n.occasions - 1]
cat("Last phi: truth", phi[n.occasions - 1], " 95% CrI", round(quantile(phi.last, c(0.025, 0.975)), 2),
    " (flat prior's 95% interval: 0.03-0.98)\n")
cat("Last p:   truth", p[n.occasions - 1], " 95% CrI", round(quantile(p.last, c(0.025, 0.975)), 2), "\n")
cat("Their product: truth", phi[n.occasions - 1] * p[n.occasions - 1], " estimate",
    round(median(phi.last * p.last), 3), " 95% CrI",
    round(quantile(phi.last * p.last, c(0.025, 0.975)), 3), "\n")

est.df <- rbind(
  data.frame(parameter = "Survival phi", t = seq_along(phi), truth = phi,
             mean = cjs$mean$phi, lo = cjs$q2.5$phi, hi = cjs$q97.5$phi),
  data.frame(parameter = "Recapture p", t = seq_along(p) + 1, truth = p,
             mean = cjs$mean$p, lo = cjs$q2.5$p, hi = cjs$q97.5$p))
p1 <- ggplot(est.df, aes(t)) +
  geom_pointrange(aes(y = mean, ymin = lo, ymax = hi, color = "Estimate, 95% CrI")) +
  geom_point(aes(y = truth, color = "Truth"), shape = 4, size = 3, stroke = 1.2) +
  scale_color_manual(NULL, values = c("Estimate, 95% CrI" = "steelblue", "Truth" = "red")) +
  facet_wrap(~ parameter, ncol = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Occasion (survival: interval starting at t; recapture: occasion t)",
       y = "Probability") +
  theme_bw(base_size = 13) + theme(legend.position = "top")
print(p1)

# 3. Goodness of fit: a posterior predictive check (BPA 7.10.2) ----
# For each posterior draw, compare the observed m-array's Freeman-Tukey
# discrepancy with that of a data set simulated from the model. If the model
# fits, simulated data look about as discrepant as the real data, and the
# Bayesian p-value is far from 0 and 1. Lesson 06 uses the same check.

cat("\nTime-varying model, Bayesian p-value:", round(mean(cjs$sims.list$fit.new > cjs$sims.list$fit), 3), "\n")

# Now the wrong model: one survival and one recapture probability for all
# years, fitted to the same time-varying data. Only the priors differ from
# cjs-mnl.jags. With this seed the check catches it. With 30 animals marked
# per year it won't always: the check's power depends on how much data
# there are (exercise 3).
cat(file = "cjs-mnl-const.jags", "
model {
# Priors and constraints
for (t in 1:(n.occasions-1)){
   phi[t] <- mean.phi
   p[t] <- mean.p
   }
mean.phi ~ dunif(0, 1)         # Prior for mean survival
mean.p ~ dunif(0, 1)           # Prior for mean recapture
# Define the multinomial likelihood
for (t in 1:(n.occasions-1)){
   marr[t,1:n.occasions] ~ dmulti(pr[t, ], rel[t])
   }
# Define the cell probabilities of the m-array
# Main diagonal
for (t in 1:(n.occasions-1)){
   q[t] <- 1-p[t]                # Probability of non-recapture
   pr[t,t] <- phi[t]*p[t]
   # Above main diagonal
   for (j in (t+1):(n.occasions-1)){
      pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t,j] <- 0
      } #j
   } #t
# Last column: probability of non-recapture
for (t in 1:(n.occasions-1)){
   pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
   } #t

# Assess model fit using Freeman-Tukey statistic
for (t in 1:(n.occasions-1)){
   for (j in 1:n.occasions){
      expmarr[t,j] <- rel[t]*pr[t,j]
      E.org[t,j] <- pow((pow(marr[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      } #j
   } #t
for (t in 1:(n.occasions-1)){
   marr.new[t,1:n.occasions] ~ dmulti(pr[t, ], rel[t])
   for (j in 1:n.occasions){
      E.new[t,j] <- pow((pow(marr.new[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      } #j
   } #t
fit <- sum(E.org[,])
fit.new <- sum(E.new[,])
}
")
cjs.const <- jags(jags.data, function() list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1)),
                  c("mean.phi", "mean.p", "fit", "fit.new"), "cjs-mnl-const.jags",
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
if (!fit_ok(cjs.const)) warning("CJS constant m-array fit: diagnostics need attention before interpreting results.")
cat("Constant model, Bayesian p-value:",
    round(mean(cjs.const$sims.list$fit.new > cjs.const$sims.list$fit), 3), "\n")

ppc.df <- rbind(
  data.frame(model = "Time-varying phi and p", observed = cjs$sims.list$fit, simulated = cjs$sims.list$fit.new),
  data.frame(model = "Constant phi and p", observed = cjs.const$sims.list$fit, simulated = cjs.const$sims.list$fit.new))
p2 <- ggplot(ppc.df, aes(observed, simulated)) +
  geom_point(alpha = 0.15, size = 0.8) +
  geom_abline(color = "red", linewidth = 1) +
  facet_wrap(~ model) +
  labs(x = "Discrepancy, actual data", y = "Discrepancy, replicate data") +
  theme_bw(base_size = 13)
print(p2)

# 4. Real data: female Leisler's bats (BPA 7.11) ----
# The book's m-array for 19 years of captures of female Leisler's bats.
# Survival now varies among years as a random effect (BPA 7.4.2): logit(phi)
# is the mean plus a yearly deviation with SD sigma, as in lesson 06.

m.leisleri <- matrix(c(4,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
0,5,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
0,0,9,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
0,0,0,10,2,0,0,0,0,0,0,0,0,0,0,0,0,0,5,
0,0,0,0,10,2,1,0,0,0,0,0,0,0,0,0,0,0,6,
0,0,0,0,0,15,0,0,0,0,0,0,0,0,0,0,0,0,6,
0,0,0,0,0,0,11,2,0,1,0,0,0,0,0,0,0,0,19,
0,0,0,0,0,0,0,12,1,1,0,0,0,0,0,0,0,0,6,
0,0,0,0,0,0,0,0,13,2,0,0,0,0,0,0,0,0,4,
0,0,0,0,0,0,0,0,0,14,0,0,0,0,0,0,0,0,6,
0,0,0,0,0,0,0,0,0,0,13,1,0,0,0,1,0,0,8,
0,0,0,0,0,0,0,0,0,0,0,15,3,1,0,0,0,0,12,
0,0,0,0,0,0,0,0,0,0,0,0,12,4,0,1,0,0,7,
0,0,0,0,0,0,0,0,0,0,0,0,0,19,2,0,0,0,3,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,28,1,0,0,4,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,22,7,2,21,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,2,21,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,14,18), ncol = 19, nrow = 18, byrow = TRUE)

cat(file = "cjs-mnl-ran.jags", "
model {
# Priors and constraints
for (t in 1:(n.occasions-1)){
   logit(phi[t]) <- mu + epsilon[t]
   epsilon[t] ~ dnorm(0, tau)
   p[t] <- mean.p
   }
mean.phi ~ dunif(0, 1)             # Prior for mean survival
mu <- log(mean.phi / (1-mean.phi)) # Logit transformation
sigma ~ dunif(0, 5)                # Prior for standard deviation
tau <- pow(sigma, -2)
sigma2 <- pow(sigma, 2)
mean.p ~ dunif(0, 1)               # Prior for mean recapture
# Define the multinomial likelihood
for (t in 1:(n.occasions-1)){
   marr[t,1:n.occasions] ~ dmulti(pr[t,], rel[t])
   }
# Define the cell probabilities of the m-array:
# Main diagonal
for (t in 1:(n.occasions-1)){
   q[t] <- 1-p[t]
   pr[t,t] <- phi[t]*p[t]
   # Above main diagonal
   for (j in (t+1):(n.occasions-1)){
      pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t,j]<-0
      } #j
   } #t
# Last column: probability of non-recapture
for (t in 1:(n.occasions-1)){
   pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
   } # t
# Assess model fit using Freeman-Tukey statistic
for (t in 1:(n.occasions-1)){
   for (j in 1:n.occasions){
      expmarr[t,j] <- rel[t]*pr[t,j]
      E.org[t,j] <- pow((pow(marr[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      }
   }
for (t in 1:(n.occasions-1)){
   marr.new[t,1:n.occasions] ~ dmulti(pr[t,], rel[t])
   for (j in 1:n.occasions){
      E.new[t,j] <- pow((pow(marr.new[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      }
   }
fit <- sum(E.org[,])
fit.new <- sum(E.new[,])
}
")

leis.data <- list(marr = m.leisleri, n.occasions = dim(m.leisleri)[2], rel = rowSums(m.leisleri))
leis.inits <- function(){list(mean.phi = runif(1, 0, 1), sigma = runif(1, 0, 5), mean.p = runif(1, 0, 1))}
# The book runs 5,000 iterations. Longer chains pass the diagnostic screen.
leis.result <- jags(leis.data, leis.inits, c("mean.p", "mean.phi", "sigma", "fit", "fit.new"),
                    "cjs-mnl-ran.jags", n.chains = nc, n.thin = 6, n.iter = 25000, n.burnin = 5000,
                    parallel = TRUE)
print(leis.result, digits = 3)
if (!fit_ok(leis.result)) warning("Bat CJS fit: diagnostics need attention before interpreting results.")
cat("Bats: mean survival", round(leis.result$mean$mean.phi, 3), " 95% CrI",
    round(quantile(leis.result$sims.list$mean.phi, c(0.025, 0.975)), 3),
    "; recapture", round(leis.result$mean$mean.p, 3),
    "; Bayesian p-value", round(mean(leis.result$sims.list$fit.new > leis.result$sims.list$fit), 3), "\n")

# Try this ----
# Predict before rerunning; reset each change before trying the next.
# 1. In section 1, set marked <- rep(10, n.occasions-1). How much wider do the
#    intervals for phi and p get?
# 2. In cjs-mnl.jags (section 2), replace p[t] ~ dunif(0, 1) with p[t] <- mean.p
#    and add the line mean.p ~ dunif(0, 1) after the loop. Is the last phi
#    identified now? What does that buy you, and what does it assume?
# 3. Change set.seed(1) to set.seed(7) at the top. Does the check in section 3
#    still catch the constant model? What does a passing check prove?

# Lessons ----
# - CJS conditions on first capture. It estimates apparent survival and
#   recapture probability, not population size.
# - Survival is a hidden state; a recapture is seen only if the animal is
#   alive. Animals seen later were alive in between, which helps the fit.
# - The m-array summarizes capture histories into multinomial counts. Same
#   model, much faster than one row per animal.
# - The last survival and recapture can't be separated without another
#   occasion, only their product. Lesson 05 needs GOJ for the same reason.
# - Real data work the same way: the bat m-array goes straight into the model.
# - A posterior predictive check compares the data with data the model could
#   have produced. A failing check is evidence of misfit; a passing check
#   only means this check couldn't tell.

# Where BPA teaches this ----
# - 7.2  The CJS model as a state-space model
# - 7.3  Models with constant parameters; 7.3.1 passing known latent states
# - 7.4  Time variation: fixed, random, and covariate effects
# - 7.5 to 7.7  Group, individual, and age effects
# - 7.9  Parameter identifiability, the example in section 2
# - 7.10  The m-array and the multinomial likelihood, with the
#         Freeman-Tukey check used in section 3
# - 7.11  Real data: survival of female Leisler's bats, section 4 here
