# Survival from dead recoveries: mark-recovery models and Swiss red kites
# (Run time: under 10 seconds. Times vary up to 2x between computers.)
# (Run 00_setup.R once first: it installs packages and checks JAGS.)
#
# Straight from Kéry and Schaub (2012), chapter 8. The simulation function,
# JAGS models, and red kite data are the book's, from the authors' JAGS
# translation (M. Schaub; jagsUI version by M. Kéry, 2025). Changes: a seed is
# set, the two red kite priors are written as two model files, and the repo's
# diagnostic checks and comments are added.
#
# In lesson 04a, marked animals were seen again alive. Here they are found
# dead and reported, like bird rings sent in by the public or tags returned
# by anglers. Two probabilities, estimated together:
#   s  survival from one year to the next
#   r  recovery: the chance a dead animal is found and reported
# A ringed bird dies in year j and is reported with probability
#   s^(j - 1) * (1 - s) * r,
# so the timing of recoveries tells survival apart from reporting.
# Section 1 uses simulated data; section 2 the book's real red kite data.
# Requires JAGS and the jagsUI package.

library(jagsUI)
source("helpers.R")
library(ggplot2)

set.seed(8)

# 1. Simulated data with constant survival and recovery (BPA 8.2.1, 8.3.1) ----

n.occasions <- 14                 # Number of release occasions
marked <- rep(50, n.occasions)    # Annual number of newly marked individuals
s <- rep(0.8, n.occasions)
r <- rep(0.2, n.occasions)

# Define matrices with survival and recovery probabilities
S <- matrix(s, ncol = n.occasions, nrow = sum(marked))
R <- matrix(r, ncol = n.occasions, nrow = sum(marked))

# Define function to simulate mark-recovery data
simul.mr <- function(S, R, marked){
   n.occasions <- dim(S)[2]
   MR <- matrix(NA, ncol = n.occasions+1, nrow = sum(marked))
   # Define a vector with the occasion of marking
   mark.occ <- rep(1:n.occasions, marked)
   # Fill the CH matrix
   for (i in 1:sum(marked)){
      MR[i, mark.occ[i]] <- 1    # Write an 1 at the release occasion
      for (t in mark.occ[i]:n.occasions){
         # Bernoulli trial: has individual survived occasion?
         sur <- rbinom(1, 1, S[i,t])
         if (sur==1) next    # If still alive, move to next occasion
         # Bernoulli trial: has dead individual been recovered?
         rp <- rbinom(1, 1, R[i,t])
         if (rp==0){
            MR[i,t+1] <- 0
            break
            }
         if (rp==1){
            MR[i,t+1] <- 1
            break
            }
      } #t
   } #i
      # Replace the NA in the file by 0
      MR[which(is.na(MR))] <- 0
      return(MR)
   }

MR <- simul.mr(S, R, marked)

# Define function to create an m-array based for mark-recovery (MR) data
marray.dead <- function(MR){
   nind <- dim(MR)[1]
   n.occasions <- dim(MR)[2]
   m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
   # Create vector with occasion of marking
   get.first <- function(x) min(which(x!=0))
   f <- apply(MR, 1, get.first)
   # Calculate the number of released individuals at each time period
   first <- as.numeric(table(f))
   for (t in 1:n.occasions){
      m.array[t,1] <- first[t]
      }
   # Fill m-array with recovered individuals
   rec.ind <- which(apply(MR, 1, sum)==2)
   rec <- numeric()
   for (i in 1:length(rec.ind)){
      d <- which(MR[rec.ind[i],(f[rec.ind[i]]+1):n.occasions]==1)
      rec[i] <- d + f[rec.ind[i]]
      m.array[f[rec.ind[i]],rec[i]] <- m.array[f[rec.ind[i]],rec[i]] + 1
      }
   # Calculate the number of individuals that are never recovered
   for (t in 1:n.occasions){
      m.array[t,n.occasions+1] <- m.array[t,1]-sum(m.array[t,2:n.occasions])
      }
   out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
   return(out)
   }

marr <- marray.dead(MR)

cat(file = "mr-mnl.jags", "
model {

# Priors and constraints
for (t in 1:n.occasions){
   s[t] <- mean.s
   r[t] <- mean.r
   }
mean.s ~ dunif(0, 1)              # Prior for mean survival
mean.r ~ dunif(0, 1)              # Prior for mean recovery

# Define the multinomial likelihood
for (t in 1:n.occasions){
   marr[t,1:(n.occasions+1)] ~ dmulti(pr[t,], rel[t])
   }
# Define the cell probabilities of the m-array
# Main diagonal
for (t in 1:n.occasions){
   pr[t,t] <- (1-s[t])*r[t]
   # Above main diagonal
   for (j in (t+1):n.occasions){
      pr[t,j] <- prod(s[t:(j-1)])*(1-s[j])*r[j]
      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t,j] <- 0
      } #j
   } #t
# Last column: probability of non-recovery
for (t in 1:n.occasions){
   pr[t,n.occasions+1] <- 1-sum(pr[t,1:n.occasions])
   } #t
}
")

jags.data <- list(marr = marr, n.occasions = dim(marr)[2]-1, rel = rowSums(marr))
inits <- function(){list(mean.s = runif(1, 0, 1), mean.r = runif(1, 0, 1))}
mr <- jags(jags.data, inits, c("mean.s", "mean.r"), "mr-mnl.jags",
           n.chains = 3, n.thin = 6, n.iter = 5000, n.burnin = 2000, parallel = TRUE)
print(mr, digits = 3)
if (!fit_ok(mr)) warning("Mark-recovery fit: diagnostics need attention before interpreting results.")
cat(sum(marked), "marked;", sum(marr[, ncol(marr)]), "never recovered\n")
cat("Survival: truth", s[1], " estimate", round(mr$mean$mean.s, 3), " 95% CrI",
    round(quantile(mr$sims.list$mean.s, c(0.025, 0.975)), 3), "\n")
cat("Recovery: truth", r[1], " estimate", round(mr$mean$mean.r, 3), " 95% CrI",
    round(quantile(mr$sims.list$mean.r, c(0.025, 0.975)), 3), "\n")

# 2. Real data: age-dependent survival in Swiss red kites (BPA 8.4) ----
# Over 50 years, 1,480 nestlings and 152 adults were ringed; 107 were found
# dead. Each m-array is one row: birds recovered in each year after ringing,
# then the number never recovered. Survival has three age classes (juvenile,
# subadult, adult) and recovery two (first year, later).
marray.juv <- c(42, 18, 5, 7, 4, 3, 2, 1, 2, 2, 1, 0, 1, 3, 0, 0, 1, 1388)
marray.ad <- c(3, 1, 1, 3, 0, 2, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 137)

# Juvenile survival is hard to estimate here, so the book fits it twice:
# A uses an informative Beta(4.2, 2.8) prior (mean 0.6, SD about 0.17), from
# earlier red kite studies; B uses a flat U(0, 1) prior. Comparing the two is
# a prior sensitivity check. Only the first prior line differs.
rk.model <- function(sjuv.prior) paste0("
model {

# Priors and constraints
sjuv ~ ", sjuv.prior, "  # Prior for juv. survival
ssub ~ dunif(0, 1)       # Prior for subad. survival
sad ~ dunif(0, 1)        # Prior for ad. survival
rjuv ~ dunif(0, 1)       # Prior for juv. recovery
rad ~ dunif(0, 1)        # Prior for ad. recovery

# Define the multinomial likelihoods
marr.j[1:(n.age+1)] ~ dmulti(pr.j[], rel.j)
marr.a[1:(n.age+1)] ~ dmulti(pr.a[], rel.a)
# Define the cell probabilities of the juvenile m-array
# First element
pr.j[1] <- (1-sjuv)*rjuv
# Second element
pr.j[2] <- sjuv*(1-ssub)*rad
# Third and further elements
for (t in 3:n.age){
   pr.j[t] <- sjuv*ssub*pow(sad,(t-3))*(1-sad)*rad
   }
# Probability of non-recovery
pr.j[n.age+1] <- 1-sum(pr.j[1:n.age])
# Define the cell probabilities of the adult m-array
# All elements
for (t in 1:n.age){
   pr.a[t] <- pow(sad,(t-1))*(1-sad)*rad
   }
# Probability of non-recovery
pr.a[n.age+1] <- 1-sum(pr.a[1:n.age])
}
")
cat(file = "mr-mnl-age3-A.jags", rk.model("dbeta(4.2, 2.8)"))
cat(file = "mr-mnl-age3-B.jags", rk.model("dunif(0, 1)"))

rk.data <- list(marr.j = marray.juv, marr.a = marray.ad, n.age = length(marray.juv)-1,
                rel.j = sum(marray.juv), rel.a = sum(marray.ad))
rk.inits <- function(){list(sjuv = runif(1, 0, 1), ssub = runif(1, 0, 1), sad = runif(1, 0, 1),
                            rjuv = runif(1, 0, 1), rad = runif(1, 0, 1))}
rk.params <- c("sjuv", "ssub", "sad", "rjuv", "rad")

rk.ageA <- jags(rk.data, rk.inits, rk.params, "mr-mnl-age3-A.jags",
                n.chains = 3, n.thin = 10, n.iter = 30000, n.burnin = 10000, parallel = TRUE)
rk.ageB <- jags(rk.data, rk.inits, rk.params, "mr-mnl-age3-B.jags",
                n.chains = 3, n.thin = 10, n.iter = 30000, n.burnin = 10000, parallel = TRUE)
print(rk.ageA, digits = 3)
print(rk.ageB, digits = 3)
if (!fit_ok(rk.ageA)) warning("Red kite fit A: diagnostics need attention before interpreting results.")
if (!fit_ok(rk.ageB)) warning("Red kite fit B: diagnostics need attention before interpreting results.")

for (fit in list(A = rk.ageA, B = rk.ageB)) {
  cat("Juvenile survival", round(median(fit$sims.list$sjuv), 3), " 95% CrI",
      round(quantile(fit$sims.list$sjuv, c(0.025, 0.975)), 3),
      "| adult survival", round(median(fit$sims.list$sad), 3), " 95% CrI",
      round(quantile(fit$sims.list$sad, c(0.025, 0.975)), 3), "\n")
}

sjuv.df <- rbind(
  data.frame(prior = "A: Beta(4.2, 2.8) prior", sjuv = rk.ageA$sims.list$sjuv),
  data.frame(prior = "B: flat U(0, 1) prior", sjuv = rk.ageB$sims.list$sjuv))
prior.df <- data.frame(sjuv = seq(0, 1, 0.005), density = dbeta(seq(0, 1, 0.005), 4.2, 2.8))
p1 <- ggplot(sjuv.df, aes(sjuv, color = prior)) +
  geom_density(linewidth = 1) +
  geom_line(data = prior.df, aes(sjuv, density), inherit.aes = FALSE, linetype = "dashed") +
  scale_color_manual(NULL, values = c("A: Beta(4.2, 2.8) prior" = "red",
                                      "B: flat U(0, 1) prior" = "black")) +
  labs(x = "Juvenile survival (dashed: the Beta(4.2, 2.8) prior)", y = "Posterior density") +
  theme_bw(base_size = 13) + theme(legend.position = "top")
print(p1)

# Try this ----
# Predict before rerunning; reset each change before trying the next.
# 1. In section 1, set marked <- rep(10, n.occasions). How much wider are the
#    intervals for s and r?
# 2. In section 2, change "dbeta(4.2, 2.8)" to "dbeta(42, 28)": the same mean,
#    0.6, with a third of the SD. Does juvenile survival move now? Which
#    is speaking, the prior or the data?

# Lessons ----
# - Dead recoveries estimate survival and reporting together. When birds are
#   found, relative to ringing, separates the two.
# - Real data go into the same m-array the simulation used.
# - Juvenile survival is poorly determined by these data: its interval is
#   wide under either prior.
# - Fit with more than one prior and compare. Here the informative prior
#   changes little, so the uncertainty is in the data, not the prior.

# Where BPA teaches this ----
# - 8.1  Introduction
# - 8.2  The mark-recovery model as a state-space model
# - 8.3  The multinomial likelihood (m-array): constant (8.3.1) and
#        age-dependent (8.3.2) parameters
# - 8.4  Real-data example: age-dependent survival in Swiss red kites, with
#        the informative and flat priors compared here
