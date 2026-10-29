# Glossary

Short definitions of the fish and statistics terms in these scripts. If a script uses a word you don't know, look here first.

## The setting in one paragraph

Young salmon and steelhead (smolts) swim down the Snake River to the ocean in spring. At Lower Granite Dam, each fish passes one of three ways: it's guided into a juvenile bypass channel, or it goes over the spillway, or it goes through the turbines. Many fish carry a small PIT tag (see PIT tag below), read automatically by antennas at detection sites. At Lower Granite, antennas read fish in the bypass (GRJ) and on the spillway (GRS). Some bypassed fish are loaded onto barges and carried downstream. Farther downstream, antennas at Little Goose Dam (GOJ) read some fish a second time. Which sites read a fish is its detection history. The route it took at Lower Granite is never seen directly; the models work it out from the detection histories.

## Fish and fisheries terms

**Barge (transport).** Some fish guided into the bypass are loaded on barges and carried past the downstream dams. A barged fish can't be read at GOJ.

**Bypass.** A channel at the dam that guides fish away from the turbines. Fish in it are read at GRJ.

**Detection history.** The set of sites that read a tagged fish, for example "GRS and GOJ". Each history is one cell in the models.

**Efficiency trial.** Marking a known number of fish, releasing them upstream of a trap, and counting how many the trap catches again. The share recaught estimates trap efficiency.

**GE (guidance efficiency, psi).** The chance a fish arriving at the dam is guided into the bypass. See also Bypass.

**GOJ.** PIT detection site in the juvenile bypass at Little Goose Dam, downstream of Lower Granite.

**GRJ.** PIT detection site in the juvenile bypass at Lower Granite Dam.

**GRS.** PIT detection array in spillbay 1 at Lower Granite. Its rows estimate read probability within that route, not the fraction of all nonguided fish using it.

**Hatchery and wild (rear type).** Fish raised in a hatchery vs fish born in the river. Wild fish are usually the target, but most tagged fish are hatchery.

**Night share (a).** The fraction of fish that pass the dam at night.

**PIT tag.** A small passive integrated transponder tag with a unique code, read by antennas without catching the fish.

**PTAGIS.** The public database of PIT tag detections in the Columbia River basin: https://www.ptagis.org/

**RSW (removable spillway weir).** A structure providing surface passage through a spillbay. At Lower Granite, it directs surface-oriented passage toward spillbay 1.

**Sample rate.** Here, the fraction of time the trap sampling gates are open. Expanding by this rate assumes sampled periods represent passage during unsampled periods.

**SCRAPI / SCRAPI2.** The smolt abundance and composition analysis workflows discussed in the talk. The teaching comparison holds GE fixed for SCRAPI-style bootstraps and carries GE posterior draws for SCRAPI2-style expansion. The examples omit several production steps.

**Screw trap.** A floating, rotating cone in a river that catches a fraction of the fish passing it.

**Smolt.** A young salmon or steelhead migrating to the ocean.

**smoltEASE.** The modeling code used for the talk's PIT-based guidance estimates, including hierarchical guidance, detection, recovery, and transport components. The scripts here are simplified teaching examples.

**Spill.** Water released over the spillway instead of through the turbines. More spill tends to send more fish over the spillway and fewer into the bypass.

**Stratum (strata).** A block of consecutive weeks grouped together.

**Trap efficiency (p).** The chance a passing fish is caught by the trap. See also Efficiency trial. (In the dam scripts, p is the chance GRS reads a fish given that it wasn't guided into the bypass, whether it then went over the spillway or through the turbines. It is not the chance of being read given spillway passage.)

## Statistics terms

**Bernoulli trial.** One yes-or-no outcome with probability p, like one fish guided or not.

**Bias.** How far an estimate lands from the truth, on average over many repeats.

**Binomial.** The number of yeses in a fixed number of independent Bernoulli trials, like recaptures out of M marked fish.

**Bootstrap.** Estimating uncertainty by recomputing an answer many times on randomly redrawn data.

**Burn-in.** The first MCMC draws, thrown away while the chains move away from their starting points.

**Chain.** One run of the MCMC sampler. Running several chains from different starts and checking they agree is how you judge convergence. See also Rhat and Trace plot.

**Complete pooling vs partial pooling.** Complete pooling forces groups to share one value. Partial pooling (a hierarchical model) lets each group have its own value but pulls it toward the group average, more strongly for groups with less data.

**Conditioning on detection.** Fish never detected anywhere aren't in the data, so the model divides each cell's probability by the chance of being detected at all.

**Convergence.** When MCMC chains have settled into the posterior. Checked with Rhat and trace plots (see both below).

**Convex.** Curving upward, like 1 / p. Averaging a convex function gives more than applying it to the average (Jensen's inequality).

**Coverage.** Over many simulated data sets, the share of 95% intervals that contain the truth. It measures calibration under the simulated conditions. Nominal 95% Bayesian intervals do not guarantee 95% coverage at every fixed parameter value.

**Credible interval.** A Bayesian interval: given the model and the data, the parameter is inside it with the stated probability. See also Posterior and Coverage.

**Hierarchical model.** A model where group-level values (each week's GE, say) are themselves drawn from a distribution with its own parameters.

**Identifiability.** Whether the data can pin down a parameter. If two parameters can trade off and fit the data equally well, they aren't separately identified.

**Jensen's inequality.** For a curved function, the average of the function isn't the function of the average. Here: the average of 1 / p is larger than 1 / (average p). See also Convex.

**Latent state.** A true quantity that isn't observed directly, like daily passage or a fish's route.

**Likelihood.** How probable the observed data are for each possible parameter value.

**Lincoln-Petersen.** A capture-recapture estimator using the overlap between marked and subsequently sampled animals: the basic form is M*C/R. Corrections depend on the sampling design. The talk's Bailey modification is M*(C+1)/(R+1).

**Logit.** log(p / (1 − p)). It maps probabilities (0 to 1) onto the whole number line, so effects can be added without pushing a probability past 0 or 1.

**MCMC (Markov chain Monte Carlo).** The simulation method JAGS uses to draw samples from the posterior when there is no exact formula.

**Multinomial.** The joint counts across mutually exclusive categories in a fixed number of independent trials with common category probabilities. Each fish has one categorical outcome; the counts across fish are multinomial.

**n.eff.** How many independent draws the correlated MCMC draws are worth.

**Posterior.** What you believe about a parameter after seeing the data: prior times likelihood, rescaled.

**Precision.** 1 / variance. JAGS's dnorm uses precision, not standard deviation: dnorm(0, 0.1) has SD about 3.2.

**Prior.** What you believe about a parameter before seeing the data.

**Rhat.** Compares variation between chains to variation within chains. Near 1 supports agreement. These lessons also check rank-normalized split Rhat <= 1.01, bulk/tail effective sample sizes >= 400, and trace plots. No single diagnostic proves convergence.

**Sampling distribution.** All the estimates you could have gotten from repeating the same study under the same truth.

**Shrinkage.** In a hierarchical model, a group's estimate pulled toward the average of all groups. See also Hierarchical model and Complete pooling vs partial pooling.

**Simpson's paradox.** An aggregate comparison reverses the comparisons within its component groups because their proportions differ. In the talk, day and night GE both fall while pooled GE rises as night passage becomes more common.

**State-space model.** A model with two parts: a process model for the true state over time, and an observation model for what you actually measure.

**Thinning.** Keeping only every k-th MCMC draw to save memory.

**Trace plot.** Each chain's draws plotted in order. Converged chains look like grassy, overlapping fuzz.
