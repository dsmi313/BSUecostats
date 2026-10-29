# Setup: run this once before the other scripts
#
# Installs any missing R packages, checks that JAGS is installed and that R
# can find it, and prints the versions you're running. If anything fails,
# the message says what to do.

# 1. R packages ----

needed <- c("jagsUI", "ggplot2", "posterior")
missing <- needed[!sapply(needed, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  cat("Installing:", missing, "\n")
  install.packages(missing, repos = "https://cloud.r-project.org")
}

# linewidth scales used by the lessons were added in ggplot2 3.4.0.
if (requireNamespace("ggplot2", quietly = TRUE) && packageVersion("ggplot2") < "3.4.0") {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("ggplot2", quietly = TRUE) || packageVersion("ggplot2") < "3.4.0")
  stop("Install ggplot2 >= 3.4.0, restart R, and rerun setup.")
if (!requireNamespace("posterior", quietly = TRUE))
  stop("Install the posterior package and rerun setup.")

# 2. JAGS itself ----
# JAGS is a separate program, not an R package. Download it from
# https://sourceforge.net/projects/mcmc-jags/, install it, then restart R. On Windows,
# install JAGS before jagsUI. jagsUI loads rjags, which fails with an error
# if JAGS can't be found.

jags.ok <- tryCatch({
  suppressPackageStartupMessages(library(jagsUI))
  TRUE
}, error = function(e) {
  cat("\nR could not load jagsUI. Usually this means JAGS isn't installed, or R\n",
      "was open when you installed it. Install JAGS from\n",
      "https://sourceforge.net/projects/mcmc-jags/, restart R, and run this script again.\n",
      "The error was:", conditionMessage(e), "\n")
  FALSE
})

# 3. A tiny test model ----
# Fits a coin-flip model: 7 heads in 10 flips. The posterior mean should be
# close to 8 / 12 = 0.67 (a flat prior plus the data).

if (jags.ok) {
  test.file <- file.path(tempdir(), "setup_test.jags")   # temporary, keeps your folder clean
  cat(file = test.file, "
model {
theta ~ dbeta(1, 1)
heads ~ dbin(theta, flips)
}
")
  test <- jags(list(heads = 7, flips = 10), NULL, "theta", test.file,
               n.chains = 3, n.iter = 2000, n.burnin = 500, n.thin = 1, verbose = FALSE)
  cat("\nTest model posterior mean of theta:", round(test$mean$theta, 2), "(expect about 0.67)\n")
  cat("JAGS is working. You're ready for 01 to 09.\n\n")
}

# 4. Record your environment ----
# Save this output with your results. It records the R, JAGS, and package
# versions, so you (or anyone else) can rebuild the same setup later.

if (jags.ok) cat("JAGS version:", as.character(rjags::jags.version()), "\n")
for (pkg in needed) {
  if (requireNamespace(pkg, quietly = TRUE)) cat(pkg, "version:", as.character(packageVersion(pkg)), "\n")
}
sessionInfo()
