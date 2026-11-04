# Run from the repository root: Rscript tests/audit_checks.R
# Focused regression checks; these do not run the long simulation studies.
source("helpers.R")
for (path in list.files(pattern = "\\.R$")) parse(path)

# Preserve chains, and reject missing, constant, nonfinite, or poorly mixed draws.
set.seed(20261029)
good <- list(samples = lapply(1:4, function(i)
  cbind(theta = rnorm(2000), wild_ge = rbeta(2000, 4, 30))))
stopifnot(fit_ok(good), "inverse_ge" %in% fit_diagnostics(good)$variable)
stopifnot(nrow(diagnostic_failures(fit_diagnostics(good))) == 0L,
          diagnostic_failures(NULL)$reason == "diagnostics unavailable")
# Report every failed criterion with its parameter, including ESS failures
# that a maximum-Rhat summary alone would miss.
d <- data.frame(variable = c("high_rhat", "low_bulk", "low_tail", "missing_mcse"),
                rhat = c(1.02, 1, 1, 1), ess_bulk = c(1000, 399, 1000, 1000),
                ess_tail = c(1000, 1000, 399, 1000), mcse_median = c(.01, .01, .01, NA))
issues <- diagnostic_failures(d)
stopifnot(identical(issues$variable, d$variable),
          identical(issues$reason, c("Rhat > 1.01", "bulk ESS < 400", "tail ESS < 400", "median MCSE unavailable")))
stopifnot(!fit_ok(NULL), !fit_ok(list(samples = list())))
bad <- good
bad$samples[[1]][1, 1] <- NA_real_
stopifnot(!fit_ok(bad))
bad <- good
bad$samples <- lapply(bad$samples, function(x) { x[, 1] <- 1; x })
stopifnot(!fit_ok(bad))
bad <- good
bad$samples[[1]][, 1] <- bad$samples[[1]][, 1] + 5
stopifnot(!fit_ok(bad))
bad <- good
bad$samples[[1]] <- bad$samples[[1]][-1, ]
stopifnot(!fit_ok(bad))

# Load only the formatter definition, not the simulations or JAGS fits.
expressions <- parse("05_six_cell.R")
for (expr in expressions) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("cell.of"))) eval(expr)
}
stopifnot(identical(cell.of(c(1,1,0,0,0,1), c(0,0,1,1,0,0),
                          c(0,1,0,1,1,0), c(0,0,0,0,0,1)),
                    c("c1","c2","c3","c4","c5","c6")))
must_error <- function(expr) inherits(tryCatch({ force(expr); NULL }, error = identity), "error")
stopifnot(must_error(cell.of(1, 1, 0, 0)),
          must_error(cell.of(1, 0, 1, 1)),
          must_error(cell.of(0, 0, 0, 1)),
          must_error(cell.of(NA, 0, 0, 0)),
          must_error(cell.of(2, 0, 0, 0)),
          must_error(cell.of(c(1, 0), 0, 0, 0)))

# A zero-success beta posterior has an infinite inverse mean but a finite
# median. The numerical median is a different estimand, not a clipped mean.
ge <- rbeta(200000, 1, 81)
stopifnot(abs(median(1 / ge) / (1 / qbeta(0.5, 1, 81)) - 1) < 0.02)
cat("Audit regression checks passed.\n")
