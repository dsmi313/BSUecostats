# Focused checks only: no JAGS fits or full lesson execution.
# Run from the repository root: Rscript tests/low_recovery.R
source("helpers.R")

load_definition <- function(path, name, envir = .GlobalEnv) {
  for (expr in parse(path)) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        identical(expr[[2]], as.name(name))) {
      eval(expr, envir = envir)
      return(invisible(NULL))
    }
  }
  stop("Definition not found: ", name)
}

reference <- dam_recovery_defaults()
stopifnot(abs(reference["phi_S"] - 123 / 3785) < 1e-12,
          abs(reference["phi_B"] - 35 / 355) < 1e-12)

# Evaluate only the four calibration assignments in each lesson, never its
# data simulation or JAGS calls. All routes must agree with the shared inputs.
for (path in c("05_six_cell.R", "06_hierarchical.R", "07_day_night.R")) {
  e <- new.env(parent = .GlobalEnv)
  for (name in c("recovery", "phi_S", "phi_B", "delta")) load_definition(path, name, e)
  stopifnot(abs(e$phi_S - reference["phi_S"]) < 1e-12,
            abs(e$phi_B - reference["phi_B"]) < 1e-12,
            abs(plogis(qlogis(e$phi_S) - e$delta) - e$phi_B) < 1e-12,
            e$delta < 0) # bypass recovery is higher in this calibration
}

load_definition("09_recover_truth.R", "draw_low_recovery_truth")
set.seed(20261104)
draws <- replicate(500, draw_low_recovery_truth())
stopifnot(all(draws["phi_S", ] >= .01 & draws["phi_S", ] <= .06),
          all(draws["phi_B", ] >= .04 & draws["phi_B", ] <= .15),
          max(abs(plogis(qlogis(draws["phi_S", ]) - draws["delta", ]) -
                    draws["phi_B", ])) < 1e-12)

# The low-recovery regime genuinely permits zero overlaps: a 500-fish
# sample at the low end should not be treated as automatically informative.
p.no.overlap <- (1 - (1 - .10) * .40 * .01)^500
stopifnot(p.no.overlap > .1)

# Exercise both prior settings by generating only the model text in a temp
# directory. Never simulate a season or invoke JAGS. The switch must change
# prior support while leaving the rest of the likelihood/model unchanged.
check_spill_models <- function() {
  expressions <- parse("06_hierarchical.R")
  writers <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("cat")) &&
                      identical(as.list(e)$file, "ge_hierarchical.jags"), expressions)
  switches <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("if")) &&
                       identical(e[[2]], quote(!restrict.spill)), expressions)
  stopifnot(length(writers) == 1L, length(switches) == 1L)
  original <- getwd()
  folder <- tempfile("spill-model-check-")
  dir.create(folder)
  on.exit({ setwd(original); unlink(folder, recursive = TRUE) })
  setwd(folder)
  model <- list()
  for (restricted in c(TRUE, FALSE)) {
    e <- new.env(parent = .GlobalEnv)
    e$restrict.spill <- restricted
    eval(writers[[1]], e)
    eval(switches[[1]], e)
    model[[as.character(restricted)]] <- readLines("ge_hierarchical.jags")
  }
  expected <- sub("beta ~ dnorm(0, 0.1) T(, 0)", "beta ~ dnorm(0, 0.1)",
                  model$`TRUE`, fixed = TRUE)
  stopifnot(identical(expected, model$`FALSE`),
            sum(grepl("T(, 0)", model$`TRUE`, fixed = TRUE)) == 1L,
            !any(grepl("T(, 0)", model$`FALSE`, fixed = TRUE)))
}
check_spill_models()
cat("Low-recovery calibration checks passed; no JAGS models were fitted.\n")
