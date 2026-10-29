# Run from the repository root: Rscript tests/talk_examples.R
# These checks use base R only and do not run any JAGS fits.
load_definition <- function(path, name) {
  for (expr in parse(path)) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        identical(expr[[2]], as.name(name))) {
      eval(expr, envir = .GlobalEnv)
      return(invisible(NULL))
    }
  }
  stop("Definition not found: ", name)
}
load_definition("02_lincoln_petersen.R", "bailey")
stopifnot(round(bailey(200, 1000, 20)) == 9533,
          round(1000 * 201 / 21) == 9571,
          200 * 1000 / 20 == 10000)
stopifnot(isTRUE(all.equal((1 - c(.1, .4)) * c(.04, .02) +
                            c(.1, .4) * c(.30, .25), c(.066, .112))))
stopifnot(abs(1 - prod(1 - c(.64, .61, .53)) - .934012) < 1e-12)
stopifnot(identical(as.numeric(round(qbinom(c(.05, .95), 4000, .0325) / .0325)),
                    c(3446, 4585)))

load_definition("06_hierarchical.R", "choose_strata")
# Strictly greater than 10, plus backward merging of an inadequate tail.
stopifnot(identical(choose_strata(c(5,5,1,2), rep(0,4), c(5,5,1,2)), rep(1L,4)),
          identical(choose_strata(c(6,6,6,6,1,1), rep(0,6), c(6,6,6,6,1,1)),
                    c(1L,1L,2L,2L,2L,2L)))
# Both conditions must pass; adequate recovery alone is insufficient.
stopifnot(identical(choose_strata(c(20,20), c(0,0), c(5,6)), c(1L,1L)))
warned <- FALSE
withCallingHandlers(choose_strata(c(1,1), c(0,0), c(1,1)),
                    warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") })
stopifnot(warned)

load_definition("05_six_cell.R", "transport.summary")
n <- c(c1=40, c2=30, c3=180, c4=180, c5=270, c6=90)
changed <- n
changed["c1"] <- n["c1"] + n["c6"]
changed["c6"] <- 0
before <- transport.summary(n)
after <- transport.summary(changed)
stopifnot(before["GE"] == after["GE"], after["transport"] == 0,
          after["bypass.recovery"] < before["bypass.recovery"])
cat("Talk example checks passed.\n")
