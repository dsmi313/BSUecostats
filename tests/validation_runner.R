# Exercise the R-only scheduler without fitting models or needing JAGS.
# Run from the repository root: Rscript tests/validation_runner.R
source("validation/run_all.R")
repo <- normalizePath(".")
fixture <- tempfile("validation runner ")
dir.create(file.path(fixture, "validation"), recursive = TRUE)
stopifnot(all(file.copy(list.files(repo, "\\.R$", full.names = TRUE), fixture)))
stopifnot(file.copy("validation/run_all.R", file.path(fixture, "validation")))
# A small child script lets this test check scheduling, patches, paths with
# spaces, logs, and nonzero exits without rerunning the statistical lessons.
writeLines(c(
  "args <- commandArgs(TRUE)",
  "setwd(args[1])",
  "cat('child ran:', args[2], '\\n')",
  "saveRDS(list(lesson = args[2]), 'result.rds')",
  "quit(status = if (args[2] == '03_bernoulli_logit.R') 7L else 0L)"
), file.path(fixture, "validation", "run_lesson.R"))
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
script <- file.path(fixture, "validation", "run_all.R")
run <- function(output, extra) suppressWarnings(system2(
  rscript, shQuote(c(script, "--output", output, extra)),
  stdout = "", stderr = ""))
out <- file.path(fixture, "parallel output")
stopifnot(run(out, c("--workers", "2", "--cases", "01_small", "02_small")) == 0L)
m <- read.csv(file.path(out, "manifest.csv"))
stopifnot(nrow(m) == 2L, all(m$exit == 0L),
          identical(m$case, c("01_small", "02_small")),
          nrow(readRDS(file.path(out, "manifest.rds"))) == 2L,
          grepl("M <- 10", read_text(file.path(out, "01_small", "01_first_bayesian_model.R")), fixed = TRUE),
          grepl("M <- 15", read_text(file.path(out, "02_small", "02_lincoln_petersen.R")), fixed = TRUE),
          grepl("child ran:", read_text(file.path(out, "01_small", "run.log")), fixed = TRUE))
# Failed children must make the scheduler fail, while preserving other cases.
stopifnot(run(out, c("--workers", "1", "--cases", "01_small", "03_default")) != 0L)
m <- read.csv(file.path(out, "manifest.csv"))
stopifnot(m$exit[m$case == "01_small"] == 0L, m$exit[m$case == "03_default"] == 7L)
# A preparation failure must clear stale results and be reported.
p <- file.path(fixture, "01_first_bayesian_model.R")
cat(gsub("M <- 50", "M <- 51", read_text(p), fixed = TRUE), file = p)
stopifnot(run(out, c("--workers", "1", "--cases", "01_small")) != 0L,
          !file.exists(file.path(out, "01_small", "result.rds")))
m <- read.csv(file.path(out, "manifest.csv"))
stopifnot(grepl("Missing exercise text", m$error, fixed = TRUE))
stopifnot(run(out, c("--cases", "nonexistent")) != 0L,
          run(out, c("--workers", "0")) != 0L)
unlink(fixture, recursive = TRUE)
cat("R validation runner checks passed.\n")
