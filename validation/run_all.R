# Run lessons and exercises in isolated folders using R only.
# Rscript validation/run_all.R --output /tmp/BSUecostats-validation --workers 2
# Requires R/JAGS and the packages installed by 00_setup.R.

read_text <- function(path) {
  paste0(paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n"), "\n")
}

replace_first <- function(text, old, new) {
  at <- regexpr(old, text, fixed = TRUE)[1]
  if (at < 0L) stop("Missing exercise text: ", substr(old, 1, 100))
  paste0(substr(text, 1, at - 1L), new, substring(text, at + nchar(old)))
}

validation_cases <- function(repo) {
  paths <- sort(list.files(repo, "^0.*\\.R$"))
  files <- setNames(paths, sub("_.*", "", paths))
  cases <- list()
  add <- function(name, key, patches = list()) {
    if (!key %in% names(files)) stop("Missing lesson: ", key)
    cases[[length(cases) + 1L]] <<- list(name = name, filename = files[[key]], patches = patches)
  }
  for (key in names(files)) add(paste0(key, "_default"), key)
  add("01_small", "01", list(list("M <- 50", "M <- 10"), list("R <- 8", "R <- 2")))
  add("01_prior", "01", list(list("M <- 50", "M <- 10"), list("R <- 8", "R <- 2"), list("print(rbind(prior.table(50, 8), prior.table(500, 80)), row.names = FALSE)", "print(prior.table(10, 2), row.names = FALSE)")))
  add("01_catch", "01", list(list("C <- 600", "C <- 1200")))
  add("02_small", "02", list(list("M <- 50", "M <- 15")))
  add("02_lowp", "02", list(list("p <- 0.15", "p <- 0.03")))
  add("02a_heterogeneity", "02a", list(list("sim.Mh(sd = 1)", "sim.Mh(sd = 2)")))
  add("02a_two", "02a", list(list("sim.Mt(T = 3)", "sim.Mt(T = 2)")))
  add("02a_notrap", "02a", list(list("sim.Mb(N = 200)", "sim.Mb(N = 200, c = 0.3)")))
  add("03_small", "03", list(list("n.fish <- 400", "n.fish <- 50")))
  add("03_zero", "03", list(list("b1 <- -0.8", "b1 <- 0")))
  add("04_outage", "04", list(list("trap.out <- 41:45", "trap.out <- 31:55")))
  add("04_small", "04", list(list("M <- rep(30,", "M <- rep(10,")))
  add("04_outage_long", "04", list(list("trap.out <- 41:45", "trap.out <- 31:55"), list("ni <- 60000; nb <- 30000", "ni <- 240000; nb <- 120000")))
  add("04_small_long", "04", list(list("M <- rep(30,", "M <- rep(10,"), list("ni <- 60000; nb <- 30000", "ni <- 480000; nb <- 240000")))
  add("04a_fewer", "04a", list(list("marked <- rep(50, n.occasions-1)", "marked <- rep(10, n.occasions-1)")))
  add("04a_constp", "04a", list(list("   p[t] ~ dunif(0, 1)           # Priors for recapture\n   }", "   p[t] <- mean.p\n   }\nmean.p ~ dunif(0, 1)")))
  add("04a_seed", "04a", list(list("set.seed(1)", "set.seed(7)")))
  add("04b_lowp", "04b", list(list("p <- rep(0.5, n.occasions)", "p <- rep(0.2, n.occasions)")))
  add("04b_fewrows", "04b", list(list("nz <- 300", "nz <- 30")))
  add("04c_fewer", "04c", list(list("marked <- rep(50, n.occasions)", "marked <- rep(10, n.occasions)")))
  add("04c_tight", "04c", list(list("rk.model(\"dbeta(4.2, 2.8)\")", "rk.model(\"dbeta(42, 28)\")")))
  add("05_facility", "05", list(list("n.facility <- 1500", "n.facility <- 5000")))
  add("05_small", "05", list(list("n.tagged <- 6000", "n.tagged <- 600")))
  source <- read_text(file.path(repo, files[["05"]]))
  at <- regexpr("# 8.", source, fixed = TRUE)[1]
  if (at < 0L) stop("Missing section 8 in lesson 05.")
  tail <- substring(source, at)
  if (!grepl("psi ~ dbeta(1, 1)", tail, fixed = TRUE)) stop("Missing no-GOJ prior.")
  add("05_prior", "05", list(list(tail, gsub("psi ~ dbeta(1, 1)", "psi ~ dbeta(1, 9)", tail, fixed = TRUE))))
  add("06_more", "06", list(list("log(8), log(150)", "log(8), log(1500)")))
  add("06_variable", "06", list(list("sigma_psi <- 0.3", "sigma_psi <- 0.8")))
  add("06_free", "06", list(list("restrict.spill <- TRUE", "restrict.spill <- FALSE")))
  add("06_positive_restricted", "06", list(list("beta <- -0.8", "beta <- 0.4")))
  add("06_positive_free", "06", list(list("beta <- -0.8", "beta <- 0.4"), list("restrict.spill <- TRUE", "restrict.spill <- FALSE")))
  add("06_repeats", "06", list(list("run.repeats <- FALSE", "run.repeats <- TRUE")))
  add("06_variable_long", "06", list(list("sigma_psi <- 0.3", "sigma_psi <- 0.8")))
  add("06_repeats_long", "06", list(list("run.repeats <- FALSE", "run.repeats <- TRUE"), list("\"ge_hierarchical_ppc.jags\",\n       n.chains = 3, n.iter = 10000, n.burnin = 5000", "\"ge_hierarchical_ppc.jags\",\n       n.chains = 3, n.iter = 40000, n.burnin = 20000")))
  add("07_equal", "07", list(list("a.true <- c(wild = 0.18", "a.true <- c(wild = 0.09")))
  add("07_zero", "07", list(list("night_eff <- 3.0", "night_eff <- 0")))
  add("07_small", "07", list(list("n.tagged <- c(wild = 80", "n.tagged <- c(wild = 20")))
  add("07_repeats", "07", list(list("run.repeats <- FALSE", "run.repeats <- TRUE")))
  add("07_small_long", "07", list(list("n.tagged <- c(wild = 80", "n.tagged <- c(wild = 20")))
  add("07a_equal", "07a", list(list("result <- run_borrow_simulation()", "wk <- borrow_inputs(); wk$aW <- wk$aH\n  result <- run_borrow_simulation(wk)")))
  add("07a_small", "07a", list(list("result <- run_borrow_simulation()", "wk <- borrow_inputs(); wk$NW <- pmax(1, round(wk$NW / 10))\n  result <- run_borrow_simulation(wk)")))
  add("07a_departure", "07a", list(list("result <- run_borrow_simulation()", "result <- run_borrow_simulation(wild.ge.multiplier = 0.8)")))
  add("08_small", "08", list(list("n.pit <- 600", "n.pit <- 100")))
  add("08_fewer", "08", list(list("\nB <- 2000", "\nB <- 200")))
  add("09_small", "09", list(list("n.tagged <- 3000", "n.tagged <- 500")))
  add("09_more", "09", list(list("n.reps <- 50", "n.reps <- 200")))
  cases
}

run_case <- function(case, repo, root, runner) {
  started <- Sys.time()
  tryCatch({
    d <- file.path(root, case$name)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    # Clear stale results before attempting a rerun.
    unlink(file.path(d, c("result.rds", "direct_fits.rds", "run.log")))
    sources <- list.files(repo, "\\.R$", full.names = TRUE)
    if (!all(file.copy(sources, d, overwrite = TRUE))) stop("Could not copy lesson files.")
    path <- file.path(d, case$filename)
    text <- read_text(path)
    for (patch in case$patches) text <- replace_first(text, patch[[1]], patch[[2]])
    if (case$name == "07_small_long")
      text <- gsub("n.iter = 6000, n.burnin = 2000",
                   "n.iter = 24000, n.burnin = 8000", text, fixed = TRUE)
    if (case$name == "06_variable_long")
      text <- gsub("n.iter = 10000, n.burnin = 5000",
                   "n.iter = 40000, n.burnin = 20000", text, fixed = TRUE)
    cat(text, file = path)
    cat("START", case$name, "\n")
    rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
    status <- system2(rscript, shQuote(c(runner, d, case$filename)),
                      stdout = file.path(d, "run.log"), stderr = file.path(d, "run.log"))
    elapsed <- round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1)
    cat("DONE", case$name, "exit =", status, "seconds =", elapsed, "\n")
    data.frame(case = case$name, exit = as.integer(status), seconds = elapsed, error = "")
  }, error = function(e) {
    data.frame(case = case$name, exit = 1L,
               seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1),
               error = conditionMessage(e))
  })
}

validation_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  output <- file.path(dirname(tempdir()), "BSUecostats-validation")
  workers <- 2L
  selected <- NULL
  i <- 1L
  while (i <= length(args)) {
    option <- args[[i]]
    if (option %in% c("--help", "-h")) {
      cat("Usage: Rscript validation/run_all.R [--output DIR] [--workers N] [--cases NAME ...]\n",
          "Run 00_setup.R first. Omit --cases to run all lessons and exercises.\n", sep = "")
      return(0L)
    }
    if (option %in% c("--output", "--workers")) {
      if (i == length(args) || startsWith(args[[i + 1L]], "--"))
        stop("Missing value for ", option)
      i <- i + 1L
      if (option == "--output") output <- args[[i]] else {
        if (!grepl("^[0-9]+$", args[[i]])) stop("--workers must be a positive integer.")
        workers <- suppressWarnings(as.integer(args[[i]]))
        if (is.na(workers) || workers < 1L) stop("--workers must be a positive integer.")
      }
    } else if (option == "--cases") {
      selected <- character()
      while (i < length(args) && !startsWith(args[[i + 1L]], "--")) {
        i <- i + 1L
        selected <- c(selected, args[[i]])
      }
      if (!length(selected)) stop("--cases needs at least one case name.")
    } else stop("Unknown option: ", option)
    i <- i + 1L
  }
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script) != 1L) stop("Run this file with Rscript.")
  script <- sub("^--file=", "", script)
  # Rscript encodes spaces in its own --file argument as ~+~.
  if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
  script <- normalizePath(script, mustWork = TRUE)
  repo <- dirname(dirname(script))
  runner <- file.path(dirname(script), "run_lesson.R")
  cases <- validation_cases(repo)
  case_names <- vapply(cases, function(x) x$name, character(1))
  if (!is.null(selected)) {
    unknown <- setdiff(selected, case_names)
    if (length(unknown)) stop("Unknown cases: ", paste(unknown, collapse = ", "))
    cases <- cases[case_names %in% selected]
  }
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  root <- normalizePath(output, mustWork = TRUE)
  results <- list()
  workers <- min(workers, length(cases))
  if (workers > 1L) {
    cluster <- parallel::makePSOCKcluster(workers, outfile = "")
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterExport(cluster, c("run_case", "read_text", "replace_first"),
                            envir = environment(validation_main))
  }
  # Checkpoint after each batch, and use PSOCK workers on Windows as well.
  for (first in seq.int(1L, length(cases), by = workers)) {
    batch <- cases[first:min(first + workers - 1L, length(cases))]
    completed <- if (workers == 1L) {
      lapply(batch, run_case, repo = repo, root = root, runner = runner)
    } else parallel::parLapply(cluster, batch, run_case, repo = repo, root = root, runner = runner)
    results <- c(results, completed)
    manifest <- do.call(rbind, results)
    write.csv(manifest, file.path(root, "manifest.csv"), row.names = FALSE)
    saveRDS(manifest, file.path(root, "manifest.rds"))
  }
  as.integer(any(manifest$exit != 0L))
}

if (sys.nframe() == 0L) {
  status <- tryCatch(validation_main(), error = function(e) {
    message(conditionMessage(e))
    1L
  })
  quit(save = "no", status = status)
}
