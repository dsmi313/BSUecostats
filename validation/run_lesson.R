args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])
options(warn = 1)
source("helpers.R")
pdf("plots.pdf", width = 11, height = 7)
e <- new.env(parent = .GlobalEnv)
records <- list()
warning.messages <- character()
e$jags <- function(...) {
  arguments <- list(...)
  f <- jagsUI::jags(...)
  d <- fit_diagnostics(f)
  records[[length(records) + 1L]] <<- list(
    model = if (!is.null(arguments$model.file)) arguments$model.file else as.character(arguments[[4]]),
    passed = fit_ok(f), diagnostics = d, issues = diagnostic_failures(d), summary = f$summary)
  saveRDS(records, "direct_fits.rds")
  f
}
started <- proc.time()[3]
error.message <- NULL
tryCatch(withCallingHandlers(source(args[2], local = e, echo = FALSE), warning = function(w) {
  warning.messages <<- c(warning.messages, conditionMessage(w))
  cat("VALIDATION WARNING:", conditionMessage(w), "\n")
  invokeRestart("muffleWarning")
}), error = function(err) {
  error.message <<- conditionMessage(err)
  cat("VALIDATION ERROR:", error.message, "\n")
})
elapsed <- unname(proc.time()[3] - started)
keep <- c("post.mean","post.ci","prior.table","N","p","C","R","M","N.simple","N.adj",
          "simple.sim","adj.sim","R.sim","p.draws","n.sims","N_seen","slide.est","flag.check",
          "truth","est","lo","hi","season","w1","w2","beta","n.tagged","sigma_psi",
          "all.weeks","tab","passed","rhats","one.week","ge.true","errors","summ","zero.guided",
          "prior.runs","wk","res","cover.fixed","cover.draws","N.week","N.day","rate","GE",
          "n.weeks","slide.expanded","slide.exact","summary.table","results","rhat.all",
          "n.kept","n.reps","overlap.counts","failed.diagnostics","replicate.summary")
values <- lapply(intersect(keep, ls(e)), function(n) get(n, e))
names(values) <- intersect(keep, ls(e))
values <- values[!vapply(values, is.function, logical(1))]
while (dev.cur() > 1) dev.off()
saveRDS(list(error = error.message, warnings = warning.messages, seconds = elapsed,
             fits = records, values = values, session = sessionInfo()), "result.rds")
cat("VALIDATION COMPLETE; seconds:", elapsed, "; fits:", length(records),
    "; failed direct fits:", sum(!vapply(records, function(x) x$passed, logical(1))),
    "; warnings:", length(warning.messages), "\n")
quit(status = if (is.null(error.message)) 0 else 1)
