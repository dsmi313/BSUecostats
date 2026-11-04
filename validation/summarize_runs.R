# Collect completed validation runs. Later directories replace earlier cases.
# Rscript validation/summarize_runs.R OUTPUT_DIRECTORY RUN_DIRECTORY [RUN_DIRECTORY ...]
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2)
dir.create(args[1], recursive = TRUE, showWarnings = FALSE)
paths <- unlist(lapply(args[-1], list.files, pattern = 'result.rds$', recursive = TRUE, full.names = TRUE))
names(paths) <- basename(dirname(paths))
paths <- paths[!duplicated(names(paths), fromLast = TRUE)]
paths <- paths[order(names(paths))]
status <- diagnostics <- recovery <- failed <- repeats <- summaries <- list()
for (case in names(paths)) {
  path <- paths[[case]]; x <- readRDS(path); v <- x$values
  source <- list.files(dirname(path), '^0.*\\.R$', full.names = TRUE)
  lesson <- source[startsWith(basename(source), paste0(sub('_.*', '', case), '_'))]
  status[[case]] <- data.frame(case = case, seconds = round(x$seconds, 1),
    error = if (is.null(x$error)) '' else x$error, warnings = length(x$warnings),
    direct.fits = length(x$fits), direct.failures = sum(!vapply(x$fits, function(f) f$passed, logical(1))),
    lesson.md5 = if (length(lesson) == 1) unname(tools::md5sum(lesson)) else NA_character_,
    helpers.md5 = unname(tools::md5sum(file.path(dirname(path), 'helpers.R'))))
  if (length(x$fits)) for (i in seq_along(x$fits)) {
    f <- x$fits[[i]]
    summaries[[paste(case, i)]] <- data.frame(case = case, fit = i, model = f$model, passed = f$passed,
      variable = rownames(f$summary), f$summary, row.names = NULL)
    if (!is.null(f$diagnostics)) diagnostics[[paste(case, i)]] <- data.frame(case = case, fit = i,
      model = f$model, passed = f$passed, as.data.frame(f$diagnostics))
  }
  if (startsWith(case, '09_') && !is.null(v$summary.table)) {
    t <- v$summary.table
    widths <- tapply(v$results$hi - v$results$lo, v$results$par, mean)
    t$mean.interval.width <- unname(widths[as.character(t$par)])
    recovery[[case]] <- data.frame(case = case, attempted = v$n.reps, passed = v$n.kept, t, row.names = NULL)
    if (nrow(v$failed.diagnostics)) failed[[case]] <- data.frame(case = case, v$failed.diagnostics)
  }
  if (grepl('_repeats', case, fixed = TRUE) && !is.null(v$passed)) {
    pass <- v$passed
    if (is.matrix(pass)) {
      for (i in seq_len(ncol(pass))) repeats[[paste(case, i)]] <- data.frame(case = case,
        replicate = seq_len(nrow(pass)), model = colnames(pass)[i], fit.passed = pass[, i], paired.passed = apply(pass, 1, all),
        max.rhat = v$rhats[, i])
    } else for (i in seq_len(ncol(v$rhats))) repeats[[paste(case, i)]] <- data.frame(case = case,
      replicate = seq_along(pass), model = colnames(v$rhats)[i], fit.passed = NA, paired.passed = pass,
      max.rhat = v$rhats[, i])
  }
}
write.table.list <- function(x, filename) {
  if (length(x)) write.csv(do.call(rbind, x), file.path(args[1], filename), row.names = FALSE)
}
write.table.list(status, 'run_status.csv')
write.table.list(diagnostics, 'fit_diagnostics.csv')
write.table.list(summaries, 'posterior_summary.csv')
write.table.list(recovery, 'recovery_summary.csv')
write.table.list(failed, 'failed_recovery_diagnostics.csv')
write.table.list(repeats, 'repeat_status.csv')
capture.output(readRDS(paths[[1]])$session, file = file.path(args[1], 'session.txt'))
print(do.call(rbind, status), row.names = FALSE)
