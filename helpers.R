# Shared checks for the JAGS lessons. Keep this file beside the scripts.
# Modern rank-normalized split Rhat and bulk/tail ESS supplement the
# traditional diagnostics printed by jagsUI. Thresholds are initial screens;
# inspect traces and Monte Carlo error for the actual summaries of interest.

fit_diagnostics <- function(fit) {
  if (is.null(fit) || is.null(fit$samples)) return(NULL)
  tryCatch({
    chains <- lapply(fit$samples, as.matrix)
    if (length(chains) < 2L || any(vapply(chains, nrow, integer(1)) < 4L))
      return(NULL)
    if (!all(vapply(chains, function(x) identical(dim(x), dim(chains[[1]])), logical(1))))
      return(NULL)
    draws <- array(NA_real_, c(nrow(chains[[1]]), length(chains), ncol(chains[[1]])),
                   dimnames = list(NULL, NULL, colnames(chains[[1]])))
    for (i in seq_along(chains)) draws[, i, ] <- chains[[i]]
    if (any(!is.finite(draws))) return(NULL)
    # Inverse GE is reported in lesson 07: check its sampling precision too.
    if ("wild_ge" %in% dimnames(draws)[[3]]) {
      inverse <- 1 / draws[, , "wild_ge"]
      if (any(!is.finite(inverse))) return(NULL)
      extended <- array(NA_real_, dim(draws) + c(0L, 0L, 1L),
                        dimnames = list(NULL, NULL, c(dimnames(draws)[[3]], "inverse_ge")))
      extended[, , seq_len(dim(draws)[3])] <- draws
      extended[, , dim(extended)[3]] <- inverse
      draws <- extended
    }
    posterior::summarise_draws(posterior::as_draws_array(draws),
                              "rhat", "ess_bulk", "ess_tail", "mcse_median")
  }, error = function(e) NULL)
}

fit_ok <- function(fit) {
  d <- fit_diagnostics(fit)
  !is.null(d) && nrow(d) > 0L &&
    all(is.finite(as.matrix(d[, c("rhat", "ess_bulk", "ess_tail", "mcse_median")]))) &&
    all(d$rhat <= 1.01 & d$ess_bulk >= 400 & d$ess_tail >= 400)
}

fit_rhat <- function(fit) {
  d <- fit_diagnostics(fit)
  if (is.null(d) || !nrow(d) || any(!is.finite(d$rhat))) return(Inf)
  max(d$rhat)
}

# Simulation loops retry NULL just as they retry poor mixing. Keep and count
# the failure instead of aborting a study or treating absent diagnostics as OK.
safe_jags <- function(...) {
  tryCatch(jagsUI::jags(...), error = function(e) {
    message("Fit failed: ", conditionMessage(e))
    NULL
  })
}
