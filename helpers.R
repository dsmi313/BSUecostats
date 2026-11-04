# Shared checks for the JAGS lessons. Keep this file beside the scripts.
# Modern rank-normalized split Rhat and bulk/tail ESS supplement the
# traditional diagnostics printed by jagsUI. Thresholds are initial screens;
# inspect traces and Monte Carlo error for the actual summaries of interest.

# Low GOJ recovery used to generate the dam teaching examples (05--07).
# MY2025 wild season ratios from the talk. These are survival times GOJ
# detection, NOT survival probabilities. They calibrate simulations, not priors.
# Bypass recovery is conditional on return to river, excluding transport.
dam_recovery_defaults <- function() {
  c(phi_S = 123 / 3785, phi_B = 35 / 355)
}

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
    # Pass the functions themselves: measure names given as strings are looked
    # up by name, so a variable called rhat in the workspace would break them.
    posterior::summarise_draws(posterior::as_draws_array(draws),
                              rhat = posterior::rhat, ess_bulk = posterior::ess_bulk,
                              ess_tail = posterior::ess_tail, mcse_median = posterior::mcse_median)
  }, error = function(e) NULL)
}

fit_ok <- function(fit) {
  d <- fit_diagnostics(fit)
  !is.null(d) && nrow(d) > 0L &&
    all(is.finite(as.matrix(d[, c("rhat", "ess_bulk", "ess_tail", "mcse_median")]))) &&
    all(d$rhat <= 1.01 & d$ess_bulk >= 400 & d$ess_tail >= 400)
}

# Explain the same screen used by fit_ok(), retaining the parameter names.
# Accept a precomputed diagnostic table so simulation loops can save it once.
diagnostic_failures <- function(diagnostics) {
  if (is.null(diagnostics) || !nrow(diagnostics)) {
    return(data.frame(variable = NA_character_, rhat = NA_real_, ess_bulk = NA_real_,
                      ess_tail = NA_real_, mcse_median = NA_real_,
                      reason = "diagnostics unavailable", stringsAsFactors = FALSE))
  }
  d <- as.data.frame(diagnostics)
  d$reason <- vapply(seq_len(nrow(d)), function(i) {
    paste(c(
      if (!is.finite(d$rhat[i])) "Rhat unavailable" else if (d$rhat[i] > 1.01) "Rhat > 1.01",
      if (!is.finite(d$ess_bulk[i])) "bulk ESS unavailable" else if (d$ess_bulk[i] < 400) "bulk ESS < 400",
      if (!is.finite(d$ess_tail[i])) "tail ESS unavailable" else if (d$ess_tail[i] < 400) "tail ESS < 400",
      if (!is.finite(d$mcse_median[i])) "median MCSE unavailable"), collapse = "; ")
  }, character(1))
  d[nzchar(d$reason), c("variable", "rhat", "ess_bulk", "ess_tail", "mcse_median", "reason"), drop = FALSE]
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
