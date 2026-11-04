# Borrowing hatchery GE: pooled vs split by day and night
# (Run time: about 5 seconds. Times vary up to 2x between computers.)
# Run from the repository root: Rscript 07a_day_night_closed_form.R
# Plain R plus ggplot2 >= 3.4.0; no JAGS needed for this example.
#
# This reproduces the talk's "What Splitting Buys" figure. MY2025 PIT counts
# set synthetic arrival numbers, timing, and hatchery day/night guidance.
# They are calibration inputs, not known biological truth. We deliberately
# assign both rears identical day GE, night GE, GRS detection, and nonguided
# GOJ recovery; only their night shares differ. Within each week, detection
# and recovery are constant across day/night and GRS-read/missed fish.
# The optional wild.ge.multiplier exercise relaxes the shared-GE assumption.
#
# Each simulation redraws detections and estimates GE, detection, and timing
# in closed form. GOJ recovery is a fixed simulation input, not estimated.
# We compare COMPLETE borrowing, not the partial pooling in lesson 06:
#   pooled: wild borrows hatchery's pooled weekly GE
#   split: wild borrows hatchery's day/night GE and uses its own night share
#
# Error = true wild GE / borrowed GE - 1, with expected catch held fixed.
# This isolates error from GE borrowing. It omits trap-catch, composition,
# and genetic uncertainty. Synthetic tagged-arrival counts also weight the
# synthetic season total; they are not measured whole-run weekly abundance.
# Points are median errors, bars the middle 90% ACROSS simulated seasons.
# They are not confidence/credible intervals from one fitted data set.
#
# Failed overlap estimates NEVER receive the known simulation probability.
# Both methods use the same valid weeks; season comparisons require all
# weeks to be valid. Failures are reported, including their causes. Summaries
# are conditional on successful estimation, which matters for sparse trials.

# 0. The talk's MY2025 season counts by day and night ----
# Slide "Day and Night GE Are Very Different": wild steelhead cells split by
# LGR passage time. PTAGIS times are PST all year; night is before sunrise or
# after sunset at the dam. c5 fish were never seen at LGR, so they have no
# passage time and only the total is known.
# Time zone matters: reading the times as daylight time instead of PST
# shifts the split by an hour around dawn and dusk and gives about 1.2% by
# day, 23% at night, and 17% passing at night.
diel <- rbind(day   = c(c1 = 54,  c2 = 15, c3 = 3213, c4 = 100, c6 = 15),
              night = c(c1 = 266, c2 = 20, c3 = 449,  c4 = 23,  c6 = 141))
c5 <- 195
stopifnot(all(colSums(diel) == c(320, 35, 3662, 123, 156)))   # 05 section 3b totals
p.season <- sum(diel[, "c4"]) / (sum(diel[, "c4"]) + c5)   # one p, day and night
guided <- diel[, "c1"] + diel[, "c2"] + diel[, "c6"]
nonguided <- (diel[, "c3"] + diel[, "c4"]) / p.season
diel.ge <- guided / (guided + nonguided)
night.share <- unname((guided + nonguided)["night"] / sum(guided + nonguided))
print(round(diel.ge, 3))
cat("Share of wild fish passing at night:", round(night.share, 3), "\n")
stopifnot(isTRUE(all.equal(round(c(diel.ge, night.share), 3), c(0.010, 0.259, 0.160),
                           check.attributes = FALSE)))

# 1. Observed calibration inputs ----
# MY2025 weekly PIT counts for weeks 15 to 22. The day/night convention used
# to split them (PST or daylight time) is not recorded here; section 0 uses PST.
borrow_inputs <- function() {
  wk <- data.frame(
    week = 15:22,
    H.Gd = c(121, 101, 137, 145, 96, 127, 6, 5),
    H.Gn = c(287, 394, 829, 498, 497, 456, 72, 48),
    H.Sd = c(4631, 9373, 17284, 13641, 8478, 3204, 1781, 1844),
    H.Sn = c(742, 665, 2673, 1159, 617, 326, 213, 137),
    H.p = c(.42391, .58363, .63397, .55235, .45013, .41516, .72059, .48571),
    W.Gd = c(8, 3, 3, 9, 2, 16, 2, 2),
    W.Gn = c(19, 40, 23, 83, 48, 83, 27, 28),
    W.Sd = c(122, 226, 303, 732, 425, 327, 418, 453),
    W.Sn = c(32, 28, 60, 132, 51, 47, 48, 32),
    W.p = c(.36364, .42857, .53659, .46154, .26154, .22857, .38889, .39048))
  for (r in c("H", "W")) {
    wk[[paste0(r, ".Nd")]] <- wk[[paste0(r, ".Gd")]] + wk[[paste0(r, ".Sd")]] / wk[[paste0(r, ".p")]]
    wk[[paste0(r, ".Nn")]] <- wk[[paste0(r, ".Gn")]] + wk[[paste0(r, ".Sn")]] / wk[[paste0(r, ".p")]]
  }
  wk$ge.d <- wk$H.Gd / wk$H.Nd
  wk$ge.n <- wk$H.Gn / wk$H.Nn
  wk$p.true <- wk$H.p
  wk$aW <- wk$W.Nn / (wk$W.Nd + wk$W.Nn)
  wk$aH <- wk$H.Nn / (wk$H.Nd + wk$H.Nn)
  wk$NW <- round(wk$W.Nd + wk$W.Nn)
  wk$NH <- round(wk$H.Nd + wk$H.Nn)
  wk
}

# 2. Closed-form estimation from observable counts only ----
# No truth parameter is accepted by this function.
estimate_borrow_rear <- function(Gd, Gn, Sd, Sn, c4, c5) {
  counts <- c(Gd, Gn, Sd, Sn, c4, c5)
  stopifnot(length(counts) == 6L, all(is.finite(counts)),
            all(counts >= 0), all(counts == floor(counts)), c4 <= Sd + Sn)
  fail <- function(reason) list(ok = FALSE, reason = reason,
                                ge.d = NA_real_, ge.n = NA_real_, a = NA_real_, pool = NA_real_)
  if (c4 + c5 == 0) return(fail("no GOJ recoveries"))
  if (c4 == 0) return(fail("no GRS-GOJ overlap"))
  p.hat <- c4 / (c4 + c5)
  Nd.hat <- Gd + Sd / p.hat
  Nn.hat <- Gn + Sn / p.hat
  if (Nd.hat <= 0 || Nn.hat <= 0) return(fail("empty day/night sample"))
  list(ok = TRUE, reason = "ok", ge.d = Gd / Nd.hat, ge.n = Gn / Nn.hat,
       a = Nn.hat / (Nd.hat + Nn.hat), pool = (Gd + Gn) / (Nd.hat + Nn.hat))
}

simulate_borrow_rear <- function(N, a, ge.d, ge.n, p, phi) {
  stopifnot(length(N) == 1L, is.finite(N), N >= 1, N == floor(N),
            all(is.finite(c(a, ge.d, ge.n, p, phi))),
            all(c(a, ge.d, ge.n, p, phi) >= 0), all(c(a, ge.d, ge.n, p, phi) <= 1))
  Nn <- rbinom(1, N, a); Nd <- N - Nn
  Gd <- rbinom(1, Nd, ge.d); Gn <- rbinom(1, Nn, ge.n)
  ngd <- Nd - Gd; ngn <- Nn - Gn
  Sd <- rbinom(1, ngd, p); Sn <- rbinom(1, ngn, p)
  # Independent recovery given nonguided passage: the identifying assumption.
  c4 <- rbinom(1, Sd + Sn, phi)
  c5 <- rbinom(1, ngd + ngn - Sd - Sn, phi)
  estimate_borrow_rear(Gd, Gn, Sd, Sn, c4, c5)
}

# 3. Paired comparisons over synthetic seasons ----
run_borrow_simulation <- function(wk = borrow_inputs(), n.sims = 2000L,
                                  seed = 2025L, phi = 123 / 3785,
                                  wild.ge.multiplier = 1) {
  stopifnot(n.sims >= 1, n.sims == floor(n.sims), nrow(wk) > 0,
            all(wk$NW > 0), all(wk$NH > 0))
  # Fix the RNG algorithms as well as the seed for reproducible figures.
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  stopifnot(length(wild.ge.multiplier) == 1L, is.finite(wild.ge.multiplier),
            wild.ge.multiplier > 0,
            all(wild.ge.multiplier * c(wk$ge.d, wk$ge.n) <= 1))
  truth <- wild.ge.multiplier * ((1 - wk$aW) * wk$ge.d + wk$aW * wk$ge.n)
  stopifnot(all(truth > 0))
  weekly <- vector("list", n.sims)
  season <- data.frame(simulation = seq_len(n.sims), pooled = NA_real_, split = NA_real_)
  for (i in seq_len(n.sims)) {
    out <- data.frame(simulation = i, week = wk$week, pooled = NA_real_, split = NA_real_,
                      valid = FALSE, hatchery.reason = "ok", wild.reason = "ok", pair.reason = "ok")
    for (j in seq_len(nrow(wk))) {
      h <- simulate_borrow_rear(wk$NH[j], wk$aH[j], wk$ge.d[j], wk$ge.n[j], wk$p.true[j], phi)
      w <- simulate_borrow_rear(wk$NW[j], wk$aW[j], wild.ge.multiplier * wk$ge.d[j],
                                wild.ge.multiplier * wk$ge.n[j], wk$p.true[j], phi)
      out$hatchery.reason[j] <- h$reason
      out$wild.reason[j] <- w$reason
      if (!h$ok || !w$ok) { out$pair.reason[j] <- "rear estimation failed"; next }
      borrowed <- c(pooled = h$pool, split = (1 - w$a) * h$ge.d + w$a * h$ge.n)
      if (any(!is.finite(borrowed)) || any(borrowed <= 0)) {
        out$pair.reason[j] <- "zero or nonfinite borrowed GE"
        next
      }
      out[j, c("pooled", "split")] <- as.list(truth[j] / borrowed - 1)
      out$valid[j] <- TRUE
    }
    weekly[[i]] <- out
    if (all(out$valid)) {
      season$pooled[i] <- sum(wk$NW * (1 + out$pooled)) / sum(wk$NW) - 1
      season$split[i] <- sum(wk$NW * (1 + out$split)) / sum(wk$NW) - 1
    }
  }
  list(inputs = wk, weekly = do.call(rbind, weekly), season = season,
       n.sims = n.sims, seed = seed, phi = phi, wild.ge.multiplier = wild.ge.multiplier)
}

summarize_borrow_simulation <- function(result) {
  weekly <- result$weekly
  rows <- list(); k <- 0L
  for (week in result$inputs$week) for (method in c("pooled", "split")) {
    x <- weekly[weekly$week == week & weekly$valid, method]
    q <- if (length(x)) quantile(x, c(.05, .5, .95), names = FALSE) else rep(NA_real_, 3)
    k <- k + 1L
    rows[[k]] <- data.frame(week = week, method = method, n = length(x),
                           lo = q[1], median = q[2], hi = q[3])
  }
  valid <- is.finite(result$season$pooled) & is.finite(result$season$split)
  seasons <- do.call(rbind, lapply(c("pooled", "split"), function(method) {
    x <- result$season[valid, method]
    q <- if (length(x)) quantile(x, c(.05, .5, .95), names = FALSE) else rep(NA_real_, 3)
    data.frame(method = method, n = length(x), lo = q[1], median = q[2], hi = q[3])
  }))
  reasons <- do.call(rbind, lapply(c("hatchery.reason", "wild.reason", "pair.reason"), function(column) {
    tab <- as.data.frame(table(week = weekly$week, reason = weekly[[column]]), stringsAsFactors = FALSE)
    tab$component <- sub(".reason", "", column, fixed = TRUE)
    tab[tab$Freq > 0, c("component", "week", "reason", "Freq")]
  }))
  list(weekly = do.call(rbind, rows), season = seasons, diagnostics = reasons,
       failed.weeks = sum(!weekly$valid), failed.seasons = sum(!valid))
}

# 4. Slide figure: exactly the same summaries as the printed table ----
plot_borrow_simulation <- function(result, summary) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || utils::packageVersion("ggplot2") < "3.4.0")
    stop("Install ggplot2 >= 3.4.0 to create the figure.")
  d <- summary$weekly
  if (any(!is.finite(as.matrix(d[, c("lo", "median", "hi")]))))
    stop("A week has no successful pairs; inspect diagnostics rather than drawing an incomplete figure.")
  d$method <- factor(d$method, levels = c("pooled", "split"),
                     labels = c("Borrow hatchery, pooled", "Borrow hatchery, split"))
  wk <- result$inputs
  labels <- sprintf("night share\nW %.0f%%  H %.0f%%\n\n%d", 100 * wk$aW, 100 * wk$aH, wk$week)
  dodge <- ggplot2::position_dodge(width = .5)
  ggplot2::ggplot(d, ggplot2::aes(factor(week), median, color = method)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = .7, color = "#222222") +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lo, ymax = hi), linewidth = 1.1, position = dodge) +
    ggplot2::geom_point(size = 3.2, position = dodge) +
    ggplot2::scale_color_manual(values = c("Borrow hatchery, pooled" = "#2A78D6",
                                          "Borrow hatchery, split" = "#1BAF7A"), name = NULL) +
    ggplot2::scale_x_discrete(labels = labels) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), breaks = c(0, .5)) +
    ggplot2::labs(x = "ISO week", y = "Error in the week's wild abundance") +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(legend.position = "top", panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 10.5, lineheight = .95),
                   plot.background = ggplot2::element_rect(fill = "white", color = NA))
}

main_borrow_simulation <- function(output.dir = "figures") {
  result <- run_borrow_simulation()
  summary <- summarize_borrow_simulation(result)
  cat("Synthetic seasons calibrated to MY2025 PIT counts; complete borrowing.\n")
  cat("Failed paired weeks:", summary$failed.weeks, "of", nrow(result$weekly), "\n")
  cat("Seasons excluded from BOTH methods:", summary$failed.seasons, "of", result$n.sims, "\n")
  cat("Summaries describe successful pairs. Bars span the middle 90% of simulated errors.\n")
  print(summary$diagnostics[summary$diagnostics$reason != "ok", ], row.names = FALSE)
  print(summary$weekly, row.names = FALSE)
  cat("\nSynthetic season total: median error (middle 90% of successful seasons)\n")
  for (i in seq_len(nrow(summary$season))) {
    s <- summary$season[i, ]
    cat(sprintf("  %-7s %+.1f%% (%+.1f%% to %+.1f%%), n = %d\n",
                s$method, 100 * s$median, 100 * s$lo, 100 * s$hi, s$n))
  }
  g <- plot_borrow_simulation(result, summary)
  print(g)
  dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)
  prefix <- file.path(output.dir, "day_night_borrowing")
  ggplot2::ggsave(paste0(prefix, ".png"), g, width = 12, height = 5.8, dpi = 300, bg = "white")
  # Base R's SVG device avoids an extra svglite dependency.
  ggplot2::ggsave(paste0(prefix, ".svg"), g, device = grDevices::svg,
                  width = 12, height = 5.8, bg = "white")
  write.csv(summary$weekly, paste0(prefix, "_weekly.csv"), row.names = FALSE)
  write.csv(summary$season, paste0(prefix, "_season.csv"), row.names = FALSE)
  write.csv(summary$diagnostics, paste0(prefix, "_diagnostics.csv"), row.names = FALSE)
  capture.output(sessionInfo(), file = paste0(prefix, "_session.txt"))
  invisible(list(result = result, summary = summary, plot = g))
}

# Try this ----
# 1. Set wk <- borrow_inputs(); wk$aW <- wk$aH, then run
#    result <- run_borrow_simulation(wk); summarize_borrow_simulation(result).
#    Does the pooled mismatch disappear when timing is the same?
# 2. Reduce wk$NW tenfold (round and keep at least one fish). Inspect failure
#    counts as well as conditional error summaries. Dropping difficult weeks
#    can make successful-only performance look more favorable.
# 3. Run run_borrow_simulation(wild.ge.multiplier = 0.8). Wild's true day/night
#    GE are now 20% below hatchery's, and wild truth is recomputed accordingly.
#    Does borrowing split GE still remove the bias when the components differ?
#
# Where BPA teaches this ----
# Chapter 7.9: identifiability; Chapter 9: detection histories and latent routes.
# This is a closed-form simulation, not a fit of the production JAGS model.

# Standard source() runs the lesson too. Tests can load definitions only.
if (!isTRUE(getOption("BSUecostats.functions_only", FALSE))) main_borrow_simulation()
