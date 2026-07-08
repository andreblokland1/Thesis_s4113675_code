############################################################
# adr_within_MUS.R
# ----------------------------------------------------------
# ADR within-project sampling — Monetary Unit Sampling
# Reference: EGESIF_16-0014-01 §6.3.1 + §7.6.3
############################################################

options(scipen = 999)

within_MUS = function(data, conf_lev, AE = 0) {
  N   = nrow(data)
  BV  = sum(data$values)
  z   = qnorm(conf_lev)
  TE  = 0.02 * BV

  # sigma_r from UVB sample items (adjusted error rates, EGESIF fn. 27)
  uvb_samp              = data[data$draw_count_uvb > 0, ]
  SI_init               = BV / max(30L, nrow(uvb_samp))
  uvb_samp$err_rate_adj = ifelse(
    uvb_samp$values > SI_init,
    uvb_samp$error / SI_init,
    uvb_samp$error / uvb_samp$values
  )
  sigma = if (nrow(uvb_samp) > 1) sd(uvb_samp$err_rate_adj) else 0

  n = max(30L, ceiling(((z * BV * sigma) / (TE - AE))^2))
  n = min(n, N)

  # Stratify into top/bottom stratum
  SI  = BV / n
  top = data[data$values > SI, ]
  bot = data[data$values <= SI, ]

  repeat {
    n_bs    = n - nrow(top)
    if (n_bs <= 0 || nrow(bot) == 0) break
    SI_bs   = sum(bot$values) / n_bs
    to_move = which(bot$values > SI_bs)
    if (length(to_move) == 0) break
    top = rbind(top, bot[to_move, ])
    bot = bot[-to_move, ]
  }

  n_ts = nrow(top)
  n_bs = n - n_ts

  if (n_bs <= 0 || nrow(bot) == 0) {
    # All items in TS; full coverage, no sampling uncertainty
    return(list(
      EE_project = sum(top$error),
      BV_project = BV,
      TER        = sum(top$error) / BV,
      n          = n_ts,
      cov_item   = n_ts / N,
      cov_exp    = 1
    ))
  }

  # Systematic PPS selection from bottom stratum
  SI_bs      = sum(bot$values) / n_bs
  bot$cumsum = cumsum(bot$values)
  start      = runif(1, 0, SI_bs)
  points     = start + (seq_len(n_bs) - 1) * SI_bs
  points     = pmin(points, max(bot$cumsum) - 1e-9)
  sel_idx    = sapply(points, function(p) which(bot$cumsum >= p)[1])
  sel_idx    = unique(sel_idx)
  bot$cumsum = NULL
  samp_bot   = bot[sel_idx, ]

  # Evaluation: §6.3.1.4
  EE_ts      = sum(top$error)
  EE_bs      = sum(samp_bot$error / samp_bot$values) * SI_bs
  EE_project = EE_ts + EE_bs

  BV_samp = sum(samp_bot$values) + sum(top$values)

  list(
    EE_project = EE_project,
    BV_project = BV,
    TER        = EE_project / BV,
    n          = n,
    cov_item   = n / N,
    cov_exp    = BV_samp / BV
  )
}
