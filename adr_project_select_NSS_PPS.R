############################################################
# adr_nss_pps.R
# ----------------------------------------------------------
# ADR Stage 1: project selection from the payment claim
# Non-statistical sampling, probability proportional to size
# Reference: EGESIF_16-0014-01 §6.3.1.3, §6.4
#
# This function performs SELECTION ONLY. Projection (TS/BS
# extrapolation) is handled by run_one_iteration.R using
# the Extrap() helper, because the same extrapolation must
# be applied to multiple quantities (true errors, ADR
# estimates, incremental corrections) for the L0–L5
# decomposition (thesis §4.3).
############################################################

options(scipen = 999)

adr_project_select_NSS_PPS = function(data, cov_item = 0.10, cov_exp = 0.10) {

  # ----------------------------------------------------------
  # Expected columns in data:
  #   id     : project identifier
  #   values : book value (declared expenditure per project)
  # ----------------------------------------------------------

  N           = nrow(data)
  BV          = sum(data$values)
  min_n       = ceiling(cov_item * N)
  min_BV_samp = cov_exp * BV

  # ==========================================================
  # 1. SELECTION — coverage loop
  #    Increase n until expenditure coverage is met
  # ==========================================================

  repeat {
    # --- 1a. Initial stratification: cut-off = BV / n -------
    SI  = BV / min_n
    top = data[data$values > SI, ]
    bot = data[data$values <= SI, ]

    # --- 1b. Iterative refinement (EGESIF 6.3.1.3) ----------
    #     Recalculate SI for bottom stratum; move items that
    #     exceed the new SI to the top stratum until stable.
    repeat {
      n_bs    = min_n - nrow(top)
      if (n_bs <= 0 || nrow(bot) == 0) break # break stops the loop and jumps to the first row after the loop
      SI_bs   = sum(bot$values) / n_bs
      to_move = which(bot$values > SI_bs)
      if (length(to_move) == 0) break
      top = rbind(top, bot[to_move, ])
      bot = bot[-to_move, ]
    }

    n_ts = nrow(top) # final sample size top stratum
    n_bs = min_n - n_ts # final sample size bottom stratum

    if (n_bs <= 0 || nrow(bot) == 0) {
      sel_idx = integer(0)  # no BS selection needed --> sel_idx is used to sample projects from the BS
      break
    }

    # --- 1c. Systematic PPS from bottom stratum --------------
    SI_bs      = sum(bot$values) / n_bs # define bs sampling interval
    bot$cumsum = cumsum(bot$values) # add col that contains the cumulative value of each projects
    start      = runif(1, min = 0, max = SI_bs) # random startingpoint between 0 and the value of the SI
    points     = start + (seq_len(n_bs) - 1) * SI_bs # select points (euros) from startingpoint with SI
    points     = pmin(points, max(bot$cumsum) - 1e-9) # ensure that we dont get outside the range of BV_BS
    sel_idx    = sapply(points, function(p) which(bot$cumsum >= p)[1]) # select the project that contains the sampled euro
    sel_idx    = unique(sel_idx) # ensure that no project can be selected more than once
    bot$cumsum = NULL # drop cumsum column

    # --- 1d. Check expenditure coverage ----------------------
    BV_samp = sum(bot$values[sel_idx]) + sum(top$values) #sum up the BV of sampled projects and check against coverage rate expenditure
    if (BV_samp >= min_BV_samp) break
    min_n = min_n + 1L # increase samplesize by 1
  }

  SI_bs = if (n_bs > 0 && nrow(bot) > 0) sum(bot$values) / n_bs else NA_real_

  # ==========================================================
  # 2. ASSEMBLE AUDITED SET
  # ==========================================================

  # Top stratum: all projects, 100% audited
  top_data = if (nrow(top) > 0) top else data[integer(0), ]

  # Bottom stratum: PPS-selected projects
  bot_selected = if (n_bs > 0 && length(sel_idx) > 0) {
    bot[sel_idx, , drop = FALSE]
  } else {
    data[integer(0), ]
  }

  # ==========================================================
  # 3. COVERAGE DIAGNOSTICS
  # ==========================================================

  BV_audited        = sum(top_data$values) + sum(bot_selected$values)
  cov_exp_realised  = BV_audited / BV
  cov_item_realised = (nrow(top_data) + nrow(bot_selected)) / N

  # ==========================================================
  # 4. OUTPUT
  # ==========================================================
  # Selection results only. Projection (TS/BS extrapolation)
  # is performed in run_one_iteration.R via extrap_TS_BS().

  list(
    # --- Population -----------------------------------------
    N                  = N,
    BV                 = BV,

    # --- Sample design --------------------------------------
    n_planned          = min_n,
    n_exhaustive       = nrow(top_data),
    n_sample           = n_bs,
    n_total_audited    = nrow(top_data) + nrow(bot_selected),
    N_bs               = nrow(bot),
    BV_exhaustive      = sum(top_data$values),
    BV_bs_total        = if (nrow(bot) > 0) sum(bot$values) else 0,
    SI_bs              = SI_bs,

    # --- Coverage -------------------------------------------
    cov_item_realised  = cov_item_realised,
    cov_exp_realised   = cov_exp_realised,

    # --- IDs for traceability -------------------------------
    exhaustive_ids     = top_data$id,
    sampled_ids        = bot_selected$id,
    all_audited_ids    = c(top_data$id, bot_selected$id)
  )
}
