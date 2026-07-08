############################################################
# run_one_iteration.R
# ----------------------------------------------------------
# Runs ONE complete iteration of the simulation pipeline:
#   Step 1: Generate payment claim (PC)
#   Step 2: UVB sampling + evaluation (all projects)
#   Step 3: ADR within-project sampling (all projects)
#           Dispatch: NSS PPS if N <= 300, MUS if N > 300
#   Step 4: Build project-level summary table
#   Step 5: ADR project selection (NSS PPS, selection only)
#   Step 6: Incremental corrections and L0–L5 decomposition
#   Step 7: Performance measures (Aims 1–3)
#
# Returns: a named list with all PC-level results needed for
#          the analysis scripts.
#
# Dependencies (must be sourced before calling):
#   Function_generate_payment_claim.R
#   Function_UVB_sampling_and_eval.R
#   adr_within_NSS_PPS.R
#   adr_within_MUS.R
#   adr_project_select_NSS_PPS.R
############################################################


# ===========================================================
# Extrap helper: TS/BS extrapolation (thesis 4.3)
# ===========================================================
# Generic implementation of the Extrap() operator defined in
# thesis 4.3. Applies the two-stratum PPS extrapolation to
# an arbitrary project-level quantity q_j. Used for L1–L4.
#
# Returns the raw extrapolated TOTAL (not divided by BV_PC).
# The caller divides by the appropriate denominator.
#
# Arguments:
#   project_table : data.frame with columns 'id' and 'values'
#   q_col         : name of the column in project_table to extrapolate (this varies based on the level of the decomposition)
#   ts_ids        : character vector of exhaustive (TS) project ids
#   bs_ids        : character vector of sampled BS project ids
#   SI_bs         : sampling interval for the BS stratum
# -----------------------------------------------------------

extrap_TS_BS = function(project_table, q_col, ts_ids, bs_ids, SI_bs) {

  # Exhaustive stratum: direct sum of q_j
  if (length(ts_ids) > 0) {
    ts_rows = project_table[match(ts_ids, project_table$id), ] # reconstruct TS based on ts_ids and project id
    EE_ts   = sum(ts_rows[[q_col]])
  } else {
    EE_ts = 0
  }

  # Non-exhaustive stratum: taint-based extrapolation
  # EE_bs = SI_bs * sum(q_j / BV_j) for sampled BS projects
  if (length(bs_ids) > 0 && !is.na(SI_bs)) {
    bs_rows = project_table[match(bs_ids, project_table$id), ]
    taints  = bs_rows[[q_col]] / bs_rows$values
    EE_bs   = SI_bs * sum(taints)
  } else {
    EE_bs = 0
  }

  EE_ts + EE_bs
}


# ===========================================================
# Main function
# ===========================================================

run_one_iteration = function(templates,
                             decile_breaks,
                             pool_by_decile,
                             n_projects,
                             meanlog,
                             sdlog,
                             k_neighbors,
                             conf_lev_uvb,
                             conf_lev_adr,
                             uvb_by,
                             uvb_max,
                             materiality,
                             adr_project_select_fn,
                             cov_item    = 0.10,
                             cov_exp     = 0.10,
                             MUS_threshold = 300L,
                             verbose     = FALSE) {

  # Arguments:
  # templates, decile_breaks, pool_by_decile : from create_dbase_and_templates.R (via readRDS)
  # n_projects      : number of projects in the payment claim
  # meanlog, sdlog  : lognormal parameters for project book values
  # k_neighbors     : k for k-NN template selection
  # conf_lev_uvb    : confidence level for UVB Beta-distribution evaluation
  # conf_lev_adr    : confidence level for ADR MUS (thesis 2.4)
  # uvb_by          : UVB sample size increment (policy; default 5)
  # uvb_max         : UVB maximum sample size (policy; default 150)
  # materiality     : materiality threshold (default 0.02)
  # adr_project_select_fn : function for ADR project selection (e.g. adr_project_select_NSS_PPS)
  # cov_item        : minimum item coverage for ADR (default 0.10)
  # cov_exp         : minimum expenditure coverage for ADR (default 0.10)
  # MUS_threshold   : number of cost items above which ADR uses MUS instead of
  #                   NSS PPS for within-project sampling (thesis 2.2.2; default 300)
  # verbose         : if TRUE, print progress to console
  #
  # Note: within-project ADR sampling uses dispatch:
  #   N_items <= MUS_threshold -> within_NSS_PPS (must be in environment)
  #   N_items >  MUS_threshold -> within_MUS     (must be in environment)

  # ============================================================
  # Step 1: Generate payment claim
  # ============================================================
  if (verbose) cat("  Step 1: Generating payment claim ...\n")

  pc = generate_payment_claim(
    templates      = templates,
    decile_breaks  = decile_breaks,
    pool_by_decile = pool_by_decile,
    n_projects     = n_projects,
    meanlog        = meanlog,
    sdlog          = sdlog,
    k_neighbors    = k_neighbors
  )

  # ============================================================
  # Step 2: UVB — audit every project in the PC
  # ============================================================
  if (verbose) cat("  Step 2: UVB sampling and evaluation ...\n")

  uvb_results = lapply(pc$projects, function(proj_df) { #apply planning and eval function to all projects
    plan = UVB_planning(
      data     = proj_df,
      conf_lev = conf_lev_uvb,
      by       = uvb_by,
      max      = uvb_max
    )
    eval = UVB_evaluation(
      data     = plan$population,
      conf_lev = conf_lev_uvb,
      usedMax  = plan$usedMax
    )
    list(
      population      = plan$population,
      EE_uvb          = eval$EE_uvb,
      EE_uvb_adj      = eval$EE_uvb_adj,
      b_taint_applied = eval$b_taint_applied,
      true_error      = eval$true_error,
      BV              = sum(proj_df$values),
      n_uvb           = plan$n,
      UB              = eval$UB
    )
  })

  # ============================================================
  # Step 3: ADR within-project sampling — ALL N projects
  #         Dispatch: NSS PPS if N <= MUS_threshold,
  #                   MUS     if N >  MUS_threshold
  #         (thesis §2.2.2 / §2.2.3)
  # ============================================================
  if (verbose) cat("  Step 3: ADR within-project sampling (all projects) ...\n")

  proj_ids = names(uvb_results)

  adr_within_results = vector("list", length(proj_ids))
  names(adr_within_results) = proj_ids

  for (k in seq_along(proj_ids)) {
    pid             = proj_ids[k]
    proj_population = uvb_results[[pid]]$population
    N_items         = nrow(proj_population)

    # Dispatch based on number of cost items (thesis 2.2.2)
    if (N_items > MUS_threshold) {
      within = within_MUS(
        data     = proj_population,
        conf_lev = conf_lev_adr
      )
    } else {
      within = within_NSS_PPS(
        data     = proj_population,
        cov_item = cov_item,
        cov_exp  = cov_exp
      )
    }

    adr_within_results[[k]] = list(
      project_id = pid,
      EE_project = within$EE_project,
      n_within   = within$n,
      BV_project = sum(proj_population$values),
      method     = if (N_items > MUS_threshold) "MUS" else "NSS_PPS"
    )
  }

  # ============================================================
  # Step 4: Build project-level summary table
  # ============================================================
  if (verbose) cat("  Step 4: Building project_table ...\n")

  project_table = data.frame(
    id              = proj_ids,
    values          = sapply(uvb_results, function(x) x$BV),
    EE_uvb          = sapply(uvb_results, function(x) x$EE_uvb),
    EE_uvb_adj      = sapply(uvb_results, function(x) x$EE_uvb_adj),
    b_taint_applied = sapply(uvb_results, function(x) x$b_taint_applied),
    true_error      = sapply(uvb_results, function(x) x$true_error),
    EE_adr          = sapply(adr_within_results, function(x) x$EE_project),
    n_uvb           = sapply(uvb_results, function(x) x$n_uvb),
    n_adr_within    = sapply(adr_within_results, function(x) x$n_within),
    UB              = sapply(uvb_results, function(x) x$UB),
    adr_method      = sapply(adr_within_results, function(x) x$method),
    stringsAsFactors = FALSE
  )
  rownames(project_table) = NULL

  # --- MLE estimates of the lognormal parameters of the PC ---------
  # Fit lognormal to the N_pc project book values in this iteration.
  # MLE: meanlog_est = mean(log(BV_j)), sdlog_est = sd(log(BV_j))
  bv_log      = log(project_table$values)
  meanlog_est = mean(bv_log)
  sdlog_est   = sd(bv_log)

  # ============================================================
  # Step 5: ADR project selection (NSS PPS — selection only)
  # ============================================================
  if (verbose) cat("  Step 5: ADR project selection ...\n")

  adr_proj = adr_project_select_fn(
    data     = project_table,
    cov_item = cov_item,
    cov_exp  = cov_exp
  )

  ts_ids   = adr_proj$exhaustive_ids
  bs_ids   = adr_proj$sampled_ids
  SI_bs    = adr_proj$SI_bs
  n_top    = length(ts_ids)
  n_bottom = length(bs_ids)

  # ============================================================
  # Step 6: Incremental corrections and L0–L5 decomposition
  #         (thesis §2.3 and §4.3)
  # ============================================================
  if (verbose) cat("  Step 6: Decomposition (L0-L5) ...\n")

  BV_pc         = sum(project_table$values)
  true_error_pc = sum(project_table$true_error)

  # --- 6a. Derived columns for decomposition ------------------
  # These are computed for ALL projects; Extrap will select the
  # relevant subset using ts_ids and bs_ids.

  # Difference between ADR and UVB estimates (can be negative)
  project_table$diff_adr_uvb = project_table$EE_adr - project_table$EE_uvb

  # Incremental correction: floored at zero (thesis §2.3)
  project_table$Delta = pmax(project_table$EE_adr - project_table$EE_uvb, 0)

  # --- 6b. L0: truth (no sampling, no extrapolation) ----------
  L0 = true_error_pc / BV_pc

  # --- 6c. L1: stage-1 extrapolation on true errors -----------
  # Isolates the effect of the first-stage project selection.
  L1_raw = extrap_TS_BS(project_table, "true_error", ts_ids, bs_ids, SI_bs)
  L1     = L1_raw / BV_pc
  
  # stage-1 extrapolation on UVB's per-project estimates
  # Needed to compute the new residual-scale L1 in the analysis.
  extrap_EE_UVB_raw = extrap_TS_BS(project_table, "EE_uvb", ts_ids, bs_ids, SI_bs)
  extrap_EE_UVB     = extrap_EE_UVB_raw / BV_pc

  # --- 6d. L2: within-project audit (= theta_PC_ADR) ---------
  # True errors replaced by ADR audit estimates.
  L2_raw = extrap_TS_BS(project_table, "EE_adr", ts_ids, bs_ids, SI_bs)
  L2     = L2_raw / BV_pc

  # --- 6e. L3: subtraction of UVB's estimate ------------------
  # Per-project difference, can be negative. No floor.
  L3_raw = extrap_TS_BS(project_table, "diff_adr_uvb", ts_ids, bs_ids, SI_bs)
  L3     = L3_raw / BV_pc

  # --- 6f. L4: incremental floor (= TER) ---------------------
  # Negative differences truncated at zero.
  L4_raw = extrap_TS_BS(project_table, "Delta", ts_ids, bs_ids, SI_bs)
  L4     = L4_raw / BV_pc      # this is the TER (thesis §2.3)

  # --- 6g. L5: denominator shrinkage (= RTER) ----------------
  # Corrections applied to both numerator and denominator for the corrections that ADR applied.
  all_selected_ids = c(ts_ids, bs_ids)
  sum_Delta_sample = sum(project_table$Delta[project_table$id %in% all_selected_ids]) # all corrections: sum of Delta_j for all selected projects

  # RTER = (EE_PC - sum_Delta_sample) / (BV_PC - sum_Delta_sample)
  # where EE_PC = L4_raw (the raw extrapolated total based on Delta_j)
  denom_RTER = BV_pc - sum_Delta_sample
  L5 = if (denom_RTER > 0) (L4_raw - sum_Delta_sample) / denom_RTER else NA_real_ # built-in a safeguard

  # ============================================================
  # Step 7: Performance measures (Aims 1–3)
  # ============================================================
  if (verbose) cat("  Step 7: Performance measures ...\n")

  # --- Aim 1: point-estimate accuracy -------------------------

  # UVB
  theta_PC_UVB = sum(project_table$EE_uvb) / BV_pc

  # ADR: stage-1 extrapolation on absolute projected errors (= L2)
  theta_PC_ADR = L2

  # --- Aim 2: operational consequences ------------------------

  # TER and RTER (already computed as L4 and L5)
  TER  = L4
  RTER = L5

  # --- Aim 3: UB coverage (project level) ---------------------
  # Proportion of projects where UVB's 95% upper bound covers
  # the true error rate (thesis 4.5, Aim 3).
  project_table$true_error_rate = project_table$true_error / project_table$values
  UB_coverage = mean(project_table$true_error_rate <= project_table$UB)

  # --- Diagnostics --------------------------------------------
  cov_exp_realised  = adr_proj$cov_exp_realised
  cov_item_realised = adr_proj$cov_item_realised

  # ============================================================
  # Return
  # ============================================================

  list(
    # --- Population -------------------------------------------
    BV_pc             = BV_pc,
    N_pc              = nrow(project_table),
    true_error_pc     = true_error_pc,
    TER_true          = L0,
    materiality       = materiality,

    # --- Aim 1: point estimates -------------------------------
    theta_PC_UVB      = theta_PC_UVB,
    theta_PC_ADR      = theta_PC_ADR,

    # --- Aim 2: operational consequences ----------------------
    TER               = TER,
    RTER              = RTER,
    sum_Delta_sample  = sum_Delta_sample,

    # --- Aim 3: UB coverage -----------------------------------
    UB_coverage       = UB_coverage,

    # --- Decomposition (thesis 4.3) --------------------------
    L0                = L0,
    L1                = L1,
    L2                = L2,
    L3                = L3,
    L4                = L4,
    L5                = L5,
    extrap_EE_UVB     = extrap_EE_UVB,

    # --- ADR selection ----------------------------------------
    n_adr_planned     = adr_proj$n_planned,
    n_adr_top         = n_top,
    n_adr_bottom      = n_bottom,
    n_adr_total       = n_top + n_bottom,
    N_bs              = adr_proj$N_bs,
    BV_exhaustive     = adr_proj$BV_exhaustive,
    BV_bs_total       = adr_proj$BV_bs_total,
    SI_bs             = SI_bs,

    # --- Coverage ---------------------------------------------
    cov_item_realised = cov_item_realised,
    cov_exp_realised  = cov_exp_realised,

    # --- Lognormal parameter estimates (PC book values) ------
    meanlog_est       = meanlog_est,
    sdlog_est         = sdlog_est,

    # --- Sample size diagnostics ------------------------------
    n_uvb_mean        = mean(project_table$n_uvb),
    n_adr_within_mean = mean(project_table$n_adr_within[project_table$id %in% all_selected_ids]),
    n_MUS_projects    = sum(project_table$adr_method == "MUS"),
    n_b_taint_applied = sum(project_table$b_taint_applied),

    # --- Project-level detail (one row per project) -----------
    project_table     = project_table
  )
}
