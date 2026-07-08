############################################################
# run_simulation.R
# ----------------------------------------------------------
# Main simulation script for the thesis:
#   "Comparing UVB and ADR audit sampling methods
#    in EU subsidy audits"
#
# Structure:
#   Section A — Parameters (edit before running)
#   Section B — Execution (load data, parallel MC loop)
#   Section C — Analysis (compute performance metrics)
#
# Prerequisites:
#   1. Run create_dbase_and_templates.R first to generate:
#      decile_breaks.rds, pool_by_decile.rds,
#      costitem_db.rds, templates.rds
#   2. All function scripts in the working directory
#
# Parallel execution:
#   Uses a PSOCK cluster (parallel package), compatible with Windows.
#   Iterations are processed in chunks; after each chunk the console
#   reports progress and an intermediate checkpoint is saved.
#   If a single iteration fails, it returns NA for all fields
#   (tryCatch) so one bad draw cannot abort the full run.
#
# Reproducing a single iteration:
#   set.seed(cfg_seed) and run one iteration to get the same
#   results as run_one_iteration_debug.R with the same seed.
#   In the MC loop, each iteration uses set.seed(cfg_seed + i - 1),
#   so iteration 1 always matches the debug script.
############################################################

library(parallel)
options(scipen = 999)

# ==============================================================
# SECTION A — PARAMETERS
# ==============================================================
# Edit these before running the simulation.

# --- Simulation control ---------------------------------------
cfg_seed      = 20260430      # base seed (iteration i uses seed = cfg_seed + i - 1)
cfg_n_sim     = 100     # number of Monte Carlo iterations
cfg_chunk_size = 50      # iterations per chunk; progress is reported after each chunk

# --- Parallel processing --------------------------------------
cfg_n_cores = max(1L, detectCores() - 1L)
# Leave one core free for the OS. Lower this manually if memory
# becomes an issue (each worker holds a full copy of all data objects).

# --- Payment claim generation ---------------------------------
cfg_n_projects  = 135   # number of projects per payment claim
cfg_meanlog     = 14.15    # lognormal mean (log scale) for project book values in bookyear 23-24
cfg_sdlog       = 1.43    # lognormal sd (log scale) for project book values in bookyear 23-24
cfg_k_neighbors = 10     # k for k-NN template selection (thesis §4.2.2)

# --- UVB parameters -------------------------------------------
cfg_conf_lev_uvb = 0.95  # confidence level for UVB Beta distribution
cfg_uvb_by       = 5L    # sample size increment (policy decision)
cfg_uvb_max      = 150L  # maximum sample size (policy decision)

# --- ADR parameters -------------------------------------------
cfg_conf_lev_adr = 0.70  # confidence level for ADR precision (EGESIF Table 3)
                          # 0.70 corresponds to system reliability category 2:
                          # "works; some improvements needed"
cfg_materiality  = 0.02  # 2% materiality threshold (EU Regulation 1303/2013)
cfg_cov_item     = 0.10  # minimum item coverage (10%, category 2)
cfg_cov_exp      = 0.10  # minimum expenditure coverage (10%)

# --- Method selection -----------------------------------------
# Within-project ADR sampling uses dispatch inside run_one_iteration:
#   N_items <= 300 -> within_NSS_PPS
#   N_items >  300 -> within_MUS
# This is not configurable; it follows thesis 2.2.2.
#
# Project selection is configurable:
#   Project selection: adr_project_select_NSS_PPS (adr_project_select_NSS_PPS.R)
#
# Interface requirements:
#   Project-selection fn(data, cov_item, cov_exp)
#     -> list with at least: n_planned, n_exhaustive, n_sample,
#        N_bs, BV_bs_total, SI_bs, exhaustive_ids, sampled_ids




# ==============================================================
# SECTION B — EXECUTION
# ==============================================================

# --- B1. Source function scripts (main process) ---------------
cat("Sourcing function scripts ...\n")

source("Function_generate_payment_claim.R")
source("Function_UVB_sampling_and_eval.R")
source("adr_within_NSS_PPS.R")
source("adr_within_MUS.R")
source("adr_project_select_NSS_PPS.R")
source("run_one_iteration.R")


adr_project_select_fn = adr_project_select_NSS_PPS

# --- B2. Load pre-computed data -------------------------------
cat("Loading pre-computed data ...\n")

decile_breaks  = readRDS("decile_breaks.rds")
pool_by_decile = readRDS("pool_by_decile.rds")
costitem_db    = readRDS("costitem_db.rds")
templates      = readRDS("templates.rds")

cat("  templates:", length(templates),
    " | deciles:", length(decile_breaks),
    " | costitem_db:", nrow(costitem_db), "rows\n\n")

# --- B3. Output fields ----------------------------------------
# create a list of all the statistics we will extract from the simulator. 
# This list will be looped over after each simulation to select the output we will store.
# there is some redundant / unused output that is stored and used for diagnostics

scalar_fields = c(
  # Population
  "BV_pc", "N_pc", "true_error_pc", "TER_true", "materiality",

  # Aim 1: point estimates
  "theta_PC_UVB", "theta_PC_ADR",

  # Aim 2: operational consequences
  "TER", "RTER", "sum_Delta_sample",

  # Aim 3: UB coverage
  "UB_coverage",

  # Decomposition (thesis §4.3)
  "L0", "L1", "L2", "L3", "L4", "L5", "extrap_EE_UVB",

  # ADR selection
  "n_adr_planned", "n_adr_top", "n_adr_bottom", "n_adr_total",
  "N_bs", "BV_exhaustive", "BV_bs_total", "SI_bs",

  # Coverage
  "cov_item_realised", "cov_exp_realised",

  # Sample size diagnostics
  "n_uvb_mean", "n_adr_within_mean", "n_MUS_projects", "n_b_taint_applied",

  # Lognormal parameter estimates (PC book values)
  "meanlog_est", "sdlog_est"
)

# --- B4. Set up PSOCK cluster ---------------------------------
cat("Starting parallel cluster on", cfg_n_cores, "cores ...\n")
cl = makeCluster(cfg_n_cores)

# Source all scripts on every worker
clusterEvalQ(cl, {
  source("Function_generate_payment_claim.R")
  source("Function_UVB_sampling_and_eval.R")
  source("adr_within_NSS_PPS.R")
  source("adr_within_MUS.R")
  source("adr_project_select_NSS_PPS.R")
  source("run_one_iteration.R")
  options(scipen = 999)
  NULL
})

# Export data objects and all config parameters to workers
clusterExport(cl, varlist = c(
  "templates", "decile_breaks", "pool_by_decile", "costitem_db",
  "cfg_seed",
  "cfg_n_projects", "cfg_meanlog", "cfg_sdlog", "cfg_k_neighbors",
  "cfg_conf_lev_uvb", "cfg_uvb_by", "cfg_uvb_max",
  "cfg_conf_lev_adr", "cfg_materiality", "cfg_cov_item", "cfg_cov_exp",
  "adr_project_select_fn",
  "scalar_fields"
))

# --- B5. Worker function --------------------------------------
# Defined here and exported so workers know it.
# Each worker receives iteration index i, sets its own seed,
# runs one iteration, and returns a named numeric vector.
# tryCatch ensures a failed iteration returns NA rather than
# aborting the entire chunk.

run_iter_worker = function(i) {
  # Per-iteration seed: iteration 1 uses cfg_seed, matching the debug script
  set.seed(cfg_seed + i - 1L) # gives each iteration a specific seed for reproducability
  result = tryCatch(
    run_one_iteration(
      templates             = templates,
      decile_breaks         = decile_breaks,
      pool_by_decile        = pool_by_decile,
      n_projects            = cfg_n_projects,
      meanlog               = cfg_meanlog,
      sdlog                 = cfg_sdlog,
      k_neighbors           = cfg_k_neighbors,
      conf_lev_uvb          = cfg_conf_lev_uvb,
      conf_lev_adr          = cfg_conf_lev_adr,
      uvb_by                = cfg_uvb_by,
      uvb_max               = cfg_uvb_max,
      materiality           = cfg_materiality,
      adr_project_select_fn = adr_project_select_fn,
      cov_item              = cfg_cov_item,
      cov_exp               = cfg_cov_exp,
      verbose               = FALSE
    ),
    error = function(e) {
      # Return a named list of NAs so the calling code can still
      # insert this row into result_matrix without crashing.
      warning(sprintf("Iteration %d failed: %s", i, conditionMessage(e)))
      setNames(
        as.list(rep(NA_real_, length(scalar_fields))),
        scalar_fields
      )
    }
  )
  result
}

clusterExport(cl, varlist = "run_iter_worker")

# --- B6. Chunk-based parallel MC loop -------------------------
cat("Starting simulation:", cfg_n_sim, "iterations on", cfg_n_cores, "cores\n")
cat("  chunk size:", cfg_chunk_size,
    "->", ceiling(cfg_n_sim / cfg_chunk_size), "chunks\n")
# cat("  Project select:", adr_project_select_fn,
#    " | within dispatch: NSS_PPS (<=300) / MUS (>300)\n")
cat("  Parameters: n_projects =", cfg_n_projects,
    " | materiality =", cfg_materiality,
    " | conf_lev_adr =", cfg_conf_lev_adr, "\n\n")

# Generate results matrix for all output (defined in scalar_fields) in the columns and 1 row per iteration
result_matrix = matrix(NA_real_, nrow = cfg_n_sim, ncol = length(scalar_fields))
colnames(result_matrix) = scalar_fields

# Build list of chunks: each element is a vector of iteration indices
chunks = split(seq_len(cfg_n_sim),
               ceiling(seq_len(cfg_n_sim) / cfg_chunk_size))

t_start      = Sys.time()
n_failed     = 0L
n_done       = 0L

# Intermediate checkpoint filename (overwritten after every chunk)
checkpoint_file = paste0("checkpoint_", cfg_n_sim, "sim.rds")

for (chunk_idx in seq_along(chunks)) {

  idx = chunks[[chunk_idx]]  # iteration indices in this chunk

  chunk_results = parLapply(cl, idx, run_iter_worker)

  # Insert chunk results into result_matrix
  for (j in seq_along(idx)) {
    i    = idx[[j]]
    iter = chunk_results[[j]]
    # Store numeric scalars
    for (f in scalar_fields) {
      result_matrix[i, f] = iter[[f]] # go to row i, column f and enter the value of list entry f
    }
    if (anyNA(unlist(iter))) n_failed = n_failed + 1L # adds 1 to the counter if there are any NA's in the current iteration
  }

  n_done   = max(idx)
  elapsed  = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  rate     = n_done / elapsed
  eta      = (cfg_n_sim - n_done) / rate

  cat(sprintf("  chunk %d/%d | iter %d/%d done | %.1f iter/s | ETA %.0fs | failed: %d\n",
              chunk_idx, length(chunks),
              n_done, cfg_n_sim,
              rate, eta,
              n_failed))

  # Save checkpoint so partial results survive a crash
  saveRDS(result_matrix, checkpoint_file)
}

stopCluster(cl)

t_end = Sys.time()
cat(sprintf(
  "\nSimulation complete: %d iterations in %.1f seconds (%.1f iter/s) | failed: %d\n\n",
  cfg_n_sim,
  as.numeric(difftime(t_end, t_start, units = "secs")),
  cfg_n_sim / as.numeric(difftime(t_end, t_start, units = "secs")),
  n_failed
))

if (n_failed > 0L)
  warning(sprintf("%d iteration(s) returned NA — inspect result_matrix rows with NA values.", n_failed))

# --- B7. Assemble pc_results ----------------------------------
pc_results           = as.data.frame(result_matrix)
pc_results$iteration = seq_len(cfg_n_sim) # add iteration number

cat("pc_results:", nrow(pc_results), "rows x", ncol(pc_results), "cols\n\n")

# --- B8. Save final results -----------------------------------
results_file = paste0("pc_results_",
                       cfg_n_sim, "sim_",
                       format(Sys.time(), "%Y%m%d_%H%M"),
                       ".rds")
saveRDS(pc_results, results_file)
cat("Results saved to:", results_file, "\n")
cat("Checkpoint file (intermediate):", checkpoint_file, "\n\n")

