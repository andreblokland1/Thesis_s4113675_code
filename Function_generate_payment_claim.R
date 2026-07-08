############################################################
# Generate a simulated payment claim (ETBR method)
# ----------------------------------------------------------
# Implements thesis §4.2: Data-Generating Mechanism
#   Stage 1: Empirical item pool      (costitem_db, prepared externally)
#   Stage 2: Project BV + template    (lognormal draw + k-NN selection)
#   Stage 3: Template-decile drawing  (draw from same deciles as template)
#   Stage 4: Population finalization  (sum of drawn items = actual BV)
#
# Input:
#   templates      — named list of project templates (from create_dbase_and_templates.R)
#                    each element: $BV_template, $n_items, $costitems (with $template_decile)
#   decile_breaks  — upper bound per decile (from create_dbase_and_templates.R)
#   pool_by_decile — item pool split by decile (from create_dbase_and_templates.R)
#
# Output:
#   Named list of dataframes, one per simulated project.
#   Each dataframe has columns: id, values, error, decile_drawn, project_id
#   Ready as input for UVB_planning / ADR functions.
############################################################

############################
# FUNCTION
############################

generate_payment_claim = function(templates,
                                  decile_breaks,
                                  pool_by_decile,
                                  n_projects,
                                  meanlog,
                                  sdlog,
                                  k_neighbors = 10) {
  # Arguments:
  # templates:      named list of template objects from build_templates()
  #                 each template$costitems must contain column template_decile
  # decile_breaks:  (retained for backward compatibility, no longer used in drawing)
  # pool_by_decile: list of dataframes, one per decile (from create_dbase_and_templates.R)
  # n_projects:     number of projects in the payment claim
  # meanlog:        mean of the log-normal distribution for project book values
  # sdlog:          sd of the log-normal distribution for project book values
  # k_neighbors:    number of nearest neighbors for template selection (thesis §4.2.2)
  
  # --------------------------------------------------------------------------
  # Preparation: template summary for k-NN
  # --------------------------------------------------------------------------
  
  # Template summary table for k-NN distance calculation
  template_ids = names(templates)
  template_bvs = vapply(templates, function(x) x$BV_template, numeric(1))
  
  # --------------------------------------------------------------------------
  # Stage 2: Draw project book values from lognormal
  # --------------------------------------------------------------------------
  
  bv_targets = rlnorm(n = n_projects, meanlog = meanlog, sdlog = sdlog)
  
  # --------------------------------------------------------------------------
  # Build each project
  # --------------------------------------------------------------------------
  
  projects     = vector("list", n_projects)
  gen_metadata = vector("list", n_projects)
  names(projects)     = paste0("P_", sprintf("%04d", seq_len(n_projects)))
  names(gen_metadata) = paste0("P_", sprintf("%04d", seq_len(n_projects)))
  
  for (p in seq_len(n_projects)) {
    
    bv_target = bv_targets[p]
    project_id = names(projects)[p]
    
    # Stage 2 (cont): k-NN template selection
    distances = abs(template_bvs - bv_target)
    k = min(k_neighbors, length(template_ids))
    nearest_idx = order(distances)[1:k]
    chosen_idx = nearest_idx[sample.int(k, 1)]
    template = templates[[chosen_idx]]
    
    # Stage 3: Draw items from the same deciles as the template
    template_deciles = template$costitems$template_decile
    n_items = length(template_deciles)
    drawn_values  = numeric(n_items)
    drawn_errors  = numeric(n_items)
    drawn_ids     = character(n_items)
    drawn_deciles = integer(n_items)
    
    for (j in seq_len(n_items)) {
      
      decile_j = template_deciles[j]
      
      # Draw a random item from that decile
      pool = pool_by_decile[[as.character(decile_j)]]
      row_idx = sample.int(nrow(pool), 1)
      
      drawn_values[j]  = pool$values[row_idx]
      drawn_errors[j]  = pool$error[row_idx]
      drawn_ids[j]     = paste0(project_id, "_", j)
      drawn_deciles[j] = decile_j
    }
    
    # Stage 4: Finalization — use actual drawn values (BV may differ from target)
    projects[[project_id]] = data.frame(
      id            = drawn_ids,
      values        = drawn_values,
      error         = drawn_errors,
      decile_drawn  = drawn_deciles,
      project_id    = rep(project_id, n_items),
      stringsAsFactors = FALSE
    )
    
    # Stage 4 (cont): store generation metadata
    decile_counts = tabulate(drawn_deciles, nbins = length(pool_by_decile))
    gen_metadata[[project_id]] = list(
      bv_target     = bv_target,
      template_id   = template_ids[[chosen_idx]],
      template_bv   = template$BV_template,
      decile_counts = decile_counts        # integer vector, length = n_deciles
    )
  }
  
  list(projects = projects, gen_metadata = gen_metadata)
}