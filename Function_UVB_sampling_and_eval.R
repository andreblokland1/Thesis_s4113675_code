options(scipen = 999)

##################################
# UVB SAMPLING PLANNING FUNCTION #
##################################
UVB_planning = function(data, idName = "id", bookValueName = "values", conf_lev = 0.95, precision = 0.02, by = 5, max = 150, start = NULL) {
  # Arguments:
  # data:          data frame with records to be sampled
  # idName:        column name representing IDs
  # bookValueName: column name representing book values
  # conf_lev:      confidence level for the sample
  # precision:     required minimum precision (difference between most likely error rate and upper bound)
  # by:            increment step for sample size (policy decision; default 5)
  # max:           maximum sample size to test (policy decision; default 150)
  # start:         optional fixed starting point for interval selection (for testing/debugging)
  #
  # Note on RNG: the interval starting point is drawn using the active RNG stream
  # of the calling environment. In parallel execution this is the worker's
  # L'Ecuyer-CMRG stream set by clusterSetRNGStream(). No internal set.seed()
  # call is made; seed management is the caller's responsibility.
  
  N.items    = nrow(data)                # Total number of items
  bookValues = data[[bookValueName]]     # Extract book values
  BV         = sum(bookValues)           # Total book value (monetary units)
  IDs        = data[[idName]]            # Extract IDs
  
  sufficient = FALSE  # Flag to indicate whether sampling conditions are met
  usedMax    = FALSE  # Flag if the max sample size was used
  
  # Try increasing sample sizes until a sufficient one is found.
  # Loop starts at 30 to ensure the minimum sample size is never below 30.
  for (n in seq(30, max, by = by)) {
    interval   = BV / n                                    # Calculate the sampling interval
    topStratum = subset(data, bookValues > interval)       # Items greater than interval form the top stratum
    m_seen     = sum(topStratum[[bookValueName]])           # Seen monetary units from top stratum
    
    # Determine interval starting point
    if (!is.null(start)) {
      intervalStartingPoint = start                        # Fixed starting point if provided
    } else {
      intervalStartingPoint = sample.int(interval - 1, size = 1)  # Pick a random interval starting point
    }
    
    # Compute interval selections for MUS
    intervalSelection = intervalStartingPoint + 0:(n - 1) * interval
    index = NULL
    for (i in 1:n) {
      index = c(index, which(intervalSelection[i] < cumsum(bookValues))[1])
    }
    
    selected_data = data[index, ]         # Sampled data rows
    selected_ids  = IDs[index]
    counts        = table(selected_ids)   # Count occurrences of selected IDs
    selected_data = unique(selected_data) # Remove duplicates
    
    bottomStratumSample = selected_data[which(selected_data[[bookValueName]] <= interval), ]
    
    m_seen            = m_seen + sum(bottomStratumSample[[bookValueName]])  # Add bottom stratum seen value
    m_seen_percentage = m_seen / BV                                          # Percentage of population seen
    
    # Beta distribution to evaluate confidence bounds
    a                   = 1 + 0:n
    b                   = 1 + n - 0:n
    v95                 = qbeta(conf_lev, a, b)             # Upper bound using beta distribution
    v                   = ((a - 1) / (a + b - 2))           # Mode of beta distribution
    relativeInaccuracy  = v95 - v
    correctedInaccuracy = precision * (1 / (1 - m_seen_percentage))  # Adjust precision for unseen part
    diff                = relativeInaccuracy - correctedInaccuracy
    
    # Check if all scenarios are within allowed error
    if (all(diff <= 0)) {
      sufficient = TRUE
      break
    }
    
    # If max reached, mark as sufficient anyway
    if (n == max) {
      sufficient = TRUE
      usedMax    = TRUE
    }
  }
  
  # Tag population with selection counts.
  # Fix: use match() to couple counts to rows by name.
  # The previous approach used which(...) %in% names(counts), which returns row
  # indices in population order, then assigned counts in table() order (alphabetical).
  # For IDs like "P_0054_86" vs "P_0054_111" these orderings differ, causing
  # draw counts to be assigned to the wrong items.
  population                = data
  population$draw_count_uvb = 0
  m = match(names(counts), population[[idName]])
  population$draw_count_uvb[m] = as.integer(counts)
  
  # Return results
  result = list(
    n          = n,
    usedMax    = usedMax,   # Indicates whether or not the maximum allowed sample size was used
    start      = intervalStartingPoint,
    population = population
  )
  return(result)
}

##################################
# UVB SAMPLE EVALUATION FUNCTION #
##################################
UVB_evaluation = function(data,
                          idName        = "id",
                          bookValueName = "values",
                          error         = "error",
                          indicatorName = "draw_count_uvb",
                          conf_lev      = 0.95,
                          usedMax       = FALSE) {
  # Arguments:
  # data:          data frame with full population (including draw_count_uvb)
  # idName:        column name representing IDs
  # bookValueName: column name representing book values
  # error:         column name representing true error in the population (known from DGM)
  # indicatorName: column name representing how many times each item is selected
  # conf_lev:      confidence level for Beta upper bound (default 0.95)
  # usedMax:       logical; TRUE if UVB_planning reached the sample size cap (from UVB_planning$usedMax)
  
  # Select sampled rows
  sample = subset(data, data[[indicatorName]] > 0)
  
  # Calculate error rate as error divided by book value (protect against division by zero)
  sample[["error_rate"]] = ifelse(sample[[bookValueName]] == 0, 0, sample[[error]] / sample[[bookValueName]])
  
  # Compute taint: weighted error rate by selection frequency
  sample[["taint"]] = sample[["error_rate"]] * sample[[indicatorName]]
  
  # Sum of known (observed) errors
  known_error = sum(sample[[error]])
  
  # Total number of monetary units selected in the sample
  n = sum(sample[[indicatorName]])
  
  # Total taint in the sample (i.e., weighted sum of error rates)
  total_taint = sum(sample[["taint"]])
  
  # Total value of the items that were not selected in the sample
  unseen_value = sum(data[[bookValueName]][data[[indicatorName]] == 0])
  
  # Total book value of the population
  BV = sum(data[[bookValueName]])
  
  # ------------------------
  # Four estimation methods:
  # ------------------------
  
  # Method 0: Extrapolated error = avg taint * unseen value + known error [STANDARD]
  extr_err    = (total_taint / n) * unseen_value
  tot_est_err = extr_err + known_error
  
  # Method 1: Total error = average taint * total book value
  tot_est_err_1 = (total_taint / n) * BV
  
  # Method 2: Conservative total error = max(method 1, known error)
  tot_est_err_2 = max(tot_est_err_1, known_error)
  
  # Method 3: Apply sample error rate to total book value
  tot_est_err_3 = (known_error / sum(sample[[bookValueName]])) * BV
  
  # True total population error
  true_error = sum(data[[error]])
  
  # Deviations from true error
  dev_est   = true_error - tot_est_err
  dev_est_1 = true_error - tot_est_err_1
  dev_est_2 = true_error - tot_est_err_2
  dev_est_3 = true_error - tot_est_err_3
  
  # Total Estimated Error Rate (TER) for each method
  TER_0 = tot_est_err   / BV
  TER_1 = tot_est_err_1 / BV
  TER_2 = tot_est_err_2 / BV
  TER_3 = tot_est_err_3 / BV
  
  # True Error Rate
  true_error_rate = true_error / BV
  
  # Deviations between estimated and true error rates
  dev_err_rate   = true_error_rate - TER_0
  dev_err_rate_1 = true_error_rate - TER_1
  dev_err_rate_2 = true_error_rate - TER_2
  dev_err_rate_3 = true_error_rate - TER_3
  
  # Bayesian upper bound for the error rate (thesis §2.5)
  # p_upper bounds the error rate in the unaudited portion;
  # the overall UB adds the known errors from the audited items.
  p_upper = qbeta(conf_lev, 1 + total_taint, 1 + n - total_taint)
  UB      = (p_upper * unseen_value + known_error) / BV

  # -------------------------------------------------------
  # Policy rule: b_taint adjustment (thesis §2.5)
  # -------------------------------------------------------
  # When all three conditions are met:
  #   (1) maximum sample size was reached (usedMax)
  #   (2) more than 50% of expenditure was audited
  #   (3) average taint s/n exceeds 5%
  # the extrapolation basis is adjusted upward so that the
  # precision condition is exactly met:
  #   b_taint = p_upper - 2% * BV / (BV - BV_sample)
  # The extrapolation then uses max(b_taint, s/n) instead
  # of s/n alone. This produces a more conservative error
  # estimate (EE_uvb_adj) when the standard estimate would
  # violate the precision requirement.
  # -------------------------------------------------------
  sample_cov = sum(sample[[bookValueName]]) / BV
  avg_taint  = total_taint / n

  b_taint_applied = usedMax && (sample_cov > 0.50) && (avg_taint > 0.05)

  if (b_taint_applied) {
    adjusted_upper_limit = 0.02 * BV / unseen_value
    b_taint_value  = p_upper - adjusted_upper_limit
    extrap_basis   = max(b_taint_value, avg_taint)
    EE_uvb_adj     = known_error + extrap_basis * unseen_value
  } else {
    b_taint_value  = NA_real_
    EE_uvb_adj     = tot_est_err   # same as standard Method 0
  }
  
  # Return results
  result = list(
    EE_uvb           = tot_est_err,       # Method 0 (standard; used in pipeline)
    EE_uvb_1         = tot_est_err_1,     # Method 1
    EE_uvb_2         = tot_est_err_2,     # Method 2
    EE_uvb_3         = tot_est_err_3,     # Method 3
    
    true_error       = true_error,
    
    dev_est          = dev_est,
    dev_est_1        = dev_est_1,
    dev_est_2        = dev_est_2,
    dev_est_3        = dev_est_3,
    
    TER_0            = TER_0,             # Method 0
    TER_1            = TER_1,             # Method 1
    TER_2            = TER_2,             # Method 2
    TER_3            = TER_3,             # Method 3
    
    true             = true_error_rate,
    
    deviance         = dev_err_rate,
    dev_err_rate_1   = dev_err_rate_1,
    dev_err_rate_2   = dev_err_rate_2,
    dev_err_rate_3   = dev_err_rate_3,
    
    sample_coverage  = sample_cov,
    total_taint      = total_taint,
    n                = n,
    BV = BV,
    UB = UB,

    # Policy rule: b_taint adjustment (thesis §2.5)
    EE_uvb_adj       = EE_uvb_adj,       # adjusted EE (= EE_uvb when rule does not trigger)
    b_taint_applied  = b_taint_applied,   # logical: TRUE if all three conditions were met
    b_taint_value    = b_taint_value      # the computed b_taint (NA when rule does not trigger)
  )
  
  return(result)
}