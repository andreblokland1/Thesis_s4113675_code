############################################################
# adr_within_NSS_PPS.R
# ----------------------------------------------------------
# ADR within-project sampling — Non-statistical, probability proportional to size
# Reference: EGESIF_16-0014-01 §6.4 + §7.6.3
############################################################

options(scipen = 999)


# guidance notes
# In summary, non-statistical sampling is considered appropriate for cases where 
# it is not possible to achieve an adequate sample size that would be required to 
# support statistical sampling. It is not possible to state the exact population 
# size below which non-statistical sampling is needed as it depends on several 
# population characteristics, but usually this threshold is somewhere between 
# 50 and 150 sampling units. The final decision should of course take into 
# consideration the balance between the cost and benefit associated with each of the methods.

# For 2014-2020, the regulation also sets criteria to be respected when non-statistical
# sampling is applied, namely to cover a minimum of 5% operations and 10% of the
# expenditure declared (Article 127(1) CPR).
# this is also linked to system audit score. In the case of the Netherlands that score
# puts us in category 2: the system works. Some improvement needed which leads to
# cov_item: 5-10% (I will assume 10% in the simulation) and cov_exp: 10% 
# (guideline_sampling_method page 148-149)

# since both coverage of items and expenditure need to be met, sampling proportional to size
# is the most efficient way to achieve both. Therefor I will use this as the default in the simulator


within_NSS_PPS = function(data, cov_item = 0.10, cov_exp = 0.10, conf_lev = NULL) {
  N           = nrow(data)
  BV          = sum(data$values)

  # If the project contains 30 or fewer cost items, all items are
  # audited exhaustively (thesis §2.2.2, Article 28(9) CDR 480/2014).
  if (N <= 30L) {
    return(list(
      EE_project = sum(data$error),
      BV_project = BV,
      TER        = sum(data$error) / BV,
      n          = N,
      cov_item   = 1,
      cov_exp    = 1
    ))
  }

  # Minimum sample size is 30 items (Article 28(9) CDR 480/2014),
  # or the coverage-based minimum if that is larger.
  min_n       = max(30L, ceiling(cov_item * N))
  min_BV_samp = cov_exp * BV

  repeat {
    SI  = BV / min_n
    top = data[data$values > SI, ] #top stratum
    bot = data[data$values <= SI, ] # bot stratum

    # Refinement --> adjusting top and bottom stratum based on intervals
    repeat {
      n_bs    = min_n - nrow(top)
      if (n_bs <= 0 || nrow(bot) == 0) break
      SI_bs   = sum(bot$values) / n_bs
      to_move = which(bot$values > SI_bs) #items in bs larger than bs interval are moved to top stratum
      if (length(to_move) == 0) break
      top = rbind(top, bot[to_move, ])
      bot = bot[-to_move, ]
    }

    n_ts = nrow(top)
    n_bs = min_n - n_ts

    if (n_bs <= 0 || nrow(bot) == 0) {
      BV_samp = sum(top$values)  # everything is in TS; full coverage
      break
    }

    SI_bs      = sum(bot$values) / n_bs
    bot$cumsum = cumsum(bot$values) #add column with cumulative values
    start      = runif(1, 0, SI_bs) # random starting point
    points     = start + (seq_len(n_bs) - 1) * SI_bs #take a sample from the start point and then add interval for the next sampled euro
    points     = pmin(points, max(bot$cumsum) - 1e-9) #vector of point value. takes the actual value of the point unless it is larger than cumsum of bott stratum. then it takes bot$cumsum - 1e-9
    sel_idx    = sapply(points, function(p) which(bot$cumsum >= p)[1]) # selection index is the first point where cumsum is larger than p
    sel_idx    = unique(sel_idx) # for stability. in theory there are no duplicates
    bot$cumsum = NULL #drop helper col
    samp_bot   = bot[sel_idx, ] # store sample in samp_bot df

    BV_samp = sum(samp_bot$values) + sum(top$values)
    if (BV_samp >= min_BV_samp) break # if sample is larger than cov_exp break
    min_n = min_n + 1L # otherwise increase sample size and start over
  }

  # Evaluation: 6.4.5.3 (page 151)
  EE_ts = sum(top$error)
  if (n_bs > 0 && nrow(bot) > 0) { #check if bs exists and is sampled
    SI_bs = sum(bot$values) / n_bs
    EE_bs = sum(samp_bot$error / samp_bot$values) * SI_bs
  } else {
    EE_bs = 0
  }
  EE_project = EE_ts + EE_bs

  list( #return list as result with error, book value of the project,  sample size, item coverage and expenditure coverage
    EE_project = EE_project,
    BV_project = BV,
    TER        = EE_project/BV,
    n          = min_n,
    cov_item   = min_n/N,
    cov_exp    = BV_samp/BV
  )
}


