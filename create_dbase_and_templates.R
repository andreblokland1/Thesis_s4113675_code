############################################################
# Build empirical item pool and project templates
# ----------------------------------------------------------
# Run once. Produces two objects used by the simulation:
#
# 1. costitem_db: the empirical item pool for decile lookup
#    Contains ONLY items that were in a UVB sample (Steekproef >= 1),
#    because only for those items the true error is known.
#    columns: project_id, id, values, error, global_decile
#
# 2. templates: structural blueprints based on ALL items within each project
#    Each template represents the  structure of a historical project
#    (including items that were not sampled), so that proportions
#    reflect the true internal composition.
#
# Input: Excel file with (amongst others:) columns Projectnummer, Code, Subtotaal, SOLL, Steekproef
############################################################

library(dplyr)
library(readxl)
library(ggplot2)

############################
# PARAMETERS
############################


EXCEL_PATH   = "BA23-25.xlsx"
SHEET_NAME   = "Sheet1"
VALUE_FLOOR  = 0.01 # safeguard against negative item values

# Decile boundaries: 9 monetary cut points that define 10 deciles.
# Intervals are right-closed: (lower, upper], so an item sits in the first
# decile whose upper bound is >= its value.
# Adjust these values based on inspection of the data:
#   run quantile(sampled_items$values, probs = seq(.1, .9, .1)) on a first pass,
#   then nudge any boundary that causes an empty decile.
DECILE_CUTS = c(1719, 1720, 1741, 1788, 3000, 9834.102, 15000, 20161.092, 29285.678)  # <<< fill in 9 values based on your data

# Helper: assign a monetary value to a decile using the fixed cut points above.
assign_to_decile = function(bv, decile_cuts = DECILE_CUTS) {
  as.integer(cut(bv,
                 breaks = c(-Inf, sort(decile_cuts), Inf),
                 labels = seq_len(length(decile_cuts) + 1),
                 right  = TRUE,
                 include.lowest = TRUE))
}



############################
# READ, PREPARE DATA AND CREATE DATABASE
############################

df = read_excel(EXCEL_PATH, sheet = SHEET_NAME)

# check the number of rows with the excel file
count_steekproef = df %>%
  filter(Steekproef > 0) %>%
  nrow()

# Base dataframe: all items, standardized column names
all_items = df |> # used for templatebuilding
  rename(
    project_id = project_id,
    id         = Code,
    values     = Subtotaal,
    soll       = SOLL,
    sampled    = Steekproef
  ) |>
  mutate(
    project_id = as.character(project_id), # enforce right format
    id         = as.character(id),
    values     = as.numeric(values),
    soll       = as.numeric(soll),
    sampled    = as.numeric(sampled),
    error      = values - soll
  ) |>
  filter(is.finite(values), values > VALUE_FLOOR)

# --------------------------------------------------------------------------
# Assign stratum (TS/BS) and audited flag to each item
# --------------------------------------------------------------------------
# The sampling interval (SI) determines the TS/BS boundary:
#   SI = BV_project / n_sample = sum(Subtotaal) / sum(Steekproef)
# Items with values >= SI are Top Stratum (TS): examined exhaustively.
# Items with values <  SI are Bottom Stratum (BS): sampled via MUS.
# The audited flag marks items that were part of the UVB sample
# (Steekproef >= 1). Only for audited items the error is known from audit;
# non-audited items have error = 0 by construction (not by observation).
# --------------------------------------------------------------------------

proj_SI = all_items |>
  group_by(project_id) |>
  summarise(
    BV_project   = sum(values, na.rm = TRUE),
    n_sample     = sum(sampled, na.rm = TRUE),
    .groups      = "drop"
  ) |>
  mutate(SI = BV_project / n_sample)

all_items = all_items |>
  inner_join(proj_SI |> select(project_id, SI), by = "project_id") |>
  mutate(
    stratum = ifelse(values >= SI, "TS", "BS"),
    audited = as.integer(sampled >= 1)
  )

# if lines are removed: reverse filter to check which lines and if it is correct that they are removed
uitgevallen_posten <- df |> 
  rename(
    project_id = project_id,
    id         = Code,
    values     = Subtotaal,
    soll       = SOLL,
    sampled    = Steekproef
  ) |>
  mutate(
    values     = as.numeric(values)
  ) |>
  
  filter(!is.finite(values) | values <= VALUE_FLOOR)

# check range of items with excel file
range(all_items$values)



# costitem_db: only sampled items (known errors), with global deciles
costitem_db = all_items |>
  filter(sampled >= 1) |>
  mutate(global_decile = assign_to_decile(values)) |>
  select(project_id, id, values, error, global_decile)

range(costitem_db$values)

# decile_breaks retained for backward compatibility with generate_payment_claim
decile_breaks = DECILE_CUTS

# pool_by_decile: costitem_db split by decile --> split the dbase into 10 parts for easy drawing of cost items
pool_by_decile = split(costitem_db, costitem_db$global_decile)

############################
# BUILD TEMPLATES
############################

build_templates = function(all_items, value_floor = 0.01) {
  
  # Project-level totals (based on ALL items)
  proj_totals = all_items |>
    group_by(project_id) |>
    summarise(
      BV_template = sum(values, na.rm = TRUE),
      n_items     = n(),
      .groups     = "drop"
    ) |>
    filter(BV_template > value_floor)
  
  # Cost items with proportions and decile assignment.
  # Uses assign_to_decile() with the same DECILE_CUTS as costitem_db, so that
  # template items and pool items with the same monetary value always get the
  # same decile label.
  db_prop = all_items |>
    inner_join(proj_totals, by = "project_id") |>
    mutate(
      proportion      = values / BV_template,
      template_decile = assign_to_decile(values)
    ) |>
    select(project_id, id, values, error, proportion, template_decile, stratum, audited)
  
  # One template per project
  templates_obj = split(db_prop, db_prop$project_id)
  
  templates_obj = lapply(names(templates_obj), function(pid) {
    
    items = templates_obj[[pid]] |> arrange(id)
    
    # -------------------------------------------------------------------
    # Compute TS and BS error rates from audited items only.
    # Audited items are those with Steekproef >= 1 in the source data.
    # Non-audited items (all in BS) have error = 0 by construction,
    # not by observation, so they must be excluded from error rate
    # calculation to avoid downward bias.
    # Error rate per stratum = sum(error) / sum(values) for audited items.
    # -------------------------------------------------------------------
    audited_items = items |> filter(audited == 1)
    ts_items = audited_items |> filter(stratum == "TS")
    bs_items = audited_items |> filter(stratum == "BS")
    
    list(
      project_id     = pid,
      BV_template    = sum(items$values, na.rm = TRUE),
      n_items        = nrow(items),
      n_audited_TS   = nrow(ts_items),
      n_audited_BS   = nrow(bs_items),
      TS_error_rate  = if (nrow(ts_items) > 0) sum(ts_items$error) / sum(ts_items$values) else NA_real_,
      BS_error_rate  = if (nrow(bs_items) > 0) sum(bs_items$error) / sum(bs_items$values) else NA_real_,
      costitems      = items |>
        transmute(
          id              = id,
          values          = values,
          error           = error,
          proportion      = proportion,
          template_decile = template_decile,
          stratum         = stratum,
          audited         = audited
        )
    )
  })
  
  names(templates_obj) = sort(unique(db_prop$project_id))
  templates_obj
}

templates = build_templates(all_items, value_floor = VALUE_FLOOR)



############################
# PLOT FOR VISUAL INSPECTION
############################

# Decile summary
df_decile_mean = costitem_db |>
  group_by(global_decile) |>
  summarise(mean_value = mean(values, na.rm = TRUE),
            .groups = "drop")

ggplot(df_decile_mean, aes(x = factor(global_decile), y = mean_value)) +
  geom_col(fill = "steelblue") +
  labs(
    title = "Mean cost item value per decile",
    x = "Global decile",
    y = "Mean item value (EUR)"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

############################
# SAVE
############################

saveRDS(costitem_db, "costitem_db.rds")
saveRDS(templates, "templates.rds")
saveRDS(decile_breaks, "decile_breaks.rds")
saveRDS(pool_by_decile, "pool_by_decile.rds")
