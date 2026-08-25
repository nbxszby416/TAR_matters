## main analysis

library(dplyr)
library(lubridate)
library(MatchIt)
library(sandwich)
library(EmpiricalCalibration)

source("source_function.R")

study_start_date <- "2022-01-01"
study_end_date   <- "2022-11-30"
entry_end_date   <- "2022-11-16"

min_age <- 12
max_age <- 21
n <- 14  # days between sequential trial starts

cohort <- read.csv("demo_cohort.csv", stringsAsFactors = FALSE) %>%
  mutate(across(ends_with("_date"), as.Date))
visit_table <- read.csv("demo_visits.csv", stringsAsFactors = FALSE) %>%
  mutate(visit_start_date = as.Date(visit_start_date))

match_vars <- c("age_trial_start", "sex_cat", "eth_cat", "pmca_index", "site")

### sequential target trial emulation
df <- data.frame()
trial_start_date <- study_start_date
i <- 1

while (as.Date(trial_start_date) <= as.Date(entry_end_date)) {
  trial_end_date <- min(as.Date(trial_start_date) + days(n), as.Date(entry_end_date) + 1)
  cat(sprintf("Trial %d: %s to %s\n", i, trial_start_date, trial_end_date))

  cohort_elig <- get_eligible_cohort(cohort, visit_table, trial_start_date, min_age, max_age)

  data_vac <- get_vaccinated_elig(cohort_elig, trial_start_date, trial_end_date)
  data_unvac <- get_unvaccinated_elig(cohort_elig, trial_start_date, trial_end_date)

  mydata <- bind_rows(data_vac, data_unvac) %>% as.data.frame()

  if (n_distinct(mydata$trt) < 2) {
    i <- i + 1
    trial_start_date <- format(as.Date(trial_end_date))
    next
  }

  ## exact matching
  formula_match <- as.formula(paste("trt ~", paste(match_vars, collapse = " + ")))

  m.out <- tryCatch(
    matchit(formula_match, data = mydata, method = "cem", k2k = TRUE),
    error = function(e) NULL
  )
  m.data <- if (is.null(m.out)) NULL else tryCatch(match.data(m.out), error = function(e) NULL)

  if (is.null(m.data) || nrow(m.data) == 0 || n_distinct(m.data$trt) < 2) {
    i <- i + 1
    trial_start_date <- format(as.Date(trial_end_date))
    next
  }

  m.data$trial_id <- i
  m.data$trial_start_date <- trial_start_date
  m.data$trial_end_date <- trial_end_date

  m.treated <- m.data %>% filter(trt == 1)
  m.untreated <- m.data %>% filter(trt == 0)

  m.treated <- m.treated %>%
    left_join(m.untreated %>% select(subclass, comparator_imm_date = imm_date_first), by = "subclass") %>%
    mutate(censor_date = comparator_imm_date) %>%
    select(-comparator_imm_date)

  m.untreated <- m.untreated %>%
    mutate(censor_date = imm_date_first)

  df <- bind_rows(df, m.treated, m.untreated)

  i <- i + 1
  trial_start_date <- format(as.Date(trial_end_date))
}

cat(sprintf("Total matched person-trials: %d\n", nrow(df)))

### outcome preprocessing
outcomes <- list(primary = list(col = "covid_visit_date",        windows = c("short", "long")),
                 cand_short_cov2 = list(col = "cand_short_sarscov2_date", windows = c("short")),
                 cand_short_resp = list(col = "cand_short_resp_date",     windows = c("short")),
                 cand_long_resp  = list(col = "cand_long_resp_date",      windows = c("long")))
trusted_nco_cols <- grep("^nco_trusted_", names(df), value = TRUE)
for (v in trusted_nco_cols){
  outcomes[[v]] <- list(col = v, windows = c("long"))}

window_bounds <- list(short = c(lo = 0, hi = 14), long = c(lo = 14, hi = Inf))

for (nm in names(outcomes)) {
  spec <- outcomes[[nm]]
  for (w in spec$windows) {
    bnd <- window_bounds[[w]]
    df <- build_window_vars(df,
                            event_col = spec$col, censor_col = "censor_date", entry_col = "cohort_entry_date",
                            study_end_date = study_end_date,
                            day_lo = unname(bnd["lo"]), day_hi = unname(bnd["hi"]),
                            prefix = paste0(nm, "_", w))
  }
}


### outcome model
results <- data.frame()
for (nm in names(outcomes)) {
  spec <- outcomes[[nm]]
  for (w in spec$windows) {
    rslt <- run_modified_poisson(df, ptime_col = paste0(nm, "_", w, "_ptime"), event_col = paste0(nm, "_", w, "_event"))
    rslt$outcome <- nm
    rslt$window <- w
    results <- bind_rows(results, rslt)
  }
}
print(results)

### hypothesis testing
trusted_rslt <- results %>% filter(outcome %in% trusted_nco_cols, window == "long")
fitnull <- fitNull(trusted_rslt$logRR, trusted_rslt$seLogRR)
cat(sprintf("Trusted-NCO null distribution: mean = %.3f, sd = %.3f\n", null_fit[1], null_fit[2]))

candidate_rslt <- results %>% filter(outcome %in% c("cand_short_cov2", "cand_short_resp", "cand_long_resp"))
p <- 2*pmin(pnorm((fitnull[1]-candidate_rslt$logRR)/sqrt(fitnull[2]^2+candidate_rslt$seLogRR^2)), 
                   pnorm((candidate_rslt$logRR-fitnull[1])/sqrt(fitnull[2]^2+candidate_rslt$seLogRR^2)))
print(p)


