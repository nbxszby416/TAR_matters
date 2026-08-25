## source functions

get_eligible_cohort <- function(cohort, visit_table, trial_start_date,
                                 min_age = 12, max_age = 21, lookback_months = 24) {
  trial_start <- as.Date(trial_start_date)

  ## no documented infection or vaccination before trial start
  cohort1 <- cohort %>%
    filter(is.na(covid_visit_date) | covid_visit_date > trial_start) %>%
    filter(is.na(imm_date_first)  | imm_date_first  >= trial_start)

  ## >= 1 visit in the lookback window, and none in the 3 days right before
  ## trial start (avoids indexing on a sick visit)
  visit_flags <- visit_table %>%
    filter(person_id %in% cohort1$person_id) %>%
    group_by(person_id) %>%
    summarise(
      has_lookback_visit = any(visit_start_date >= trial_start %m-% months(lookback_months) &
                                  visit_start_date < trial_start),
      has_recent_visit   = any(visit_start_date >= trial_start - days(3) &
                                  visit_start_date < trial_start),
      .groups = "drop"
    )

  cohort2 <- cohort1 %>%
    inner_join(visit_flags, by = "person_id") %>%
    filter(has_lookback_visit, !has_recent_visit) %>%
    select(-has_lookback_visit, -has_recent_visit)

  ## age constraint
  cohort2 %>%
    mutate(age_trial_start = as.numeric(difftime(trial_start, birth_date, units = "days")) / 365.25) %>%
    filter(age_trial_start >= min_age, age_trial_start < max_age)
}


get_vaccinated_elig <- function(cohort, trial_start_date, trial_end_date) {
  trial_start <- as.Date(trial_start_date)
  trial_end   <- as.Date(trial_end_date)

  cohort %>%
    filter(imm_date_first >= trial_start, imm_date_first < trial_end) %>%
    filter(is.na(covid_visit_date) | covid_visit_date > imm_date_first) %>%
    mutate(cohort_entry_date = as.Date(imm_date_first), trt = 1)
}


get_unvaccinated_elig <- function(cohort, trial_start_date, trial_end_date) {
  trial_start <- as.Date(trial_start_date)
  trial_end   <- as.Date(trial_end_date)

  cohort %>%
    filter(is.na(imm_date_first) | imm_date_first >= trial_end) %>%
    filter(is.na(covid_visit_date) | covid_visit_date > trial_start) %>%
    mutate(cohort_entry_date = as.Date(trial_start), trt = 0)
}



build_window_vars <- function(data, event_col, censor_col, entry_col,
                              study_end_date, day_lo, day_hi, prefix) {
  entry     <- as.Date(data[[entry_col]])
  event     <- as.Date(data[[event_col]])
  censor    <- as.Date(data[[censor_col]])
  study_end <- as.Date(study_end_date)

  window_start <- entry + day_lo
  window_end   <- if (is.infinite(day_hi)) study_end else entry + day_hi
  admin_end    <- pmin(study_end, window_end)

  ## earliest of: administrative end of the window, or comparator crossing
  ## over to vaccinated
  followup_end <- pmin(admin_end, censor, na.rm = TRUE)

  ptime    <- pmax(0, as.numeric(followup_end - window_start))
  occurred <- !is.na(event) & event > window_start & event <= followup_end

  out <- data.frame(ptime = ptime, event = as.integer(occurred))
  names(out) <- paste0(prefix, c("_ptime", "_event"))
  cbind(data, out)
}


run_modified_poisson <- function(data, ptime_col, event_col, trt_col = "trt") {
  d <- data.frame(event = data[[event_col]], ptime = data[[ptime_col]], trt = data[[trt_col]])
  d <- d[!is.na(d$ptime) & d$ptime > 0, ]

  fit <- glm(event ~ trt, offset = log(ptime), family = poisson(link = "log"), data = d)
  robust_se <- sqrt(diag(sandwich::vcovHC(fit, type = "HC0")))

  logRR   <- unname(coef(fit)["trt"])
  seLogRR <- unname(robust_se["trt"])

  data.frame(logRR = logRR, seLogRR = seLogRR,
             RR = exp(logRR), RR_lower = exp(logRR - 1.96 * seLogRR), RR_upper = exp(logRR + 1.96 * seLogRR))
}

