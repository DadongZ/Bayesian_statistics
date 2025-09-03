# Simulate two-arm exponential survival with uniform accrual, admin censoring,
# and dropout specified as an overall percentage of subjects (not a hazard).
simulate_survival_data <- function(
    n_control,
    n_treatment,
    hr_true,
    median_survival_control = 12,    # months
    accrual_duration = 12,           # months (uniform entry in [0, accrual_duration])
    followup_duration = 12,          # months after accrual ends
    dropout_overall = 0,             # overall dropout % across the study (e.g., 20 or 0.20)
    # calibration controls
    calib_n_per_arm = 20000,         # Monte Carlo size per arm for hazard calibration
    tol = 0.005,                     # acceptable abs error in dropout proportion (e.g., 0.5 percentage point)
    max_iter = 30,                   # max iterations for root-finding
    seed = NULL                      # optional seed for reproducibility of calibration
) {
  # --- Helpers ---
  # Interpret dropout_overall as percentage (20 or 0.20 are both OK)
  if (!is.numeric(dropout_overall) || length(dropout_overall) != 1 || dropout_overall < 0) {
    stop("dropout_overall must be a single non-negative numeric value.")
  }
  target_dropout <- if (dropout_overall > 1) dropout_overall / 100 else dropout_overall
  if (target_dropout > 0.999) stop("dropout_overall too large; must be < 100%.")
  
  # Hazards
  lambda_c <- log(2) / median_survival_control
  lambda_t <- lambda_c * hr_true
  study_end <- accrual_duration + followup_duration
  
  if (!is.null(seed)) set.seed(seed)
  
  sim_arm <- function(n, lambda_event, lambda_drop) {
    entry <- runif(n, 0, accrual_duration)
    admin_time <- pmax(study_end - entry, 0)
    t_event <- rexp(n, rate = lambda_event)
    t_drop  <- if (lambda_drop > 0) rexp(n, rate = lambda_drop) else rep(Inf, n)
    time_obs <- pmin(t_event, t_drop, admin_time)
    status   <- as.integer(t_event <= time_obs & is.finite(t_event))
    cause <- ifelse(status == 1, "event",
                    ifelse(t_drop <= admin_time, "dropout", "administrative"))
    data.frame(
      time = time_obs,
      event = status,
      entry_time = entry,
      exit_time  = entry + time_obs,
      cause = cause
    )
  }
  
  simulate_given_dropout_hazard <- function(drop_h, n_c, n_t, for_calibration = FALSE) {
    d_c <- sim_arm(n_c, lambda_c, drop_h)
    d_t <- sim_arm(n_t, lambda_t, drop_h)
    d <- rbind(
      transform(d_c, treatment = 0, arm = "Control"),
      transform(d_t, treatment = 1, arm = "Treatment")
    )
    rownames(d) <- NULL
    if (for_calibration) {
      mean(d$cause == "dropout")
    } else {
      d
    }
  }
  
  # If target_dropout == 0, skip calibration
  if (target_dropout == 0) {
    drop_hazard <- 0
  } else {
    # ---- Calibrate dropout hazard to hit target overall dropout proportion ----
    # We search for h in [0, h_upper] such that observed_dropout(h) ~= target_dropout.
    # Increase h_upper until dropout proportion is well above target, then do bisection.
    # Start with a sensible upper bound: a few multiples of average event hazard + admin censoring effect.
    avg_event_h <- (lambda_c + lambda_t) / 2
    h_upper <- max(5 * avg_event_h, 0.5) # start upper bound
    obs_upper <- simulate_given_dropout_hazard(h_upper, calib_n_per_arm, calib_n_per_arm, TRUE)
    
    tries <- 0
    while (obs_upper < min(0.98, target_dropout + 0.1) && tries < 10) {
      h_upper <- h_upper * 2
      obs_upper <- simulate_given_dropout_hazard(h_upper, calib_n_per_arm, calib_n_per_arm, TRUE)
      tries <- tries + 1
    }
    if (obs_upper < target_dropout) {
      warning("Could not reach the target dropout proportion even with a very high hazard. ",
              "The achieved dropout will be below the target.")
      drop_hazard <- h_upper
    } else {
      # Bisection search
      h_low <- 0
      h_high <- h_upper
      obs_low <- simulate_given_dropout_hazard(h_low,  calib_n_per_arm, calib_n_per_arm, TRUE) # ~0
      obs_high <- obs_upper
      
      iter <- 0
      while (iter < max_iter) {
        h_mid <- 0.5 * (h_low + h_high)
        obs_mid <- simulate_given_dropout_hazard(h_mid, calib_n_per_arm, calib_n_per_arm, TRUE)
        if (abs(obs_mid - target_dropout) <= tol) {
          drop_hazard <- h_mid
          break
        }
        if (obs_mid < target_dropout) {
          h_low  <- h_mid
          obs_low <- obs_mid
        } else {
          h_high <- h_mid
          obs_high <- obs_mid
        }
        iter <- iter + 1
        if (iter == max_iter) drop_hazard <- h_mid
      }
    }
  }
  
  # ---- Final simulation with calibrated hazard ----
  out <- simulate_given_dropout_hazard(drop_hazard, n_control, n_treatment, FALSE)
  attr(out, "calibrated_dropout_hazard") <- drop_hazard
  attr(out, "target_dropout_overall") <- target_dropout
  attr(out, "settings") <- list(
    median_survival_control = median_survival_control,
    hr_true = hr_true,
    accrual_duration = accrual_duration,
    followup_duration = followup_duration,
    calib_n_per_arm = calib_n_per_arm,
    tol = tol,
    max_iter = max_iter
  )
  out
}
