# sim_setup.R

source("fission_eval_helpers.R")
source("method_eval_helpers.R")

library(patchwork)
library(purrr)
library(tibble)
library(dplyr)

run_one_replicate <- function(seed, counter, d, n, p, sigma, eps_vec) {
  set.seed(seed)
  
  d_vals <- paste(d, collapse = ",")
  
  dat <- generate_dataset_with_intercept(n = n, p = p, d = d, sigma = sigma) #sd
  X <- dat$X
  true_r <- dat$true_r
  true_mu <- dat$mu
  
  X_raw <- dat$X_raw
  true_mu_raw <- dat$mu_raw
  
  # --- BCV ---
  bicv <- map(2:3, ~ bicv_eval(X_raw, nfold = .x, true_mu_raw))
  
  bcv_tbl <- tibble(
    Counter = counter,
    Seed = seed,
    d_vals = d_vals,
    True_r = true_r,
    Best_r = map_int(bicv, ~ which.min(.x$err_true) - 1),
    Method = "BCV",
    k = c(2L, 3L),
    eps = NA_real_,
    sigtil = NA_real_,
    sigma = sigma,
    r_hat = map_int(bicv, ~ which.min(.x$err) - 1)
  )
  
  # --- Data thinning ---
  thin_out <- map(eps_vec, function(eps_i) {

    res <- datathin(
      X_raw,
      arg = sigma^2,
      family = "gaussian",
      epsilon = c(eps_i, 1 - eps_i)
    )

    X1_raw <- res[, , 1]
    X2_raw <- res[, , 2]

    X1_scaled <- scale(X1_raw, center = TRUE, scale = TRUE) #was true before

    center_vec <- attr(X1_scaled, "scaled:center")
    scale_vec  <- attr(X1_scaled, "scaled:scale")

    X2_scaled <- scale(X2_raw,
                       center = center_vec,
                       scale  = scale_vec)

    mu2_raw <- (1 - eps_i) * dat$mu_raw
    mu2_scaled <- scale(mu2_raw,
                        center = center_vec,
                        scale  = scale_vec)

    thin_eval(
      X1_scaled,
      X2_scaled,
      eps_i,
      mu2_scaled
    )
  })

  thin_tbl <- tibble(
    Counter = counter,
    Seed = seed,
    d_vals = d_vals,
    True_r = true_r,
    Best_r = map_int(thin_out, ~ which.min(.x$err_true) - 1),
    Method = "Thin",
    k = NA_integer_,
    eps = eps_vec,
    sigtil = NA_real_,
    sigma = sigma,
    r_hat = map_int(thin_out, ~ which.min(.x$err) - 1)
  )
  
  # --- K-fold thinning (multifold) ---
  K_vec <- 2:5

  mthin_out <- map(K_vec, function(K) {
    multithin_eval(
      X_raw = X_raw,
      K = K,
      mu_raw = dat$mu_raw,
      sigma = sigma
    )
  })

  mthin_tbl <- tibble(
    Counter = counter,
    Seed = seed,
    d_vals = d_vals,
    True_r = true_r,
    Best_r = map_int(mthin_out, ~ which.min(.x$err_true) - 1),
    Method = paste0(K_vec, "fold_Thin"),
    k = K_vec,
    eps = round((K_vec - 1) / K_vec, 2),
    sigtil = NA_real_,
    sigma = sigma,
    r_hat = map_int(mthin_out, ~ which.min(.x$err) - 1)
  )

  # --- Data fission ---
  sigma_tilde <- c(sqrt(0.0001), sqrt(0.01), sqrt(0.7), sqrt(1), sqrt(1.3), sqrt(2))

  
  fiss_out <- map(sigma_tilde, function(sig_t) {
    res <- datathin(X, arg = sig_t^2, family = "gaussian")
    fiss_eval(
      res[, , 1],
      res[, , 2],
      sig_true  = sigma,
      sig_tilde = sig_t,
      eps = 0.5,
      true_mu
    )
  })

  fiss_tbl <- tibble(
      Counter = counter,
      Seed    = seed,
      d_vals  = d_vals,
      True_r  = true_r,
    Best_r = map_int(fiss_out, ~ which.min(.x$err_true) - 1),
    Method = "Fission",
    k = NA_integer_,
    eps = 0.5,
    sigtil = sigma_tilde,
    sigma = sigma,
    r_hat = map_int(fiss_out, ~ which.max(.x$loglike) - 1)
  )

#  --- K-fold fission ---
  K_vec <- 5:5

  mfiss_tbl <- map_dfr(K_vec, function(K) {

    mfiss_out <- map(sigma_tilde, function(sig_t) {
      multifiss_eval(
        X = X_raw,
        K = K,
        sig_true  = sigma,
        sig_tilde = sig_t,
        mu = true_mu_raw
      )
    })

    tibble(
      Counter = counter,
      Seed = seed,
      d_vals = d_vals,
      True_r = true_r,
      Best_r = map_int(mfiss_out, ~ which.min(.x$err_true) - 1),
      Method = paste0(K, "fold_Fission"),
      k = K,
      eps = (K - 1) / K,
      sigtil = sigma_tilde,
      sigma = sigma,
      r_hat = map_int(mfiss_out, ~ which.max(.x$loglike) - 1)
    )
  })
  
  bind_rows(bcv_tbl, thin_tbl, mthin_tbl, fiss_tbl, mfiss_tbl)
}
