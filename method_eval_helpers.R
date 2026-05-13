# method_eval_helpers.R

library(datathin)

rand_orth <- function(n, rank) {
  Z <- matrix(rnorm(n * rank), n, rank)
  qr.Q(qr(Z))
}

generate_dataset_with_intercept <- function(n, p, d_vals, sigma) {
  rank <- length(d_vals)
  U <- rand_orth(n, rank) # n x rank
  V <- rand_orth(p, rank) # p x rank
  D <- diag(d_vals, nrow = rank, ncol = rank)
  
  mu_centered <- U %*% D %*% t(V)
  mu_centered <- scale(mu_centered, center = TRUE, scale = FALSE)
  
  mu_colmeans <- rnorm(p, mean = 0, sd = 1)
  mu0 <- matrix(rep(mu_colmeans, each = n), nrow = n, ncol = p)
  mu_raw <- mu0 + mu_centered
  
  X_raw <- mu_raw + matrix(rnorm(n * p, sd = sigma), n, p)
  
  X_scaled <- scale(X_raw, center = TRUE, scale = TRUE)
  
  mu_scaled <- scale(mu_raw, 
                     center = attr(X_scaled, "scaled:center"),
                     scale  = attr(X_scaled, "scaled:scale"))
  
  X <- X_scaled
  mu <- mu_scaled
  
  pc <- princomp(X, cor = FALSE)
  
  list(
    U = U, V = V, D = D,
    mu0 = mu0, mu_centered = mu_centered,
    mu = mu, X = X,
    mu_raw = mu_raw, X_raw = X_raw,
    true_r = rank,
    pc = pc
  )
}

bicv_eval <- function(X, nfold, mu) {
  ranks <- 0:floor((ncol(X) * (nfold - 1)) / nfold)
  n <- nrow(X)
  p <- ncol(X)
  
  row_folds <- sample(rep(1:nfold, length.out = n))
  col_folds <- sample(rep(1:nfold, length.out = p))
  
  errs <- lapply(ranks, function(k) numeric(0))
  errs_true <- lapply(ranks, function(k) numeric(0))
  names(errs) <- as.character(ranks)
  names(errs_true) <- as.character(ranks)
  
  for (i in 1:nfold) {
    for (j in 1:nfold) {
      I <- which(row_folds == i)
      J <- which(col_folds == j)
      if (length(I) == 0 || length(J) == 0) next
      
      IJ   <- X[I, J, drop = FALSE]
      IJc  <- X[I, -J, drop = FALSE]
      IcJ  <- X[-I, J, drop = FALSE]
      IcJc <- X[-I, -J, drop = FALSE]
      
      train_means_all  <- colMeans(X[-I, , drop = FALSE])
      train_means_train <- train_means_all[-J]
      train_means_hold  <- train_means_all[J]
      
      IcJc_centered <- scale(IcJc, center = train_means_train, scale = FALSE)
      IJc_centered  <- scale(IJc,  center = train_means_train, scale = FALSE)
      IcJ_centered  <- scale(IcJ,  center = train_means_hold,  scale = FALSE)
      
      svd_fit <- svd(IcJc_centered)
      
      for (k in ranks) {
        if (k == 0) {
          pred <- matrix(train_means_hold, nrow = length(I), ncol = length(J), byrow = TRUE)
        } else {
          Uk <- svd_fit$u[, 1:k, drop = FALSE]
          sk <- svd_fit$d[1:k]
          Vk <- svd_fit$v[, 1:k, drop = FALSE]
          Sk_inv <- diag(1 / sk, nrow = k)
          
          pred_centered <- IJc_centered %*% Vk %*% Sk_inv %*% t(Uk) %*% IcJ_centered
          pred <- pred_centered + matrix(train_means_hold, nrow = length(I), ncol = length(J), byrow = TRUE)
        }
        
        err <- mean((IJ - pred)^2)
        errs[[as.character(k)]] <- c(errs[[as.character(k)]], err)
        err_true <- mean((mu[I,J] - pred)^2)
        errs_true[[as.character(k)]] <- c(errs_true[[as.character(k)]], err_true)
      }
    }
  }
  list(err = sapply(errs, mean), err_true = sapply(errs_true, mean))
}

thin_eval <- function(X1, X2, eps, mu) {
  n <- nrow(X1)
  p <- ncol(X1)
  
  mu_hat <- matrix(colMeans(X1), nrow = n, ncol = p, byrow = TRUE) #scale of X1
  
  X1_centered <- scale(X1, center = TRUE, scale = FALSE)
  sv <- svd(X1_centered)
  
  ranks <- 0:p
  errs <- lapply(ranks, function(k) numeric(0))
  names(errs) <- as.character(ranks)
  errs_true <- lapply(ranks, function(k) numeric(0))
  names(errs_true) <- as.character(ranks)
  
  for (k in ranks) {
    if (k == 0) {
      #pred_X2_1 <- mu_hat #scale of X1
      pred_X2_1 <- matrix(colMeans(X1), nrow=n, ncol=p, byrow=TRUE)  # rank-0
    } else {
      X1_k <- sv$u[, 1:k, drop = FALSE] %*%
        diag(sv$d[1:k], nrow = k) %*%
        t(sv$v[, 1:k, drop = FALSE])
      pred_X2_1 <- matrix(colMeans(X1), nrow=n, ncol=p, byrow=TRUE) + X1_k
    }
    pred_X2 <- ((1-eps)/(eps)) * pred_X2_1 #scale of X2
    
    err <- mean((X2 - pred_X2)^2)
    errs[[as.character(k)]] <- c(errs[[as.character(k)]], err)
    err_true <- mean(((mu) - pred_X2)^2)
    errs_true[[as.character(k)]] <- c(errs_true[[as.character(k)]], err_true)
  }
  list(err = sapply(errs, mean), err_true = sapply(errs_true, mean))
}

cond_lik_gaus <- function(X1, X2, mu, sigma_sq_hat, sig_tilde, eps) {
  
  mu_1 <- eps * mu
  mu_2 <- (1 - eps) * mu
  
  sig_1 <- eps^2 * sigma_sq_hat + eps * (1 - eps) * sig_tilde^2
  sig_2 <- (1 - eps)^2 * sigma_sq_hat + eps * (1 - eps) * sig_tilde^2
  
  sig_1 <- max(sig_1, 1e-8)
  sig_2 <- max(sig_2, 1e-8)
  
  rho <- eps * (1 - eps) * (sigma_sq_hat - sig_tilde^2) / sqrt(sig_1 * sig_2)
  rho <- max(min(rho, 0.999), -0.999)  # clamp
  
  cond_mean <- mu_2 + rho * sqrt(sig_2 / sig_1) * (X1 - mu_1)
  cond_var  <- sig_2 * (1 - rho^2)
  cond_var  <- max(cond_var, 1e-8)
  
  loglik <- dnorm(X2, mean = cond_mean, sd = sqrt(cond_var), log = TRUE)
  
  loglik
}

fiss_eval <- function(X1, X2, sig_true, sig_tilde, eps, mu) {
  n <- nrow(X1)
  p <- ncol(X1)
  
  mu0_hat <- matrix(colMeans(X1), nrow = n, ncol = p, byrow = TRUE)
  
  X1_centered <- scale(X1, center = TRUE, scale = FALSE)
  sv <- svd(X1_centered)
  
  ranks <- 0:p
  ll_by_rank <- vector("list", length = length(ranks))
  SSE <- numeric(length(ranks))
  SST <- sum((X1 - mu0_hat)^2)
  err_true <- numeric(length(ranks))
  
  for (h in ranks) {
    
    if (h == 0) {
      mu_hat <- mu0_hat
    } else {
      mu_centered_hat <- sv$u[, 1:h, drop = FALSE] %*%
        diag(sv$d[1:h], nrow = h) %*%
        t(sv$v[, 1:h, drop = FALSE])
      mu_hat <- mu0_hat + mu_centered_hat
    }
    
    mu_1_full <- eps * mu_hat
    v1_hat <- mean((X1 - mu_1_full)^2)
    
    sigma_sq_hat <- (v1_hat - eps * (1 - eps) * sig_tilde^2) / eps^2
    sigma_sq_hat <- max(sigma_sq_hat, 1e-8)
    
    ll_mat <- matrix(0, n, p)
    
    for (i in 1:n) {
      for (j in 1:p) {
        ll_mat[i, j] <- cond_lik_gaus(
          X1[i, j],
          X2[i, j],
          mu_hat[i, j],
          sigma_sq_hat,
          sig_tilde,
          eps
        )
      }
    }
    
    ll_by_rank[[h + 1]] <- ll_mat
    
    pred_X2 <- ((1 - eps) / eps) * mu_hat
    err_true[h + 1] <- sum(((1 - eps) * mu - pred_X2)^2)
  }
  
  means_by_rank <- sapply(ll_by_rank, mean)
  
  if (all(is.infinite(means_by_rank))) {
    warning("All log-likelihoods are -Inf")
  }
  
  R2_by_rank <- 1 - SSE / SST
  
  list(
    loglike = means_by_rank,
    SSE = SSE,
    SST = SST,
    R2 = R2_by_rank,
    err_true = err_true
  )
}

multithin_eval <- function(X_raw, K, mu_raw, sigma) {
  
  res <- datathin(
    X_raw,
    family = "gaussian",
    K = K,
    arg = sigma^2
  )
  
  eps_train <- (K - 1) / K
  p <- ncol(X_raw)
  ranks <- 0:p
  
  fold_err <- matrix(0, nrow = K, ncol = length(ranks))
  fold_err_true <- matrix(0, nrow = K, ncol = length(ranks))
  
  for (f in seq_len(K)) {
    
    X2_raw <- res[, , f]
    X1_raw <- apply(res[, , -f, drop = FALSE], c(1,2), sum)
    
    X1_scaled <- scale(X1_raw, center = TRUE, scale = TRUE)
    center_vec <- attr(X1_scaled, "scaled:center")
    scale_vec  <- attr(X1_scaled, "scaled:scale")
    
    X2_scaled <- scale(X2_raw,
                       center = center_vec,
                       scale  = scale_vec)
    
    mu2_raw <- (1 - eps_train) * mu_raw
    mu2_scaled <- scale(mu2_raw,
                        center = center_vec,
                        scale  = scale_vec)
    
    out <- thin_eval(
      X1_scaled,
      X2_scaled,
      eps_train,
      mu2_scaled
    )
    
    fold_err[f, ] <- out$err
    fold_err_true[f, ] <- out$err_true
  }
  
  list(
    err = colMeans(fold_err),
    err_true = colMeans(fold_err_true)
  )
}

multifiss_eval <- function(
    X, 
    K,
    sig_true,
    sig_tilde,
    mu,
    family = "gaussian"
) {
  res <- datathin(
    X,
    arg = sig_tilde^2,
    family = family,
    K = K
  )
  
  eps_train <- (K - 1) / K
  p <- ncol(X)
  ranks <- 0:p
  
  fold_ll <- matrix(0, nrow = K, ncol = length(ranks))
  fold_err_true <- matrix(0, nrow = K, ncol = length(ranks))
  
  for (f in seq_len(K)) {
    X2 <- res[, , f]
    X1 <- apply(res[, , -f, drop = FALSE], c(1, 2), sum)
    
    out <- fiss_eval(
      X1 = X1,
      X2 = X2,
      sig_true = sig_true,
      sig_tilde = sig_tilde,
      eps = eps_train,
      mu = mu
    )
    
    fold_ll[f, ] <- out$loglike
    fold_err_true[f, ] <- out$err_true
  }
  
  list(
    loglike = colMeans(fold_ll),
    err_true = colMeans(fold_err_true)
  )
}
