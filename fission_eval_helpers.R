library(datathin)

## We need this more complicated data generating mechanism if we want sctransform to work like at all for esitmating overdispersions. 
## But for them moment I am not using this. 
makeData_tau <- function(n,p, K,clusterMeans, tau) {
  mean_mat <- matrix(0, nrow=n, ncol=p)
  true_clust <- sample(1:K, size=n, replace=T)
  if (K==1) {
    mean_mat <- matrix(clusterMeans, nrow=n, ncol=p, byrow=T)
  } else {
    for (k in 1:K) {
      mean_mat[true_clust==k,] <- matrix(clusterMeans[k,], nrow=sum(true_clust==k), ncol=p, byrow=T)
    }
  }
  
  lambda_bars <- colMeans(mean_mat)
  overdisps <- lambda_bars/tau
  
  X <- sapply(1:p, function(u) rnbinom(n, mu=mean_mat[,u], size=overdisps[u]))
  return(list(dat=X, clusters=true_clust, overdisps=overdisps))
}

## Generate simple data. 
## clusterMeans must have dimension K*p. 
## Overdisp should just be a vector of length p, for gene-specific overdispersion.
makeData <- function(n,p, K,clusterMeans, overdisps) {
  mean_mat <- matrix(0, nrow=n, ncol=p)
  true_clust <- sample(1:K, size=n, replace=T)
  if (K==1) {
    mean_mat <- matrix(clusterMeans, nrow=n, ncol=p, byrow=T)
  } else {
    for (k in 1:K) {
      mean_mat[true_clust==k,] <- matrix(clusterMeans[k,], nrow=sum(true_clust==k), ncol=p, byrow=T)
    }
  }
  X <- sapply(1:p, function(u) rnbinom(n, mu=mean_mat[,u], size=overdisps[u]))
  return(list(dat=X, clusters=true_clust))
}

### Evaluate the log likelihood of entire dataset given estimated clusters. No train/test
### split. 
naive.lik <- function(X, clusterlabs) {
  
  p <- NCOL(X)
  mean.hats <- matrix(NA, nrow=n, ncol=p)
  disp.hats <- matrix(NA, nrow=n, ncol=p)
  
  if (length(unique(clusterlabs))==1) {
    for (j in 1:NCOL(X)) {
      mod <- MASS::glm.nb(X[,j]~1)
      mean.hats[,j] <- predict(mod, type="response")
      disp.hats[,j] <- mod$theta
    }
  } else {
    for (j in 1:NCOL(X)) {
      mod <- MASS::glm.nb(X[,j]~as.factor(clusterlabs))
      mean.hats[,j] <- predict(mod, type="response")
      disp.hats[,j] <- mod$theta
    }}
  
  log.liks <- sapply(1:p, function(u) dnbinom(X[,u], mu=mean.hats[,u], size=disp.hats[,u], log=T))
  return(sum(-log.liks))
}


### Evaluate the log likelihood of a test set using parameters estimated on a training set.
### Assumes independence between training set and test set.
thin.lik <- function(dat.train, dat.test, overdisps, clusterlabs, epsilon) {
  p <- NCOL(dat.train)
  mean.hats.train <- matrix(NA, nrow=n, ncol=p)
  disp.hats.train <- rep(NA,p)
  if (length(unique(clusterlabs))==1) {
    for (j in 1:NCOL(dat.train)) {
      mod <- MASS::glm.nb(dat.train[,j]~1)
      mean.hats.train[,j] <- predict(mod, type="response")
      disp.hats.train[j] <- mod$theta
    }
  } else {
    for (j in 1:NCOL(dat.train)) {
      mod <- MASS::glm.nb(dat.train[,j]~as.factor(clusterlabs))
      
      mean.hats.train[,j] <- predict(mod, type="response")
      disp.hats.train[j] <- mod$theta
    }
  }
  overdisps.hat <- disp.hats.train/epsilon
  mean.hats <- mean.hats.train/epsilon
  
  log.liks <- sapply(1:p, function(u) dnbinom(dat.test[,u], 
                                              mu=(1-epsilon)*mean.hats[,u], size=overdisps.hat[u]*(1-epsilon), log=T))
  return(sum(-log.liks))
}


## Evaluate the conditional log likelihood of dataset X2 given X1 and the
## estimated clusters. 
## Obtain estimates of mean and overdispersion from the training set.
## Evaluate conditional log-lik on test set.
## Pretty straightforward.  
conditional.lik <- function(X1, X2, clusterlabs, eps) {
  
  p <- NCOL(X1)
  mean.hats.train <- matrix(NA, nrow=n, ncol=p)
  disp.hats <- rep(NA, p)
  
  if (length(unique(clusterlabs))==1) {
    for (j in 1:NCOL(X1)) {
      mod <- MASS::glm.nb(X1[,j]~1)
      mean.hats.train[,j] <- predict(mod, type="response")
      disp.hats[j] <- mod$theta
    }
  } else {
    for (j in 1:p) {
      mod <- MASS::glm.nb(X1[,j]~as.factor(clusterlabs))
      mean.hats.train[,j] <- predict(mod, type="response")
      disp.hats[j] <- mod$theta
    }
  }
  
  theta.hats <- sapply(1:p, function(u) 1/(1+mean.hats.train[,u]/(eps*disp.hats[u])))
  
  log.liks <- sapply(1:p, function(u) dnbinom(X2[,u], prob=theta.hats[,u]+eps-theta.hats[,u]*eps, size=disp.hats[u]+X1[,u], log=T))
  return(sum(-log.liks))
}
