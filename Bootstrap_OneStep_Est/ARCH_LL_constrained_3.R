# Constrained Gaussian LL Objective (on alpha2)
log_likelihood_arch3_constrained <- function(par) {
  alpha1 <- par[1]
  alpha2 <- par[2]
  alpha3 <- alpha3_null
  omega <- par[3]
  sigma_squared <- rep(0, n-3)
  log_likelihood <- 0
  eps <- 1e-8
  for (t in 4:n) {
    sigma_squared[t-3] <- omega + alpha1 * X[t - 1]^2 + alpha2 * X[t - 2]^2 + alpha3 * X[t - 3]^2
    log_likelihood <- log_likelihood - 0.5 * (log(sigma_squared[t-3] + eps) + (X[t]^2 /  sigma_squared[t-3]))
  }
  return(-log_likelihood)  
}




