log_likelihood_arch2 <- function(theta) {
  alpha1 <- theta[1]
  alpha2 <- theta[2]
  omega <- theta[3]
  sigma_squared <- rep(0, n-2)
  log_likelihood = 0
  eps <- 1e-8
  for (t in 3:n) {
    sigma_squared[t-2] <- omega + alpha1 * X[t - 1]^2 + alpha2 * X[t - 2]^2 
    log_likelihood <- log_likelihood - 0.5 * (log(sigma_squared[t-2] + eps) + (X[t]^2 / sigma_squared[t-2]))
  }
  return(-log_likelihood) 
}
