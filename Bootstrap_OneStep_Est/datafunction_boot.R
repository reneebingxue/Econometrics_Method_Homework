# Generating data in bootstrap
simulate_arch2_boot <- function(n, z_boot, theta_hat) {
  theta_hat <- theta_hat
  X_b <- numeric(n)  
  # initial variance assuming stationary process
  X_b[1] <- 0
  X_b[2] <- 0
  
  for (t in 3:n) {
    X_b[t] <- sqrt(theta_hat[3] + theta_hat[1] * X_b[t - 1]^2 + theta_hat[2] * X_b[t - 2]^2) * z_boot[t]
  }
  
  return(X_b)
}
