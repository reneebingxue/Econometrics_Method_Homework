# Function to simulate ARCH(3) process
simulate_arch3 <- function(n, theta) {
  z <- rnorm(n)  # generating i.i.d. standard normal variables
  sigma_squared <- numeric(n)  # variance will be stored here
  X <- numeric(n)  # ARCH(3) process will be stored here
  
  # initial values of X assuming mean zero
  X[1] <- 0
  X[2] <- 0
  X[3] <- 0
  
  for (t in 4:n) {
    sigma_squared[t] <- theta[4] + theta[1] * X[t - 1]^2 + theta[2] * X[t - 2]^2 + theta[3] * X[t - 3]^2
    X[t] <- sqrt(sigma_squared[t]) * z[t]
  }
  
  return(list(X = X, sigma_squared = sigma_squared, z = z))
}
