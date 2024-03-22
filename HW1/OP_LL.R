OP_LL <- function(y,x,par) {
  
  n = length(y) 
  alpha1 = par[1]
  alpha2 = par[2]
  beta = par[3]
  
  #construct conditional probability mass
  phi0 <- pnorm(alpha1-beta * x, mean=0, sd=1) #conditional probability of y_i=0 | x_i
  phi2 <- pnorm(beta * x - alpha2, mean=0, sd=1) #conditional probability of y_i=2 | x_i
  phi1 <- 1- phi0 - phi2 #conditional probability of y_i=1 | x_i
  
  eps <- 1e-10  # Small positive constant to avoid problem with logarithm calculation
  phi0 <- phi0 + eps
  phi1 <- phi1 + eps
  phi2 <- phi2 + eps
  
  # Calculate the log-likelihood components
  log_likelihood_y0 <- sum(log(phi0[y == 0]))
  log_likelihood_y1 <- sum(log(phi1[y == 1]))
  log_likelihood_y2 <- sum(log(phi2[y == 2]))
  
  # Calculate the overall negative log-likelihood
  f = -(1/n)* (log_likelihood_y0 + log_likelihood_y1 + log_likelihood_y2)
  
	return(f)
}
