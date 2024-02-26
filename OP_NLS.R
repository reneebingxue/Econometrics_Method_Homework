OP_NLS <- function(y,x,par) {
  
  n = length(y) 
  alpha1 = par[1]
  alpha2 = par[2]
  beta = par[3]
  
  #construct conditional probabilities
  phi0 <- pnorm(alpha1-beta * x, mean=0, sd=1) #conditional probability of y_i=0 | x_i
  phi2 <- pnorm(beta * x - alpha2, mean=0, sd=1) #conditional probability of y_i=2 | x_i
  phi1 <- 1 - phi0 - phi2 #conditional probability of y_i=1 | x_i
  
  #construct m(x_i) function 
  m_xt <- (0)*phi0 + (1)*phi1 + (2)*phi2 #m(x_i) function designed based on nonlinear least square
  #m_xt <- pnorm(beta*x-alpha1, mean=0, sd=1)+pnorm(beta*x-alpha2, mean=0, sd=1)

  f = (1/n)*sum((y-m_xt)^2)

	return(f)
}
