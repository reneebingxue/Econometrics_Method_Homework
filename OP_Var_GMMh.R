OP_Var_GMMh <- function(y,x,par) {
  
  n = length(y) 
  alpha1 = par[1]
  alpha2 = par[2]
  beta = par[3]
  
  #construct conditional probabilities
  phi0 <- pnorm(alpha1-beta * x, mean=0, sd=1) #conditional probability of y_i=0 | x_i
  phi2 <- pnorm(beta * x - alpha2, mean=0, sd=1) #conditional probability of y_i=2 | x_i
  phi1 <- 1- phi0 - phi2 #conditional probability of y_i=1 | x_i
  
  #construct indicator functions
  y0<- ifelse(y==0,1,0)
  y1<- ifelse(y==1,1,0)
  
  #generate moment conditions
  g <- cbind(1, x)
  prem1 <- y0 - phi0
  prem2 <- y1 - phi1
  m0 <- (prem1*g)
  m1 <- (prem2*g)
  m <- cbind(m0,m1)
  
  #Define G in theorem 7: jacobian of gmm moment condition matrix: 4 by 3
  r1 <- colSums(cbind(-dnorm(alpha1-beta*x),0,dnorm(alpha1-beta*x)*x))*(1/n)
  r2 <- colSums(cbind(-dnorm(alpha1-beta*x)*x,0,dnorm(alpha1-beta*x)*x^2))*(1/n)
  r3 <- colSums(cbind(dnorm(alpha1-beta*x),-dnorm(alpha2-beta*x),(dnorm(alpha2-beta*x)-dnorm(alpha1-beta*x))*x))*(1/n)
  r4 <- colSums(cbind(dnorm(alpha1-beta*x)*x,-dnorm(alpha2-beta*x)*x,(dnorm(alpha2-beta*x)-dnorm(alpha1-beta*x))*x^2))*(1/n)
  G = rbind(r1,r2,r3,r4)
  
  #Define Om in theorem 7: (optimal weight matrix) variance-covariance matrix of gmm moment condition: 4 by 4 
  O = (1/n)* (t(m)%*%m)
  
  #Define GMM Var-Cov matrix
  VC = solve(t(G)%*%G)%*%(t(G)%*%O%*%G)%*%solve(t(G)%*%G)
	return(VC)
}
