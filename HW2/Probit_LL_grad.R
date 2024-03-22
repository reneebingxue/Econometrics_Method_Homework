Probit_LL_grad<- function (y,x,par) {
  
  n = length(y) 
  alpha1 = par[1]
  alpha2 = par[2]
  beta = par[3]
  
  #construct conditional probability mass
  phi0 <- dnorm(alpha1-beta * x, mean=0, sd=1) #conditional probability of y_i=0 | x_i
  phi2 <- dnorm(beta * x - alpha2, mean=0, sd=1) #conditional probability of y_i=2 | x_i
  phi1 <- dnorm(alpha2 - beta * x, mean=0, sd=1) - dnorm(alpha1 - beta * x, mean=0, sd=1)
 # phic <-   1- phi0 - phi2 #conditional probability of y_i=1 | x_i
  Phi0 <- pnorm(alpha1 - beta * x, mean=0, sd=1)
  Phi2 <- pnorm(beta * x - alpha2, mean=0, sd=1) #conditional probability of y_i=0 | x_i
  Phi1 <- pnorm(alpha2-beta * x, mean=0, sd=1)-pnorm(alpha1-beta * x, mean=0, sd=1) #conditional probability of y_i=0 | x_i

  # Computing the gradient
  g1 = (y==0)*(phi0/Phi0) - (y==1)*(phi0/Phi1)
  g2 = - (y==2)*(phi2/Phi2) + (y==1)*(phi2/Phi1)
  g3 = - (y==0)*(phi0/Phi0)*x - (y==1)*(phi1/Phi1)*x + (y==2)*(phi2/Phi2)*x
  
  G = colMeans(cbind(g1,g2,g3))
  
  return(G) 
}