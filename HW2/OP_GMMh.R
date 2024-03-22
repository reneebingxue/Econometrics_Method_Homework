OP_GMMh <- function(y,x,par) {
  
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
  
  #generate the premultiply matrix
  prem1 <- y0 - phi0
  prem2 <- y1 - phi1
  prem_matrix <- rbind(t(prem1),t(prem2))
  
  #generate moment conditions
  G <- cbind(1, x) 
  check = (prem_matrix) %*%  G
  H = (prem_matrix %*%  G)/n #H has the 4 elements of moment conditions in 2*2 matrix form 
  H_squared <- H^2 #since we consider identity weight matrix
  f= sum(H_squared)
  
  if (FALSE){
  #generate the premultiply matrix
  g <- cbind(1, x)
  prem1 <- y0 - phi0
  prem2 <- y1 - phi1
  m0 <- (prem1*g)
  m1 <- (prem2*g)
  m <- cbind(m0,m1)
  
  h <- colSums(m)*(1/n)
  
  f= t(h) %*% h}
  
  return(f)
}
