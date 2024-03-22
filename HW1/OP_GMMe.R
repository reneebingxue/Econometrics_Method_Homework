OP_GMMe <- function(y,x,par) {
  
  n = length(y) 
  alpha1 = par[1]
  alpha2 = par[2]
  beta = par[3]
  
  #construct conditional probabilities
  phi0 <- pnorm(alpha1-beta * x, mean=0, sd=1) #conditional probability of y_i=0 | x_i
  phi2 <- pnorm(beta * x - alpha2, mean=0, sd=1) #conditional probability of y_i=2 | x_i
  phi1 <- 1 - phi0 - phi2 #conditional probability of y_i=1 | x_i
  
  #construct m(x_i) function 
  m_xi <- (0)*phi0 + (1)*phi1 + (2)*phi2 #m(x_i) function designed based on nonlinear least square
  
  #generate the moment condition function
  g <- cbind(1, x, x^2)
  G <- colSums((y-m_xi)*g)/n #the final moment condition with dimension 3
  
  # Define the objective function 
  f= t(G) %*% G  #this is actually redundant: under just ID, just sum and optimize give you the MoM fast approach. 
  return(f)
  if (FALSE) {
    OP_GMMe <- function(y,x,par) {
      
      n = length(y) 
      alpha1 = par[1]
      alpha2 = par[2]
      beta = par[3]
      
      #generate the moment condition function
      G <- cbind(1, x^1, x^2) 
      prem <- y - pnorm(beta * x - alpha1,mean=0,sd=1) - pnorm(beta * x - alpha2,mean=0,sd=1)
      rows <- lapply(prem, rep, times = 3)
      prem_matrix <- do.call(rbind, rows)
      H = G  *  prem_matrix
      # Calculate the column sums of H
      HS <- colSums(H)/n
      
      # Define the objective function 
      f= t(HS) %*% HS
      return(f)
    }}
}





