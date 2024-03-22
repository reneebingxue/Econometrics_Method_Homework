OP_LL <- function(y,x,par) {
  
  n = length(y) 
  alpha1 = par[1]
  alpha2 = par[2]
  beta = par[3]
	
  Phi1 = pnorm(alpha1 - x*beta)
  Phi2 = pnorm(alpha2 - x*beta)
  
  f = mean((y==0)*log(Phi1) + (y==1)*log(Phi2-Phi1) + (y==2)*log(1-Phi2))

	return(-f)
}
