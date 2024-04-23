score_gradient<- function (par) {
  T <- length(X)
  alpha1 = par[1]
  alpha2 = par[2]
  omega = par[3]
  g1=0
  g2=0
  g3=0
  for (t in 3:T) {
    
    # g1=g1-(1/2)*((1/(omega+alpha1*X[t-1]^2+alpha2*X[t-2]^2) -(X[t])^2/((omega+alpha1*X[t-1]^2+alpha2*X[t-2]^2)^2))*X[t-1]^2)
    # g2=g2-(1/2)*((1/(omega+alpha1*X[t-1]^2+alpha2*X[t-2]^2) -(X[t])^2/((omega+alpha1*X[t-1]^2+alpha2*X[t-2]^2)^2))*X[t-2]^2)
    # g3=g3-(1/2)*(1/(omega+alpha1*X[t-1]^2+alpha2*X[t-2]^2) -(X[t])^2/((omega+alpha1*X[t-1]^2+alpha2*X[t-2]^2)^2))
    
    com = omega+alpha1*X[t-1]^2+alpha2*X[t-2]^2
    g1=g1-(1/2)*((1/(com) -(X[t])^2/((com)^2))*X[t-1]^2)
    g2=g2-(1/2)*((1/(com) -(X[t])^2/((com)^2))*X[t-2]^2)
    g3=g3-(1/2)*(1/(com) -(X[t])^2/(com)^2)
  }
  g1=g1/(T-2)
  g2=g2/(T-2)
  g3=g3/(T-2)
  
  G = cbind(g1,g2,g3)
  return(G) 
}
