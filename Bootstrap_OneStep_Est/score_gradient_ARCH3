score_gradient_ARCH3 <- function (par) {
  
  T <- length(X)
  alpha1 = par[1]
  alpha2 = par[2]
  alpha3 = par[3]
  omega = par[4]
  g1=0
  g2=0
  g3=0
  g4=0
  for (t in 4:T) {
    com = omega+alpha1*X[t-1]^2+alpha2*X[t-2]^2+alpha3*X[t-3]^2
    g1=g1-(1/2)*((1/(com) -(X[t])^2/((com)^2))*X[t-1]^2)
    g2=g2-(1/2)*((1/(com) -(X[t])^2/((com)^2))*X[t-2]^2)
    g3=g3-(1/2)*((1/(com) -(X[t])^2/((com)^2))*X[t-3]^2)
    g4=g4-(1/2)*(1/(com) -(X[t])^2/(com)^2)
  }
  g1=g1/(T-3)
  g2=g2/(T-3)
  g3=g3/(T-3)
  g4=g4/(T-3)
  G = cbind(g1,g2,g3,g4)
  return(G) 
}
