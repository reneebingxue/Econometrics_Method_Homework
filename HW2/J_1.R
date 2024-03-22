J_1 <- function(y,x,par) {
  
  n = length(y) 
  alpha1 = par[1]
  alpha2 = par[2]
  beta = par[3]
  
  Phi1 = pnorm(alpha1 - x*beta)
  Phi2 = pnorm(alpha2 - x*beta)
  
  phi1 = dnorm(alpha1 - x*beta)
  phi2 = dnorm(alpha2 - x*beta)
  
  J_1 = matrix(0,k,k)
  
  for (i in 1:n) {
    f1 = (y[i]==0)*phi1[i]/Phi1[i] - (y[i]==1)*phi1[i]/(Phi2[i]-Phi1[i])
    f2 = (y[i]==1)*phi2[i]/(Phi2[i]-Phi1[i]) - (y[i]==2)*phi2[i]/(1-Phi2[i])
    f3 = -(y[i]==0)*phi1[i]*x[i]/Phi1[i] + (y[i]==1)*(phi1[i]-phi2[i])*x[i]/(Phi2[i]-Phi1[i]) + (y[i]==2)*phi2[i]*x[i]/(1-Phi2[i])
    s = rbind(f1,f2,f3)
    j_1 = s %*% t(s)
    J_1 = J_1 + j_1
  }
  
  J = J_1/n
  
  return(J)
}
