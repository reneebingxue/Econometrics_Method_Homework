J_hessian <- function (par) {
  
  T <- length(X)
  
  alpha1 = par[1]
  alpha2 = par[2]
  omega = par[3]
  
  h11=0
  h22=0
  h33=0
  h12=0
  h13=0
  h23=0
  
  H=matrix(rep(0,9),3)
  
  for (t in 3:T) {
    
    com = (omega+alpha1*X[t-1]^2+alpha2*X[t-2]^2)
    
    h11=h11-(1/2)*(-1/(com)^2 +(2*(X[t])^2)/(com^3))*X[t-1]^4
    h22=h22-(1/2)*(-1/(com)^2 +(2*(X[t])^2)/(com^3))*X[t-2]^4
    h33=h33-(1/2)*(-1/(com)^2 +(2*(X[t])^2)/(com^3))
    
    h13=h13-(1/2)*(-1/(com)^2 +(2*(X[t])^2)/(com^3))*X[t-1]^2
    h12=h12-(1/2)*(-1/(com)^2 +(2*(X[t])^2)/(com^3))*X[t-1]^2*X[t-2]^2
    h23=h23-(1/2)*(-1/(com)^2 +(2*(X[t])^2)/(com^3))*X[t-2]^2
  }
  
  h1=cbind(h11,h12,h13)
  h2=cbind(h12,h22,h23)
  h3=cbind(h13,h23,h33)
  H = rbind(h1,h2,h3)
  
  H=H/(T-2)
  
  return(H) 
}
