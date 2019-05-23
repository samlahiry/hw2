algo_leverage=function(X,y,r,type)## function to estimate beta with sampling probability pi and subsample size r
{
  if(missing(type))
    type = "unif"
  n=length(y)
  if(type=="unif")
    pi=rep(1/n,n)
  if(type=="lev")
    pi=diag(X%*%solve(t(X)%*%X)%*%t(X))
  samp=sample(1:n,size=r,replace=FALSE,prob=pi)
  Phi=diag(1/sqrt(pi[samp]))
  fit=lm((Phi%*%y[samp])~0+(Phi%*%X[samp,]) ) ##Weighted Leveraging
  return(as.vector(fit$coeff))
}
