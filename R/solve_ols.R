#' Solving linear systems of the form Ax=b
#'
#' @name solve_ols
#' @param A matrix of the system to be solved
#' @param b vector of the system to be solved
#' @param method Gauss, Jacobi or Jacobi parallel
#' @param tol tolerance level of relative errors (default 1e-10)
#' @param maxiter maximum number of iterations (default 10000)
#' @param core number of core used for parallel computing (only valid for Jacobi Parallel)
#'
#' @return The solution of the system
#' @export

solve_ols=function(A,b,method,tol,maxiter,core)
{
  if(missing(tol))
    tol=1e-10
  if(missing(maxiter))
    maxiter=10000
  if(missing(core))
    core=3
  norm=function(x)sqrt(sum(x^2))
  n=dim(A)[1]
  v_0=numeric(n)
  D=diag(diag(A))
  L=matrix(0,n,n)
  U=matrix(0,n,n)
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      L[i,j]=A[i,j]*(i>j)
      U[i,j]=A[i,j]*(i<j)
    }
  }
  err=100000
  iter=0
  if(method=="Gauss")
  {
  v=v_0
  ptm=proc.time()
  while(err>tol & iter<maxiter)
  {
   gs=solve(L+D)%*%(b-U%*%v)
   err=norm(gs-v)/norm(v)
   iter=iter+1
   v=gs
  }
  print(proc.time()-ptm)
  return(gs)
  }
  if(method=="Jacobi")
  {
    v=v_0
    ptm=proc.time()
    while(err>tol & iter<maxiter)
    {
      jseq=diag(diag(A)^(-1))%*%(b-(L+U)%*%v)
      err=norm(jseq-v)/norm(v)
      iter=iter+1
      v=jseq
    }
    print(proc.time()-ptm)
    return(jseq)
  }
  if(method=="Jacobi-parallel")
  {
    v=v_0
    ptm=proc.time()
    require(doParallel)
    cl= makeCluster(core)
    registerDoParallel(cl)
    ptm=proc.time()
    while(err>tol & iter<maxiter)
    {

        vec1=foreach(l=1:n) %dopar%
        {(L+U)[l,]%*%v}
        jpar=diag(diag(A)^(-1))%*%(b-unlist(vec1))
        err=norm(jpar-v)/norm(v) #relative error of Jacobi(parallel)
        iter=iter+1
        v=jpar
    }
    print(proc.time()-ptm)
    return(jpar)
  }
}
