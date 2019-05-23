#' Solving elastic net via coordinate descent
#'
#' @name  elnet_coord 
#' @param y The response vector
#' @param X Design matrix
#' @param Lambda Sequence of penalizing weights 
#' @param alpha weightage on L_1 norm
#' @param beta_init initial coefficient vector
#' @param iterlength no. of iteration of the coordinate descent algorithm
#'
#' @return None
#' @export






elnet_coord=function(y,X,Lambda,alpha,beta_init,iterlength) ## coordinate descent elastic net
{
  st=function(x,lambda)## soft thresholding operator
  {
    return((x>=0)*(abs(x)-lambda)*(abs(x)>lambda)-(x<0)*(abs(x)-lambda)*(abs(x)>lambda))
  }
  p=ncol(X)
  n=nrow(X)
  solpath=matrix(0,nrow=length(Lambda),ncol=p)
  for(i in 1:length(Lambda))
  {
    beta=beta_init
    for(t in 1:iterlength)
    {
      for(j in 1:p)
      {
        z_j=y-X[,-j]%*%beta[-j]
        x_j=X[,j]
        beta[j]=st(t(x_j)%*%z_j/n,Lambda[i]*alpha)/(t(x_j)%*%x_j/n+Lambda[i]*(1-alpha))
      }
    }
    solpath[i,]=beta
  }
  return(solpath)
}



#solution paths
#plt=function(n,Lambda,alpha,beta_init,iterlength)## plots solution paths for different values of n and alpha
#{
  #Sol=elnet_coord(y,X,Lambda,alpha,beta_init,iterlength)
  #Sol1=cbind(Sol,apply(abs(Sol),1,sum))
  #dat=data.frame(number=as.character(1:10),t(Sol[,1:10]))
  #dat_melt=melt(dat, id.vars = 'number')
  #dat1=data.frame(dat_melt,norm=rep(Sol1[,21],each=10))
  #ggplot(dat1, aes(x = norm, y = value)) + geom_line(aes(color = number,group=number),size=1)
  
#}