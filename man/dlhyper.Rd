\name{dlhyper}
\alias{dlhyper}
\title{Function to tune the hyperparameter in the prior distribtuion}
\usage{
dlhyper<-function(x,y)
}

\arguments{

\item{y}{Response variable, a n*1 vector. Same as the argument in dlmain}
\item{x}{input matrix, each row is an observation vector, dimension n*p. Same as the argument in dlmain}
}

\value{
\item{hyper}{A value that can use in the following posterior computation}
}


\description{This function is to tune the value of hyperparameter in the prior, which can be [1/max(n,p),1/2]. We use the method proposed by Zhang et al. (2018). This method tune the hyperparameter by incorporating a prior on R^2. And they give a direct way to minimize KL directed divergence for special condition.}
\examples{
  rho=0.5
  p=1000
  n=100
  #set up correlation matrix
  m<-matrix(NA,p,p)
  for(i in 1:p){
    for(j in i:p)
      m[i,j]=rho^(j-i)}
  m[lower.tri(m)]<-t(m)[lower.tri(m)]
  #generate x
  library("mvtnorm")
  x=rmvnorm(n,mean=rep(0,p),sigma=m)
  #generate beta
  beta=c(rep(0,10),runif(n=5,min=-1,max=1),rep(0,20),runif(n=5,min=-1,max=1),rep(0,p-40))
  #generate y
  y=x\%*\%beta+rnorm(n)
  hyper=dlhyper(x,y)
}
