\name{dlvs}
\alias{dlvs}
\title{Function to do Bayesian variable selection via penalized credible region}
\usage{
dlvs(dlresult)
}

\arguments{
\item{dlresult}{Posterior samples of beta. A large matrix (nmc/thin)*p}
}

\value{
\item{betatil}{Variable selection result of beta, a p*1 vector. Most of the values shrinks to 0}
}



\description{
This is a function using the algorithm doing variable selection via penalized credible interval proposed by Bondell et al. (2012). The computation of the proposed sequence is doing matrix computing and using existing LASSO software.
}
\examples{
\dontrun{
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
  dlresult=dl(x,y,hyper=hyper)
  dlvs(dlresult)}
}
