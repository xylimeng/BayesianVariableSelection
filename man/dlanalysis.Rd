\name{dlanalysis}
\alias{dlanalysis}
\title{Function to analyse dirichlet laplace posterior by using gibbs sampling data}
\usage{
dlanalysis(dlresult,alpha=0.05)
}

\arguments{

\item{dlresult}{Posterior samples of beta. A large matrix (nmc/thin)*p}
\item{alpha}{Level for the credible intervals. For example,the default is alpha = 0.05 means 95\% credible intervals.}
}

\value{
\item{betamean}{Posterior mean of beta, a p*1 vector.}
\item{LeftCI}{The left bounds of the credible intervals.}
\item{RightCI}{The right bounds of the credible intervals.}
\item{betamedian}{Posterior median of Beta, a p*1 vector.}
}


\description{This is a function that analyse the MCMC sampling result by computing the posterior mean, median and credible intervals.}

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
  dlresult=dl(x,y,hyper=hyper)
  da=dlanalysis(dlresult,alpha=0.05)
  da$betamean
  da$betamedian
  da$LeftCI
  da$RightCI
}
