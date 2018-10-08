# dlbayes

The goal of dlbayes is to implement the Dirichlet Laplace shrinkage prior in Bayesian linear regression and variable selection, featuring: 
- utility functions in implementing Dirichlet-Lapace priors such as visualization; 
- scalability in Bayesian linear regression; 
- penalized credible regions for variable selection. 

## Installation

You can install the released version of dlbayes from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("dlbayes")
```

## Example

``` r
## example
  rho=0.5
  p=1000
  n=100
  #set up correlation matrix
  m<-matrix(NA,p,p)
  for(i in 1:p){
    for(j in i:p){
      m[i,j]=rho^(j-i)
    }
  }
  m[lower.tri(m)]<-t(m)[lower.tri(m)]
  #generate x
  library("mvtnorm")
  x=rmvnorm(n,mean=rep(0,p),sigma=m)
  #generate beta
  beta=c(rep(0,10),runif(n=5,min=-1,max=1),rep(0,20),runif(n=5,min=-1,max=1),rep(0,p-40))
  #generate y
  y=x%*%beta+rnorm(n)
  hyper=dlhyper(x,y)
  dlresult=dl(x,y,hyper=hyper)
  dlvs(dlresult)
  da=dlanalysis(dlresult,alpha=0.05)
  da$betamean
  da$betamedian
  da$LeftCI
  da$RightCI
  theta=dlprior(hyper=1/2,p=10000000,plt=TRUE,min=-5,max=5,sigma=1)
```

## Reference 



