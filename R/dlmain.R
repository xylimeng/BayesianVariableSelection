#' Function to implement the Dirichlet Laplace shrinkage prior in Bayesian linear regression
#'
#' This function is the baysian linear regression version of the algorithm proposed in
#' Bhattacharya et al. (2015). The function is fast because we use fast sampling method
#' compute posterior samples. The method proposed in Bhattacharya et al. (2015) is used
#' in the second step perfectly solving the large p problem. The local shrinkage
#' controlling parameter psi_j are updated via a slice sampling scheme given by
#' Polson et al. (2014). And the parameters phi_j have various inverse gaussian
#' distribution. We generate variates with transformation into multiple roots
#' by Michael et al. (1976).
#'
#' @param x input matrix, each row is an observation vector, dimension n*p.
#' @param y Response variable, a n*1 vector.
#' @param burn Number of burn-in MCMC samples. Default is 5000.
#' @param nmc Number of posterior draws to be saved. Default is 5000.
#' @param thin Thinning parameter of the chain. Default is 1 means no thinning.
#' @param hyper The value of hyperparameter in the prior, can be [1/max(n,p),1/2].
#' It controls local shrinkage scales through psi. Small values of hyperparameter would
#' lead most of the result close to zero; while large values allow small singularity at
#' zero. We give a method and a function to tuning this parameter. See the function
#' called "dlhyper" for details.
#'
#'
#' @return \item{betamatrix}{Posterior samples of beta. A large matrix (nmc/thin)*p}
#'
#' @examples \dontrun{
#' library("mvtnorm")
#' rho=0.5
#' p=1000
#' n=100
#' #set up correlation matrix
#' m<-matrix(NA,p,p)
#' for(i in 1:p){
#'   for(j in i:p)
#'     m[i,j]=rho^(j-i)}
#' m[lower.tri(m)]<-t(m)[lower.tri(m)]
#' #generate x
#' library("mvtnorm")
#' x=rmvnorm(n,mean=rep(0,p),sigma=m)
#' #generate beta
#' beta=c(rep(0,10),runif(n=5,min=-1,max=1),rep(0,20),runif(n=5,min=-1,max=1),rep(0,p-40))
#' #generate y
#' y=x%*%beta+rnorm(n)
#' hyper=dlhyper(x,y)
#' dlresult=dl(x,y,hyper=hyper)}
#'
#'
#' @export
#'


dl<-function(x,y,burn=5000,nmc=5000,thin=1,hyper=1/2){
  p=ncol(x)
  n=nrow(x)
  #calculate hyperparameter
  xtx=t(x)%*%x
  #initial parameters
  a=rep(hyper,p)
  psi=stats::rexp(p,rate=1/2)
  psi1=rep(0,p)
  phi=LaplacesDemon::rdirichlet(n=1,alpha=a)
  phi[phi <= (1e-40)]<-(1e-40)
  tau=stats::rgamma(n=1,shape=p*a,rate=1/2)
  Ti=rep(1,p)
  beta=rep(0,p)
  hi=rep(1,p)
  ti=1
  betamatrix<-matrix(rep(NA,nmc*p),nrow=5000)

  #Niter iterations
  for(i in 1:(burn+nmc)){
    #step1:sample sigma^2
    s=c(psi*phi^2*tau^2)
    E_1=max(t(y-x%*%beta)%*%(y-x%*%beta),1e-8)
    E_2=max(sum(beta^2*s),1e-8)
    sigma2=1/stats::rgamma(1,(n+p)/2,rate=(E_1+E_2)/2)
    sigma1=sqrt(sigma2)
    if(sigma1>1e20) print("Please choose a better hyperparameter, it is too big")

    #step2:sample beta
    u=stats::rnorm(p)*sqrt(s)
    delta=stats::rnorm(n)
    v=x%*%u+delta
    stx=as.numeric(s)*t(x)
    w=ginv(x%*%stx+diag(n))%*%(y/sigma1-v)
    beta=(u+(stx%*%w))*sigma1

    mix1=abs(beta)/sigma1
    mix2=mix1/c(phi)
    #step3:sample psi
    mu=tau/mix2
    pv=(stats::rnorm(p))^2
    pu=stats::runif(p)
    temp2=mu*pv
    temp3=sqrt(4*temp2+temp2^2)
    temp4=mu+0.5*(pv*(mu^2))-0.5*(mu*temp3)
    locs=(pu<=mu/(mu+temp4))
    psi1=locs*temp4+(1-locs)*(mu^2/temp4)
    psi=1/psi1

    #step4:sample tau
    tau=GIGrvg::rgig(n=1,lambda=p*a-p,psi=1,chi=2*sum(mix2))

    #step5:sample phi
    hu=stats::runif(p,0,exp(-1/(2*hi)))
    hl=1/(2*log(1/hu))
    hf=stats::pgamma(hl,shape=1-a,rate=mix1)
    hr=pmin(runif(p,hf,1),(1-(1e-20)))
    hi=stats::qgamma(hr,shape=1-a,rate=mix1)
    Ti=1/hi
    phi=Ti/sum(Ti)
    phi[phi<=(1e-20)]=(1e-20)

    #beta output
    if(i>burn&&i%%thin==0) betamatrix[(i-burn)/thin,]=beta
    if(i%%1000==0) print(i)
  }
  return(betamatrix)
}
