dlprior<-function(hyper=1/2,p=10000000,plt=TRUE,min=-5,max=5,sigma=1){
  #dirichlet-laplace
  #prior for beta: beta[j] ~ N(0, sigma^2*psi[j]*phi[j]^2*tau^2)
  #phi ~ Dir(a,..., a)
  #tau ~ gamma(p*a, 1/2)
  #psi[j] iid~ exp(1/2)
  a=rep(hyper,p)
  psi=stats::rexp(p,1/2)
  phi=LaplacesDemon::rdirichlet(1,a)
  tau=stats::rgamma(1,shape=p*a,rate=1/2)
  beta=rep(1,p)
  for(j in 1:p){beta[j]=rnorm(1,sd=sigma*sqrt(psi[j]*(phi[j]^2)*(tau^2)))}
  if(plt==TRUE){
  d1=density(beta,from=min,to=max)
  plot(d1)
  }
  return(beta)
}







