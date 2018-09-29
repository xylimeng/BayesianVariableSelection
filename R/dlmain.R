
dl<-function(x,y,burn=5000,nmc=5000,thin=1,hyper=1/2){
  p=ncol(x)
  n=nrow(x)
  #calculate hyperparameter
  xtx=t(x)%*%x
  #initial parameters
  a=rep(hyper,p)
  psi=rexp(p,rate=1/2)
  psi1=rep(0,p)
  phi=rdirichlet(n=1,alpha=a)
  phi[phi <= (1e-40)]<-(1e-40)
  tau=rgamma(n=1,shape=p*a,rate=1/2)
  Ti=rep(1,p)
  beta=rep(0,p)
  hi=rep(1,p)
  betamatrix<-matrix(rep(NA,nmc*p),nrow=5000)
  
  #Niter iterations
  for(i in 1:(burn+nmc)){
    #step1:sample sigma^2
    s=c(psi*phi^2*tau^2)
    E_1=max(t(y-x%*%beta)%*%(y-x%*%beta),1e-8)
    E_2=max(sum(beta^2*s),1e-8)
    sigma2=1/stats::rgamma(1,(n+p)/2,rate=(E_1+E_2)/2)
    sigma1=sqrt(sigma2)
    
    #step2:sample beta
    u=rnorm(p)*sqrt(s)
    delta=rnorm(n)
    v=x%*%u+delta
    stx=as.numeric(s)*t(x)
    w=ginv(x%*%stx+diag(n))%*%(y/sigma1-v)
    beta=(u+(stx%*%w))*sigma1
    
    mix1=abs(beta)/sigma1
    mix2=mix1/c(phi)
    #step3:sample psi
    mu=tau/mix2
    pv=(rnorm(p))^2
    pu=runif(p)
    temp2=mu*pv
    temp3=sqrt(4*temp2+temp2^2)
    temp4=mu+0.5*(pv*(mu^2))-0.5*(mu*temp3)
    locs=(pu<=mu/(mu+temp4))
    psi1=locs*temp4+(1-locs)*(mu^2/temp4)
    psi=1/psi1
    
    #step4:sample tau
    tau=GIGrvg::rgig(n=1,lambda=p*a-p,psi=1,chi=2*sum(mix2))
    
    #step5:sample phi
    hu=runif(p,0,exp(-1/(2*hi)))
    hl=1/(2*log(1/hu))
    hf=pgamma(hl,shape=1-a,rate=mix1)
    hr=pmin(runif(p,hf,1),(1-(1e-20)))
    hi=qgamma(hr,shape=1-a,rate=mix1)
    Ti=1/hi
    phi=Ti/sum(Ti)
    phi[phi<=(1e-20)]=(1e-20)
    
    #beta output
    if(i>burn&&i%%thin==0) betamatrix[(i-burn)/thin,]=beta
    if(i%%1000==0) print(i)
  }
  return(betamatrix)
}