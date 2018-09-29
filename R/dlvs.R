
dlvs<-function(dlresult){
  p=ncol(dlresult)
  betamean<-apply(dlresult,2,mean)
  betacov<-cov(dlresult)
  D<-betamean^2
  cov1=expm::sqrtm(betacov)
  covinv<-ginv(cov1)
  xstar=t(as.numeric(D)*covinv)
  ystar=covinv%*%betamean
  xstar<-scale(xstar, scale=F)
  ystar<-ystar-mean(ystar)
  model<-glmnet(xstar,ystar,standardize=FALSE,alpha=1,family="gaussian")
  lam=cv.glmnet(xstar,ystar,standardize=FALSE,alpha=1,family="gaussian")$lambda.1se
  betastar=coef(model,s=lam)[2:(p+1)]
  betatil=D*betastar
  return(betatil)
}
