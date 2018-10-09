dlvs<-function(dlresult){
  #use penalized credible region to do variable selection
  #matrix computation to make the solutions be accomplished by LASSO
  p=ncol(dlresult)
  betamean<-apply(dlresult,2,mean)
  betacov<-cov(dlresult)
  D<-betamean^2
  cov1=expm::sqrtm(betacov)
  covinv<-MASS::ginv(cov1)
  xstar=t(as.numeric(D)*covinv)
  ystar=covinv%*%betamean
  #scale
  xstar<-scale(xstar, scale=F)
  ystar<-ystar-mean(ystar)
  #solve LASSO problem by glmnet package
  model<-glmnet::glmnet(xstar,ystar,standardize=FALSE,alpha=1,family="gaussian")
  lam=glmnet::cv.glmnet(xstar,ystar,standardize=FALSE,alpha=1,family="gaussian")$lambda.1se
  betastar=coef(model,s=lam)[2:(p+1)]
  betatil=D*betastar
  return(betatil)
}
