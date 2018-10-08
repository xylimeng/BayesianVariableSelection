dlanalysis<-function(dlresult,alpha=0.05){
  betamean=apply(dlresult,2,mean)
  betamedian=apply(dlresult,2,median)
  min=apply(dlresult,2,quantile,prob=alpha/2)
  max=apply(dlresult,2,quantile,prob=1-alpha/2)
  result=list("betamean"=betamean,"betamedian"=betamedian,"LeftCI"=min,"RightCI"=max)
  return(result)
}
