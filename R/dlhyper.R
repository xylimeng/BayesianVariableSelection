dlhyper<-function(x,y){
  p=ncol(x)
  n=nrow(x)
  #calculate hyperparameter
  xtx=t(x)%*%x
  d=eigen(xtx/n)$values
  P=sum(d)
  Q=4*sum(d^2)-sum(d)^2
  R=-sum(d)^3
  C=P^2/9-Q/3
  A=P*Q/6-P^3/27-R/2
  B=A^2-C^3
  hyper=sqrt(2/((A+sqrt(B))^(1/3)+sign(A-sqrt(B))*abs(A-sqrt(B))^(1/3)-P/3))
  hyper[hyper<1/p]=1/p
  return(hyper)
}



