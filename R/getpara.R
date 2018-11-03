
getpara <- function(n0,n1,SE,beta,pval,pi0Hat,sigma02){
  n <- n0+n1
  m <- length(pval)
  c <- 2*sqrt(n0*n1/n)
  lfdr <- rep(0,m)
  sigmaHat2 <- rep(0,m)
  f1 <- rep(0,m)
  f0 <- rep(0,m)
  # pval <- rep(0,m)
  z <- rep(0,m)
  k <- n1/n
  muHat <- rep(0,m)
  for(i in 1:m){
    #root <- get_root(k=k,beta=beta[i],ff=f[i])
    # f1[i] <- root$root1
    # f0[i] <- root$root0
    # x0 <- rbinom(n0,2,f0Hat[i])
    # x1 <- rbinom(n1,2,f1Hat[i])
    sigmaHat2[i]<- 4/c^2/(SE[i])^2
    if(beta[i]==0){
      z[i] <- 0
    }
    if(beta[i]!=0){

      z[i] <- -qnorm(pval[i]/2)*beta[i]/abs(beta[i])
    }

    #c*(f1[i]-f0[i])/sqrt(sigmaHat2[i])
    # pval[i] <- 2*pnorm(-abs(z[i]))
    # if((length(root$root1)>1)||(length(root$root0)>1)){
    #     print(i)
    # }

  }
  # if(missing(pi0Hat)){
  #   pi0Hat <- snpNullS(pval)$pi0
  #   if(pi0Hat==1){
  #     print("warning: pi0Hat=1")
  #     pi0Hat=0.99
  #   }
  # }
  sigma0Hat2 <- sigma02/c^2
  for(j in 1:m ){
    lfdr[j] <- pi0Hat*dnorm(z[j])/(pi0Hat*dnorm(z[j])+
                                     (1-pi0Hat)*dnorm(z[j]/sqrt(1+c^2*sigma0Hat2)))
    muHat[j] <- (1-lfdr[j])*c*sigma0Hat2*z[j]/(1+c^2*sigma0Hat2)
  }

  return(list(lfdr=lfdr,muHat=muHat,pi0Hat=pi0Hat,sigma0Hat2=sigma0Hat2,sigmaHat2=sigmaHat2))

}

