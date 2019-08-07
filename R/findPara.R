#' Derive the parameters
#'
#' @param train train set after processed by 'extractInfo'
#' @param bed bed file after processed by 'extractInfo'
#' @param bim bim file after processed by 'extractInfo'
#' @param fam fam file after processed by 'extractInfo'
#' @param N1 case number
#' @param N0 control number
#' @return A list including
#' estimated mu (muHat)
#' estimated sigma2 (sigmaHat2)
#' estimated proportion of non-associated SNPs (pi0)
#' estimated variance of effect sizes of associated SNPs (sigma02)
#' @description
#' All the input files can be generated from 'extractInfo'
#' @seealso
#' \code{\link{extractInfo}}
#' @references
#' Song, S., Jiang, W., Hou, L. and Zhao, H. Leveraging effect size distributions to improve polygenic risk scores derived from genome-wide association studies. \emph{Submitted}.
#' @author
#' Shuang Song, Wei Jiang, Lin Hou and Hongyu Zhao
#' @importFrom utils write.table
#' @export


findPara <- function(train,bed,bim,fam,N1,N0){
  #colnames(train)[1:7] <- c("chr","snpid","a1","a2","pos","or","p")
  print("Coordinating the ref alleles...")
  colnames(train)[which(colnames(train)=="a1")] <- "A1"
  colnames(train)[which(colnames(train)=="a2")] <- "A2"
  colnames(train)[which(colnames(train)=="or")] <- "OR"
  colnames(train)[which(colnames(train)=="p")] <- "P"
  a=train$A1
  b=bim[,5]
  sig=agtc(a,b)
  train$OR=(train$OR)^sig

  z <- -qnorm(train$P/2)*log(train$OR)/abs(log(train$OR))
  se <- log(train$OR)/z
  z[which(is.na(z))] <- 0
  se[which(is.na(se))] <- 1

  #res <- findpara_RD1(n0=N0,n1=N1,SE=se,beta=log(train1$or),pval=train1$p)
  print("Utilizing EM algorithm to derive pi0 and sigma02")
  temp <- snpEM(z)
  if(temp$pi0==1){
    print("Warning: pi0=1")
    temp$pi0 <- 0.999
  }
  print("Generating the parameters")
  res <- getpara(n0=N0,n1=N1,SE=se,beta=log(train$OR),
                                pval=train$P,pi0Hat=temp$pi0,sigma02=temp$sigma2)


  #res <- findpara_RD1_pi0(n0=N0,n1=N1,SE=se,beta=log(train1$or),pval=train1$p,pi0=0.9984)
  result=data.frame(muHat=res$muHat,sigmaHat2=res$sigmaHat2)
  write.table(result,"res_para.txt")
  print("Completed.")
  return(list(muHat=res$muHat,sigmaHat2=res$sigmaHat2,pi0=temp$pi0,sigma02=temp$sigma02))

}
