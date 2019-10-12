#' Generate effect size from training set
#'
#' @param train train dataset
#' @param N1 case number
#' @param N0 control number
#' @return The effect sizes for each SNP
#' @description
#' Only the training set is needed. This function is designed for the condition
#' that testing data is very large and hard to be loaded into R. The effect size
#' will be generated directly. Users can calculate scores in plink with the generated
#' effect size.
#'
#' An example training dataset can be acquired using data("traindat")
#' @details
#' The raw training data should be a file with
#' 8 columns including CHROM, POS, A1, A2, OR, P, SNP, N in order.
#' The CHROM column and the SNP column is used for indexing.
#'
#' @examples
#' data("traindat")
#'  \dontrun{
#'  getEffectSize(traindat, N1=364, N0=2063)
#'  }
#'
#' @references
#' Song, S., Jiang, W., Hou, L. and Zhao, H. Leveraging effect size distributions to improve polygenic risk scores derived from genome-wide association studies. \emph{Submitted}.
#' @author
#' Shuang Song, Wei Jiang, Lin Hou and Hongyu Zhao
#' @export

getEffectSize <- function(train,N1,N0){

  #train <- fread(trainpath,header=T)
  #colnames(train)[1:7] <- c("chr","snpid","a1","a2","pos","or","p")
  print("Coordinating the ref alleles...")
  #a=train$a1
  # b=bim[,5]
  # sig=agtc(a,b)
  #train$or=(train$or)^sig
  colnames(train)[which(colnames(train)=="a1")] <- "A1"
  colnames(train)[which(colnames(train)=="a2")] <- "A2"
  colnames(train)[which(colnames(train)=="or")] <- "OR"
  colnames(train)[which(colnames(train)=="p")] <- "P"
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
  #result=data.frame(muHat=res$muHat,sigmaHat2=res$sigmaHat2)
  muHatnew <- res$muHat/sqrt(res$sigmaHat2)
  result <- data.frame(train,effectsize=muHatnew)

  write.table(result,"res_para.txt")
  print("Completed.")
  return(result)
}
