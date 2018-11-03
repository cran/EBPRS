#' Extract information from raw data
#'
#' @param trainpath train dataset path
#' @param test test dataset(list) including fam, bed, bim(generated from plink files, plink2R::read_plink is recommended)
#' @return A list including processed training data (train) and testing data (bed, bim, fam)
#' @description
#' The first step of the algorithm, to clean the dataset and extract information from raw data.
#' (Please notice that there are some requirements for the training and testing datasets.)
#' @import data.table
#' @details
#' The raw training data should be a file with
#' 8 columns including CHROM, POS, A1, A2, BETA, P, SNP, N in order.
#' The CHROM column should only be a number from 1 to 22. The SNP column
#' is the rsid number.
#'
#' "test" file can be generated from read_plink("test_plink_file")
#' The raw testing data could be the files transformed from plink2R (using plink bfiles).
#'
#' test is a list including fam (6 columns with information on samples), bim (6 columns with information on SNPs), bed (genotypes 0, 1, 2)
#' @references
#' Song, S., Jiang, W., Hou, L. and Zhao, H. Leveraging effect size distributions to improve polygenic risk scores derived from genome-wide association studies. \emph{Submitted}.
#' @seealso
#' \url{https://github.com/gabraham/plink2R}
#' @author
#' Shuang Song, Wei Jiang, Lin Hou and Hongyu Zhao
#' @export



extractInfo <- function(trainpath,test){
  # if (!is.installed("plink2R")){
  #   options(unzip = "internal")
  #   devtools::install_github("gabraham/plink2R/plink2R")
  #   #devtools::use_package("plink2R")
  # }
  #require("plink2R")
    #source("read_plink.R")
  print1 <- paste("Reading the training data from",trainpath)
  print(print1)

  train <- fread(trainpath,header=T)
  #print2 <- paste("Reading the testing data from",testpath)
  #print(print2)
  #devtools::install_github("gabraham/plink2R/plink2R")
  #test <- read_plink(testpath)
  print("Merging the data...")
  trainindex <- data.frame(snp=train[,2],indextrain=1:nrow(train))
  testindex <- data.frame(snp=test$bim[,2],indextest=1:nrow(test$bim))
  colnames(trainindex)[1] <- "snp"
  colnames(testindex)[1] <- "snp"
  index <- merge(trainindex,testindex,by="snp")
  index=cbind(index,1:nrow(index))
  #write.table(index,"index.txt")

  train1 <- train[index[,2],]
  bed1 <- test$bed[,index[,3]]
  bim1 <- test$bim[index[,3],]
  fam1 <- test$fam
  #rm(train)
  #rm(test)
  print("Processing NAs in bed file...")
  bed1 <- bedNA0(bed1)
  snpNum <- dim(train1)[1]
  peoNum <- dim(fam1)[1]

  Info <- list(train=train1,bed=bed1,bim=bim1,fam=fam1)
  print("Completed.")
  temp1 <- paste(snpNum,"SNPs are included in calculating PRS.")
  print(temp1)
  temp2 <- paste(peoNum,"individuals in total.")
  print(temp2)
  return(Info)

}
