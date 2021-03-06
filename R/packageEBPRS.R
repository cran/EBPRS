#' Description of the package
#'
#' @description
#' Description of the package. This is the 2.0.3 version.
#' @details
#' EB-PRS is a novel method that leverages information for effect sizes across all the markers to improve the prediction accuracy.  No parameter tuning is needed in the method, and no external information is needed. This R-package provides the calculation of polygenic risk scores from the given training summary statistics and test data. We can use EB-PRS to extract main information, estimate Empirical Bayes parameters, derive polygenic risk scores for  each individual in test data, and evaluate the PRS according to AUC and predictive r2.
#' \tabular{ll}{
#' Package: \tab EBPRS\cr
#' Type: \tab Package\cr
#' Date: \tab 2019-12\cr
#' Version: \tab 2.1.0\cr
#' }
#'
#' The package contains three main functions for users,\code{read_plink}, \code{EBPRS}, and \code{validate}.
#'
#' 1. \code{read_plink}. Thie function is used to read plink bfiles into R and reformat to suit the input of function \code{EBPRS()}.
#'
#' 2. \code{EBPRS}. This function integrate three parts: (1) merge the train and test (if have)
#' data, (2) estimate effectsize (3) generate polygenic risk scores (if test data provided.)
#'
#'
#' There is a strict requirement for the format of imput, which is
#' detailedly illustrated in details in function \code{EBPRS()}. The training summary statistics are necessary.
#' The test data can either
#' be included in the input or not. If test data are provided. The function will first
#' merge the data, as well as generate scores for each person in the result.
#' Users could first use the function \code{read_plink()} implemented in our package to read plink files into R.
#'
#'
#' 3. \code{validate}. We use this to validate the performance of the PRS.
#'
#'
#' 4. \code{data("traindat")} for the example training dataset.
#'
#'
#' A complete pipeline can be:
#'
#' train <- fread('trainpath')  (pay attention to the format, detailed in \code{EBPRS()})
#'
#' test <- read_plink('testpath')   (path to the plink bfile without extensions)
#'
#' result <- EBPRS(train=traindat, test=plinkfile, N1, N0)
#'
#' validate(result$S, truey)
#'
#' or
#'
#' train <- fread('trainpath')  (pay attention to the format)
#'
#' result <- EBPRS(train=traindat, N1, N0)  (will only provide estimated effect sizes)
#'
#'
#' @references
#' Song S, Jiang W, Hou L, Zhao H (2020) Leveraging effect size distributions to improve polygenic risk scores derived from summary statistics of genome-wide association studies. PLoS Comput Biol 16(2): e1007565. https://doi.org/10.1371/journal.pcbi.1007565
#' @seealso
#' \code{\link{EBPRS}},
#' \code{\link{validate}},
#'
#' @author
#' Shuang Song, Wei Jiang, Lin Hou and Hongyu Zhao
#' @export

EBPRSpackage <- function(){}
