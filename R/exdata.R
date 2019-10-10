#' Example data for training set
#'
#' Summary statistics simulated in the manuscript Leveraging effect size distributions to improve polygenic risk scores derived from genome-wide association studies. Data from a QTL experiment on gravitropism in
#'
#'
#' @docType data
#'
#' @usage data("traindat")
#'
#' @format An object of class \code{"cross"}; see \code{\link[qtl]{read.cross}}.
#'
#' @keywords datasets
#'
#' @references
#' Song, S., Jiang, W., Hou, L. and Zhao, H. Leveraging effect size distributions to improve polygenic risk scores derived from genome-wide association studies. \emph{Submitted}.
#'
#'
#' @examples
#' data("traindat")
#'  \dontrun{
#'  getEffectSize(traindat, N1=364, N0=2063)
#'  }
"traindat"