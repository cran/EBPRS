#' Read plink bfiles to R and reformat
#' @param path path to the test files (plink bfiles, without extension)
#' @description
#' To read plink files into R and transfer the files to the format we use in the EBPRS() function.
#' @details
#' The input should not include the extension. For example, the test files are AA.bed,
#' AA.bim and AA.fam, then the input should be 'AA' instead of 'AA.bed'.
#' @references
#' Song S, Jiang W, Hou L, Zhao H (2020) Leveraging effect size distributions to improve polygenic risk scores derived from summary statistics of genome-wide association studies. PLoS Comput Biol 16(2): e1007565. https://doi.org/10.1371/journal.pcbi.1007565
#' @seealso
#' \code{\link{EBPRS}}
#' @author
#' Shuang Song, Wei Jiang, Lin Hou and Hongyu Zhao
#' @import data.table BEDMatrix
#' @export

read_plink <- function(path){
  file <- list()
  fam <- fread(paste0(path,'.fam'))
  bim <- fread(paste0(path,'.bim'))
  bed <- as.matrix(BEDMatrix(path))
  file$fam <- fam
  file$bim <- bim
  file$bed <- bed
  return(file)
}


