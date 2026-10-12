library(tidyverse)

##you will need to change this file name
SNPF <-read.table("SNP.accuracy.freqXrayA.txt", header=T, fill = TRUE)


GrannySNP <- function(inputfile, founderfile, outdir){ 
  yy <- read.table(file = founderfile, header = FALSE, sep = "\t")
  founderroots <- yy[,1] 
  
  #load the SNP data
  if (is.data.frame(inputfile)) {
    xx <- inputfile
  } else if (is.character(inputfile)) {
    xx <- read.table(file = inputfile, header = TRUE)
  } else {
    stop("inputfile should be a data frame OR a file")
  }

  names(xx)[1] <- "CHROM"
  names(xx)[2] <- "POS"

  ##load the allele frequencies in the founders, based on Illumina data  
  founders <- match(paste("freq_", founderroots, sep=''), names(xx))
  if (any(is.na(founders))) {
    stop("Some founder columns are not found in the input data.")
  }
 
  Nf <- length(founders)
 #count the number of founders where this SNP cannot be called (it is NA)
  myna <- apply(xx[,founders], 1, function(x) sum(is.na(x)))
#keep only the SNPs that can be called in every founder
  xxx <- xx[myna == 0, ]
#estimate heterzygosity in each founder
  conthet <- xxx[,founders] * (1.0 - xxx[,founders])
#calculate mean heterozygosity across the founders
  mcon <- apply(conthet, 1, mean)
#exclude SNPs where heterozygosity is greater than 1% on average
  xxxx <- xxx[mcon < 0.01, ]     
  filename <- paste(outdir, "SNP.accuracy.freqXrayA1.txt", sep = '')
#keep the good SNPs that even your grandmother would believe.
  write.table(xxxx, filename, row.names = FALSE)
}

GrannySNP(SNPF1, "ff1.txt", "OUT")
