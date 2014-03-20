library(Rcpp, quietly=TRUE)
library(locfdr, quietly=TRUE)
library(weights, quietly=TRUE)
library(seqPatch, quietly=TRUE)
library(pbh5, quietly=TRUE)

Args <- commandArgs(trailingOnly = TRUE)

cmpH5.file <- as.character(Args[1])
out.dir <- as.character(Args[2])
reagent <- as.character(Args[3])
mapQV.cutoff <- as.numeric(Args[4])

cat('start to pileup reads\n')
cat('cmpH5file is: ',cmpH5.file,'\n')
cat('outdir is: ', out.dir, '\n')
cat('reagent is: ', reagent, '\n')
cat('mapQVthreshold is: ', mapQV.cutoff, '\n')

loadCmpH5toGenomeF(cmpH5.file, out.dir, reagent=reagent, mapQV.cutoff=mapQV.cutoff)


