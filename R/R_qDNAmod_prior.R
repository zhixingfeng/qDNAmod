library(Rcpp, quietly=TRUE)
library(locfdr, quietly=TRUE)
library(weights, quietly=TRUE)
library(seqPatch, quietly=TRUE)

Args <- commandArgs(trailingOnly = TRUE)

native.pileup.file <- as.character(Args[1])
wga.pileup.file <- as.character(Args[2])
out.dir <- as.character(Args[3])

cat('start to calculate prior\n')
cat('native pileup file is in: ',native.pileup.file,'\n')
cat('wga pileup file is in: ',wga.pileup.file,'\n')
cat('output directory is: ', out.dir, '\n')

if (!file.exists(paste(native.pileup.file,'/genomeF.Rdata',sep='')))
	stop(paste(native.pileup.file,'/genomeF.Rdata',sep=''),' does not exist')

if (!file.exists(paste(wga.pileup.file,'/genomeF.Rdata',sep='')))
	stop(paste(wga.pileup.file,'/genomeF.Rdata',sep=''),' does not exist')

cat('load ',paste(native.pileup.file,'/genomeF.Rdata',sep=''),'\n')
load(paste(native.pileup.file,'/genomeF.Rdata',sep=''))
genomeF.native <- genomeF
rm(genomeF);gc()

cat('load ',paste(wga.pileup.file,'/genomeF.Rdata',sep=''),'\n')
load(paste(wga.pileup.file,'/genomeF.Rdata',sep=''))
genomeF.wga <- genomeF
rm(genomeF);gc()


f1 <- estimate.f1(genomeF.native, genomeF.wga, out.dir=out.dir);
cat('save f1\n')
save(f1,file= paste(out.dir,'/f1.Rdata',sep=''))


