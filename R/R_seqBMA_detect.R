library(Rcpp, quietly=TRUE)
library(locfdr, quietly=TRUE)
library(weights, quietly=TRUE)
library(seqPatch, quietly=TRUE)

Args <- commandArgs(trailingOnly = TRUE)

native.pileup.file <- as.character(Args[1])
wga.pileup.file <- as.character(Args[2])
prior.file <- as.character(Args[3])
out.dir <- as.character(Args[4])

cat('start to calculate modification proportion\n')
cat('native pileup file is in: ',native.pileup.file,'\n')
cat('wga pileup file is in: ',wga.pileup.file,'\n')
cat('prior file is in: ', prior.file, '\n')
cat('output directory is: ', out.dir, '\n')

if (!file.exists(paste(native.pileup.file,'/genomeF.Rdata',sep='')))
        stop(paste(native.pileup.file,'/genomeF.Rdata',sep=''),' does not exist')

if (!file.exists(paste(wga.pileup.file,'/genomeF.Rdata',sep='')))
        stop(paste(wga.pileup.file,'/genomeF.Rdata',sep=''),' does not exist')

if (!file.exists(paste(prior.file,'/f1.Rdata',sep='')))
	stop(paste(prior.file,'/f1.Rdata',sep=''), ' does not exist')


cat('load ',paste(native.pileup.file,'/genomeF.Rdata',sep=''),'\n')
load(paste(native.pileup.file,'/genomeF.Rdata',sep=''))
genomeF.native <- genomeF
rm(genomeF);gc()

cat('load ',paste(wga.pileup.file,'/genomeF.Rdata',sep=''),'\n')
load(paste(wga.pileup.file,'/genomeF.Rdata',sep=''))
genomeF.wga <- genomeF
rm(genomeF);gc()

load(paste(prior.file,'/f1.Rdata',sep=''))

rl.detect <- detectModPropEB(genomeF.native, genomeF.wga, f1$x, f1$y, 50)

if (!file.exists(out.dir))
	dir.create(out.dir,recursive=TRUE)

save(rl.detect, file=paste(out.dir,'/detect.Rdata',sep=''))





