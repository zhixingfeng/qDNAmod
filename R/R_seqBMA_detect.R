library(Rcpp, quietly=TRUE)
library(locfdr, quietly=TRUE)
library(weights, quietly=TRUE)
library(seqPatch, quietly=TRUE)

Args <- commandArgs(trailingOnly = TRUE)

native.pileup.file <- as.character(Args[1])
wga.pileup.file <- as.character(Args[2])
prior.file <- as.character(Args[3])
out.dir <- as.character(Args[4])
n.iter <- as.integer(Args[5])

cat('start to calculate modification proportion\n')
cat('native pileup file is in: ',native.pileup.file,'\n')
cat('wga pileup file is in: ',wga.pileup.file,'\n')
cat('prior file is in: ', prior.file, '\n')
cat('output directory is: ', out.dir, '\n')
cat('number of iteration is: ', n.iter, '\n')

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

if (n.iter >= 2){
	detect.iter <- detectModPropEB.iter(genomeF.native, genomeF.wga, f1$x, f1$y, is_f1_var=FALSE, 100, is.EM=FALSE, 0.005, n.iter)
	f1$x <- detect.iter$f1.x
	f1$y <- detect.iter$f1.y
}
rl.detect <- detectModPropEB(genomeF.native, genomeF.wga, f1$x, f1$y, is_f1_var=TRUE, 100, is.EM=FALSE)

if (!file.exists(out.dir))
	dir.create(out.dir,recursive=TRUE)
if (n.iter >= 2){
	save(rl.detect, detect.iter, file=paste(out.dir,'/detect.Rdata',sep=''))
}else{
	save(rl.detect, file=paste(out.dir,'/detect.Rdata',sep=''))
}



