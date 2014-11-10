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
is.two.round <- as.logical(Args[6])


cat('start to calculate modification proportion\n')
cat('native pileup file is in: ',native.pileup.file,'\n')
cat('wga pileup file is in: ',wga.pileup.file,'\n')
cat('prior file is in: ', prior.file, '\n')
cat('output directory is: ', out.dir, '\n')
cat('number of iteration is: ', n.iter, '\n')
cat('using two round estimation? ', is.two.round, '\n')


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
	detect.iter <- detectModPropEB.iter(genomeF.native, genomeF.wga, f1$x, f1$y, is_f1_var=FALSE, 50, is.EM=FALSE, 0.005, n.iter)
	f1$x <- detect.iter$f1.x
	f1$y <- detect.iter$f1.y
}

# use mixture model with fixed f1 to remove loci with low modification proportion
if (is.two.round == TRUE){
	rl.detect.1 <- detectModPropEB(genomeF.native, genomeF.wga, f1$x, f1$y, is_f1_var=FALSE, 25, is.EM=FALSE)
	# remove loci with prop < 0.01	
	for (i in 1:length(rl.detect.1$pos)){
		cur.idx <- which(rl.detect.1$pos[[i]]$prop < 0.01)
		cur.name <- names(rl.detect.1$pos[i])
		for (j in cur.idx){
			genomeF.native$features$ipd_pos[[cur.name]][[j]] <- numeric(0)
			genomeF.native$features$moleculeID_pos[[cur.name]][[j]] <- numeric(0)
		}
	}	
	for (i in 1:length(rl.detect.1$neg)){
                cur.idx <- which(rl.detect.1$neg[[i]]$prop < 0.01)
                cur.name <- names(rl.detect.1$neg[i])
                for (j in cur.idx){
                        genomeF.native$features$ipd_neg[[cur.name]][[j]] <- numeric(0)
                        genomeF.native$features$moleculeID_neg[[cur.name]][[j]] <- numeric(0)
                }
        }
}

# use mixture model with varible f1 to reduce bias
rl.detect.2 <- detectModPropEB(genomeF.native, genomeF.wga, f1$x, f1$y, is_f1_var=TRUE, 25, is.EM=FALSE)

# combine rl.detect.1 and rl.detect.2
for (i in 1:length(rl.detect.2$pos)){
	idx.0 <- which(is.na(rl.detect.2$pos[[i]]$prop) & rl.detect.2$pos[[i]]$n_mol==0)	
	rl.detect.2$pos[[i]]$prop[idx.0] <- 0	
}

for (i in 1:length(rl.detect.2$neg)){
        idx.0 <- which(is.na(rl.detect.2$neg[[i]]$prop) & rl.detect.2$neg[[i]]$n_mol==0)
        rl.detect.2$neg[[i]]$prop[idx.0] <- 0
}

rl.detect <- rl.detect.2

if (!file.exists(out.dir))
	dir.create(out.dir,recursive=TRUE)
if (n.iter >= 2){
	save(rl.detect, detect.iter, file=paste(out.dir,'/detect.Rdata',sep=''))
}else{
	save(rl.detect, file=paste(out.dir,'/detect.Rdata',sep=''))
	rl.detect.report(rl.detect, out.file.prefix = paste(out.dir,'/detect',sep=''))
}



