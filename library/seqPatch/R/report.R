rl.detect.report <- function(rl.detect, out.file.prefix = './detect')
{
	if (any(names(rl.detect$pos) != names(rl.detect$neg)) | any(names(rl.detect$genome.start.pos) != names(rl.detect$genome.start.neg)))
		stop('incompatible forward strand and backward strand names')
	for (i in 1:length(rl.detect$pos)){
		out.file <- paste(out.file.prefix, '_', names(rl.detect$pos[i]),'.txt', sep='')
		cur.locus <- c(rl.detect$genome.start.pos[[i]] + 1:length(rl.detect$pos[[i]]$prop) - 1,
				rl.detect$genome.start.neg[[i]] + 1:length(rl.detect$neg[[i]]$prop) - 1)
		cur.prop <- c(rl.detect$pos[[i]]$prop, rl.detect$neg[[i]]$prop)
		cur.strand <- c(rep(0,length(rl.detect$pos[[i]]$prop)), rep(1,length(rl.detect$neg[[i]]$prop)))
		cur.N_1 <- c(rl.detect$pos[[i]]$N_1, rl.detect$neg[[i]]$N_1)
		cur.N_0 <- c(rl.detect$pos[[i]]$N_0, rl.detect$neg[[i]]$N_0)	
		cur.avg_n <- c(rl.detect$pos[[i]]$avg_n, rl.detect$neg[[i]]$avg_n)
		cur.cvg_wga <- c(rl.detect$pos[[i]]$cvg_wga, rl.detect$neg[[i]]$cvg_wga)
		
		rl.report <- cbind(cur.locus, cur.strand, cur.prop, cur.N_1, cur.N_0, cur.avg_n, cur.cvg_wga) 
		rl.report <- as.data.frame(rl.report)
		names(rl.report) <- c('locus', 'strand', 'prop', 'N_1', 'N_0', 'avg_n', 'cvg_wga')	
		
		rl.report <- rl.report[order(rl.report$locus),]
		write.table(rl.report, file = out.file, row.names = FALSE, quote = FALSE, sep = '\t')
	}
}


