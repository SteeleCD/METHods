#refFile="/media/DATA/Analysis_tools/GenomeAnalysisTK-3.3/resources/human_g1k_v37/human_g1k_v37.fasta"
#bamFile1="/media/DATA/ONB_project_exome/chris/docker/sampleM.bam"
#bamFile2="/media/DATA/ONB_project_exome/chris/docker/sampleN.bam"
runSequenza = function(refFile,bamFileM,bamFileN,outDir,sampleName="Test")
	{
	sequenzaFile=system.file("exec", "sequenza-utils.py", package="sequenza")
	dataDir = paste(outDir,"data",sep="/")
	plotsDir = paste(outDir,"plots",sep="/")
	infoDir = paste(outDir,"info",sep="/")
	seqOutDir = paste(outDir,"seqOut",sep="/")
	system(paste0("mkdir ",dataDir))
	system(paste0("mkdir ",plotsDir))
	system(paste0("mkdir ",seqOutDir))
	system(paste0("mkdir ",infoDir))
	# get pileup
	system(paste0("samtools mpileup -f ",refFile," -Q 20 ",bamFileM,"|gzip>",dataDir,"/mutant.pileup.gz")) 
	system(paste0("samtools mpileup -f ",refFile," -Q 20 ",bamFileN,"|gzip>",dataDir,"/normal.pileup.gz")) 
	system(paste0(sequenzaFile," GC-windows -w 50 ",refFile," | gzip > ",dataDir,"/human_g1k_v37.gc50Base.txt.gz"))
	system(paste0(sequenzaFile," pileup2seqz -gc ",dataDir,"/human_g1k_v37.gc50Base.txt.gz -n ",dataDir,"/normal.pileup.gz -t ",dataDir,"/mutant.pileup.gz | gzip > ",dataDir,"/out.seqz.gz"))
	system(paste0(sequenzaFile," seqz-binning -w 50 -s ",dataDir,"/out.seqz.gz | gzip > ",dataDir,"/out_small.seqz.gz"))
	# normalization of depth ratio
	library(sequenza)
	data.file=paste0(dataDir,"/out_small.seqz.gz")
	seqz.data=read.seqz(data.file)
	gc.stats = gc.norm(x = seqz.data$depth.ratio, gc = seqz.data$GC.percent)
	gc.vect = setNames(gc.stats$raw.mean, gc.stats$gc.values)
	seqz.data$adjusted.ratio = seqz.data$depth.ratio / gc.vect[as.character(seqz.data$GC.percent)]
	# plot adjustment
	pdf(paste0(plotsDir,"/plotAdjustment.pdf"))
	par(mfrow = c(1,2), cex = 1, las = 1, bty = 'l')
	matplot(gc.stats$gc.values, gc.stats$raw,
	type = 'b', col = 1, pch = c(1, 19, 1), lty = c(2, 1, 2),
	xlab = 'GC content (%)', ylab = 'Uncorrected depth ratio')
	legend('topright', legend = colnames(gc.stats$raw), pch = c(1, 19, 1))
	hist2(seqz.data$depth.ratio, seqz.data$adjusted.ratio,
	breaks = prettyLog, key = vkey, panel.first = abline(0, 1, lty = 2),
	xlab = 'Uncorrected depth ratio', ylab = 'GC-adjusted depth ratio')
	dev.off()
	# extract info from seqz
	seq.ex <- sequenza.extract(data.file,chromosome.list=c(1:22,'X','Y'))
	pdf(paste0(plotsDir,"/chromView.pdf"))
		for(i in 1:length(seq.ex$chromosomes))
			{
			chromosome.view(mut.tab = seq.ex$mutations[[i]], baf.windows = seq.ex$BAF[[i]],
			ratio.windows = seq.ex$ratio[[i]], min.N.ratio = 1,
			segments = seq.ex$segments[[i]], main = seq.ex$chromosomes[i])
			}
	dev.off()
	# inference of cellularity and ploidy
	CP.example <- sequenza.fit(seq.ex)
	sequenza.results(sequenza.extract = seq.ex, cp.table = CP.example,
	sample.id = sampleName, out.dir=seqOutDir)
	# confidence intervals
	cint <- get.ci(CP.example)
	pdf(paste0(plotsDir,"/confidenceIntervals.pdf"))
		cp.plot(CP.example)
		cp.plot.contours(CP.example, add = TRUE, likThresh = c(0.95))
		par(mfrow = c(2,2))
		cp.plot(CP.example)
		cp.plot.contours(CP.example, add = TRUE)
		plot(cint$values.cellularity, ylab = "Cellularity",
		xlab = "posterior probability", type = "n")
		select <- cint$confint.cellularity[1] <= cint$values.cellularity[,2] &
		cint$values.cellularity[,2] <= cint$confint.cellularity[2]
		#polygon(y = c(cint$confint.cellularity[1], cint$values.cellularity[select, 2], cint$confix = c(0, cint$values.cellularity[select, 1], 0), col='red', border=NA)
		lines(cint$values.cellularity)
		abline(h = cint$max.cellularity, lty = 2, lwd = 0.5)
		plot(cint$values.ploidy, xlab = "Ploidy",
		ylab = "posterior probability", type = "n")
		select <- cint$confint.ploidy[1] <= cint$values.ploidy[,1] &
		cint$values.ploidy[,1] <= cint$confint.ploidy[2]
		#polygon(x = c(cint$confint.ploidy[1], cint$values.ploidy[select, 1], cint$confint.ploidy[+ y = c(0, cint$values.ploidy[select, 2], 0), col='red', border=NA)
		lines(cint$values.ploidy)
		abline(v = cint$max.ploidy, lty = 2, lwd = 0.5)
	dev.off()
	# call CNVs and mutations
	cellularity <- cint$max.cellularity
	ploidy <- cint$max.ploidy
	avg.depth.ratio <- mean(seq.ex$gc$adj[, 2])
	# mutations
	mut.tab <- na.exclude(do.call(rbind, seq.ex$mutations))
	mut.alleles <- mufreq.bayes(mufreq = mut.tab$F,
	depth.ratio = mut.tab$adjusted.ratio,
	cellularity = cellularity, ploidy = ploidy,
	avg.depth.ratio = avg.depth.ratio)
	outTable = cbind(mut.tab[,c("chromosome","position","F","adjusted.ratio", "mutation")],
mut.alleles)
	save(list="outTable",file=paste0(infoDir,"/mutTable.Rdata"))
	# cnvs
	seg.tab <- na.exclude(do.call(rbind, seq.ex$segments))
	cn.alleles <- baf.bayes(Bf = seg.tab$Bf, depth.ratio = seg.tab$depth.ratio,
	cellularity = cellularity, ploidy = ploidy,
	avg.depth.ratio = avg.depth.ratio)
	seg.tab <- cbind(seg.tab, cn.alleles)
	save(list="seg.tab",file=paste0(infoDir,"/cnvTable.Rdata"))
	# visualize
	pdf(paste0(plotsDir,"/visualization.pdf"))
		chromosome.view(mut.tab = seq.ex$mutations[[3]], baf.windows = seq.ex$BAF[[3]],
		ratio.windows = seq.ex$ratio[[3]], min.N.ratio = 1,
segments = seg.tab[seg.tab$chromosome == seq.ex$chromosomes[3],],
		main = seq.ex$chromosomes[3],
		cellularity = cellularity, ploidy = ploidy,
		avg.depth.ratio = avg.depth.ratio)
		genome.view(seg.cn = seg.tab, info.type = "CNt")
		legend("bottomright", bty="n", c("Tumor copy number"),col = c("red"),
		inset = c(0, -0.4), pch=15, xpd = TRUE)
		genome.view(seg.cn = seg.tab, info.type = "AB")
		legend("bottomright", bty = "n", c("A-allele","B-allele"), col= c("red", "blue"),
		inset = c(0, -0.45), pch = 15, xpd = TRUE)
	dev.off()
	}
