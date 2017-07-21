# ================================================================
# Functions for annotating methylation DMRs
# ================================================================


# load up transcripts etc
preDMR = function() 
	{
	# get gene info
	library(Homo.sapiens)
	genesRanges = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
	# get gene names
	library(org.Hs.eg.db)
	unmapped = org.Hs.eg.db::org.Hs.egSYMBOL
	mapped = mappedkeys(unmapped)
	genes = unlist(as.list(unmapped[mapped]))
	# get transcripts
	library(bumphunter)
	transcripts = annotateTranscripts(txdb=TxDb.Hsapiens.UCSC.hg19.knownGene,by="gene")
	return(list(genesRanges=genesRanges,
			genes=genes,
			transcripts=transcripts))
	}

# cross reference DMRs with gene info
DMRinfo = function(dmrFile,# DMR file from bumphunter
		fwerThresh=0.1,
		pThresh = 0.05,
		doP = TRUE,
		outDir,
		doGO=TRUE,
		doKEGG=FALSE) 
	{
	preDMRinfo = preDMR()
	# set threshold
	if(doP)
		{
		thresh = pThresh
		threshType = "p.value"
		} else {
		thresh = fwerThresh
		threshType = "fwer"
		}
	# load DMRs
	fileEnding = rev(strsplit(dmrFile,"[.]")[[1]])[1]
	if(grepl("R",fileEnding))
		{
		load(dmrFile)
		dmr = dmr$table
		} else {
		dmr = read.table(dmrFile,sep="\t",as.is=TRUE)
		}
	# what chromosome encoding
	chromBool = all(grepl("chr",dmr$chr))
	# make bed files for homer
	forBed = dmr
	if(!chromBool) 
		{
		forBed[,"chr"] = paste0("chr",forBed[,"chr"])
		} else {
		dmr$chr = gsub("chr","",dmr$chr)
		}
	rtracklayer::export.bed(makeGRangesFromDataFrame(forBed[which(forBed[,threshType]<thresh),]),paste0(outDir,"/allSig.bed"))
	indexDownBed = which(forBed[,threshType]<thresh&forBed[,"value"]<0)
	indexUpBed = which(forBed[,threshType]<thresh&forBed[,"value"]>0)
	if(length(indexDownBed)>0)
		{
		rtracklayer::export.bed(makeGRangesFromDataFrame(forBed[indexDownBed,]),paste0(outDir,"/downSig.bed"))
		}
	if(length(indexUpBed)>0)
		{
		rtracklayer::export.bed(makeGRangesFromDataFrame(forBed[indexUpBed,]),paste0(outDir,"/upSig.bed"))
		}
	# volcano plot
	zeroBool = dmr[,threshType]==0
  if(any(zeroBool))
  {
	dmr[which(zeroBool),threshType] = min(dmr[-which(zeroBool),threshType])
  }
	pdf(paste0(outDir,"/volcano.pdf"))
	colours = ifelse(dmr[,threshType]<thresh,"green4","black")
	colours[which(zeroBool)] = "red"
	PCH = ifelse(zeroBool,2,1)
	plot(dmr[,"value"],-log(dmr[,threshType]),col=colours,xlab="Beta difference",ylab=paste0("-log(,",threshType,")"),pch=PCH)
	dev.off()
	# overlapping genes all
	sigIndex = which(dmr[,threshType]<thresh)
	infoAll = paste0("chr",dmr[sigIndex,"chr"],":",dmr[sigIndex,"start"],"-",dmr[sigIndex,"end"])
	gRangesAll = as(infoAll,"GRanges")
	# overlaps between genes and DMR
	overlapAll = subsetByOverlaps(preDMRinfo$genesRanges,gRangesAll)
	DMRgenesAll = preDMRinfo$genes[which(names(preDMRinfo$genes)%in%overlapAll$gene_id)]
	write.table(unique(DMRgenesAll),file=paste0(outDir,"/allGenes.txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
	# overlapping genes upreg
	indexUp = which(dmr[,"value"]>0&dmr[,threshType]<thresh)
	if(length(indexUp)>0)
		{
		infoUp = paste0("chr",dmr[indexUp,"chr"],":",dmr[indexUp,"start"],"-",dmr[indexUp,"end"])
		gRangesUp = as(infoUp,"GRanges")
		# overlaps between genes and DMR
		overlapUp = subsetByOverlaps(preDMRinfo$genesRanges,gRangesUp)
		DMRgenesUp = preDMRinfo$genes[which(names(preDMRinfo$genes)%in%overlapUp$gene_id)]
		write.table(unique(DMRgenesUp),file=paste0(outDir,"/upGenes.txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
		}	
	# overlapping genes downreg
	indexDown = which(dmr[,"value"]<0&dmr[,threshType]<thresh)
	if(length(indexDown)>0)
		{	
		infoDown = paste0("chr",dmr[indexDown,"chr"],":",dmr[indexDown,"start"],"-",dmr[indexDown,"end"])
		gRangesDown = as(infoDown,"GRanges")
		# overlaps between genes and DMR
		overlapDown = subsetByOverlaps(preDMRinfo$genesRanges,gRangesDown)
		DMRgenesDown = preDMRinfo$genes[which(names(preDMRinfo$genes)%in%overlapDown$gene_id)]
		write.table(unique(DMRgenesDown),file=paste0(outDir,"/downGenes.txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
		}	
	# feature DMRs e.g. promoter DMRs
	index = which(dmr[,threshType]<thresh)
	info = paste0("chr",dmr[index,"chr"],":",dmr[index,"start"],"-",dmr[index,"end"])
	colours = colours[index]
	PCH = PCH[index]
	gRanges = as(info,"GRanges")
	matched = matchGenes(x=gRanges,subject=preDMRinfo$transcripts)
	sigDMRs = dmr[index,]	
	# subset to different features
	# 3' UTR
	threePrimeIndex = which(matched$UTR%in%c("3'UTR","overlaps 3'UTR"))
	if(length(threePrimeIndex)>0) {
		pdf(paste0(outDir,"/3prime-volcano.pdf"))
		plot(sigDMRs[threePrimeIndex,"value"],-log(sigDMRs[threePrimeIndex,threshType]),pch=PCH[threePrimeIndex],col=colours[threePrimeIndex],xlab="Beta difference",ylab=paste0("-log(",threshType,")"))
		dev.off()}
	# 5' UTR
	fivePrimeIndex = which(matched$UTR%in%c("5'UTR","overlaps 5'UTR"))
	if(length(fivePrimeIndex)>0) {
		pdf(paste0(outDir,"/5prime-volcano.pdf"))
		plot(sigDMRs[fivePrimeIndex,"value"],-log(sigDMRs[fivePrimeIndex,threshType]),pch=PCH[fivePrimeIndex],col=colours[fivePrimeIndex],xlab="Beta difference",ylab="-log(",threshType,")")
		dev.off()}
	# promoter
	promoterIndex = which(matched$region=="promoter")
	write.table(matched$name[promoterIndex],file=paste0(outDir,"/promAllGenes.txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
	upIndex = which(sigDMRs[,"value"]>0)
	write.table(matched$name[upIndex[which(upIndex%in%promoterIndex)]],file=paste0(outDir,"/promUpGenes.txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
	downIndex = which(sigDMRs[,"value"]<0)
	write.table(matched$name[downIndex[which(downIndex%in%promoterIndex)]],file=paste0(outDir,"/promDownGenes.txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
	if(length(promoterIndex)>0) 
		{
		pdf(paste0(outDir,"/prom-volcano.pdf"))
		plot(sigDMRs[promoterIndex,"value"],-log(sigDMRs[promoterIndex,threshType]),pch=PCH[promoterIndex],col=colours[promoterIndex],xlab="Beta difference",ylab=paste0("-log(",threshType,")"))
		dev.off()
		}
	# upstream
	upstreamIndex = which(matched$region=="upstream")
	if(length(upstreamIndex)>0) {
		pdf(paste0(outDir,"/upstream-volcano.pdf"))
		plot(sigDMRs[upstreamIndex,"value"],-log(sigDMRs[upstreamIndex,threshType]),pch=PCH[upstreamIndex],col=colours[upstreamIndex],xlab="Beta difference",ylab=paste0("-log(",threshType,")"))
		dev.off()}
	# downstream
	downstreamIndex = which(matched$region=="downstream")
	if(length(downstreamIndex)>0) {
		pdf(paste0(outDir,"/downstream-volcano.pdf"))
		plot(sigDMRs[downstreamIndex,"value"],-log(sigDMRs[downstreamIndex,threshType]),pch=PCH[downstreamIndex],col=colours[downstreamIndex],xlab="Beta difference",ylab=paste0("-log(",threshType,")"))
		dev.off()}
	# intronic
	intronicIndex = which(matched$description=="inside intron")
	if(length(intronicIndex)>0) {
		pdf(paste0(outDir,"/intronic-volcano.pdf"))
		plot(sigDMRs[intronicIndex,"value"],-log(sigDMRs[intronicIndex,threshType]),pch=PCH[intronicIndex],col=colours[intronicIndex],xlab="Beta difference",ylab=paste0("-log(",threshType,")"))
		dev.off()}
	# exonic
	exonicIndex = which(matched$description%in%c("inside exon","covers exon(s)","overlaps exon upstream","overlaps exon downstream","overlaps two exons","covers"))
	if(length(exonicIndex)>0) {
		pdf(paste0(outDir,"/exonic-volcano.pdf"))
		plot(sigDMRs[exonicIndex,"value"],-log(sigDMRs[exonicIndex,threshType]),pch=PCH[exonicIndex],col=colours[exonicIndex],xlab="Beta difference",ylab=paste0("-log(",threshType,")"))
		dev.off()}
	write.csv(cbind(sigDMRs,matched),paste0(outDir,"/annotatedDMRs.csv"))
	# GO enrichment
	if(doGO)
	{
	system(paste0("mkdir ",outDir,"/GO"))
	system(paste0("mkdir ",outDir,"/GO/up"))
	system(paste0("mkdir ",outDir,"/GO/down"))
	GOenrichmentMethy(geneFile=paste0(outDir,"/promUpGenes.txt"),
		outDir=paste0(outDir,"/GO/up"),
		array="EPIC")
	GOenrichmentMethy(geneFile=paste0(outDir,"/promDownGenes.txt"),
		outDir=paste0(outDir,"/GO/down"),
		array="EPIC")		
	}
	# KEGG enrichment
	if(doKEGG)
	{
	enrichKEGG(dataDir=outDir,
	           upFile="promUpGenes.txt",
	           downFile="promDownGenes.txt")
	}
	return(cbind(sigDMRs,matched))
	}

# cross reference gene list with COSMIC

# cross reference gene list with epigenetic genes 
