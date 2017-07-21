# ================================================================
# Functions for plotting DMRs
# ================================================================

# function to get genes that overlap a region
getOverlapGenes = function(chr,start,end)
	{
	# get gRanges from this DMR
	ranges = gsub(" ","",paste0(chr,":",start,":",end))
	library(dplyr)
	gRanges = sapply(ranges, function (x) {res=strsplit(x, ':')}) %>%
		unlist %>%
		as.numeric %>%
		matrix(ncol=3, byrow=T) %>%
		as.data.frame %>%
		dplyr::select(chrom=V1, start=V2, end=V3) %>%
		mutate(chrom=paste0('chr', chrom)) %>%
		makeGRangesFromDataFrame
	# get Hsapiens genes
	library(Homo.sapiens)
	genesRanges = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
	# overlaps between genes and DMR
	overlap = subsetByOverlaps(genesRanges,gRanges); overlap
	geneStart = overlap@ranges@start
	geneEnd = overlap@ranges@start+overlap@ranges@width
	OVERLAP <<- overlap
	# get gene names
	library(org.Hs.eg.db)
	INFO <<- org.Hs.egSYMBOL
	print(str(org.Hs.egSYMBOL))
	unmapped = org.Hs.eg.db::org.Hs.egSYMBOL
	mapped = mappedkeys(unmapped)
	genes = unlist(as.list(unmapped[mapped]))
	GENES <<- genes
	DMRgenes = genes[which(names(genes)%in%overlap$gene_id)]
	DMRGENES <<- DMRgenes
	# return
	return(list(geneStart=geneStart,geneEnd=geneEnd,geneName=DMRgenes))
	}

# function to get beta values of all probes within a region, and any genes that overlap the region
getRegion = function(betas,chr,start,end,manifest,flank=10000)
	{
	region = sort(c(start,end))
	region = region[c(1,length(region))]
	region[1] = region[1]-flank 
	region[2] = region[2]+flank
	# down to chromosome
	manifest = manifest[which(manifest$CHR==chr),,drop=FALSE]
	# down to region
	index = which(manifest$MAPINFO>region[1]&manifest$MAPINFO<region[2])
	manifest = manifest[index,,drop=FALSE]
	# subset betas
	betas = betas[which(rownames(betas)%in%manifest$IlmnID),,drop=FALSE]
	# order
	betas = betas[order(rownames(betas)),,drop=FALSE]
	manifest = manifest[order(manifest$IlmnID),,drop=FALSE]
	# get info ready for output
	info = manifest$MAPINFO
	names(info) = manifest$IlmnID
	info = info[which(names(info)%in%rownames(betas))]
	# get overlapping genes
	overlapGenes = steeleLib:::getOverlapGenes(chr,region[1],region[2])
	# return
	return(list(betas=betas,pos=info,overlaps=overlapGenes,manRegion=manifest))
	}

fillColours = function(colour,n,doSep=FALSE)
	{
	if(doSep)
		{
		out = rainbow(n)
		} else {

		}
	}

# function to plot DMRs with overlapping genes displayed (if any)
plotDMR = function(betas,dmrs,index,manifest,flank=10000,
		groupIndices,doInvLogit=TRUE,xOffset=10,
		plotRatio=TRUE,plotAll=FALSE,allSepCols=FALSE,colours=NULL,
		customTitle=NULL,doLegend=FALSE,legloc="topleft",doSmooth=TRUE)
	{
	# colours
	if(is.null(colours)) colours = 1:length(groupIndices)
	# plotAll
	if(length(plotAll)==1) 
		{
		plotAll = rep(plotAll,length(groupIndices))		
		}
	if(length(allSepCols)==1)
		{
		allSepCols = rep(allSepCols,length(groupIndices))	
		}
	nGroups = sapply(groupIndices,length)
	allCols = vector(length=length(unlist(groupIndices)))
	for(j in 1:length(groupIndices)) allCols[groupIndices[[j]]] = colours[j]
	#allCols = unlist(sapply(1:length(groupIndices),FUN=function(x) rep(colours[x],times=length(groupIndices[[x]]))))
	if(any(allSepCols)) 
		{
		sepCols = rainbow(sum(nGroups[which(allSepCols)]))
		allCols[unlist(groupIndices[which(allSepCols)])] = sepCols
		}
	# get values in region
	region = steeleLib:::getRegion(betas,dmrs[index,"chr"],dmrs[index,"start"],dmrs[index,"end"],manifest,flank)
	# get positions
	pos1 = matrix(region$pos,ncol=length(groupIndices[[1]]),nrow=nrow(region$betas))
	pos2 = matrix(region$pos,ncol=length(groupIndices[[2]]),nrow=nrow(region$betas))
	# inverse logit for Ms
	if(doInvLogit) region$betas = invlogit(region$betas)
	# layout
	if(length(region$overlaps$geneStart)>0) 
		{
		layout(matrix(1:2,ncol=1,nrow=2),widths=1,heights=c(3,1))
		par(mar=c(0,4,4,2))
		XLAB = NA
		XAXT = 'n'
		} else {
		par(mar=c(5,4,4,2))
		layout(matrix(1,ncol=1,nrow=1),widths=1,heights=c(1))
		XLAB = "Position"
		XAXT = NULL
		}
	# plot DMR
	toOrder = order(pos1[,1])
	TITLE = ifelse(is.null(customTitle),paste0("Chromosome ",dmrs[index,"chr"]),customTitle)
	# plot betas
	positions = pos1[toOrder,1]
	if(plotRatio)
		{
		# plot ratio of betas in two groups		
		ratios = log(rowMeans(region$betas[,groupIndices[[1]],drop=FALSE])[toOrder]/rowMeans(region$betas[,groupIndices[[2]],drop=FALSE])[toOrder])
		if(doSmooth) ratios = smooth(ratios) 
		plot(positions,ratios,type="l",col="black",lwd=2,xlab=XLAB,ylab="log(ratio)",main=TITLE,xaxt=XAXT)
		abline(h=0,lty=2)
		} else {
		# plot separate betas
		plot(NA,xlab=XLAB,ylab="Beta",main=TITLE,xaxt=XAXT,ylim=c(0,1),xlim=range(positions))
		# add DMR limits
		mapply(FUN=function(a,b) polygon(x=c(a,b,b,a),y=c(9000,9000,-9000,-9000),col=rgb(0.5,0.5,0.5,0.5),lty=2),a=dmrs[index,"end"],b=dmrs[index,"start"])
		for(i in 1:length(groupIndices))
			{
			# plot group means
			#lines(positions,smooth(rowMeans(region$betas[,groupIndices[[i]],drop=FALSE])[toOrder]),lty=1,lwd=2,col=colours[i])
			if(doSmooth) 
				{
				plotBetas = apply(region$betas[toOrder,groupIndices[[i]],drop=FALSE],MARGIN=2,smooth)
				} else {
				plotBetas = region$betas[toOrder,groupIndices[[i]],drop=FALSE]
				}
			if(!is.matrix(plotBetas)) plotBetas = t(plotBetas)
			lines(positions,rowMeans(plotBetas),lty=1,lwd=2,col=colours[i])
			# plot individual betas
			if(plotAll[i])
				{
				#matplot(x=positions,y=plotBetas,lty=2,col=allCols[groupIndices[[i]]],add=TRUE,type="l")
				sapply(groupIndices[[i]],FUN=function(x) 
					{
					if(doSmooth) 
						{
						toPlot = smooth(region$betas[,x][toOrder])
						} else {
						toPlot = region$betas[,x][toOrder]
						}
					lines(positions,toPlot,lty=2,lwd=1,col=allCols[x])
					})
				}
			}
		}
	# legend
	if(doLegend)
		{
		# lines
		if(any(plotAll))
			{
			linesLeg = c(1,2)
			linesCols = c("gray","gray")
			linesNames = c("Mean","Separate")
			linesLWD = c(2,1)
			} else {
			linesLeg = NULL
			linesCols = NULL
			linesNames = NULL
			linesLWD = NULL
			}
		if(any(allSepCols))
			{
			indices = unlist(groupIndices[which(allSepCols)])
			linesLeg = c(linesLeg,rep(2,length(indices)))
			linesCols = c(linesCols,allCols[indices])
			linesNames = c(linesNames,colnames(betas)[indices])
			linesLWD = c(linesLWD,rep(1,length(indices)))
			}
		# colours
		colsLeg = rep(15,times=length(groupIndices))
		colsCols = colours
		colsNames = names(groupIndices)
		# legend
		legend(legloc,
			legend=c(colsNames,linesNames),
			col=c(colsCols,linesCols),
			lty=c(rep(NA,times=length(colsLeg)),linesLeg),
			lwd=c(rep(NA,times=length(colsLeg)),linesLWD),
			pch=c(colsLeg,rep(NA,times=length(linesLeg))))
		}
	# plot genes
	if(length(region$overlaps$geneStart)>0)
		{
		par(mar=c(5,4,0,2))
		yVals = 1:length(region$overlaps$geneStart)
		YLIM = range(yVals)
		YLIM[1] = YLIM[1]-1
		yVals = yVals-1
		if(length(xOffset)==1) xOffset = rep(xOffset,times=length(region$overlaps$geneStart)) 
		plot(NA,xlim=range(pos1),ylim=YLIM,yaxt="n",ylab=NA,xlab="Position")
		for(i in 1:length(yVals))
			{
			xRange = range(positions)
			polygon(c(region$overlaps$geneStart[i],region$overlaps$geneStart[i],region$overlaps$geneEnd[i],region$overlaps$geneEnd[i]),c(yVals[i]+0.35,yVals[i]+0.65,yVals[i]+0.65,yVals[i]+0.35),col="black")
			if(region$overlaps$geneStart[i]<xRange[1]&region$overlaps$geneEnd[i]>xRange[2])
				{
				xText = mean(xRange)
				yText = yVals[i]+0.75
				} else if(region$overlaps$geneStart[i]<xRange[1]) {
				xText = region$overlaps$geneEnd[i]+abs(diff(xRange))/xOffset[i]
				yText = yVals[i]+0.5
				} else if(region$overlaps$geneEnd[i]>xRange[2]) {
				xText = region$overlaps$geneStart[i]-abs(diff(xRange))/xOffset[i]
				yText = yVals[i]+0.5
				} else {
				xText = region$overlaps$geneEnd[i]+abs(diff(xRange))/xOffset[i]
				yText = yVals[i]+0.5
				}
			text(x=xText,y=yText,region$overlaps$geneName[i])
			}
		}
	}

# volcano plot
plotVolcano = function(data,manifestProbes=NULL,
		file=NULL,PAR=list(mfrow=c(1,1),
		mar=c(5,4,4,2)),title="test",
		folder=getwd(),
		MAIN="test")
	{
	base = paste0(folder,title,"/")
	if(!is.null(file)) pdf(paste0(base,title,file))
	if(!is.null(PAR)) do.call(par,PAR)
	index = which(data$adj.P.Val<0.05)
	Xvals = data$logFC[index]
	Yvals = -log(data$adj.P.Val[index])
	probes = paste0(data$probeID[index])
	XLIM = range(Xvals)
	YLIM = range(Yvals)
	if(!is.null(manifestProbes))
		{
		manifestProbes = paste0(manifestProbes)
		Xvals = Xvals[which(probes%in%manifestProbes)]
		Yvals = Yvals[which(probes%in%manifestProbes)]
		probes = probes[which(probes%in%manifestProbes)]
		}
	plot(Xvals,Yvals,col=ifelse(abs(Xvals)<0.3,"black","red"),ylab="-log(P)",xlab="log(fold change)",xlim=XLIM,ylim=YLIM,main=MAIN)
	if(!is.null(file)) dev.off()
	return(list(hyper=probes[which(Xvals>0.3)],hypo=probes[which(Xvals<c(-0.3))]))
	}
