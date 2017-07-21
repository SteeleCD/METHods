# ==================================================================
# plotting copy number segments
# ==================================================================

check = function(toCheck)
	{
	isNull = is.null(toCheck)
	if(isNull) return(TRUE)
	isEmpty = length(toCheck)==0
	if(isEmpty) return(TRUE)
	isNA = is.na(toCheck)
	if(isNA) return(TRUE)
	return(FALSE)
	}

# combined segments with low numbers of markers
combineSegs = function(chromData,nMin=15)
	{
	toSkip = NULL
	index = which(chromData[,5]<nMin)
	if(length(index)==0) return(chromData)
	for(i in 1:length(index))
		{
		distanceFromStart = abs(chromData[index[i],3]-chromData[index[i]-1,4])
		print(distanceFromStart)
		DS <<- distanceFromStart
		distanceFromEnd = abs(chromData[index[i],4]-chromData[index[i]+1,3])
		print(distanceFromEnd)
		DE <<- distanceFromEnd
		closest = which.min(c(distanceFromStart,distanceFromEnd))
		if(check(distanceFromStart)) closest=2
		if(check(distanceFromEnd)) closest=1
		if(check(distanceFromStart)&check(distanceFromEnd)) 
			{
			toSkip = c(toSkip,i)
			next
			}
		print(closest)
		CLOSEST  <<- closest
		if(closest==1) 
			{
			nprobes1 = chromData[index[i]-1,5]
			nprobes2 = chromData[index[i],5]
			logR1 = chromData[index[i]-1,6]
			logR2 = chromData[index[i],6]
			} else if(closest==2) {
			nprobes1 = chromData[index[i]+1,5]
			nprobes2 = chromData[index[i],5]
			logR1 = chromData[index[i]+1,6]
			logR2 = chromData[index[i],6]
			}
			chromData[index[i]-1,6]=(nprobes1*logR1+nprobes2*logR2)/(nprobes1+nprobes2)
			chromData[index[i]-1,5]=nprobes1+nprobes2
		}
	if(!is.null(toSkip)) index = index[-toSkip]
	chromData[-index,]
	}

# plot a single chromosome
plotFun = function(seg.mean,num.mark=NULL,seg.start,seg.end,
		YLIM,colours,highlightStarts,highlightEnds,fusionPos,chrom)
	{
	#scatterplot3d(z=seg.mean,y=num.mark,x=seg.position)
	plot(NA,ylim=range(seg.mean),xlim=range(c(seg.start,seg.end)),col=colours[num.mark],xaxt="n",main=chrom)
	abline(h=0)
	if(!is.null(num.mark))
		{
		sapply(1:length(seg.mean),FUN=function(x) lines(x=c(seg.start[x],seg.end[x]),y=rep(seg.mean[x],2),col=colours[num.mark][x],lwd=3))
		} else {
		sapply(1:length(seg.mean),FUN=function(x) lines(x=c(seg.start[x],seg.end[x]),y=rep(seg.mean[x],2),col=colours,lwd=3))
		}
	# higlight regions
	if(length(highlightStarts)>0)
		{
		for(i in 1:length(highlightStarts))
			{
			polygon(x=c(highlightStarts[i],
					highlightEnds[i],
					highlightEnds[i],
					highlightStarts[i]),
				y=c(YLIM[2]+1,YLIM[2]+1,YLIM[1]-1,YLIM[1]-1),
				col=rgb(0.1,0.1,0.1,0.1))
			}
		}
	# plot lines at fusions
	if(length(fusionPos)>0)
		{
		abline(v=fusionPos,lty=2)
		}
	}

# plot all chromosomes
plotByChrom = function(segFile=NULL,segObj=NULL,dataDir,
		outDir=NULL,fileName=NULL,combineSegs=FALSE,
		sampleCol=2,chromCol=1,startCol=3,endCol=4,nMarkCol=NULL,segMeanCol=6,
		plotAber=FALSE,highlightChroms=NULL,highlightStarts=NULL,highlightEnds=NULL,
		fusionChroms=NULL,fusionPos=NULL)
	{
	library(copynumber)
	if(!is.null(segFile)) data = read.table(paste0(dataDir,"/",segFile),head=TRUE,sep="\t")
	if(!is.null(segObj)) data = segObj
	# plot by chms
	chms = unique(data[,chromCol])
	individual = unique(data[,sampleCol])
	# combine small segs
	if(combineSegs)
		{	
		data = do.call(rbind,sapply(individual,FUN=function(x) do.call(rbind,sapply(chms,FUN=function(y) combineSegs(data[which(data[,sampleCol]==x&data[,chromCol]==y),]),simplify=FALSE)),simplify=FALSE))
		} else {
		data=seg
		}
	# colours
	if(!is.null(nMarkCol)) 
		{
		colours = colorRampPalette(c("red","green"))(max(data[,nMarkCol]))
		} else {
		colours="black"
		}
	library(scatterplot3d)
	# plot
	if(!is.null(fileName)) pdf(paste0(outDir,'/',fileName))
	sapply(individual,FUN=function(x) {
		par(mfrow=c(4,6),mar=c(1,2,2,0))
		sapply(chms,FUN=function(y) {
			print(paste0(x,":",y))
			# any highlighted regions in this chms?
			if(!is.null(highlightChroms))
				{
				highlightIndex = which(highlightChroms==y)
				} else {
				highlightIndex = c()
				} 
			# any fusions in this chms?
			if(!is.null(fusionChroms))
				{
				fusionIndex = which(fusionChroms==y)
				} else {
				fusionIndex = c()
				}
			index = which(data[,sampleCol]==x&data[,chromCol]==y)
			if(!is.null(nMarkCol))
				{
				# coloured by number of markers
				plotFun(data[index,segMeanCol],
					data[index,nMarkCol],
					data[index,startCol],
					data[index,endCol],
					range(data[,segMeanCol]),
					colours,
					highlightStarts=highlightStarts[highlightIndex],
					highlightEnds=highlightEnds[highlightIndex],
					fusionPos=fusionPos[fusionIndex],
					chrom=y)
				} else {
				# coloured black
				plotFun(seg.mean=data[index,segMeanCol],
					seg.start=data[index,startCol],
					seg.end=data[index,endCol],
					YLIM=range(data[,segMeanCol]),
					colours=colours,
					highlightStarts=highlightStarts[highlightIndex],
					highlightEnds=highlightEnds[highlightIndex],
					fusionPos=fusionPos[fusionIndex],
					chrom=y)
				}
			})
		})
	if(!is.null(fileName)) dev.off()
	# abberation plot
	if(plotAber)
		{
		data = cbind(data[,c(sampleCol,chromCol)],rep("q",times=nrow(data)),data[,c(startCol,endCol,nMarkCol,segMeanCol)])
		colnames(data) = c("sampleID","chrom","arm","start.pos","end.pos","n.probes","logR.mean")
		if(!is.null(fileName)) pdf(paste0(outDir,'/copynumber.pdf'))
	 	plotAberration(data,thres.gain=0.2)
		if(!is.null(fileName)) dev.off()
		}
	}
	
# wrapper for plotting CN heatmap
plotAbber = function(segFile=NULL,segObj=NULL,dataDir,
		outDir=NULL,fileName="abberationPlot.pdf",
		HEAD=TRUE,logCN=TRUE,thresh=0.2,combine=TRUE,
		chrom=NULL)
	{
	library(copynumber)
	if(!is.null(segFile)) data = read.table(paste0(dataDir,"/",segFile),head=HEAD,sep="\t")
	if(!is.null(segObj)) data = segObj
	if(!logCN) data[,6] = log2(data[,6])-1
	# plot by chms
	chms = unique(data[,2])
	individual = unique(data[,1])
	# combine small segs
	if(combine)
		{
		data = do.call(rbind,sapply(individual,FUN=function(x) do.call(rbind,sapply(chms,FUN=function(y) combineSegs(data[which(data[,1]==x&data[,2]==y),]),simplify=FALSE)),simplify=FALSE))	
		}
	data = cbind(data[,1:2],rep("q",times=nrow(data)),data[,3:ncol(data)])
	colnames(data) = c("sampleID","chrom","arm","start.pos","end.pos","n.probes","logR.mean")
	if(!is.null(outDir)) pdf(paste0(outDir,'/',fileName))
	plotAberration(data,thres.gain=thresh,chrom=chrom)
	if(!is.null(outDir)) dev.off()
	}

