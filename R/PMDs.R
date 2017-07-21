# ====================================================================
# Functions to define partially methylated domains - WIP
# ====================================================================

# split betas by chromosome
splitBetas = function(chmrs,betas=betasSarc,intergenic=TRUE,mani)
	{
	# probes that match chms
	index = which(mani$chr==chmrs)
	subManifest = mani[index,,drop=FALSE]
	# intergenic
	if(intergenic)
		{
		index = which(subManifest$UCSC_RefGene_Name=="")
		subManifest = subManifest[index,]
		}
	# probes in array
	index = which(rownames(subManifest)%in%rownames(betas))
	subManifest = subManifest[index,]
	# subset array
	subBetas = betas[rownames(subManifest),,drop=FALSE]
	# order
	subBetas = subBetas[order(subManifest$pos),,drop=FALSE]
	rownames(subBetas) = subManifest$pos[order(subManifest$pos)]
	return(subBetas)
	}

# get info for a window
getMeanSDn = function(start,end,pos,betas)
	{
	# which probes are in window
	index = which(pos>=start&pos<end)
	betas = betas[index]
	# mean, SD, n
	return(list(mean=mean(betas),sd=sd(betas),n=length(betas)))
	}

# moving window mean and SD
windowMeanSD = function(pos,betas,windowSize=50000)
	{
	# make windows
	pos = as.numeric(pos)
	starts = seq(from=min(pos),to=max(pos),by=windowSize)
	ends = starts+windowSize
	# get window mean etc
	means = mapply(FUN=getMeanSDn,start=starts,end=ends,MoreArgs=list(pos=pos,betas=betas))
	colnames(means) = starts
	return(means)
	}

# get CPs 
# previous defaults
# methodCP="PELT",pen="MBIC"
definePMDs = function(betas,changeFilter=c(-0.1),colour="blue",
	windowSize=30000,pen="MBIC",methodCP="PELT",penVal=0,
	stat="Normal",plotGen=FALSE,plotCP=FALSE,split=TRUE,
	plotLines=TRUE,plotPoints=FALSE,plotAll=TRUE,nLimit=1,
	method="changepoint",segLimit=0,retAll=FALSE)
	{
	library(changepoint)
	# split into p and q
	if(split)
		{
		centreIndex = which.max(diff(as.numeric(names(betas))))
		betas = list(betas[1:centreIndex],betas[(centreIndex+1):length(betas)])
		} else {
		betas = list(betas)
		}
	outAll = NULL
	for(i in 1:length(betas))
		{
		subBetas = betas[[i]]
		# get windowed means
		info = windowMeanSD(pos=names(subBetas),betas=subBetas,windowSize=windowSize)
		# remove NAs
		naIndex = !is.na(info["mean",])
		nIndex = info["n",]>=nLimit
		keepIndex = which(naIndex&nIndex)
		# get changepoints
		if(method=="changepoint")
			{
			# changepoints library
			CP = cpt.meanvar(as.vector(logit(as.numeric(info["mean",keepIndex]))),
				method=methodCP,penalty=pen,
				test.stat=stat,pen.value=penVal)
			CP@param.est$mean = as.vector(invlogit(CP@param.est$mean))
			# previous = MBIC and not logit
			# filters
			#meanB = mean(subBetas)
			indexStart = which(diff(CP@param.est$mean)<changeFilter)# reduction > 0.1 
			indexEnd = which(diff(CP@param.est$mean)>c(-changeFilter))
			# start and end points of CP
			firstTest = which(indexEnd==1)
			if(length(firstTest)>0)
			{
			  starts = c(CP@cpts[indexStart],CP@cpts[indexEnd[-firstTest]-1],1)
			  vals = c(CP@param.est$mean[indexStart],
			           CP@param.est$mean[indexEnd[-firstTest]-1],
			           CP@param.est$mean[1])
			  ends = c(CP@cpts[indexStart+1],CP@cpts[indexEnd[-firstTest]],CP@cpts[indexEnd[firstTest]])
			} else {
			  starts = c(CP@cpts[indexStart],CP@cpts[indexEnd-1])
			  vals = c(CP@param.est$mean[indexStart],CP@param.est$mean[indexEnd-1])
			  ends = c(CP@cpts[indexStart+1],CP@cpts[indexEnd])
			}
      DMPinfo = data.frame(starts,ends,vals)
      DMPinfo = unique(DMPinfo)
			# plot changepoints
			if(plotCP) 
				{
				plot(CP)
				abline(v=c(DMPinfo$starts,DMPinfo$ends),col="red",lty=2)
				}
			} else {
			# segment with DNAcopy
			seg = DNAcopy::segment(DNAcopy::CNA(logit(as.numeric(info["mean",keepIndex])),rep(1,times=length(keepIndex)),1:length(keepIndex)))
			seg$output$seg.mean = as.vector(invlogit(seg$output$seg.mean))
			# filter
			indexStart = which(diff(seg$output$seg.mean)<changeFilter)# reduction > 0.1 
			indexEnd = which(diff(seg$output$seg.mean)>c(-changeFilter))
			# start and end points of CP
			  starts = c(seg$output$loc.start[indexStart+1],
					seg$output$loc.start[indexEnd])
			  ends = c(seg$output$loc.end[indexStart+1],
				seg$output$loc.end[indexEnd])
			firstTest = which(indexEnd==1)
			if(length(firstTest)>0)
			{	
			vals = c(seg$output$seg.mean[indexStart],
			           seg$output$seg.mean[indexEnd[-firstTest]-1],
			           seg$output$seg.mean[1])
			} else {
			vals = c(seg$output$seg.mean[indexStart],seg$output$seg.mean[indexEnd-1])
			}
      			DMPinfo = data.frame(starts,ends,vals)
      			DMPinfo = unique(DMPinfo)
			
			# start and end points of CP
			#index = which(seg$output$seg.mean<segLimit)
			#starts = seg$output[index,"loc.start"]
			#ends = seg$output[index,"loc.end"]
			#vals = invlogit(seg$output[index,"seg.mean"])
			# plot changepoints
			if(plotCP) 
				{
				plot(seg)
				abline(v=c(starts,ends),col="red",lty=2)
				}
			}
		# genome scale start and end points
		DMPinfo$starts=as.numeric(colnames(info)[keepIndex][DMPinfo$starts])
		DMPinfo$ends=as.numeric(colnames(info)[keepIndex][DMPinfo$ends])
		# plot changepoints
		if(plotGen)
			{
			plot(NA,xlim=range(as.numeric(names(subBetas))),ylim=0:1)
			mapply(FUN=function(a,b){polygon(x=c(a,a,b,b),y=c(-2,2,2,-1),col="gray")},
				a=DMPinfo$starts,b=DMPinfo$ends)
			points(names(subBetas),subBetas)
			lines(names(subBetas),smooth(subBetas),col=colour)
			}
		# return
		outAll = rbind(outAll,DMPinfo)
		}
	# plot combined
	if(plotAll)
		{
		betas = do.call(c,betas)
		plot(NA,xlim=range(as.numeric(names(betas))),ylim=0:1,ylab="Beta",xlab="Position")
		mapply(FUN=function(a,b){polygon(x=c(a,a,b,b),y=c(-2,2,2,-1),col="gray")},
			a=outAll[,1],b=outAll[,2])
		if(plotPoints) points(names(betas),betas)
		if(plotLines) lines(names(betas),smooth(betas),col=colour)
		}
	return(outAll)
	}



