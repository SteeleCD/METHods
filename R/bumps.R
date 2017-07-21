# ===========================================================
# Functions for defining methylation bumps in density plots
# ===========================================================

# get bump limits using a sliding window and change of gradient
# this is the preferred method
getBumpLimsWindow = function(betas,windowSize=50,negCount=c(-20),gradLimit=10,sevLimit=c(-100),retLim=TRUE,combine=TRUE,plotDiag=FALSE,doAUC=FALSE)
{
if(plotDiag)
	{
	par(mfrow=c(3,1),mar=c(3,4,0,0))
	} else {
	par(mfrow=c(1,1),mar=c(3,4,0,0))
	}
if(doAUC) library(flux)
  # GET DENSITY
  dens = density(betas)
  # get gradients of density
  index1 = 2:length(dens$y)
  gradient = (dens$y[index1]-dens$y[index1-1])/(dens$x[index1]-dens$x[index1-1])
  # get change in gradients of density
  index2 = 2:length(gradient)
  d2 = (gradient[index2]-gradient[index2-1])/(dens$x[index2+1]-dens$x[index2])
  # new method - sliding window of change in gradient
  negIndex = sapply(1:(length(d2)-windowSize),FUN=function(x) sum(sign(d2[x:(x+windowSize)])))
  negDegree = sapply(1:(length(d2)-windowSize),FUN=function(x) sum(d2[x:(x+windowSize)]))
  # check that enough of window is negative
  test = which(negIndex<negCount)
  if(length(test)>0)
  {
    # runs of -ve change in gradient
    runs = rle(negIndex<negCount)
    runIndex = which(runs$values)
    ends = cumsum(runs$lengths)[runIndex]
    lengths = runs$lengths[runIndex]
    starts = ends-lengths
    # window limits
    tmp = cbind(starts,ends)
    severity = apply(tmp,MARGIN=1,FUN=function(x) sum(negDegree[x[1]:x[2]])/abs(diff(x)))
    # add window size to limits
    tmp[,2] = tmp[,2]+windowSize
    # remove those with low severity
    sevIndex = which(severity>sevLimit)
    if(length(sevIndex)>0)
	{
	tmp = tmp[-sevIndex,,drop=FALSE]
	severity = severity[-sevIndex]
	}
    # remove those that have too large a gradient
    largeGrad = apply(tmp,MARGIN=1,FUN=function(x) any(abs(gradient[x[1]:x[2]])>gradLimit))
    tmp = tmp[which(!largeGrad),,drop=FALSE]
    severity = severity[which(!largeGrad)]
    # remove those with too extreme limits
    lims = matrix(dens$x[tmp],ncol=2)
    extremeIndex = which(apply(lims,MARGIN=1,FUN=function(x) any(x<0.1|x>0.9)))
    if(length(extremeIndex)>0)
	{
    	tmp = tmp[-extremeIndex,,drop=FALSE]
    	lims = lims[-extremeIndex,,drop=FALSE]
	severity = severity[-extremeIndex]
	}
    if(nrow(tmp)==0) 
	{
	tmp = NA
	lims = NA	
	}
  } else {
    # no bump
    tmp = NA
    lims = NA
  }
  # combine if multiple bumps
if(!any(is.na(tmp)))
	{
  if(combine)
	{
	if(nrow(tmp)>1)
		{
		tmp = matrix(range(tmp),ncol=2)
		lims = matrix(range(lims),ncol=2)
		severity = sum(severity)
		}
	}
}
  # plot
  # plot density
  plot(dens,main=NA)
  abline(v=lims,lty=2)
  if(plotDiag)
	{
	# plot gradient
  	plot(dens$x[index1],gradient,col=as.factor(sign(gradient)))
  	abline(h=0,lty=2,col="gray")
  	abline(v=lims,lty=2)
  	# plot change in gradient
  	plot(dens$x[index2+1],d2,col=as.factor(sign(d2)))
  	abline(h=0,lty=2,col="gray")
  	abline(v=lims,lty=2)
	}
  # return bump limits
  if(retLim)
  {
  if(!any(is.na(tmp)))
  {
if(doAUC)
{
   AUC = sum(apply(tmp,MARGIN=1,FUN=function(i) auc(x=dens$x[i[1]:i[2]],y=dens$y[i[1]:i[2]])))
score = AUC*severity
} else {
AUC = list(x=dens$x[tmp[1,1]:tmp[1,2]],y=dens$y[tmp[1,1]:tmp[1,2]])
score = NA
}
   out = list(lims=lims,
	height=max(dens$y[min(tmp[,1]):max(tmp[,2])]),
	relHeight=max(dens$y[min(tmp[,1]):max(tmp[,2])])/max(dens$y),
	auc=AUC,
	severity=severity,
	score=score)
  } else {out=NA}
  } else {
    # TRUE/FALSE output
    if(!any(is.na(tmp)))
    {
      out = TRUE
    } else {out=FALSE}
  }
  return(out)
}


# get bump limits using change in sign of gradient
getBumpLims = function(betas)
{
  # GET DENSITY
  dens = density(betas)
  index = 2:length(dens$y)
  # get gradients of density
  gradient = (dens$y[index]-dens$y[index-1])/(dens$x[index]-dens$x[index-1])
  # signs of gradient
  signs = sign(gradient)
  # minimums of density
  mins = which(signs[2:length(signs)]-signs[1:(length(signs)-1)]==2)+1
  plot(dens)
  abline(v=dens$x[mins],lty=2)
  # check if have bump
  check = length(mins)>1
  if(check) 
  {
    # get maximums
    maxs = which(signs[2:length(signs)]-signs[1:(length(signs)-1)]==c(-2))+1
    maxValue = dens$y[maxs[2]]
    minValues = dens$x[mins]
    minValue = min(minValues)
    bumpHeight = maxValue-minValue
    relHeight = bumpHeight/max(dens$y)
    return(list(lims=dens$x[mins],height=bumpHeight,relHeight=relHeight))
  } else {
    return(NA)
  }
}
