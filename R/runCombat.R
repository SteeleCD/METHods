# ==========================================================================
# ComBAT normalisation
# ==========================================================================


# logit function
logit = function(p,offset=1e-9) 
	{
  p = as.matrix(p)
	maxIndex = which(p==1)
	minIndex = which(p==0)
	if(length(maxIndex)>0) p[maxIndex] = 1-offset
	if(length(minIndex)>0) p[minIndex] = 0+offset
	log(p/(1-p))
	}
	
# inverse logit function
invlogit = function(x) 
	{
  x = as.matrix(x)
	maxIndex = which(is.infinite(x)&x>0)
	minIndex = which(is.infinite(x)&x<0)
	inv = exp(x)/(1+exp(x))
	if(length(maxIndex)>0) inv[maxIndex] = 1
	if(length(minIndex)>0) inv[minIndex] = 0
	return(inv)
	}
	
# combat without controls
combatNoControl = function(normFile,dataFile,removeSlides=NULL,saveIntermediates=FALSE,outDir=getwd())
	{
	# without controls
	load(normFile)
	load(dataFile)
	pd = data$pd
	rm(data)
	betaNorm = norm$beta
	if(!is.null(removeIndex))
		{
		index = which(pd$Slide%in%removeSlides)
		pd = pd[-index,]
		betaNorm = betaNorm[,-index]			
		}
	if(saveIntermediates)
		{
		save(list="betaNorm",file=paste0(outDir,"/normForCombat.Rdata"))
		save(list="pd",file=paste0(outDir,"pdForCombat.Rdata"))
		}
	library(sva)
	batch = pd$Slide
	modcombat = model.matrix(~1, data=pd)
	combat_edata = ComBat(dat=logit(betaNorm), batch=batch, mod=modcombat, par.prior=TRUE)
	save(list="combat_edata",file=paste0(outDir,"/combat.Rdata"))				
	}


# with controls
combatControls = function(normFile,dataFile,outDir=getwd())
	{
	library(ChAMP)
	load(normFile)
	load(datafile)
	combat = champ.runCombat(beta.c=norm$beta,pd=data$pd)	
	save(list="combat",file=paste0(outDir,"/combat.Rdata"))
	}
