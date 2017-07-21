# ==========================================================================
# pipeline for running minfi inc. normalisation etc
# ==========================================================================

minfiPipeline = function(dataDir,outDir=getwd(),nCells=5,legendloc="topleft",plotcolours=c("green","red","black","blue","cyan"),refactor=NULL,combat=NULL,toIgnore=NULL,badSampleCutoff=10.5)
	{
	# load in sample sheet
	sheet = read.metharray.sheet(dataDir)
	if(!is.null(toIgnore)) sheet = sheet[-which(sheet$Sample_Name%in%toIgnore),] 
	# read idats based on sample sheet
	RGset = read.metharray.exp(targets=sheet)
	# QC 
	phenoData = pData(RGset)
	Mset = preprocessRaw(RGset)
	pdf(paste0(outDir,"/QCplot.pdf"))
	plotQC(getQC(Mset),badSampleCutoff=1000)
	abline(badSampleCutoff *2 , -1, lty = 2)
	dev.off()
	pdf(paste0(outDir,"/densityPlot.pdf"))
	densityPlot(Mset,sampGroups=phenoData$Sample_Group)
	dev.off()
	pdf(paste0(outDir,"/densityBeanPlot.pdf"))
	densityBeanPlot(Mset,sampGroups=phenoData$Sample_Group)
	dev.off()
	pdf(paste0(outDir,"/controlProbes.pdf"))
	controlStripPlot(RGset)
	dev.off()
	# sex
	ratioSet = ratioConvert(Mset,what="both",keepCN=TRUE)
	Gset = mapToGenome(ratioSet)
	pdf(paste0(outDir,"/sexPlot.pdf"))
	plotSex(getSex(Gset,cutoff=-2))
	dev.off()
	rm(phenoData)
	rm(Mset)
	rm(ratioSet)
	rm(Gset)
	# detection p values
	qualCut = 0.01
	detP = detectionP(RGset)
	lowQual = detP>qualCut
	# fraction of failed positions per sample
	colMeans(lowQual) 
	probesToKeep = rownames(detP)[which(rowSums(detP>=qualCut)<=0)]
	# functional normalisation
	funnorm = preprocessFunnorm(RGset)
	# drop SNP probes
	funnorm = dropLociWithSnps(funnorm)
	save(list="funnorm",file=paste0(outDir,"/funnorm-droppedSnps.Rdata"))
	# get betas
	betas = getBeta(funnorm)
	annotation = getAnnotation(funnorm)
	# remove sex probes
	sexProbes = rownames(annotation)[which(annotation[,"chr"]%in%c("chrX","chrY"))]
	betas = betas[-which(rownames(betas)%in%sexProbes),]
	# keep high quality probes
	betas = betas[which(rownames(betas)%in%probesToKeep),]
	# save betas
	colnames(betas) = sheet$Sample_Name
	save(list=c("sheet","betas"),file=paste0(outDir,"/betas-funnorm-dropSex-dropLowQ.Rdata"))
	# refactor estimation
	if(!is.null(refactor))
		{
		# create refactor file
		tmp = cbind(rownames(betas),betas)
		colnames(tmp)[1] = "ID"
		write.table(tmp,file=paste0(outDir,"refactor_input.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE)
		rm(tmp)
		# run refactor to estimate cell proportions in each sample
		refactorOutput <- refactor(paste0(outDir,"/refactor_input.txt"), nCells, out = paste0(outDir,"/refactor_out.Rdata"))
		# make factors for combat
		covModel = model.matrix(~refactorOutput$refactor_components)
		}
	if(!is.null(combat))
		{
		if(is.null(covModel)) covModel = model.matrix(~1,data=sheet)
		# batch
		batch = as.factor(sheet[,combat])
		# run combat to normalise based on batch and cell proportions
		combatMs = ComBat(dat=logit(betas),batch=batch,mod=covModel)
		# save combat Ms
		save(list="combatMs",file=paste0(outDir,"/combatMs.Rdata"))
		forPCA = combatMs
		} else {
		forPCA = betas
		}
	# pca
	pca = runPCA(data=forPCA,fileOut=paste0(outDir,"/pca.Rdata"))
	plotPCA(pca,fileName="pca.pdf",outDir=outDir,groups=sheet$Pool_ID,arrows=FALSE,legendloc=legendloc,colours=plotcolours)
	# hclust
	groups = sapply(levels(as.factor(sheet$Pool_ID)),FUN=function(x) sheet$Sample_Name[which(sheet$Pool_ID==x)])
	hClust = runHclust(pca$x,groups=groups,file=paste0(outDir,"/hclust.pdf"),plotcolours=colours)
	save(list="hClust",file=paste0(outDir,"/hclust.Rdata"))
	}
