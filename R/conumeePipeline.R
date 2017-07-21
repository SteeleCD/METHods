# ================================================================
# Conumee pipeline for copy number analysis from methylation array
# ================================================================

# function to get idat file names
getFiles = function(base,dataDir)
	{
	basenames = list.files(file.path(base,dataDir))
	basenames = basenames[grep(".idat",basenames)]
	basenames = gsub("_Grn.idat","",basenames)
	basenames = gsub("_Red.idat","",basenames)
	basenames = unique(basenames)
	basenames = file.path(base,dataDir,basenames)
	return(basenames)
	}

# conumee pipeline
conumeePipeline = function(baseDir,dataDirs,outDir)
	{
	# get file names
	basenames = unlist(sapply(dataDirs,FUN=function(x) getFiles(baseDir,x)))
	# read in data
	rgset = read.metharray(basenames,verbose=TRUE) 
	save(list="rgset",file=paste0(outDir,"/rgset.Rdata"))
	# data from rgset to mset
	mset = preprocessIllumina(rgset)
	rm(rgset)
	save(list="mset",file=paste0(outDir,"/mset.Rdata"))
	# controls
	library("CopyNumber450kData")
	data(RGcontrolSetEx)
	MsetControls = preprocessIllumina(RGcontrolSetEx)
	# annotation
	data(exclude_regions)
	data(detail_regions)
	# load	data
	data = CNV.load(mset)
	# only exclude HLA
	anno = CNV.create_anno(exclude_regions = exclude_regions[1], detail_regions = detail_regions)
	save(list="anno",file=paste0(outDir,"/anno.Rdata"))
	# load controls data
	rm(mset)
	controlData = CNV.load(MsetControls) 
	rm(MsetControls)
	# CNV
	sampleNames = names(data@intensity)
	for(i in 1:length(sampleNames))
		{
		x = CNV.fit(data[sampleNames[i]],controlData,anno)
		x = CNV.bin(x)
		x = CNV.detail(x)
		x = CNV.segment(x)
		pdf(paste0(outDir,"/",sampleNames[i],"-plot.pdf"))
		CNV.genomeplot(x)
		dev.off()
		save(list="x",file=paste0(outDir,"/",sampleNames[i],"-res.Rdata"))
		}
 	}

