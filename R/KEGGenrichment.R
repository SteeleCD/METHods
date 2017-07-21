# ================================================================
# KEGG enrichment from methylation data
# ================================================================


enrichKEGG = function(dataDir,upFile="promUpGenes.txt",downFile="promDownGenes.txt",dmrFile="annotatedDMRs.csv",array="EPIC",type="geneList")
	{
	library(clusterProfiler)
	library(ReactomePA)
	# read in manifest
	manifest = steeleLib:::getManifestOld(array)
	allGenes = sapply(manifest$UCSC_RefGene_Name,FUN=function(x) unique(strsplit(x,split=";")[[1]]))
	allGenes = table(unlist(allGenes))
	# convert all genes to entrez id
	entrezAll = bitr(names(allGenes), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
	# KEGG ENRICHMENT
	# organise folders
	outDir = paste0(dataDir,"KEGG/")
	system(paste0("mkdir ",outDir))
	# read in DM genes
	if(type=="geneList")
		{
		DMGu = read.table(paste0(dataDir,"/",upFile),as.is=TRUE,col.names="gene")[,1]
		DMGd = read.table(paste0(dataDir,"/",downFile),as.is=TRUE,col.names="gene")[,1]
		} else {
		dmrs = read.csv(paste0(dataDir,"/",dmrFile),as.is=TRUE)
		DMGu = unique(dmrs$name[which(dmrs$value>0)])
		DMGd = unique(dmrs$name[which(dmrs$value<0)])
		}
	# convert to entrez
	if(length(DMGu)>0)
		{
		entrezU = bitr(DMGu, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
		} else {
		entrezU = data.frame(SYMBOL=character(0),ENTREZID=integer(0))
		}
	if(length(DMGd)>0)
		{
		entrezD = bitr(DMGd, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
		} else {
		entrezD = data.frame(SYMBOL=character(0),ENTREZID=integer(0))
		}
	# compare reactome clusters
	reactClusts <- compareCluster(geneCluster = list(up=entrezU[,"ENTREZID"],down=entrezD[,"ENTREZID"]), fun = "enrichPathway")
	# plot ckr
	pdf(paste0(outDir,"reactome-dotplot.pdf"))
	dotplot(reactClusts)
	dev.off()
	# save results
	save(list="reactClusts",file=paste0(outDir,"reactome-results.Rdata"))
	return(reactClusts)
	}

