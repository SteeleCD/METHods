library(minfi)

runMinfi = function(baseDir,saveIntermediates=TRUE,outDir=getwd(),normMethod="funnorm",nCores=8)
	{
	# read sheet
	sheet = read.metharray.sheet(baseDir)	
	# read idats based on sample sheet
	RGset = read.metharray.exp(targets=sheet)
	rm(sheet)
	if(saveintermediates)
		{
		save(list="RGset",file=paste0(outDir,"/RGset.Rdata"))
		}
	# get methylset
	Mset = preprocessRaw(RGset)
	if(saveintermediates)
		{
		save(list="Mset",file=paste0(outDir,"/Mset.Rdata"))
		}
	rm(Mset)
	# map to genome
	GRset = mapToGenome(RGset)
	save(list="GRset",file=paste0(outDir,"/GRset.Rdata"))
	# probe locations as GenomicRanges
	if(saveIntermediates)
		{
		gr = granges(GRset)
		save(list="gr",file=paste0(outDir,"/granges.Rdata"))
		rm(gr)
		}
	# QC plot
	qcReport(RGset,pdf=paste0(outDir,"/QC.pdf"))
	# choose normalisation	
	# not implemented yet
	if(normMethod=="funnorm")
		{
		norm = preprocessFunnorm(RGset)	
		} else {
		norm = GRset	
		}
	if(saveIntermediates)
		{
		save(list="norm",file=paste0(outDir,"/normGR-",normMethod,".Rdata"))	
		}
	# convert normalised to M
	Mvals = getM(norm)
	if(saveIntermediates)
		{
		save(list="Mvals",file=paste0(outDir,"/normM-",normMethod,".Rdata"))	
		}
	# remove SNPs
	GRset <- dropLociWithSnps(GRset, snps=c("SBE","CpG"), maf=0)
	save(list="GRset",file=paste0(,"/GRset-droppedSNPs.Rdata"))
	# run DMP
	beta <- getBeta(norm)
	CT  <- pData(norm)$Sample_Group
	dmp <- dmpFinder(beta, pheno = CT  , type = "categorical")
	save(list="dmp",file=paste0(outDir,"/DMP.Rdata"))
	# run DMR
	runDMR(norm,outDir)
	}
	
runDMR = function(norm,outDir=getwd(),nCores=8)
	{
	# DMR
	library(doParallel)
	registerDoParallel(cores = 8)
	pheno <- pData(norm)$Sample_Group
	designMatrix <- model.matrix(~ pheno)
	# run with low cutoff
	cutoff = 0.2
	dmrs <- bumphunter(norm, design = designMatrix, 
             cutoff = cutoff, B=0, type="Beta")
	# if many find a new cutoff that is more stringent
	while(nrow(dmrs)>30000)
		{
		cutoff = cutoff+0.05
		dmrs <- bumphunter(norm, design = designMatrix, 
            cutoff = cutoff, B=0, type="Beta")
		}
	# run again with new cutoff, but with many permutations
	dmrs <- bumphunter(norm, design = designMatrix, 
             cutoff = cutoff, B=1000, type="Beta")
	save(list="dmrs",file=paste0(outDir,"/DMR.Rdata"))	
	# sva
	library(sva)
	mval = getM(GRset)
	pheno = pData(GRset)
	mod = model.matrix(~as.factor(Sample_Group),data=pheno)
	mod0 = model.matrix(~1,data=pheno)
	sva.results = sva(mval,mod,mod0)
	save(list="sva.results",file=paste0(outDir,"/SVA.Rdata"))
	}

