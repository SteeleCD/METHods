# ================================================================
# utility functions for getting GISTIC inputs from various tools
# ================================================================


# from champ output to GISTIC seg file
champToGistic=function(dir)
	{
	files = list.files(dir)
	files = files[grep(".txt",files)]
	out = NULL
	for(i in 1:length(files))
		{
		out = rbind(out,read.table(paste0(dir,"/",files[i]),head=TRUE))
		}
	out = out[,-ncol(out)]
	colnames(out) = c('Sample','Chromosome','Start Position','End Position','Num markers','Seg.CN')
	write.table(out,file=paste0(dir,'/seg.txt'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
	}

# from sequenza to GISIC seg file
sequenzaToGistic = function(segFiles,sampleNames=NULL,outDir,outFile=NULL)
	{
	if(is.null(sampleNames)) sampleNames = paste0("S",1:length(segFiles))
	out = NULL
	for(i in 1:length(segFiles))
		{
		data = read.table(segFiles[i],head=TRUE,sep="\t",as.is=TRUE)
		data = data[,c("chromosome","start.pos","end.pos","N.BAF","CNt")]
		data$CNt[which(data$CNt==0)] = 1e-9
		data[,"CNt"] = log2(data[,"CNt"])-1
		data = cbind(rep(sampleNames[i],times=nrow(data)),data)
		colnames(data) = c("Sample","Chromosome","Start Position","End Position","Num markers","Seg.CN")
		out = rbind(out,data)
		}
	if(!is.null(outFile))
		{
		write.table(out,paste0(outDir,outFile),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
		} else {
		return(out)
		}
	}

# from conumee to GISTIC seg file
conumeeToGISTIC = function(dir,ignore=NULL,outFile="seg.txt")
	{
	files = list.files(dir)
	files = files[grep("-res",files)]
	if(!is.null(ignore))
		{
		index = unique(unlist(sapply(ignore,FUN=function(x) grep(x,files))))
		if(length(index)>0) files = files[-index]
		}
	out = NULL
	for(i in 1:length(files))
		{
		load(paste0(dir,"/",files[i]))
		tmp = x@seg$summary
		tmp[,"seg.median"] = tmp[,"seg.median"]-x@bin$shift
		tmp[,1] = strsplit(files[i],split="-")[[1]][1]
		tmp[,2] = gsub("chr","",tmp[,2])
		tmp = tmp[,c(1:5,8)]
		out = rbind(out,tmp)
		}
	colnames(out) = c("Sample","Chromosome","Start Position","End Position","Num markers","Seg.CN")
	write.table(out,file=paste0(dir,"/",outFile),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
	}

# create pseudo markers file from seg file
segToMarkers=function(file,dir,manifestMarkersFile = "/media/grp11/CI_Pathology_Steele/methylation/manifests/EPICmarkers.txt")
	{
	columnNames = c("Marker Name", "Chromosome", "Marker Position")
	# read segment file 
	data = read.table(paste0(dir,'/',file),head=TRUE,sep='\t',as.is=TRUE)
	tmp1 = data[,c(2,3)]
	tmp2 = data[,c(2,4)]
	colnames(tmp2) = colnames(tmp1)
	out = rbind(tmp1,tmp2)
	out = out[order(out[,1],out[,2]),]
	out = out[!duplicated(out),]
	out = cbind(paste0('marker',1:nrow(out)),out)
	colnames(out) = columnNames
	# read manifest files
	manifest = read.table(manifestMarkersFile,sep="\t",head=TRUE,as.is=TRUE)
	colnames(manifest) = columnNames
	#dup = apply(out,MARGIN=1,FUN=function(x) any(manifest[,2]==x[2]&manifest[,3]==x[3]))
	#if(sum(dup)>0)
	#	{
	#	out = out[-which(dup),]
	#	}
	out = rbind(out,manifest)
	write.table(out,file=paste0(dir,'/markers.txt'),col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')
	}
