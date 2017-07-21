getMarkers = function(dir,manifestFile,outFile)
	{
	data = read.csv(paste0(dir,"/",manifestFile))
	data = data[,c("IlmnID","CHR","MAPINFO")]
	colnames(data) = c("Marker Name","Chromosome","Marker Position")
	index = which(data[,1]==""|is.na(data[,1]))
	if(length(index)!=0) data = data[-index,]
	index = which(data[,2]==""|is.na(data[,2]))
	if(length(index)!=0) data = data[-index,]
	index = which(data[,3]==""|is.na(data[,3]))
	if(length(index)!=0) data = data[-index,]
	write.table(data,file=paste0(dir,"/",outFile),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
	}


#getMarkers(dir=getwd(),manifestFile="HumanMethylation450_15017482_v1-2.csv",outFile="450markers.txt")
