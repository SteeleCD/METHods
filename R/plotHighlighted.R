# ==================================================================
# plot densities, highlighting certain samples
# ==================================================================

highlightInds = function(data,inds,colour="red",MAIN=NULL)
	{
	#load(file)
	index = which(colnames(data)%in%inds)
	densities = apply(data,MARGIN=2,density)
	rbinded = do.call(rbind,densities)
	plot(NA,xlim=range(unlist(rbinded[,"x"])),ylim=range(unlist(rbinded[,"y"])),ylab="Density",xlab="Normalised beta",main=MAIN)
	sapply((1:length(densities))[-index],FUN=function(x) lines(densities[[x]],col="gray"))
	sapply((1:length(densities))[index],FUN=function(x) lines(densities[[x]],col=colour,lwd=2))
	}
