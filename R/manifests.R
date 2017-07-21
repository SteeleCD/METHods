# ================================================================
# functions for working with methylation manifests
# ================================================================

# function to load manifests
getManifestOld = function(array="EPIC")
{
  if(array=="EPIC")
  {
  library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  manifest = cbind(IlluminaHumanMethylationEPICanno.ilm10b2.hg19@data$Locations, 
                   IlluminaHumanMethylationEPICanno.ilm10b2.hg19@data$Other)
  } else {
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)  
  manifest = cbind(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations, 
                   IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other)
  }
  return(manifest)
}


# function to load manifests
getManifest = function(array="EPIC")
{
  if(array=="EPIC")
  {
    library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
    data(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
    data(Locations)
    data(Other)
    manifest = cbind(Locations, 
                     Other)
  } else {
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)  
    data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    data(Locations)
    data(Other)
    manifest = cbind(Locations, 
                     Other)
  }
  return(manifest)
}
