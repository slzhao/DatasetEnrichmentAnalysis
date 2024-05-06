dataRankTest=function(selectedGenes,referenceData,universGenes=NULL,nRep=1000,returnExpRank=FALSE,missingValue=0) {
  referenceDataInPackage<-data(package="DatasetEnrichmentAnalysis")$results[,"Item"]
  referenceDataInPackage=gsub("^referenceDataRank \\(","",referenceDataInPackage)
  referenceDataInPackage=gsub("\\)$","",referenceDataInPackage)

  if (is.character(referenceData)) { #referenceDataInPackage
    if (referenceData %in% referenceDataInPackage) {
      data(list=referenceData,envir=environment()) #get referenceDataRank
      #referenceData<-get(referenceDataName,envir=environment())
    } else {
      stop(paste0("The referenceData name ",referenceData," is not in the dataset of DatasetEnrichmentAnalysis package"))
    }
  } else {
    if (!is(referenceDataName, "matrix") & (!is(referenceData, "data.frame"))) {
      stop(paste0("The referenceData is not a matrix or data.frame"))
    }
    referenceDataRank=referenceData
  }

  if (is.null(universGenes)) {
    universGenes=row.names(referenceDataRank)
  }
  selectedGenes=intersect(selectedGenes,universGenes)
  expdata=referenceDataRank[intersect(row.names(referenceDataRank),universGenes),]

  print(paste0("Number of selected genes: ",length(selectedGenes)))
  print(paste0("Number of univers/experiment genes: ",nrow(expdata)))

  #expdata_norm<-abs(expdata)/(apply(expdata,1,function(x){sqrt(sum(x^2))}))
  #expdata_rank<-apply(expdata_norm,2,rank)
  expdata_rank=expdata

  topn_rank<-expdata_rank[rownames(expdata_rank) %in% selectedGenes,]
  topn_rank_sum<-apply(topn_rank,2,sum)

  if (returnExpRank) {
    return(list(expdata_rank,topn_rank))
  }

  perm_rank<-NULL
  for (i in 1:nRep) {
    perm_rank<-rbind(perm_rank,apply(expdata_rank[sample(1:nrow(expdata_rank),nrow(topn_rank)),],2,sum))
  }

  pval<-NULL

  for (i  in 1:ncol(topn_rank)) {
    pval<-c(pval,sum(perm_rank[,i]>topn_rank_sum[i]))
  }

  pval<-pval/nRep
  names(pval)<-colnames(topn_rank)

  #return(pval)

  expdataCount=apply(expdata,2,function(x) length(which(x!=missingValue)))
  selectedGenesCount=apply(expdata[selectedGenes,],2,function(x) length(which(x!=missingValue)))

  return(data.frame(Total=expdataCount,Selected=selectedGenesCount,p=pval,pAdj=p.adjust(pval,method="BH")))
}


# dataToRank=function(referenceData,universGenes=NULL) {
#   if (is.null(universGenes)) {
#     universGenes=row.names(referenceData)
#   }
#   referenceDataNormlized=abs(referenceData)/(apply(referenceData,1,function(x){sqrt(sum(x^2))}))
#   referenceDataRank=apply(referenceDataNormlized,2,rank)
#   return(referenceDataRank)
# }
