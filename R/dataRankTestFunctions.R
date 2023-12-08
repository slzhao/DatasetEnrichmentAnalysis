dataRankTest=function(selectedGenes,referenceData,universGenes=NULL,nRep=1000,returnExpRank=FALSE,missingValue=0) {
  if (is.null(universGenes)) {
    universGenes=row.names(referenceData)
  }
  selectedGenes=intersect(selectedGenes,universGenes)
  expdata=referenceData[intersect(row.names(referenceData),universGenes),]

  print(paste0("Number of selected genes: ",length(selectedGenes)))
  print(paste0("Number of univers/experiment genes: ",nrow(expdata)))

  expdata_norm<-abs(expdata)/(apply(expdata,1,function(x){sqrt(sum(x^2))}))
  expdata_rank<-apply(expdata_norm,2,rank)

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
