#!/usr/bin/env Rscript

args=(commandArgs(TRUE))


if(length(args) < 2 ){
  print("Missing input arguments :" )  
  print("args[1] Path to edgeR output")
  print("args[2] cutoff for significant genes")
  print("args[3] selection of cutoff [logFC|logCPM|LR|PValue|FDR] OPTIONAL " )
  print("Example start line :")
  print("./getGeneLists2TopGO.R malpighiantubules_vs_other_tissues.txt  0.01")
  
  print("")
}
fileName = args[1] #FileName of ath to gene annotation genes
cutoff = args[2] #Outputfrom [blast2GO/argot2]
Determinant = "FDR"
if(length(args) > 2){
  Determinant = args[3] #Genes that will be considered as background 
}



GeneList <- read.table(fileName, header=T,sep="\t", stringsAsFactors=F)
GOI = unique(rownames(GeneList[GeneList[Determinant] < cutoff, ]))
Background = unique(rownames(GeneList))

write.table(x=GOI,file=paste(fileName,cutoff,"txt", sep = "."),quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(x=Background,file=paste(fileName,"background","txt", sep = "."),quote=FALSE,row.names=FALSE,col.names=FALSE)

