
#source("http://bioconductor.org/biocLite.R") 
#biocLite("topGO")
#install.packages("memisc") 
#biocLite("Rgraphviz")



getGOtermGeneRelationship  <- function(ModelGeneGOrelationshipFile,ModelGeneGOrelationshipFileFormat = "blast2GO"){
  
  if(ModelGeneGOrelationshipFileFormat == "blast2GO"){
    # SeqName Hit-Desc        GO-Group        GO-ID   Term
    # HMEL025033-PA   homeobox protein abdominal-a homolog    F       GO:0003700      sequence-specific DNA binding transcription factor activity
    ModelGeneGOrelationship <- read.table(ModelGeneGOrelationshipFile, header=T, fill=T, sep="\t", quote="", stringsAsFactors=F)
    colnames(ModelGeneGOrelationship) <- c("Locus", "Name", "Function", "GOterm", "Description")
    
  }else if(ModelGeneGOrelationshipFileFormat == "argot2"){
    # Sequence        Aspect  GO ID   Name    Total Score     Internal Confidence     Information Content
    # HMEL015695-PA   F       GO:0004222      metalloendopeptidase activity   6.728758473824388       0.5     7.3129653461121995
    
    ModelGeneGOrelationship <- read.table(ModelGeneGOrelationshipFile, header=T, fill=T, sep="\t", quote="", stringsAsFactors=F)
    colnames(ModelGeneGOrelationship) <- c("Locus", "Function", "GOterm", "Name", "Total Score", "Internal Confidence","Information Content")
  }else{
    ##Asumes that there are three colums with 1.GeneName 2.Function  3. Goterm and it containts a header
    
    # Locus        Function  GOterm  
    # HMEL015695-PA   F       GO:0004222
    ModelGeneGOrelationship <-  read.table(ModelGeneGOrelationshipFile, header=T, fill=T, sep="\t", quote="", stringsAsFactors=F)
    colnames(ModelGeneGOrelationship) <- c("Locus", "Function", "GOterm")
    
  }
  return(ModelGeneGOrelationship)
}  


getAssociatedGenesEdgeR <- function(edgeRFile, cutoff){
  
  GeneList <- read.table(edgeRFile, header=T,sep="\t", stringsAsFactors=F)
  GOIList = unique(rownames(GeneList[GeneList[Determinant] < cutoff, ]))
  BackgroundList = unique(rownames(GeneList))
  print("")
  print(paste("Defining set of locuses that will serve as background from this file",edgeRFile, sep=" : ") )
  print(paste("Number of unique locuses ",length(BackgroundList), sep=" : ") )
  print("")
  print("")
  print(paste("Defining set of locuses that will serve as significant locuses with",Determinant,"less than",cutoff,  sep=" ") )
  print(paste("Number of unique locuses ",length(GOIList), sep=" : ") )
  print("")
  print("")
  
}
  
