#!/usr/bin/env Rscript

#source("http://bioconductor.org/biocLite.R") 
#biocLite("topGO")
#install.packages("memisc") 
#biocLite("Rgraphviz")



args=(commandArgs(TRUE))

if(length(args) < 2){
  print("Missing input arguments :" )  
  print("args[1] Full path to gene-Goterm relationship file")
  print("args[2] Path to edgeR output")
  print("args[3] GOtermFile that contains go terms that should be extracted")
  print("args[4] cutoff")
  
  
  
  print("Example line:")
  print("./extractGenes.R Hmel1.21.cdhit.blast2go.GOannotation.txt Fat_body_core_vs_extended.txt  GO0003735_with_children.txt 0.05 FDR")
  
  print("")
}
ModelGeneGOrelationshipFile = args[1]          #Full path to GOterm annotated genes
expressionFile = args[2]                            #Full path to candidate genes File
GOtermsFile= args[3]
cutoff = args[4]
cutoffParameter = args[5]

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

getGenesWithAssociatedGOTerm <- function(ModelGeneGOrelationshipFile,expressionFile,GOtermsFile, cutoff=0.05 , cutoffParameter = "FDR"){
  ## read in files
  GOtermGeneRelationship = getGOtermGeneRelationship(ModelGeneGOrelationshipFile)
  expression = read.table(expressionFile,header=TRUE,sep = "\t", stringsAsFactors=F)
  GOterms = read.table(GOtermsFile,header=FALSE,sep = "\t")
  
  
  
  #Filtering of Genes
  GOgenes = GOtermGeneRelationship[ GOtermGeneRelationship$GOterm %in%  GOterms$V1,1]
  #Filtering in expression table
  expression_GOresults =expression[rownames(expression) %in% GOgenes,]
  dim(expression_GOresults)
  # Filter to remove all above cutoff
  expression_GOresults_underCutoff = expression_GOresults[expression_GOresults[cutoffParameter]< cutoff, ] 
  return(expression_GOresults_underCutoff)
}



expressionFile = "Asterids_vs_Ranunculales.txt"
GOtermsFile = "GO0003735_with_children.txt"

# Get genes associated with goterms in goterms file
expression_GOresults_underCutoff = getGenesWithAssociatedGOTerm(ModelGeneGOrelationshipFile,expressionFile,GOtermsFile, cutoff=cutoff , cutoffParameter = cutoffParameter)

# print to file 
expression_Goresults_File = paste(expressionFile,GOtermsFile ,cutoff, "genes.txt", sep = "_")
write.table(expression_GOresults_underCutoff,expression_Goresults_File,quote = FALSE, row.names = TRUE, col.names =TRUE)  


expressionFile = "Asterids_vs_Ranunculales.txt"
GOtermsFile = "GO0003735_with_children.txt"

# print to file 
expression_Goresults_File = paste(expressionFile,GOtermsFile ,cutoff, "genes.txt", sep = "_")
write.table(expression_GOresults_underCutoff,expression_Goresults_File,quote = FALSE, row.names = TRUE, col.names =TRUE)  

