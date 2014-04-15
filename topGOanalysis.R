#!/usr/bin/env Rscript

#source("http://bioconductor.org/biocLite.R") 
#biocLite("topGO")

library(topGO)
library(memisc)
library(Rgraphviz)

args=(commandArgs(TRUE))

if(length(args) < 5 ){
  print("Missing input arguments :" )  
  print("args[1] Full path to gene-Goterm relationship file")
  print("args[2] Output type of gene-Goterm relationship file [blast2GO/argot2/other]")
  print("args[3] Type of function. This is dependent on the naming in the gene-Goterm relationship file but in most cases it is one of [F/P/C]")
  print("args[4] Full path to candidate genes, e.g DE expressed genes.")
  print("args[5] Prefix for results for filename of GO analysis")
  print("args[6] Full path to background genes [Optional if not entire gene set should be used]")
  print("args[7] Ortholog gene relationship [Optional if go terms are associated with another species genes]")
  print("")
  print("Example start line :")
  
  print("")
}
ModelGeneGOrelationshipFile = args[1] #Full path to gene annotation genes
ModelGeneGOrelationshipFileFormat = args[2] #Outputfrom [blast2GO/argot2]
GOtermFunction = args[3] #Output type of gene-Goterm relationship file [blast2GO/argot2/other]
significantGenesFile = args[4] #Full path to candidate genes File
outName = args[5] # Prefix for results for filename of GO analysis
if(length(args) > 5){
  backgroundGenesFile = args[6] #Genes that will be considered as background 
}
if(length(args) > 6){
  ModelGeneOrthologsFile = args[7] #Genes that will be considered as background 
}

getTopGOdata <- function(ModelGeneGOrelationshipFile,ModelGeneGOrelationshipFileFormat,GOtermFunction,significantGenesFile,
                         backgroundGenesFile = "Not present",ModelGeneOrthologsFile="Not present"){
  #############################################################
  # Load the ortholog information  START 
  ##############################################################
  #Data structure specific for genomes found in phytozome.
  if(ModelGeneOrthologsFile != "Not present"){
    ModelGeneOrthologs <- read.table(ModelGeneOrthologsFile, header=F, fill=T, sep="\t", quote="", stringsAsFactors=F)
    #Add names to make it functional for this script 
    names(ModelGeneOrthologs) <-c("Locus","Ortholog")
  }else{
    ModelGeneOrthologs <- "Not present"
  }
  
  #############################################################
  # Load the ortholog information  STOP
  ##############################################################
  
  #############################################################
  # Load the background list and Gene of interest information START
  ##############################################################
  
  #If named "Not present" then all Locuses will be considered 
  if(backgroundGenesFile != "Not present"){
    backgroundList <- read.table(backgroundGenesFile, header=F)[[1]]
  }else{
    BackgroundList <- "Not present"
  }
  
  #Genes of interest, should be a subset of BackgroundList  
  GOIList <- read.table(significantGenesFile, header=F)[[1]]
  
  ############################################################
  # Load the background list and Gene of interest information STOP
  ##############################################################
  
  #############################################################
  # Load the Gene GOTerm relationship table START
  ##############################################################
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
    ModelGeneGOrelationship <- read.argot2file(ModelGeneGOrelationshipFile)
    colnames(ModelGeneGOrelationship) <- c("Locus", "Function", "GOterm")
    
  }
  
  #############################################################
  # Load the Gene GOTerm relationship table  STOP
  ##############################################################
  
  #############################################################
  # Reformat the ortholog information  START
  ##############################################################
  
  if(ModelGeneOrthologs != "Not present"){
    ModelGeneOrthologsRightFormat = getRightOrthoInfo(ModelGeneOrthologs,mRNA2Gene = TRUE)
    printOrthoInfo(ModelGeneOrthologsRightFormat)
  }else{
    ModelGeneOrthologsRightFormat = data.frame(Locus=unique(ModelGeneGOrelationship$Locus),Ortholog=unique(ModelGeneGOrelationship$Locus))
  }
  
  #############################################################
  # Reformat the ortholog information  STOP
  ##############################################################
  
  
  #############################################################
  # Filter for background list of genes and
  # and mark genes of interest (GOI) START
  ##############################################################
  
  if(BackgroundList == "Not present"){
    BackgroundList = ModelGeneOrthologsRightFormat$Locus
  }
  ModelGeneOrthologsRightFormat = filterGeneList(ModelGeneOrthologsRightFormat,BackgroundList,GOIList)
  #printOrthoInfoAfterFiltering(ModelGeneOrthologsRightFormat)
  
  # Run function which filters the files based on BackgroundList and marks based on  GOIList
  
  #############################################################
  # Filter for background list of genes and
  # and mark genes of interest (GOI) STOP
  ##############################################################
  
  
  #############################################################
  # Change the GOterm table into right format
  # and filter for only those present in List START
  ##############################################################
  
  ModelGeneGOrelationship <- getGOtermFormat(ModelGeneGOrelationship,Function=GOtermFunction,LocusList =ModelGeneOrthologsRightFormat$Ortholog)
  printGOtableInfo(ModelGeneGOrelationship)
  
  #############################################################
  # Change the GOterm table into right format
  # and filter for only those present in List START
  ##############################################################
  
  
  #############################################################
  # Move the format into the way TopGO wants it. 
  # and load TopGO START
  ##############################################################
  
  geneList = getTopGOGenelistFormat(ModelGeneGOrelationship$Locus,ModelGeneOrthologsRightFormat$Ortholog[ModelGeneOrthologsRightFormat$GOI==1])
  gene2GOlist =getTopGOGene2GOrelationshipFormat(ModelGeneGOrelationship)
  
  topGOdata <- new("topGOdata", description=outName, ontology="BP", allGenes=geneList, annot = annFUN.gene2GO, gene2GO = gene2GOlist, nodeSize=0)
  
  
  #############################################################
  # Move the format into the way TopGO wants it. 
  # and load TopGO STOP
  ##############################################################
  
  return (topGOdata)
}

#statistic tests
printTopGOresults <- function(topGOdata,outName){
  test.statFis <- new("classicCount", testStatistic=GOFisherTest, name="Fisher test")
  resultFis <- getSigGroups(topGOdata, test.statFis)
  resultFisadj <- resultFis
  Fisadj <- p.adjust(score(resultFis),method="BH")
  score(resultFisadj) <- Fisadj
  
  test.statKS <- new("classicScore", testStatistic=GOKSTest, name="KS tests")
  resultKS <- getSigGroups(topGOdata, test.statKS)
  resultKSadj <- resultKS
  KSadj <- p.adjust(score(resultKS),method="BH")
  score(resultKSadj) <- KSadj
  
  test.statElim <- new("elimCount", testStatistic=GOFisherTest, name="Fisher test", cutOff=0.01)
  resultElim <- getSigGroups(topGOdata, test.statElim)
  
  test.statWei <- new("weightCount", testStatistic=GOFisherTest, name="Fisher test", sigRatio="ratio")
  resultWei <- getSigGroups(topGOdata, test.statWei)
  
  list <- list(classic=score(resultFis), KS=score(resultKS), elim=score(resultElim), weight=score(resultWei))
  
  listadj <- list(classicadj=score(resultFisadj), KSadj=score(resultKSadj), elim=score(resultElim), weight=score(resultWei))
  
  allRes <- GenTable(topGOdata, classicFisher = resultFis, classicKS = resultKS, elimFisher = resultElim, WeightedFisher = resultWei, orderBy="WeightedFisher", ranksOf="classicFisher", topNodes=30)
  
  allResadj <- GenTable(topGOdata, classicFisheradj = resultFisadj, classicKSadj = resultKSadj, elimFisher = resultElim, WeightedFisher = resultWei, orderBy="WeightedFisher", ranksOf="classicFisheradj", topNodes=30)
  
  #output of data
  #saves everything in Results order within the folder one is at
  if (!file.exists("results")){
    dir.create(file.path("results"))
  }
  
  write.table(allResadj, file=paste("results/",outName,"_stat_adj.txt",sep=""), row.names=FALSE, col.names=TRUE)
  
  
  pdf(paste("results/P_value_distribution_",outName,"_adj.pdf",sep=""))
  par(mfrow = c(2,2))
  for (nn in names(listadj)){
    p.val <- listadj[[nn]]
    hist(p.val[p.val < 1], br = 50, xlab = "p values", 
         main = paste("Histogram for method:", nn, sep = " "))
  }
  dev.off()
  
  pdf(paste("results/GO_graphFis_",outName,"_15nodes_adj.pdf",sep=""))
  showSigOfNodes(topGOdata, score(resultFisadj), firstSigNodes=15, useInfo="all")
  dev.off()
  
  write.table(allRes, file=paste("results/",outName,"_stat.txt",sep=""), row.names=FALSE, col.names=TRUE)
  
  pdf(paste("results/P_value_distribution_",outName,".pdf",sep=""))
  par(mfrow = c(2,2))
  for (nn in names(list)){
    p.val <- list[[nn]]
    hist(p.val[p.val < 1], br = 50, xlab = "p values", 
         main = paste("Histogram for method:", nn, sep = " "))
  }
  dev.off()
  
  pdf(paste("results/GO_graphFis_",outName,"_15nodes.pdf",sep=""))
  showSigOfNodes(topGOdata, score(resultFis), firstSigNodes=15, useInfo="all")
  dev.off()
  
  #pdf(paste("results/GO_graphFis_TEST_NEW_3_15nodes.pdf",sep=""))
  #showSigOfNodes(topGOdata, score(resultFis), firstSigNodes=15, useInfo="all")
  #dev.off()
  
  pdf(paste("results/GO_graphWeiFis_",outName,"_15nodes.pdf",sep=""))
  showSigOfNodes(topGOdata, score(resultWei), firstSigNodes=15, useInfo="all")
  dev.off()
  
}

getGOtermFormat <- function(ModelGeneGOrelationship, Function="P",LocusList=ModelGeneGOrelationship$Locus,GOtermList=ModelGeneGOrelationship$GOterm){
  
  
  if(checkFormat(ModelGeneGOrelationship,c("Locus","GOterm","Function"))){
    cat("GOterm table format is correct. \n")
    #Keep only columns that are interesting for this purpose
    ModelGeneGOrelationship = ModelGeneGOrelationship[,c("Locus","GOterm","Function")]
    
    #Filter based on function 
    ModelGeneGOrelationship = ModelGeneGOrelationship[which(ModelGeneGOrelationship$Function == Function),]
    
    #Filter based on Locus 
    ModelGeneGOrelationship = ModelGeneGOrelationship[ModelGeneGOrelationship$Locus %in% LocusList,]
    
    #Filter based on GOterm 
    ModelGeneGOrelationship = ModelGeneGOrelationship[ModelGeneGOrelationship$GOterm %in% GOtermList,]
    return (ModelGeneGOrelationship)
  }else{
    stop("Something wrong with the GOterm table column names")
  }
  return (NULL)
  
}







printOrthoInfo <- function(ModelGeneOrthologsRightFormat){
  cat("This is the informtion about the orthologs\n")
  print(head(ModelGeneOrthologsRightFormat))
  cat("Nr of unique Locuses:\n")
  cat(length(unique(ModelGeneOrthologsRightFormat$Locus)))
  cat("\nNr of unique orthologs:\n")
  cat(length(unique(ModelGeneOrthologsRightFormat$Ortholog)))
  cat("\n\n")
}



printOrthoInfoAdvanced <- function(ModelGeneOrthologsRightFormat){
  cat("This is the informtion about the orthologs\n")
  print(head(ModelGeneOrthologsRightFormat))
  cat("Nr of unique Locuses:\n")
  cat(length(unique(ModelGeneOrthologsRightFormat$Locus)))
  cat("\nNr of unique orthologs:\n")
  cat(length(ModelGeneOrthologsRightFormat$GOI==1))
  cat("\n\n")
}


printGOtableInfo <- function(ModelGeneOrthologsRightFormat){
  cat("This is the informtion about the GOtable\n")
  print(head(ModelGeneOrthologsRightFormat))
  cat("\nNr of unique Functions:\n")
  cat(length(unique(ModelGeneOrthologsRightFormat$Function)))
  cat("\nNr of unique Locuses:\n")
  cat(length(unique(ModelGeneOrthologsRightFormat$Locus)))
  cat("\nNr of unique GOterms:\n")
  cat(length(unique(ModelGeneOrthologsRightFormat$GOterm)))
  cat("\n\n")
}



checkFormat <- function(dat,Mandatory){
  found = Mandatory %in% colnames(dat)
  return (all(found))
}


getRightOrthoInfo <- function(ModelGeneOrthologs,mRNA2Gene=FALSE){
  
  if(checkFormat(ModelGeneOrthologs,c("Locus","Ortholog"))){
    cat("Ortholog format is correct. Getting right Ortho info\n")
    
    #Keep only interesting columns
    ModelGeneOrthologs <- ModelGeneOrthologs[ , c("Locus","Ortholog")]
    
    #Remove rows where there is no ortholog
    ModelGeneOrthologs = ModelGeneOrthologs[ModelGeneOrthologs$Ortholog>0 , ]
    
    if(mRNA2Gene){
      ##Going from isoform to gene not always applicable
      ModelGeneOrthologs$Ortholog =getLocusFromGeneName(ModelGeneOrthologs$Ortholog)
    }
    
    return (ModelGeneOrthologs)
  }else{
    stop("Something wrong with the Ortholog column names")
  }
  return (NULL)
  
}

getLocusFromGeneName <-function(GeneName){
  #This assumes geneName to be AT1GXXXXX.1
  Raw = unlist(strsplit(GeneName,"[.]") )
  Locus = Raw[seq(1,length(Raw)-1,by=2)]
  return (Locus)
}





filterGeneList <- function(GenesWithOrthologs, SubsetList=GenesWithOrthologs$Locus, GOIList){
  #GenesWithOrthologs should be a two row data.frame  with names "Locus" "Ortholog"
  
  #First remove all GeneNames that is not in the SubsetList
  FilteredGenesWithOrthologs <- GenesWithOrthologs[GenesWithOrthologs$Locus %in% SubsetList,]
  #Add column if Gene with ortholog is interesting
  FilteredGenesWithOrthologs$GOI  <- as.integer(FilteredGenesWithOrthologs$Locus %in% GOIList) 
  return(FilteredGenesWithOrthologs)
  
}

getTopGOGenelistFormat <- function(Genes,GenesOfInterest){
  #First step
  Genes = unique(Genes)
  GenesGOI <- factor(as.integer(Genes %in% GenesOfInterest))
  names(GenesGOI) <- Genes
  return (GenesGOI)
}

getTopGOGene2GOrelationshipFormat <- function(ModelGeneGOrelationship){
  gene2GOlist_sb <- split(ModelGeneGOrelationship$GOterm, ModelGeneGOrelationship$Locus)
  gene2GOlist_sb <- lapply(gene2GOlist_sb, unique) #remove duplicates
  return (gene2GOlist_sb)
}





topGOdata = getTopGOdata(ModelGeneGOrelationshipFile,ModelGeneGOrelationshipFileFormat,GOtermFunction,significantGenesFile)
printTopGOresults(topGOdata,outName)

# set workddir for the RNAs 




#nodeSize: prunes GO hierarchy for less than X genes for one GO term

#statistic tests
#Fis & KS adjusted, others not adjusted
