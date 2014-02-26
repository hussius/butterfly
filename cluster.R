library(pheatmap)
library(pvclust)

counts <- read.delim("read_counts_with_proper_headers.csv",sep="\t",row.names=1)

meta  <- read.csv("rna-seqecoevowabi_relational_table.csv")
meta <- meta[-which(meta$Customer_ID=="B1MT"),]

labelled_heatmap <- function(data, meta){
	corrs <- cor(data)
	display.names <- paste(meta$Host_Plant,meta$Tissue,meta$Host.Plant.use,sep=".")
	names(display.names) <- meta$Customer_ID
	colnames(corrs) <- display.names[colnames(corrs)]
	rownames(corrs) <- colnames(corrs)
	pheatmap(corrs)
}

labelled_heatmap(counts,meta)

# Normalize
source("/Users/mikaelhuss1/Desktop/allt/devel/RNASeqTools/RNASeqTools.R")

tmm <- cpm.tmm(counts)
labelled_heatmap(tmm,meta)

log.cpm.tmm <- normalize.voom(tmm)
labelled_heatmap(log.cpm.tmm,meta)

# Colors
tissue <- meta$Tissue
names(tissue) <- meta$Customer_ID
hostplant <- meta$Host_Plant
names(hostplant) <- meta$Customer_ID
pgroup <- meta$Phylogeny_group
names(pgroup) <- meta$Customer_ID
type <- meta$Host.Plant.use
names(type) <- meta$Customer_ID

# ToDo: some PCAs to try to identify PC:s with various factors
p <- prcomp(t(log.cpm.tmm))

comp1 <- 4
comp2 <- 5

par(mfrow=c(2,2))
# Color by tissue
plot(p$x[,c(comp1,comp2)], col=as.numeric(tissue[colnames(counts)]),pch=20,main=paste0("PCs ", comp1, ",", comp2))
legend("topright", legend=c("Fat body","Gut","Labial gland","Malpighian tubules"),pch=20,col=1:4)
# Color by host plant
plot(p$x[,c(comp1,comp2)], col=as.numeric(hostplant[colnames(counts)]),pch=20,main=paste0("PCs ", comp1, ",", comp2))
#legend("topright", legend=unique(hostplant),pch=20,col=1:length(unique(hostplant)))
# Double-check coloring
#textxy(p$x[,1],p$x[,2],labs=colnames(log.cpm.tmm))
# Color by phylogeny group
plot(p$x[,c(comp1,comp2)], col=as.numeric(pgroup[colnames(counts)]),pch=20,main=paste0("PCs ", comp1, ",", comp2))
legend("topright", legend=unique(pgroup),pch=20,col=1:length(unique(pgroup)))
# Color by core vs extended
plot(p$x[,c(comp1,comp2)], col=as.numeric(type[colnames(counts)]),pch=20,main=paste0("PCs ", comp1, ",", comp2))
legend("topright", legend=unique(type),pch=20,col=1:length(unique(type)))

# Loop over combinations of PCs and color by family (etc) 
par(mfrow=c(4,4))
for (comp1 in 1:4){
	for (comp2 in 1:4){
		if (comp1 != comp2){
	# Color by family
	family = meta$Family
	names(family) <- meta$Customer_ID
	plot(p$x[,c(comp1,comp2)], 	col=as.numeric(family[colnames(counts)]),pch=20,main=paste0("PCs ", comp1, ",", comp2))
#legend("topright", 	legend=unique(family),pch=20,col=1:length(unique(family)))
	}
}
}

# Try to filter out everything below CPM 1 
less.than.1 <- which(rowMeans(tmm)<1)
x <- tmm[-less.than.1,]

p <- prcomp(t(x))

comp1 <- 2
comp2 <- 3

par(mfrow=c(4,4))
for (comp1 in 1:4){
	for (comp2 in 1:4){
		if (comp1 != comp2){
	plot(p$x[,c(comp1,comp2)], 	col=as.numeric(pgroup[colnames(counts)]),pch=20,main=paste0("PCs ", comp1, ",", comp2))
#legend("topright", 	legend=unique(family),pch=20,col=1:length(unique(family)))
	}
}
}

# Look at gut samples only.
subset <- meta[which(meta$Tissue=="Gut"),"Customer_ID"]

columns <- intersect(subset, colnames(tmm))
x <- tmm[,columns]

x.meta <- meta[which(meta$Tissue=="Labial gland"),]
x.log <- normalize.voom(x)

#p <- prcomp(t(x.log[which(rowMeans(gut)>1),]))
p <- prcomp(t(x.log))

par(mfrow=c(5,4))
for (comp1 in 1:5){
	for (comp2 in 1:5){
		if (comp1 != comp2){
	plot(p$x[,c(comp1,comp2)], 	col=as.numeric(type[colnames(x)]),pch=20,main=paste0("PCs ", comp1, ",", comp2))
#legend("topright", 	legend=unique(family),pch=20,col=1:length(unique(family)))
	}
}
}

# PCA with labels
comp1 <- 1
comp2 <- 3
plot(p$x[,c(comp1,comp2)], 	col=as.numeric(type[colnames(gut)]),pch=20,main=paste0("PCs ", comp1, ",", comp2))
legend("topright", 	legend=unique(type),pch=20,col=1:length(unique(type)))

display.names <- paste(meta$Host_Plant,meta$Phylogeny_group,meta$Host.Plant.use,sep=".")
names(display.names) <- meta$Customer_ID
display.names.gut <- display.names[colnames(gut)]

text(p$x[,comp1],p$x[,comp2],labels=display.names.gut,cex=0.6)

# Clustering with bootstrapping

# Using contigs with CPM > 1
res <- pvclust(x.log[which(rowMeans(x)>1),],nboot=100,method.hclust="complete")

# Using all contigs
res <- pvclust(x.log,nboot=100,method.hclust="complete")

# 1000 bootstrap samples for all contigs, complete linkage
res <- pvclust(x.log,nboot=1000,method.hclust="complete")

save(res, file="gut_pvclust_complete_1000.Robj")

pdf("gut_pvclust_complete_1000.pdf")
plot(res)
dev.off()

# DESeq2
subset <- meta[which(meta$Tissue=="Gut"),"Customer_ID"]
columns <- intersect(subset, colnames(counts))
meta.gut <- meta[which(meta$Tissue=="Gut"),]

dds <- DESeqDataSetFromMatrix(countData = counts[,columns], colData = meta.gut[,c("Customer_ID","Host.Plant.use","Phylogeny_group","Host_Plant")], design = ~Phylogeny_group)

dds <- DESeq(dds, betaPrior=FALSE)
res <- results(dds)
sig <- res[which(res$padj<0.001),]
sig.o <- sig[order(sig$padj),]
sig.top <- sig.o[1:100,]


