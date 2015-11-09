Cluster analysis of butterfly transcriptomes
========================================================

Cluster & analyze only those mapping to an Hme ortholog.


```r
library(pheatmap)
library(pvclust)
```

```
## Warning: package 'pvclust' was built under R version 3.1.2
```

```r
library(DESeq2)
```

```
## Warning: package 'DESeq2' was built under R version 3.1.2
```

```
## Loading required package: S4Vectors
## Loading required package: stats4
## Loading required package: BiocGenerics
```

```
## Warning: package 'BiocGenerics' was built under R version 3.1.2
```

```
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist, unsplit
## 
## Loading required package: IRanges
```

```
## Warning: package 'IRanges' was built under R version 3.1.2
```

```
## Loading required package: GenomicRanges
```

```
## Warning: package 'GenomicRanges' was built under R version 3.1.2
```

```
## Loading required package: GenomeInfoDb
```

```
## Warning: package 'GenomeInfoDb' was built under R version 3.1.2
```

```
## Loading required package: Rcpp
## Loading required package: RcppArmadillo
```

```
## Warning: package 'RcppArmadillo' was built under R version 3.1.2
```
Let's define a function for plotting a heatmap with more descriptive names.


```r
labelled_heatmap <- function(data, meta){
  corrs <- cor(data)
	display.names <- paste(meta$Host_Plant,meta$Tissue,meta$Host.Plant.use,sep=".")
	names(display.names) <- meta$Customer_ID
	colnames(corrs) <- display.names[colnames(corrs)]
	rownames(corrs) <- colnames(corrs)
	pheatmap(corrs, cex=0.8)
}
```
... and some functions for normalization and pre-processing of counts.


```r
normalize.voom <- function(counts){
  require(limma)
	return(voom(counts)$E)
}
 
cpm.tmm <- function(counts, groups=NA){
	require(edgeR)
	if(is.na(groups)){
		d<-DGEList(counts=counts)
	}
	else{
		d<-DGEList(counts=counts, group=groups)
	}
	d <- calcNormFactors(d, method="TMM") 
	return(cpm(d, normalized.lib.sizes=TRUE))
}

common.disp <- function(counts_, group_){ 
d.edge <- DGEList(counts=counts_,group=group_)
d.edge <- calcNormFactors(d.edge)
d.edge <- estimateCommonDisp(d.edge)
return(d.edge$common.dispersion)
}

common.disp.nogrp <- function(counts_){ 
d.edge <- DGEList(counts=counts_)
d.edge <- calcNormFactors(d.edge)
d.edge <- estimateCommonDisp(d.edge)
return(d.edge$common.dispersion)
}
```
Read count table and metadata (information about samples). We remove the failed sample (B1MT) at once to keep the metadata table "in sync" with the count table - it's more convenient if they inculde exactly the same samples.


```r
counts <- read.delim("genewise_readcount2.txt",sep="\t",row.names=1)
meta  <- read.delim("rna-seqecoevowabi_relational_table2.txt",colClasses=c(rep("factor",7),"numeric"))
#meta <- meta[-which(meta$Customer_ID=="B1MT"),]
```
Define some vectors that will be useful for coloring plots later.


```r
tissue <- meta$Tissue
names(tissue) <- meta$Customer_ID
hostplant <- meta$Host_Plant
names(hostplant) <- meta$Customer_ID
pgroup <- meta$Phylogeny_group
names(pgroup) <- meta$Customer_ID
type <- meta$Host.Plant.use
names(type) <- meta$Customer_ID
```

Normalize and do a principal component analysis.


```r
tmm <- cpm.tmm(counts)
```

```
## Loading required package: edgeR
```

```
## Warning: package 'edgeR' was built under R version 3.1.2
```

```
## Loading required package: limma
## 
## Attaching package: 'limma'
## 
## The following object is masked from 'package:DESeq2':
## 
##     plotMA
## 
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
```

```r
labelled_heatmap(tmm,meta)
```

![plot of chunk :normalize_PCA](figure/:normalize_PCA-1.png) 

```r
log.cpm.tmm <- normalize.voom(tmm)

columns <- intersect(meta$Customer_ID, colnames(counts))
norm.counts <- log.cpm.tmm[,columns]

temp <- rbind(HostPlantUse=as.character(meta$Host.Plant.use), norm.counts)
temp2 <- rbind(HostPlant=as.character(meta$Host_Plant), temp)
temp3 <- rbind(Tissue=as.character(meta$Tissue), temp2)
temp4 <- rbind(Family=as.character(meta$Family), temp3)
simca <- rbind(PhylogenyGroup=as.character(meta$Phylogeny_group), temp4)

write.table(simca, file="normalized_ortho_counts.txt",sep="\t",quote=F)

labelled_heatmap(log.cpm.tmm,meta)
```

![plot of chunk :normalize_PCA](figure/:normalize_PCA-2.png) 

```r
p <- prcomp(t(log.cpm.tmm))
```

Try some Anova stuff suggested on BioStar.

```r
library(reshape)
```

```
## 
## Attaching package: 'reshape'
## 
## The following object is masked from 'package:IRanges':
## 
##     expand
## 
## The following object is masked from 'package:S4Vectors':
## 
##     rename
```

```r
m <- melt(norm.counts)
colnames(m) <- c("gene_ID","sample_ID","log-CPM")
tissue <- rep(meta$Tissue, each=nrow(norm.counts))
hostplantuse <- rep(meta$Host.Plant.use, each=nrow(norm.counts))
family <- rep(meta$Family, each=nrow(norm.counts))
phylogroup <- rep(meta$Phylogeny_group, each=nrow(norm.counts))
data <- data.frame(m, tissue=tissue, host_plant_use=hostplantuse, phylogeny_group=phylogroup, family=family)
subset <- data[sample(1:nrow(data), 50000),]
#fit <- lm(log.CPM ~ tissue + host_plant_use + phylogeny_group + family, data=data)
#fit <- lm(log.CPM ~ tissue + phylogeny_group + host_plant_use + family, data=data)

pdf("anova_comparison.pdf")
par(mfrow=c(3,1))
#barplot(a$"F value",names.arg=rownames(a), main="Anova on log-CPM",cex.names=0.7)
#dev.off()
# The above is dubious because host plant use and phylogeny group are strongly collinear, so drop each of these in turn
# Switch order
fit <- lm(log.CPM ~ phylogeny_group + host_plant_use + tissue, data=data)
a <- anova(fit)
barplot(a$"F value",names.arg=rownames(a), main="Anova on log-CPM",cex.names=0.7)

fit <- lm(log.CPM ~ host_plant_use + phylogeny_group + tissue, data=data)
a <- anova(fit)
#pdf("anova_2.pdf")
barplot(a$"F value",names.arg=rownames(a), main="Anova on log-CPM",cex.names=0.7)
#dev.off()
# Switch order
fit <- lm(log.CPM ~ tissue + phylogeny_group + host_plant_use, data=data)
a <- anova(fit)
#pdf("anova_3.pdf")
barplot(a$"F value",names.arg=rownames(a), main="Anova on log-CPM",cex.names=0.7)
dev.off()
```

```
## pdf 
##   2
```

```r
# The above is dubious because host plant use and phylogeny group are strongly collinear, so drop each of these in turn
# Without host plant use
fit <- lm(log.CPM ~ tissue + phylogeny_group + family, data=data)
a <- anova(fit)
pdf("anova_no_hostplantuse.pdf")
barplot(a$"F value",names.arg=rownames(a), main="Anova on log-CPM",cex.names=0.7)
dev.off()
```

```
## pdf 
##   2
```

```r
# Without phylogroup
fit <- lm(log.CPM ~ tissue + host_plant_use + family, data=data)
a <- anova(fit)
pdf("anova_no_phylogroup.pdf")
barplot(a$"F value",names.arg=rownames(a), main="Anova on log-CPM",cex.names=0.7)
dev.off()
```

```
## pdf 
##   2
```

```r
# Try with interactions
pdf("anova_interactions.pdf")
# (a) panel
par(mfrow=c(3,1))
fit <- lm(log.CPM ~  phylogeny_group + host_plant_use + tissue + tissue:phylogeny_group  + tissue:host_plant_use, data=data)
a <- anova(fit)
barplot(a$"F value"[-6],names.arg=rownames(a)[-6], main="Anova on log-CPM",cex.names=0.7,ylim=c(0,1000))
# (b) panel
fit <- lm(log.CPM ~ host_plant_use + phylogeny_group + tissue + tissue:phylogeny_group + tissue:host_plant_use, data=data)
a <- anova(fit)
barplot(a$"F value"[-6],names.arg=rownames(a)[-6], main="Anova on log-CPM",cex.names=0.7,ylim=c(0,1000))
# (c) panel
fit <- lm(log.CPM ~ tissue + phylogeny_group + host_plant_use + tissue:phylogeny_group + tissue:host_plant_use, data=data)
a <- anova(fit)
barplot(a$"F value"[-6],names.arg=rownames(a)[-6], main="Anova on log-CPM",cex.names=0.7,ylim=c(0,1000))
dev.off()
```

```
## pdf 
##   2
```

Visualize two of the components, colored by tissue.


```r
comp1 <- 1
comp2 <- 2
plot(p$x[,c(comp1,comp2)], col=as.numeric(tissue[colnames(counts)]),pch=20,main=paste0("PCs ", comp1, ",", comp2))
legend("topright", legend=c("Fat body","Gut","Labial gland","Malpighian tubules"),pch=20,col=1:4)
```

![plot of chunk :pca_tissue](figure/:pca_tissue-1.png) 



```r
comp1 <- 1
comp2 <- 2
par(mfrow=c(2,2))
# Color by tissue
plot(p$x[,c(comp1,comp2)], col=as.numeric(tissue[colnames(counts)]),pch=20,main=paste0("Tissue, PCs ", comp1, ",", comp2))
legend("topright", legend=c("Fat body","Gut","Labial gland","Malpighian tubules"),pch=20,col=1:4)
# Color by host plant
plot(p$x[,c(comp1,comp2)], col=as.numeric(hostplant[colnames(counts)]),pch=20,main=paste0("Host plant, PCs ", comp1, ",", comp2))
legend("topright", legend=unique(hostplant),pch=20,col=1:length(unique(hostplant)))
# Color by phylogeny group
plot(p$x[,c(comp1,comp2)], col=as.numeric(pgroup[colnames(counts)]),pch=20,main=paste0("Phylogroup, PCs ", comp1, ",", comp2))
legend("topright", legend=unique(pgroup),pch=20,col=1:length(unique(pgroup)))
# Color by core vs extended
plot(p$x[,c(comp1,comp2)], col=as.numeric(type[colnames(counts)]),pch=20,main=paste0("CoreVsExtended, PCs ", comp1, ",", comp2))
legend("topright", legend=unique(type),pch=20,col=1:length(unique(type)))
```

![plot of chunk :pca_four_colorings](figure/:pca_four_colorings-1.png) 

Loop over combinations of PCs and color by family. 

```r
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
```

![plot of chunk :pca_loop](figure/:pca_loop-1.png) 

Let's look at gut samples only, colored by core vs extended.


```r
subset <- meta[which(meta$Tissue=="Gut"),"Customer_ID"]
columns <- intersect(subset, colnames(tmm))
x <- tmm[,columns]

x.meta <- meta[which(meta$Tissue=="Labial gland"),]
x.log <- normalize.voom(x)

#p <- prcomp(t(x.log[which(rowMeans(x.log)>1),]))
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
```

![plot of chunk :pca_gut](figure/:pca_gut-1.png) 

PCA of gut samples with more informative text labels, colored by core vs extended.


```r
par(mfrow=c(1,1))
comp1 <- 1
comp2 <- 2
plot(p$x[,c(comp1,comp2)], 	col=as.numeric(type[colnames(x.log)]),pch=20,main=paste0("PCs ", comp1, ",", comp2))
legend("topright", 	legend=unique(type),pch=20,col=1:length(unique(type)))
display.names <- paste(meta$Host_Plant,meta$Phylogeny_group,meta$Host.Plant.use,sep=".")
names(display.names) <- meta$Customer_ID
display.names.gut <- display.names[colnames(x.log)]
text(p$x[,comp1],p$x[,comp2],labels=display.names.gut,cex=0.6)
```

![plot of chunk :pca_gut_textlabels](figure/:pca_gut_textlabels-1.png) 

Try clustering with bootstrapping (pvclust) on gut samples.
We have commented out the actual commands here because they take a long time to run. Instead, we read a previously saved version of the pvclust output.


```r
# Using contigs with CPM > 1
# res <- pvclust(x.log[which(rowMeans(x)>1),],nboot=100,method.hclust="complete")
# Using all contigs
# res <- pvclust(x.log,nboot=100,method.hclust="complete")
# 1000 bootstrap samples for all contigs, complete linkage
# res <- pvclust(x.log,nboot=1000,method.hclust="complete")
#save(res, file="gut_pvclust_complete_1000.Robj")
#pdf("gut_ortho_pvclust_complete_100.pdf")
#plot(res)
#dev.off()
```

Core vs extended differential gene expression in gut
====================================================
Using only gut samples (n=18).


```r
# DESeq2
subset <- meta[which(meta$Tissue=="Gut"),"Customer_ID"]
columns <- intersect(subset, colnames(counts))
meta.gut <- meta[which(meta$Tissue=="Gut"),]

dds <- DESeqDataSetFromMatrix(countData = counts[,columns], colData = meta.gut[,c("Customer_ID","Host.Plant.use","Phylogeny_group","Host_Plant","Family")], design = ~Family+Host.Plant.use)

dds <- DESeq(dds, betaPrior=FALSE)
```

```
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
```

```r
res <- results(dds)
sig <- res[which(res$padj<0.01),]
sig.o <- sig[order(sig$padj),]
print(dim(sig.o))
```

```
## [1] 224   6
```

```r
write.table(sig.o,file="sig_gut_corevsextended_DESeq2_0.01.txt",quote=F)
# 214 
```

Labial gland
-------------

```r
subset <- meta[which(meta$Tissue=="Labial gland"),"Customer_ID"]
columns <- intersect(subset, colnames(counts))
meta.lab <- meta[which(meta$Tissue=="Labial gland"),]

dds <- DESeqDataSetFromMatrix(countData = counts[,columns], colData = meta.lab[,c("Customer_ID","Host.Plant.use","Phylogeny_group","Host_Plant","Family")], design = ~Family+Host.Plant.use)

dds <- DESeq(dds, betaPrior=FALSE)
```

```
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
```

```r
res <- results(dds)
sig <- res[which(res$padj<0.01),]
sig.o <- sig[order(sig$padj),]
print(dim(sig.o))
```

```
## [1] 289   6
```

```r
write.table(sig.o,file="sig_labialgland_corevsextended_DESeq2_0.01.txt",quote=F)
```

Fat body
-------------

```r
subset <- meta[which(meta$Tissue=="Fat body"),"Customer_ID"]
columns <- intersect(subset, colnames(counts))
meta.fat <- meta[which(meta$Tissue=="Fat body"),]

dds <- DESeqDataSetFromMatrix(countData = counts[,columns], colData = meta.fat[,c("Customer_ID","Host.Plant.use","Phylogeny_group","Host_Plant","Family")], design = ~Family+Host.Plant.use)

dds <- DESeq(dds, betaPrior=FALSE)
```

```
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
```

```r
res <- results(dds)
sig <- res[which(res$padj<0.01),]
sig.o <- sig[order(sig$padj),]
print(dim(sig.o))
```

```
## [1] 692   6
```

```r
write.table(sig.o,file="sig_fatbody_corevsextended_DESeq2_0.01.txt",quote=F)
```

Malpighian tubules
-------------

```r
subset <- meta[which(meta$Tissue=="Malpighian tubules"),"Customer_ID"]
columns <- intersect(subset, colnames(counts))
meta.mal <- meta[which(meta$Tissue=="Malpighian tubules"),]
q
```

```
## function (save = "default", status = 0, runLast = TRUE) 
## .Internal(quit(save, status, runLast))
## <bytecode: 0x10dfaa2b0>
## <environment: namespace:base>
```

```r
dds <- DESeqDataSetFromMatrix(countData = counts[,columns], colData = meta.mal[,c("Customer_ID","Host.Plant.use","Phylogeny_group","Host_Plant","Family")], design = ~Family+Host.Plant.use)

dds <- DESeq(dds, betaPrior=FALSE)
```

```
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
```

```r
res <- results(dds)
sig <- res[which(res$padj<0.01),]
sig.o <- sig[order(sig$padj),]
print(dim(sig.o))
```

```
## [1] 81  6
```

```r
write.table(sig.o,file="sig_malpi_corevsextended_DESeq2_0.01.txt",quote=F)
```

Use edgeR to obtain sets of DE orthologs for a variety of comparisons.
==========================================================
Fit a model using several covariates. 


```r
library(edgeR)
columns <- intersect(meta$Customer_ID, colnames(counts))
counts_ <- counts[,columns]
design=model.matrix(~0+Tissue+as.factor(Family)+Host.Plant.use+Phylogeny_group,data=meta)
y <- DGEList(counts_, group=meta$Host.Plant.use)
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
```

```
## Loading required package: splines
```

```r
y <- estimateGLMTagwiseDisp(y,design)
edger.fit <- glmFit(y,design)
```
Extract DE orthocontigs for all the four tissues wrt the other tissues.

```r
head(edger.fit$coefficient)
```

```
##               TissueFat body  TissueGut TissueLabial gland
## HMEL000002-PA     -11.125685 -11.411361         -10.848125
## HMEL000003-PA     -14.545179 -14.952617         -13.978157
## HMEL000004-PA      -9.933157  -9.431300          -9.465831
## HMEL000006-PA     -12.005264  -8.646108          -9.402437
## HMEL000007-PA     -19.821256 -18.745180         -19.208857
## HMEL000008-PA     -11.257788 -11.406644         -10.915072
##               TissueMalpighian tubules as.factor(Family)2
## HMEL000002-PA                -9.218108        -0.04247922
## HMEL000003-PA               -14.043592         0.01370774
## HMEL000004-PA                -8.540449        -0.01675743
## HMEL000006-PA                -8.375376        -0.05433070
## HMEL000007-PA               -19.676228        -0.12131947
## HMEL000008-PA               -11.231059        -0.05092760
##               as.factor(Family)3 Host.Plant.useExtended
## HMEL000002-PA       -0.102952091            -0.05167402
## HMEL000003-PA        0.087635349             0.43358056
## HMEL000004-PA       -0.002908794            -0.05721833
## HMEL000006-PA       -0.102051527            -0.03826289
## HMEL000007-PA        0.430874760             1.27850089
## HMEL000008-PA       -0.195917674            -0.22302530
##               Phylogeny_groupRanunculales Phylogeny_groupRosids
## HMEL000002-PA                   0.2272732           -0.01041724
## HMEL000003-PA                  -0.4564330            0.13438742
## HMEL000004-PA                  -0.2026609            0.02535253
## HMEL000006-PA                  -0.7330188           -0.17288203
## HMEL000007-PA                  -0.1375566            1.03962847
## HMEL000008-PA                   0.2156424           -0.08926501
```

```r
library(lattice)
# Fat body vs. other tissues
lrt <- glmLRT(edger.fit,contrast=c(1,-1/3,-1/3,-1/3,0,0,0,0,0))
temp <- topTags(lrt,n=100000)
sig.fdr <- which(as.data.frame(temp)[,"FDR"] < 0.01)
sig.fb <- as.data.frame(temp)[sig.fdr,]
for(gene in head(rownames(sig.fb))){
  expr <- t(counts_[gene,]/colSums(counts_))
  temp <- data.frame(expr=expr, tissue=meta$Tissue, core_vs_ext=meta$Host.Plant.use)
  print(bwplot(expr~tissue, data=temp, main=gene))
  readline()
}
```

![plot of chunk :edger-tissues](figure/:edger-tissues-1.png) 
![plot of chunk :edger-tissues](figure/:edger-tissues-2.png) 
![plot of chunk :edger-tissues](figure/:edger-tissues-3.png) 
![plot of chunk :edger-tissues](figure/:edger-tissues-4.png) 
![plot of chunk :edger-tissues](figure/:edger-tissues-5.png) 
![plot of chunk :edger-tissues](figure/:edger-tissues-6.png) 

```r
write.table(as.data.frame(topTags(lrt,n=100000)),file="fatbody_vs_other_tissues.txt",quote=F,sep="\t")
# Gut vs other tissues
lrt <- glmLRT(edger.fit,contrast=c(-1/3,1,-1/3,-1/3,0,0,0,0,0))
temp <- topTags(lrt,n=100000)
sig.fdr <- which(as.data.frame(temp)[,"FDR"] < 0.01)
sig.gut <- as.data.frame(temp)[sig.fdr,]
for(gene in head(rownames(sig.gut))){
  expr <- t(counts_[gene,]/colSums(counts_))
  temp <- data.frame(expr=expr, tissue=meta$Tissue, core_vs_ext=meta$Host.Plant.use)
  print(bwplot(expr~tissue, data=temp, main=gene))
  readline()
}
```

![plot of chunk :edger-tissues](figure/:edger-tissues-7.png) 
![plot of chunk :edger-tissues](figure/:edger-tissues-8.png) 
![plot of chunk :edger-tissues](figure/:edger-tissues-9.png) 
![plot of chunk :edger-tissues](figure/:edger-tissues-10.png) 
![plot of chunk :edger-tissues](figure/:edger-tissues-11.png) 
![plot of chunk :edger-tissues](figure/:edger-tissues-12.png) 

```r
write.table(as.data.frame(topTags(lrt,n=100000)),file="gut_vs_other_tissues.txt",quote=F,sep="\t")
# Labial gland vs other tissues
lrt <- glmLRT(edger.fit,contrast=c(-1/3,-1/3,1,-1/3,0,0,0,0,0))
temp <- topTags(lrt,n=100000)
sig.fdr <- which(as.data.frame(temp)[,"FDR"] < 0.01)
sig.lg <- as.data.frame(temp)[sig.fdr,]
for(gene in head(rownames(sig.lg))){
  expr <- t(counts_[gene,]/colSums(counts_))
  temp <- data.frame(expr=expr, tissue=meta$Tissue, core_vs_ext=meta$Host.Plant.use)
  print(bwplot(expr~tissue, data=temp, main=gene))
  readline()
}
```

![plot of chunk :edger-tissues](figure/:edger-tissues-13.png) 
![plot of chunk :edger-tissues](figure/:edger-tissues-14.png) 
![plot of chunk :edger-tissues](figure/:edger-tissues-15.png) 
![plot of chunk :edger-tissues](figure/:edger-tissues-16.png) 
![plot of chunk :edger-tissues](figure/:edger-tissues-17.png) 
![plot of chunk :edger-tissues](figure/:edger-tissues-18.png) 

```r
write.table(as.data.frame(topTags(lrt,n=100000)),file="labialgland_vs_other_tissues.txt",quote=F,sep="\t")

# Malpighian tubules vs other tissues
lrt <- glmLRT(edger.fit,contrast=c(-1/3,-1/3,-1/3,1,0,0,0,0,0))
temp <- topTags(lrt,n=100000)
sig.fdr <- which(as.data.frame(temp)[,"FDR"] < 0.01)
sig.mt <- as.data.frame(temp)[sig.fdr,]
for(gene in head(rownames(sig.mt))){
  expr <- t(counts_[gene,]/colSums(counts_))
  temp <- data.frame(expr=expr, tissue=meta$Tissue, core_vs_ext=meta$Host.Plant.use)
  print(bwplot(expr~tissue, data=temp, main=gene))
  readline()
}
```

![plot of chunk :edger-tissues](figure/:edger-tissues-19.png) 
![plot of chunk :edger-tissues](figure/:edger-tissues-20.png) 
![plot of chunk :edger-tissues](figure/:edger-tissues-21.png) 
![plot of chunk :edger-tissues](figure/:edger-tissues-22.png) 
![plot of chunk :edger-tissues](figure/:edger-tissues-23.png) 
![plot of chunk :edger-tissues](figure/:edger-tissues-24.png) 

```r
write.table(as.data.frame(topTags(lrt,n=100000)),file="malpighiantubules_vs_other_tissues.txt",quote=F,sep="\t")
```

Now test plant groups against each other. (Doesn't work well to compare one against 2)

```r
design=model.matrix(~0+Phylogeny_group+Tissue+as.factor(Family)+Host.Plant.use,data=meta)
y <- DGEList(counts_, group=meta$Host.Plant.use)
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)
edger.fit <- glmFit(y,design)

# Asterids vs Ranunculales.
lrt <- glmLRT(edger.fit,contrast=c(1,-1,0,0,0,0,0,0,0))
temp <- topTags(lrt,n=100000)
sig.fdr <- which(as.data.frame(temp)[,"FDR"] < 0.01)
sig.ast <- as.data.frame(temp)[sig.fdr,]
for(gene in head(rownames(sig.ast))){
  expr <- t(counts_[gene,]/colSums(counts_))
  temp <- data.frame(expr=expr, phylogroup=meta$Phylogeny_group, core_vs_ext=meta$Host.Plant.use)
  print(bwplot(expr~phylogroup, data=temp, main=gene))
  readline()
}
```

![plot of chunk :edger-phylo](figure/:edger-phylo-1.png) 
![plot of chunk :edger-phylo](figure/:edger-phylo-2.png) 
![plot of chunk :edger-phylo](figure/:edger-phylo-3.png) 
![plot of chunk :edger-phylo](figure/:edger-phylo-4.png) 
![plot of chunk :edger-phylo](figure/:edger-phylo-5.png) 
![plot of chunk :edger-phylo](figure/:edger-phylo-6.png) 

```r
write.table(as.data.frame(topTags(lrt,n=100000)),file="Asterids_vs_Ranunculales.txt",quote=F,sep="\t")

# Ranunculales vs Rosids.
lrt <- glmLRT(edger.fit,contrast=c(0,1,-1,0,0,0,0,0,0))
temp <- topTags(lrt,n=100000)
sig.fdr <- which(as.data.frame(temp)[,"FDR"] < 0.01)
sig.ran <- as.data.frame(temp)[sig.fdr,]
for(gene in head(rownames(sig.ran))){
  expr <- t(counts_[gene,]/colSums(counts_))
  temp <- data.frame(expr=expr, phylogroup=meta$Phylogeny_group, core_vs_ext=meta$Host.Plant.use)
  print(bwplot(expr~phylogroup, data=temp, main=gene))
  readline()
}
```

![plot of chunk :edger-phylo](figure/:edger-phylo-7.png) 
![plot of chunk :edger-phylo](figure/:edger-phylo-8.png) 
![plot of chunk :edger-phylo](figure/:edger-phylo-9.png) 
![plot of chunk :edger-phylo](figure/:edger-phylo-10.png) 
![plot of chunk :edger-phylo](figure/:edger-phylo-11.png) 
![plot of chunk :edger-phylo](figure/:edger-phylo-12.png) 

```r
write.table(as.data.frame(topTags(lrt,n=100000)),file="Ranunculales_vs_Rosids.txt",quote=F,sep="\t")

# Rosids vs Asterids.
lrt <- glmLRT(edger.fit,contrast=c(-1,0,1,0,0,0,0,0,0))
temp <- topTags(lrt,n=100000)
sig.fdr <- which(as.data.frame(temp)[,"FDR"] < 0.01)
sig.ros <- as.data.frame(temp)[sig.fdr,]
for(gene in head(rownames(sig.ros))){
  expr <- t(counts_[gene,]/colSums(counts_))
  temp <- data.frame(expr=expr, phylogroup=meta$Phylogeny_group, core_vs_ext=meta$Host.Plant.use)
  print(bwplot(expr~phylogroup, data=temp, main=gene))
  readline()
}
```

![plot of chunk :edger-phylo](figure/:edger-phylo-13.png) 
![plot of chunk :edger-phylo](figure/:edger-phylo-14.png) 
![plot of chunk :edger-phylo](figure/:edger-phylo-15.png) 
![plot of chunk :edger-phylo](figure/:edger-phylo-16.png) 
![plot of chunk :edger-phylo](figure/:edger-phylo-17.png) 
![plot of chunk :edger-phylo](figure/:edger-phylo-18.png) 

```r
write.table(as.data.frame(topTags(lrt,n=100000)),file="Rosids_vs_Asterids.txt",quote=F,sep="\t")
```


Now test core vs extended within each tissue. (already done above with DESeq2 but let's try edgeR as well)


```r
# Fat body
subset <- meta[which(meta$Tissue=="Fat body"),"Customer_ID"]
columns <- intersect(subset, colnames(counts))
meta.fat <- meta[which(meta$Tissue=="Fat body"),]

design=model.matrix(~0+as.factor(Family)+Phylogeny_group+Host.Plant.use,data=meta.fat)
y <- DGEList(counts[,columns], group=meta.fat$Host.Plant.use)
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)

edger.fit <- glmFit(y, design)
lrt <- glmLRT(edger.fit,coef="Host.Plant.useExtended")

temp <- topTags(lrt,n=100000)
sig.fdr <- which(as.data.frame(temp)[,"FDR"] < 0.01)
sig.ce.fb <- as.data.frame(temp)[sig.fdr,]
for(gene in head(rownames(sig.ce.fb))){
  expr <- t(counts_[gene,]/colSums(counts_))
  temp <- data.frame(expr=expr, tissue=meta$Tissue, core_vs_ext=meta$Host.Plant.use)
  print(bwplot(expr~core_vs_ext|tissue, data=temp, main=gene))
  readline()
}
```

![plot of chunk :edger-core-ext](figure/:edger-core-ext-1.png) 
![plot of chunk :edger-core-ext](figure/:edger-core-ext-2.png) 
![plot of chunk :edger-core-ext](figure/:edger-core-ext-3.png) 
![plot of chunk :edger-core-ext](figure/:edger-core-ext-4.png) 
![plot of chunk :edger-core-ext](figure/:edger-core-ext-5.png) 
![plot of chunk :edger-core-ext](figure/:edger-core-ext-6.png) 

```r
write.table(as.data.frame(topTags(lrt,n=100000)),file="Fat_body_core_vs_extended.txt",quote=F,sep="\t")

# Gut
subset <- meta[which(meta$Tissue=="Gut"),"Customer_ID"]
columns <- intersect(subset, colnames(counts))
meta.gut <- meta[which(meta$Tissue=="Gut"),]

design=model.matrix(~0+as.factor(Family)+Phylogeny_group+Host.Plant.use,data=meta.gut)
y <- DGEList(counts[,columns], group=meta.gut$Host.Plant.use)
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)

edger.fit <- glmFit(y, design)
lrt <- glmLRT(edger.fit,coef="Host.Plant.useExtended")

temp <- topTags(lrt,n=100000)
sig.fdr <- which(as.data.frame(temp)[,"FDR"] < 0.01)
sig.ce.gut <- as.data.frame(temp)[sig.fdr,]
for(gene in head(rownames(sig.ce.gut))){
  expr <- t(counts_[gene,]/colSums(counts_))
  temp <- data.frame(expr=expr, tissue=meta$Tissue, core_vs_ext=meta$Host.Plant.use)
  print(bwplot(expr~core_vs_ext|tissue, data=temp, main=gene))
  readline()
}
```

![plot of chunk :edger-core-ext](figure/:edger-core-ext-7.png) 
![plot of chunk :edger-core-ext](figure/:edger-core-ext-8.png) 
![plot of chunk :edger-core-ext](figure/:edger-core-ext-9.png) 
![plot of chunk :edger-core-ext](figure/:edger-core-ext-10.png) 
![plot of chunk :edger-core-ext](figure/:edger-core-ext-11.png) 
![plot of chunk :edger-core-ext](figure/:edger-core-ext-12.png) 

```r
write.table(as.data.frame(topTags(lrt,n=100000)),file="Gut_core_vs_extended.txt",quote=F,sep="\t")


# Labial gland
subset <- meta[which(meta$Tissue=="Labial gland"),"Customer_ID"]
columns <- intersect(subset, colnames(counts))
meta.lab <- meta[which(meta$Tissue=="Labial gland"),]

design=model.matrix(~0+as.factor(Family)+Phylogeny_group+Host.Plant.use,data=meta.lab)
y <- DGEList(counts[,columns], group=meta.lab$Host.Plant.use)
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)

edger.fit <- glmFit(y, design)
lrt <- glmLRT(edger.fit,coef="Host.Plant.useExtended")

temp <- topTags(lrt,n=100000)
sig.fdr <- which(as.data.frame(temp)[,"FDR"] < 0.01)
sig.ce.lg <- as.data.frame(temp)[sig.fdr,]
for(gene in head(rownames(sig.ce.lg))){
  expr <- t(counts_[gene,]/colSums(counts_))
  temp <- data.frame(expr=expr, tissue=meta$Tissue, core_vs_ext=meta$Host.Plant.use)
  print(bwplot(expr~core_vs_ext|tissue, data=temp, main=gene))
  readline()
}
```

![plot of chunk :edger-core-ext](figure/:edger-core-ext-13.png) 
![plot of chunk :edger-core-ext](figure/:edger-core-ext-14.png) 
![plot of chunk :edger-core-ext](figure/:edger-core-ext-15.png) 
![plot of chunk :edger-core-ext](figure/:edger-core-ext-16.png) 
![plot of chunk :edger-core-ext](figure/:edger-core-ext-17.png) 
![plot of chunk :edger-core-ext](figure/:edger-core-ext-18.png) 

```r
write.table(as.data.frame(topTags(lrt,n=100000)),file="Labial_gland_core_vs_extended.txt",quote=F,sep="\t")


# Malpighian tubules
subset <- meta[which(meta$Tissue=="Malpighian tubules"),"Customer_ID"]
columns <- intersect(subset, colnames(counts))
meta.mal <- meta[which(meta$Tissue=="Malpighian tubules"),]

design=model.matrix(~0+as.factor(Family)+Phylogeny_group+Host.Plant.use,data=meta.mal)
y <- DGEList(counts[,columns], group=meta.mal$Host.Plant.use)
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)

edger.fit <- glmFit(y, design)
lrt <- glmLRT(edger.fit,coef="Host.Plant.useExtended")

temp <- topTags(lrt,n=100000)
sig.fdr <- which(as.data.frame(temp)[,"FDR"] < 0.01)
sig.ce.mt <- as.data.frame(temp)[sig.fdr,]
for(gene in head(rownames(sig.ce.mt))){
  expr <- t(counts[gene,columns]/colSums(counts[,columns]))
  temp <- data.frame(expr=expr, core_vs_ext=meta.mal$Host.Plant.use)
  print(bwplot(expr~core_vs_ext, data=temp, main=gene))
  readline()
}
```

![plot of chunk :edger-core-ext](figure/:edger-core-ext-19.png) 
![plot of chunk :edger-core-ext](figure/:edger-core-ext-20.png) 
![plot of chunk :edger-core-ext](figure/:edger-core-ext-21.png) 
![plot of chunk :edger-core-ext](figure/:edger-core-ext-22.png) 
![plot of chunk :edger-core-ext](figure/:edger-core-ext-23.png) 
![plot of chunk :edger-core-ext](figure/:edger-core-ext-24.png) 

```r
write.table(as.data.frame(topTags(lrt,n=100000)),file="Malpighian_tubules_core_vs_extended.txt",quote=F,sep="\t")
```

Test for differences between Family 1 and Family 3.
First test regardless of tissue. 

```{edger-family}
design=model.matrix(~0+as.factor(Family)+Phylogeny_group+Host.Plant.use+Tissue,data=meta)
y <- DGEList(counts, group=meta$Family)
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)

edger.fit <- glmFit(y, design)
lrt <- glmLRT(edger.fit,contrast=c(1,0,-1,0,0,0,0,0,0))

temp <- topTags(lrt,n=100000)
sig.fdr <- which(as.data.frame(temp)[,"FDR"] < 0.01)
sig.fam <- as.data.frame(temp)[sig.fdr,]
for(gene in head(rownames(sig.fam))){
  expr <- t(counts[gene,]/colSums(counts[,]))
  temp <- data.frame(expr=expr, fam=meta$Family)
  print(bwplot(expr~fam, data=temp, main=gene))
  readline()
}

write.table(as.data.frame(topTags(lrt,n=100000)),file="Malpighian_tubules_core_vs_extended.txt",quote=F,sep="\t")

```

In nettle.


```r
# Nettle
subset <- meta[which(meta$Host_Plant=="Nettle"),"Customer_ID"]
columns <- intersect(subset, colnames(counts))
meta.net <- meta[which(meta$Host_Plant=="Nettle"),]

design=model.matrix(~0+as.factor(Family)+Tissue,data=meta.net)
y <- DGEList(counts[,columns], group=meta.net$Family)
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)

edger.fit <- glmFit(y, design)
lrt <- glmLRT(edger.fit,contrast=c(1,0,-1,0,0,0))

temp <- topTags(lrt,n=100000)
sig.fdr <- which(as.data.frame(temp)[,"FDR"] < 0.05)
sig.fam.net <- as.data.frame(temp)[sig.fdr,]
for(gene in head(rownames(sig.fam.net))){
  expr <- t(counts[gene,columns]/colSums(counts[,columns]))
  temp <- data.frame(expr=expr, fam=meta.net$Family)
  print(bwplot(expr~fam, data=temp, main=gene))
  readline()
}
```

![plot of chunk :nettle](figure/:nettle-1.png) 
![plot of chunk :nettle](figure/:nettle-2.png) 
![plot of chunk :nettle](figure/:nettle-3.png) 
![plot of chunk :nettle](figure/:nettle-4.png) 
![plot of chunk :nettle](figure/:nettle-5.png) 
![plot of chunk :nettle](figure/:nettle-6.png) 

```r
write.table(as.data.frame(topTags(lrt,n=100000)),file="Fam1vsFam3_nettle.txt",quote=F,sep="\t")
```

In thistle.


```r
# Malpighian tubules
subset <- meta[which(meta$Host_Plant=="Thistle"),"Customer_ID"]
columns <- intersect(subset, colnames(counts))
meta.this <- meta[which(meta$Host_Plant=="Thistle"),]

design=model.matrix(~0+as.factor(Family)+Tissue,data=meta.this)
y <- DGEList(counts[,columns], group=meta.this$Family)
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)

edger.fit <- glmFit(y, design)
lrt <- glmLRT(edger.fit,contrast=c(1,0,-1,0,0,0))

temp <- topTags(lrt,n=100000)
sig.fdr <- which(as.data.frame(temp)[,"FDR"] < 0.05)
sig.fam.this <- as.data.frame(temp)[sig.fdr,]
for(gene in head(rownames(sig.fam.this))){
  expr <- t(counts[gene,columns]/colSums(counts[,columns]))
  temp <- data.frame(expr=expr, fam=meta.this$Family)
  print(bwplot(expr~fam, data=temp, main=gene))
  readline()
}
```

![plot of chunk :thistle](figure/:thistle-1.png) 
![plot of chunk :thistle](figure/:thistle-2.png) 
![plot of chunk :thistle](figure/:thistle-3.png) 
![plot of chunk :thistle](figure/:thistle-4.png) 
![plot of chunk :thistle](figure/:thistle-5.png) 
![plot of chunk :thistle](figure/:thistle-6.png) 

```r
write.table(as.data.frame(topTags(lrt,n=100000)),file="Fam1vsFam3_thistle.txt",quote=F,sep="\t")
```

Test the host plant use/tissue interaction. "Differences due to host plant use after accounting for tissue differences"


```r
design=model.matrix(~Tissue*Host.Plant.use,data=meta)
y <- DGEList(counts_, group=meta$Host.Plant.use)
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)
edger.fit <- glmFit(y,design)
lrt <- glmLRT(edger.fit,coef="Host.Plant.useExtended") 
temp <- topTags(lrt,n=10000)
sig.fdr <- which(as.data.frame(temp)[,"FDR"] < 0.01)
sig.edger <- as.data.frame(temp)[sig.fdr,]
```


```r
nf <- calcNormFactors(counts_)
design=model.matrix(~as.factor(Family)+Phylogeny_group+Tissue+Host.Plant.use,data=meta) # The same as for edgeR
voom.data <- voom(counts_, design, lib.size = colSums(counts_) * nf)

fit <- lmFit(voom.data, design)
fit <- eBayes(fit) # To calculate p values

fit2 <- contrasts.fit(fit, contrast.matrix)
```

```
## Error in contrasts.fit(fit, contrast.matrix): object 'contrast.matrix' not found
```

```r
fit2 <- eBayes(fit2)
```

```
## Error in ebayes(fit = fit, proportion = proportion, stdev.coef.lim = stdev.coef.lim, : object 'fit2' not found
```

```r
pvalues <- fit2$p.value
```

```
## Error in eval(expr, envir, enclos): object 'fit2' not found
```

```r
pdf("sig_ortho_count.pdf")
barplot(colSums(fit$p.value[,-1] < 0.05),las=2,names.arg=c("Fam2","Fam3","Ranunc","Rosids","Gut","Labial","Malpigh","Extended"),main="Orthologs with p<0.05 for coefficient")
dev.off()
```

```
## pdf 
##   2
```

```r
pdf("pvalue_dist_ortho.pdf")
boxplot(fit$p.value[,-1],las=2,names=c("Fam2","Fam3","Ranunc","Rosids","Gut","Labial","Malpigh","Extended"),main="p value distribution for coefficients")
dev.off()
```

```
## pdf 
##   2
```

```r
adj.p <- p.adjust(pvalues, method="BH")
```

```
## Error in p.adjust(pvalues, method = "BH"): object 'pvalues' not found
```

```r
alpha = 0.05
sig <- adj.p[adj.p<alpha]
```

```
## Error in eval(expr, envir, enclos): object 'adj.p' not found
```

```r
sig.o <- sig[order(sig)]
```

```
## Error in (function (classes, fdef, mtable) : unable to find an inherited method for function 'NSBS' for signature '"SimpleIntegerList"'
```

```r
print(paste(length(sig), "differentially expressed genes"))
```

```
## [1] "6 differentially expressed genes"
```

```r
# Ordinary t statistics
ord.t.stat <- fit2$coef/fit2$stdev.unscaled/fit2$sigma
```

```
## Error in eval(expr, envir, enclos): object 'fit2' not found
```

```r
ord.t.stat <- fit$coef/fit$stdev.unscaled/fit$sigma
```

Plot some examples of genes with low p values. These are differentially expressed in core vs extended overall, that is, considering all tissues.


```r
library(lattice)

for (gene in head(rownames(sig.edger))){
  expr <- t(counts_[gene,]/colSums(counts_))
  temp <- data.frame(expr=expr, tissue=meta$Tissue, core_vs_ext=meta$Host.Plant.use)
  print(bwplot(expr~core_vs_ext|tissue, data=temp, main=gene))
  readline()
}
```

![plot of chunk plot-genes](figure/plot-genes-1.png) 
![plot of chunk plot-genes](figure/plot-genes-2.png) 
![plot of chunk plot-genes](figure/plot-genes-3.png) 
![plot of chunk plot-genes](figure/plot-genes-4.png) 
![plot of chunk plot-genes](figure/plot-genes-5.png) 
![plot of chunk plot-genes](figure/plot-genes-6.png) 

Try to look for "specific"/"binary" contigs
===========================================


```r
plotcore <- function(genes, counts){
for (gene in genes){
  expr <- t(counts[gene,]/colSums(counts))
  temp <- data.frame(expr=expr, tissue=meta$Tissue, core=meta$Host.Plant.use)
  print(bwplot(expr~core|tissue, data=temp, main=gene))
  readline()
}
}
```


```r
columns <- intersect(meta$Customer_ID, colnames(counts))
counts_ <- counts[,columns]
core <- meta[which(meta$Host.Plant.use=="Core"),"Customer_ID"]
extended <- meta[which(meta$Host.Plant.use=="Extended"),"Customer_ID"]
columns.core <- intersect(core, colnames(counts_))
columns.extended <- intersect(extended, colnames(counts_))
count.core <- counts_[,columns.core]
count.extended <- counts_[,columns.extended]

core.mean <- rowMeans(count.core)
ext.mean <- rowMeans(count.extended)
length(which(core.mean==0)) 
```

```
## [1] 187
```

```r
length(which(ext.mean==0))
```

```
## [1] 8
```

```r
zero.in.core <- which(core.mean==0 & ext.mean != 0)
plotcore(names(zero.in.core),counts_)
```

![plot of chunk :specific-contigs](figure/:specific-contigs-1.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-2.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-3.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-4.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-5.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-6.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-7.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-8.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-9.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-10.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-11.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-12.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-13.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-14.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-15.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-16.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-17.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-18.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-19.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-20.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-21.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-22.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-23.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-24.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-25.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-26.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-27.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-28.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-29.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-30.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-31.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-32.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-33.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-34.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-35.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-36.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-37.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-38.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-39.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-40.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-41.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-42.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-43.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-44.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-45.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-46.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-47.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-48.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-49.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-50.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-51.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-52.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-53.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-54.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-55.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-56.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-57.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-58.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-59.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-60.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-61.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-62.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-63.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-64.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-65.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-66.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-67.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-68.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-69.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-70.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-71.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-72.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-73.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-74.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-75.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-76.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-77.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-78.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-79.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-80.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-81.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-82.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-83.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-84.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-85.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-86.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-87.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-88.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-89.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-90.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-91.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-92.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-93.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-94.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-95.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-96.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-97.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-98.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-99.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-100.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-101.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-102.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-103.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-104.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-105.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-106.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-107.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-108.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-109.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-110.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-111.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-112.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-113.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-114.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-115.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-116.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-117.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-118.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-119.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-120.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-121.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-122.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-123.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-124.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-125.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-126.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-127.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-128.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-129.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-130.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-131.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-132.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-133.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-134.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-135.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-136.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-137.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-138.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-139.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-140.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-141.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-142.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-143.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-144.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-145.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-146.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-147.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-148.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-149.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-150.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-151.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-152.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-153.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-154.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-155.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-156.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-157.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-158.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-159.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-160.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-161.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-162.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-163.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-164.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-165.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-166.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-167.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-168.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-169.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-170.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-171.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-172.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-173.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-174.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-175.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-176.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-177.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-178.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-179.png) 
![plot of chunk :specific-contigs](figure/:specific-contigs-180.png) 

```r
which(core.mean < 1 & ext.mean > 5)
```

```
## HMEL004016-PA HMEL007289-PA HMEL008728-PA HMEL011034-PA HMEL011306-PA 
##          1188          3143          4029          5466          5625 
## HMEL025050-PB 
##         10040
```

Differential expression between phylogenetic groups
===================================================

Using all samples, with family and tissue as factors.


```r
plotphylo <- function(sig, counts){
for (gene in head(rownames(sig))){
  expr <- t(counts[gene,]/colSums(counts))
  temp <- data.frame(expr=expr, tissue=meta$Tissue, phylo=meta$Phylogeny_group)
  print(bwplot(expr~phylo|tissue, data=temp, main=gene))
}
}
```

Rosids vs. Ranunculales
-----------------------


```r
columns <- intersect(meta$Customer_ID, colnames(counts))
counts_ <- counts[,columns]

dds <- DESeqDataSetFromMatrix(countData = counts_, colData = meta[,c("Customer_ID","Host.Plant.use","Phylogeny_group","Host_Plant","Tissue","Family")], design = ~Tissue+Family+Phylogeny_group)

dds <- DESeq(dds, betaPrior=FALSE, fitType="local")
```

```
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
```

```r
#res <- results(dds, contrast=c("Phylogeny_group","Rosids","Asterids"))
res <- results(dds, contrast=c("Phylogeny_group","Rosids","Ranunculales"))
#res <- results(dds, contrast=c("Phylogeny_group","Ranunculales","Asterids"))
sig <- res[which(res$padj<0.01),]
sig.o <- sig[order(sig$padj),]
head(sig.o)
```

```
## log2 fold change (MLE): Phylogeny_group Rosids vs Ranunculales 
## Wald test p-value: Phylogeny_group Rosids vs Ranunculales 
## DataFrame with 6 rows and 6 columns
##                   baseMean log2FoldChange     lfcSE      stat       pvalue
##                  <numeric>      <numeric> <numeric> <numeric>    <numeric>
## HMEL009626-PA  692.9861747       3.084385 0.1600844 19.267237 1.011956e-82
## HMEL008644-PA   82.0671144       3.694179 0.3184821 11.599331 4.153115e-31
## HMEL007175-PA    3.9051899      18.759953 1.8842936  9.955961 2.375173e-23
## HMEL005123-PA    0.1887582     -34.083737 3.5729178 -9.539469 1.435661e-21
## HMEL002164-PA    0.1455243     -20.455552 2.2332143 -9.159691 5.204625e-20
## HMEL012417-PA 1573.4694606       1.499306 0.1690040  8.871422 7.221744e-19
##                       padj
##                  <numeric>
## HMEL009626-PA 1.001128e-78
## HMEL008644-PA 2.054338e-27
## HMEL007175-PA 7.832530e-20
## HMEL005123-PA 3.550750e-18
## HMEL002164-PA 1.029787e-16
## HMEL012417-PA 1.190745e-15
```

```r
print(dim(sig.o))
```

```
## [1] 655   6
```

```r
plotphylo(sig.o, counts_)
```

![plot of chunk :deseq-phylogroups-1](figure/:deseq-phylogroups-1-1.png) ![plot of chunk :deseq-phylogroups-1](figure/:deseq-phylogroups-1-2.png) ![plot of chunk :deseq-phylogroups-1](figure/:deseq-phylogroups-1-3.png) ![plot of chunk :deseq-phylogroups-1](figure/:deseq-phylogroups-1-4.png) ![plot of chunk :deseq-phylogroups-1](figure/:deseq-phylogroups-1-5.png) ![plot of chunk :deseq-phylogroups-1](figure/:deseq-phylogroups-1-6.png) 

```r
write.table(sig.o,file="sig_rosids_ranunc_alltissues_DESeq2_0.01.txt",quote=F)
```


Rosids vs. Asterids
-----------------------


```r
res <- results(dds, contrast=c("Phylogeny_group","Rosids","Asterids"))
sig <- res[which(res$padj<0.01),]
sig.o <- sig[order(sig$padj),]
head(sig.o)
```

```
## log2 fold change (MLE): Phylogeny_group Rosids vs Asterids 
## Wald test p-value: Phylogeny group Rosids vs Asterids 
## DataFrame with 6 rows and 6 columns
##                  baseMean log2FoldChange      lfcSE      stat       pvalue
##                 <numeric>      <numeric>  <numeric> <numeric>    <numeric>
## HMEL004610-PA  1284.70123      4.1567793 0.49498770  8.397743 4.551436e-17
## HMEL003790-PA   348.23171     -0.5371077 0.07370858 -7.286908 3.171487e-13
## HMEL009525-PA   201.73237     -0.6158440 0.08568747 -7.187095 6.618434e-13
## HMEL011845-PA   142.78271     -1.2443760 0.17631647 -7.057628 1.693693e-12
## HMEL012134-PA 22542.25251      1.0330760 0.14716083  7.020047 2.217931e-12
## HMEL015874-PA    97.43393     -0.9655630 0.13851152 -6.970994 3.147092e-12
##                       padj
##                  <numeric>
## HMEL004610-PA 4.057605e-13
## HMEL003790-PA 1.413690e-09
## HMEL009525-PA 1.966778e-09
## HMEL011845-PA 3.774818e-09
## HMEL012134-PA 3.954571e-09
## HMEL015874-PA 4.676055e-09
```

```r
print(dim(sig.o))
```

```
## [1] 1547    6
```

```r
plotphylo(sig.o, counts_)
```

![plot of chunk :deseq-phylogroups-2](figure/:deseq-phylogroups-2-1.png) ![plot of chunk :deseq-phylogroups-2](figure/:deseq-phylogroups-2-2.png) ![plot of chunk :deseq-phylogroups-2](figure/:deseq-phylogroups-2-3.png) ![plot of chunk :deseq-phylogroups-2](figure/:deseq-phylogroups-2-4.png) ![plot of chunk :deseq-phylogroups-2](figure/:deseq-phylogroups-2-5.png) ![plot of chunk :deseq-phylogroups-2](figure/:deseq-phylogroups-2-6.png) 

```r
write.table(sig.o,file="sig_rosids_asterids_alltissues_DESeq2_0.01.txt",quote=F)
```

Ranunculales vs. Asterids
-----------------------


```r
res <- results(dds, contrast=c("Phylogeny_group","Ranunculales","Asterids"))
sig <- res[which(res$padj<0.01),]
sig.o <- sig[order(sig$padj),]
head(sig.o)
```

```
## log2 fold change (MLE): Phylogeny_group Ranunculales vs Asterids 
## Wald test p-value: Phylogeny group Ranunculales vs Asterids 
## DataFrame with 6 rows and 6 columns
##                  baseMean log2FoldChange     lfcSE      stat        pvalue
##                 <numeric>      <numeric> <numeric> <numeric>     <numeric>
## HMEL009626-PA   692.98617      -3.241173 0.1524731 -21.25735 2.819183e-100
## HMEL008644-PA    82.06711      -3.845777 0.3064722 -12.54853  4.049390e-36
## HMEL012417-PA  1573.46946      -1.771187 0.1600256 -11.06815  1.790558e-28
## HMEL006882-PA  3306.66033      -1.753259 0.1624629 -10.79175  3.765454e-27
## HMEL005487-PA  1693.12783       1.078915 0.1041595  10.35830  3.837185e-25
## HMEL002661-PA 20650.09445      -3.605240 0.3560531 -10.12557  4.255000e-24
##                       padj
##                  <numeric>
## HMEL009626-PA 2.513302e-96
## HMEL008644-PA 1.805015e-32
## HMEL012417-PA 5.320943e-25
## HMEL006882-PA 8.392256e-24
## HMEL005487-PA 6.841700e-22
## HMEL002661-PA 6.322221e-21
```

```r
print(dim(sig.o))
```

```
## [1] 2131    6
```

```r
plotphylo(sig.o, counts_)
```

![plot of chunk :deseq-phylogroups-3](figure/:deseq-phylogroups-3-1.png) ![plot of chunk :deseq-phylogroups-3](figure/:deseq-phylogroups-3-2.png) ![plot of chunk :deseq-phylogroups-3](figure/:deseq-phylogroups-3-3.png) ![plot of chunk :deseq-phylogroups-3](figure/:deseq-phylogroups-3-4.png) ![plot of chunk :deseq-phylogroups-3](figure/:deseq-phylogroups-3-5.png) ![plot of chunk :deseq-phylogroups-3](figure/:deseq-phylogroups-3-6.png) 

```r
write.table(sig.o,file="sig_ranunc_asterids_alltissues_DESeq2_0.01.txt",quote=F)
```

Correlations of PC scores to experimental variables
---------------------------------------------------

```r
p <- prcomp(t(log.cpm.tmm))
pdf("PC_correlations.pdf")
par(mfrow=c(2,2))
for (i in 1:4){
scores <- p$x[,i]
scores.o <- scores[meta$Customer_ID]

pvals <- vector()
pvals <- c(pvals, kruskal.test(scores.o, meta$Tissue)$p.value)
pvals <- c(pvals, kruskal.test(scores.o, meta$Phylogeny_group)$p.value)
pvals <- c(pvals, kruskal.test(scores.o, meta$Host.Plant.use)$p.value)
pvals <- c(pvals, kruskal.test(scores.o, meta$Host_Plant)$p.value)
pvals <- c(pvals, kruskal.test(scores.o, meta$Family)$p.value)
#pvals <- c(pvals, cor.test(scores.o, meta$Million_reads)$p.value)
names(pvals) <- c("Tissue","PhylogenyGroup","HostPlantUse","HostPlant","Family")
barplot(-log(pvals),las=2,ylab="Significance (-log(p))",main=paste("Principal component",i))
abline(2, 0, col="red")
}
dev.off()
```

```
## pdf 
##   2
```

Try SVA to assess importance of different factors
-------------------------------------------------
Overkill?


```r
library(sva)
```

```
## Loading required package: mgcv
```

```
## Warning: package 'mgcv' was built under R version 3.1.2
```

```
## Loading required package: nlme
## 
## Attaching package: 'nlme'
## 
## The following object is masked from 'package:IRanges':
## 
##     collapse
## 
## This is mgcv 1.8-4. For overview type 'help("mgcv-package")'.
## Loading required package: genefilter
## 
## Attaching package: 'genefilter'
## 
## The following object is masked from 'package:base':
## 
##     anyNA
```

```r
sampleinfo <- meta[,c("Host.Plant.use","Host_Plant", "Tissue", "Family", "Phylogeny_group")]
rownames(sampleinfo) <- meta$Customer_ID
data <- log.cpm.tmm[,rownames(sampleinfo)]
mod <- model.matrix(~as.factor(Tissue), data=sampleinfo)
mod0 <- model.matrix(~1, data=sampleinfo)
n.sv <- num.sv(data,mod,method="leek")
svobj <- sva(as.matrix(data),mod,mod0,n.sv=n.sv)
```

```
## Number of significant surrogate variables is:  3 
## Iteration (out of 5 ):1  2  3  4  5
```

```r
surr <- svobj$sv
```

What is in the first two surrogate variable?

```r
pdf("SVA.pdf")
par(mfrow=c(2,1))
for (i in 1:2){
pvals <- vector()
pvals <- c(pvals, cor.test(surr[,i], as.numeric(sampleinfo$Host.Plant.use))$p.value)
pvals <- c(pvals, cor.test(surr[,i], as.numeric(sampleinfo$Host_Plant))$p.value)
pvals <- c(pvals, cor.test(surr[,i], as.numeric(sampleinfo$Phylogeny_group))$p.value)
pvals <- c(pvals, cor.test(surr[,i], as.numeric(sampleinfo$Family))$p.value)
names(pvals) <- c("HostPlantUse","HostPlant","PhylogenyGroup","Family")
barplot(-log(pvals),las=2,ylab="Significance (-log(p))",main=paste("PC", i))
abline(2,0,col="red")
}
dev.off()
```

```
## pdf 
##   2
```
