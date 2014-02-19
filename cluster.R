library(pheatmap)

counts <- read.delim("read_counts_with_proper_headers.csv",sep="\t",row.names=1)

meta  <- read.csv("rna-seqecoevowabi_relational_table.csv")
# meta <- meta[-which(meta$Customer_ID=="B1MT"),]

corrs <- cor(counts)

display.names <- paste(meta$Host_Plant,meta$Tissue,meta$Host.Plant.use,sep=".")
names(display.names) <- meta$Customer_ID

colnames(corrs) <- display.names[colnames(corrs)]
rownames(corrs) <- colnames(corrs)

pheatmap(corrs)