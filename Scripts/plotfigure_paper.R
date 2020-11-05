# Load libraries
library(edgeR)
library(DESeq2)
library(stringr)
library(ggplot2)
library(matrixStats)
library(pheatmap)
library(edgeR)
## 1. Read the output from featurecounts. 
counts=read.csv("Zf_hypoxia_raw_counts.txt", sep="", head=T, skip=1, row.names = "Geneid")
## Create an object for storing names of samples and condition of the samples ##(Male Female as a first look) 
colnames(counts)[6:dim(counts)[2]]
#assign names to samples
colnames(counts)[6:dim(counts)[2]]=  c("C1","C2","C3","H1","H2","H3")
samples=cbind(colnames(counts)[6:dim(counts)[2]],rep(c("control","hypoxia"),each=3))
rownames(samples)=samples[,1]
samples=as.data.frame(samples[,-1])
colnames(samples)="treatment"
dds=DESeqDataSetFromMatrix(countData = counts[,6:dim(counts)[2]],colData = samples,design = ~ treatment)
keep <- rowSums(cpm(counts(dds)))     >= 6 # cpm average 1 per sample
dds <- dds[keep,]
dds <- DESeq(dds)

res <- results(dds, cooksCutoff=TRUE, independentFiltering=F,contrast=c("treatment","hypoxia","control"))
DE_genes<-res[which(res$padj<=0.05),]
# remove the one where there are some individuals  outliers
dim(DE_genes)
mat<-assay(rld)[which(rownames(assay(rld))%in%rownames(DE_genes)),]
mat <- mat - rowMeans(mat) 
apply(mat,1,scale)
colnames(mat)<-c( "C1", "C2", "C3","H1",  "H2","H3") #colnames(mat)<-c( "Control (s1)", "Control (s3)", "Control (s4)","Hypoxia (s6)",  "Hypoxia (s9)","Hypoxia (s10)")
pdf("heatmap.pdf")
pheatmap(t(scale(t(mat))),show_rownames = F,cellwidth=20,cellheight=2,border_color="lightgrey",cluster_row=T,cluster_col=F,fontsize=10,las=2) # easy to visualize rownames 
dev.off()
