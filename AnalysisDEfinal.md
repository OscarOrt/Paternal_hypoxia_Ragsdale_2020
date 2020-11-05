AnalysisDEfinal
================
Ludovic Dutoit
8/20/2020

Part 1: Samples check
---------------------

``` r
library(edgeR)
```

    ## Warning: package 'edgeR' was built under R version 3.5.2

    ## Loading required package: limma

    ## Warning: package 'limma' was built under R version 3.5.1

``` r
library(DESeq2)
```

    ## Warning: package 'DESeq2' was built under R version 3.5.2

    ## Loading required package: S4Vectors

    ## Warning: package 'S4Vectors' was built under R version 3.5.1

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Warning: package 'BiocGenerics' was built under R version 3.5.1

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following object is masked from 'package:limma':
    ## 
    ##     plotMA

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind,
    ##     colMeans, colnames, colSums, dirname, do.call, duplicated,
    ##     eval, evalq, Filter, Find, get, grep, grepl, intersect,
    ##     is.unsorted, lapply, lengths, Map, mapply, match, mget, order,
    ##     paste, pmax, pmax.int, pmin, pmin.int, Position, rank, rbind,
    ##     Reduce, rowMeans, rownames, rowSums, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which, which.max,
    ##     which.min

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## Warning: package 'IRanges' was built under R version 3.5.1

    ## Loading required package: GenomicRanges

    ## Warning: package 'GenomicRanges' was built under R version 3.5.1

    ## Loading required package: GenomeInfoDb

    ## Warning: package 'GenomeInfoDb' was built under R version 3.5.2

    ## Loading required package: SummarizedExperiment

    ## Warning: package 'SummarizedExperiment' was built under R version 3.5.1

    ## Loading required package: Biobase

    ## Warning: package 'Biobase' was built under R version 3.5.1

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: DelayedArray

    ## Warning: package 'DelayedArray' was built under R version 3.5.1

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## Loading required package: BiocParallel

    ## Warning: package 'BiocParallel' was built under R version 3.5.2

    ## 
    ## Attaching package: 'DelayedArray'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

    ## The following objects are masked from 'package:base':
    ## 
    ##     aperm, apply

``` r
library(stringr)
```

    ## Warning: package 'stringr' was built under R version 3.5.2

``` r
library(ggplot2)
```

    ## Warning: package 'ggplot2' was built under R version 3.5.2

``` r
library(matrixStats)
library(pheatmap)
```

    ## Warning: package 'pheatmap' was built under R version 3.5.2

``` r
library(edgeR)
```

``` r
## 1. Read the output from featurecounts. 
counts=read.csv("../Data/Zf_hypoxia_raw_counts.txt", sep="", head=T, skip=1, row.names = "Geneid")
## Create an object for storing names of samples and condition of the samples ##(Male Female as a first look) 

colnames(counts)[6:dim(counts)[2]]
```

    ## [1] "X05_bam.C1_5_sorted.bam" "X05_bam.C2_6_sorted.bam"
    ## [3] "X05_bam.C4_7_sorted.bam" "X05_bam.H1_5_sorted.bam"
    ## [5] "X05_bam.H3_8_sorted.bam" "X05_bam.H5_7_sorted.bam"

``` r
#assign names to samples
colnames(counts)[6:dim(counts)[2]]=  c("C1","C2","C3","H1","H2","H3")
samples=cbind(colnames(counts)[6:dim(counts)[2]],rep(c("control","hypoxia"),each=3))
rownames(samples)=samples[,1]
samples=as.data.frame(samples[,-1])
colnames(samples)="treatment"
samples
```

    ##    treatment
    ## C1   control
    ## C2   control
    ## C3   control
    ## H1   hypoxia
    ## H2   hypoxia
    ## H3   hypoxia

``` r
#

# esimtate new counts and variance
dds=DESeqDataSetFromMatrix(countData = counts[,6:dim(counts)[2]],colData = samples,design = ~ treatment)
keep <- rowSums(cpm(counts(dds)))     >= 6 # cpm average 1 per sample
dds <- dds[keep,]
#Look at the PCA
plotPCA(rlog(dds,blind=T), intgroup="treatment")+theme_bw()
```

![](AnalysisDEfinal_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
pcaData <- plotPCA(rlog(dds,blind=T), intgroup=c("treatment"), returnData=TRUE)
plot(pcaData$PC1, pcaData$PC2,col="white")
text(pcaData$PC1, pcaData$PC2,colnames(dds),cex=1)
```

![](AnalysisDEfinal_files/figure-markdown_github/unnamed-chunk-6-1.png)

Samples are separated but they need both axes for that, it is not optimal, but it is not terrible :)

Differential expression analysis
--------------------------------

I now look at identifying differentially expressed genes

``` r
dds <- DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

I first make sure the dispersion fits decently:

``` r
## Visualize the dispersion
plotDispEsts( dds, ylim = c(1e-6, 1e1) )
```

![](AnalysisDEfinal_files/figure-markdown_github/unnamed-chunk-8-1.png)

This plot is very important as it allows us to see that the dispersion fits the data.

in blue are dispersion values for the newly estimated counts, in red is the final dispersion estimate, The circled genes have "extreme" variance and will not be fitted the red dispersion line.

``` r
res <- results(dds, cooksCutoff=TRUE, independentFiltering=F,contrast=c("treatment","hypoxia","control"))
DE_genes<-res[which(res$padj<=0.05),]
# remove the one where there are some individuals  outliers
dim(DE_genes)
```

    ## [1] 91  6

More visualisation
------------------

Now that we have the DE genes, we will do a bit more visualisation

``` r
#visualization of the DE genes this is good too.
rld <- rlog(dds, blind=F)###regularized logarithm or rlog normalisation which is suited for data visualisation ( unlike the one we used for actually used for testing) It incorporates a prior on the sample differences (Love, Huber, and Anders 2014). , data on the log2 scale which has been normalized with respect to library size or other normalization factors. Difference to the average expression of every differentially expressed gene. the counts are transformed using regularized logarithm (Love, Huber, and Anders 2014) .


##Heatmap of significant genes

mat<-assay(rld)[which(rownames(assay(rld))%in%rownames(DE_genes)),]
mat <- mat - rowMeans(mat) 

colnames(mat)<-c( "Control (s1)", "Control (s3)", "Control (s4)","Hypoxia (s6)",  "Hypoxia (s9)","Hypoxia (s10)")



pheatmap(mat,show_rownames = T,cellwidth=20,cellheight=5,border_color="lightgrey",cluster_row=T,cluster_col=F,fontsize=6) # easy to visualize rownames 
```

![](AnalysisDEfinal_files/figure-markdown_github/unnamed-chunk-10-1.png)

That is a heatmap of all the DE genes.

We then do an MA plotand a volcano plot which show us logfold change across expression highlighting significant DE genes. That helps us to see what kind of power we have to detect DE genes.

``` r
#MA plot

plotMA(res)
```

![](AnalysisDEfinal_files/figure-markdown_github/unnamed-chunk-11-1.png)

``` r
plot(res$log2FoldChange, -log10(res$padj),pch=19,cex=0.5,frame=F,cex.lab=1.7,cex.axis=1.7,xlim=c(-10,10),xlab=)
points(DE_genes$log2FoldChange, -log10(DE_genes$padj),pch=19,cex=0.5,col="red")
```

![](AnalysisDEfinal_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
pdf("Figure_volcano.pdf")
plot(res$log2FoldChange, -log10(res$pvalue),pch=19,cex=0.5,frame=F,cex.lab=1.7,cex.axis=1.7,xlim=c(-10,10),col="black")
points(DE_genes$log2FoldChange, -log10(DE_genes$pvalue),pch=19,cex=0.5,col="red")
dev.off()
```

    ## quartz_off_screen 
    ##                 2

We can now which gene is what using the script addgenes.py that use an annotation database.

### GO analysis

In order to run a GO analysis downstream, I prepare to file, a "target" file with overexpressed genes and a "background"" file with all genes but the overexpressed ones

``` r
write(row.names(DE_genes),file="target.txt",sep="\n")
'%!in%' <- function(x,y)!('%in%'(x,y))
background<-levels(as.factor(row.names(dds)))[levels(as.factor(row.names(dds)))%!in%row.names(DE_genes)]
write(background,file="background.txt",sep="\n")
```

One can just put those two lists in [Gorilla online](http://cbl-gorilla.cs.technion.ac.il/)

The results have been downloaded and saved in excel files GOFUNCTION.xls GOFUNCTION.xls and GOPROCESS.xlsGOPROCESS.xls.

### Output results

``` r
write.table(DE_genes,"de_genes.txt",sep="\t",quote=F)
write.table(res,"table_results_allgenes.txt",sep="\t",quote=F)
```

### Compare to rashid et al. (hypoxia in fish database hrgfish )

``` r
rashid_genes=read.table("dataRashidetahrgfish.txt",h=F)
length(which(rashid_genes[,1]%in%rownames(res))) # we got 46 of them
```

    ## [1] 46

``` r
ordered_rashid_genes<-rashid_genes[order(rashid_genes[,1]),]
ordered_res<-res[order(res[,1]),]
write.table(ordered_res[which(rownames(ordered_res)%in%as.character(rashid_genes[,1])),],"Table_S2.txt",sep="\t")
ordered_res[which(rownames(ordered_res)%in%as.character(rashid_genes[,1])),]
```

    ## log2 fold change (MLE): treatment hypoxia vs control 
    ## Wald test p-value: treatment hypoxia vs control 
    ## DataFrame with 46 rows and 6 columns
    ##                            baseMean      log2FoldChange             lfcSE
    ##                           <numeric>           <numeric>         <numeric>
    ## ENSDARG00000054323 18.6412319344216 -0.0656098291985259 0.469084926987221
    ## ENSDARG00000061047 21.7682035585696   0.213707277461826 0.483793820536981
    ## ENSDARG00000078452 32.7505361713857   0.397918088833641 0.390578214928286
    ## ENSDARG00000006181 34.2697833555559  -0.223164200597423 0.369478094094763
    ## ENSDARG00000032553 40.8126491902295    0.12575189682411 0.330925492840637
    ## ...                             ...                 ...               ...
    ## ENSDARG00000019644 6156.12854233729  -0.209598581493719 0.165483214938788
    ## ENSDARG00000011665 6869.24216488502  -0.376156441343802  0.26051902959518
    ## ENSDARG00000022456 9822.34792361457  -0.122644523699183 0.155486466352072
    ## ENSDARG00000015551 14927.7750025268  -0.346535018839536 0.130869972243247
    ## ENSDARG00000016771  64945.401621324  -0.394968635409977 0.241904208214175
    ##                                  stat              pvalue
    ##                             <numeric>           <numeric>
    ## ENSDARG00000054323 -0.139867698627445   0.888764523131052
    ## ENSDARG00000061047  0.441732135446923   0.658683052558703
    ## ENSDARG00000078452   1.01879232795076   0.308301568779146
    ## ENSDARG00000006181 -0.603998462058178   0.545844664667007
    ## ENSDARG00000032553  0.380000633207994   0.703944945116524
    ## ...                               ...                 ...
    ## ENSDARG00000019644  -1.26658514321981   0.205303667479211
    ## ENSDARG00000011665  -1.44387318626325   0.148774650390903
    ## ENSDARG00000022456 -0.788779413260805   0.430240943821706
    ## ENSDARG00000015551   -2.6479337689125 0.00809853896403004
    ## ENSDARG00000016771  -1.63274809614012   0.102521987416485
    ##                                 padj
    ##                            <numeric>
    ## ENSDARG00000054323 0.993149737464813
    ## ENSDARG00000061047 0.971003972759044
    ## ENSDARG00000078452 0.896093930788175
    ## ENSDARG00000006181 0.953033939899338
    ## ENSDARG00000032553 0.978818427145694
    ## ...                              ...
    ## ENSDARG00000019644 0.851996539207929
    ## ENSDARG00000011665 0.804973519245826
    ## ENSDARG00000022456 0.933716465176177
    ## ENSDARG00000015551 0.307147636446072
    ## ENSDARG00000016771 0.741525361286179

### Grab the hemoglobin genes frpm the count for the supplementary table

    ENSDARG00000079078 hbz
    ENSDARG00000097011 hemoglobin, alpha adult 1
    ENSDARG00000097238 hbaa1
    ENSDARG00000089087 hbba1

``` r
hb_de_genes<-c("ENSDARG00000079078","ENSDARG00000097011","ENSDARG00000097238","ENSDARG00000089087")
hb_genes<-counts[which(rownames(counts)%in%hb_de_genes),6:11]
hb_genes
```

    ##                    C1 C2 C3  H1   H2 H3
    ## ENSDARG00000097238  0  0  2   9  303  1
    ## ENSDARG00000097011  2  1 17 671 1427 14
    ## ENSDARG00000089087  0  0  8 341  974  8
    ## ENSDARG00000079078  1  3 13 227 1169 21

``` r
rownames(hb_genes)<-c("hbz","hemoglobin, alpha adult 1","hbaa1","hbba1")
colnames(hb_genes)<-sub("_sorted.bam","",colnames(hb_genes))
colnames(hb_genes)<-sub("X05_bam.","",colnames(hb_genes))
write.table(hb_genes,"hb_table.txt",sep="\t",quote=F)
```

the heatmaps and the volcano plots are generated in [plotfigurepaper.R](plotfigurepaper.R)
