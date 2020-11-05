# Load libraries
library(GenomicRanges)
library(Gviz)
library(RColorBrewer)
library(ggplot2)

# Change working directory
setwd("C:/Users/OscarJavier/Documents/GitHub/Paternal_hypoxia_Ragsdale_2020/Data")

# Figure 5 panel A

# RNA Hb
data_Hb_RNA <- read.csv("RNA_Hb.txt",sep= "\t")

df_Hb_RNA <- data.frame("chr"=data_Hb_RNA$Chromosome,"start"=data_Hb_RNA$Start,
                        "end"=data_Hb_RNA$End,"strand"=rep("*",nrow(data_Hb_RNA)))


Gr_Hb_RNA_hypoxia <- makeGRangesFromDataFrame(df_Hb_RNA)  
values(Gr_Hb_RNA_hypoxia) <- DataFrame(Hypoxia = data_Hb_RNA$Hypoxia_RNA_seq)

Gr_Hb_RNA_normoxia <- makeGRangesFromDataFrame(df_Hb_RNA)  
values(Gr_Hb_RNA_normoxia) <- DataFrame(Normoxia = data_Hb_RNA$Control_RNA_seq)

Hb_RNA_track_hypoxia <- AnnotationTrack(Gr_Hb_RNA_hypoxia, name = "Hypoxia")
Hb_RNA_track_normoxia <- AnnotationTrack(Gr_Hb_RNA_normoxia, name = "Normoxia")

plotTracks(list(DataTrack(Gr_Hb_RNA_hypoxia, name ="Hypoxia",ylim= c(0,20)),
                DataTrack(Gr_Hb_RNA_normoxia, name ="Normoxia",ylim = c(0,20))),
           type = c("g","histogram"), background.title = "darkblue",
           from = 55097000, to = 55106999)

# 800 x 400

# RNA Acot21
data_Acot21_RNA <- read.csv("RNA_Acot21.txt",sep= "\t")

df_Acot21_RNA <- data.frame("chr"=data_Acot21_RNA$Chromosome,"start"=data_Acot21_RNA$Start,
                        "end"=data_Acot21_RNA$End,"strand"=rep("*",nrow(data_Acot21_RNA)))

Gr_Acot21_RNA_hypoxia <- makeGRangesFromDataFrame(df_Acot21_RNA)  
values(Gr_Acot21_RNA_hypoxia) <- DataFrame(hypoxia = data_Acot21_RNA$Hypoxia_RNA_seq)

Gr_Acot21_RNA_normoxia <- makeGRangesFromDataFrame(df_Acot21_RNA)  
values(Gr_Acot21_RNA_normoxia) <- DataFrame(male = data_Acot21_RNA$Control_RNA_seq)

Hb_Acot21_track_hypoxia <- AnnotationTrack(Gr_Acot21_RNA_hypoxia, name = "Hypoxia")
Hb_Acot21_track_normoxia <- AnnotationTrack(Gr_Acot21_RNA_normoxia, name = "Normoxia")

plotTracks(list(DataTrack(Gr_Acot21_RNA_hypoxia, name ="Hypoxia",ylim= c(0,1)),
                DataTrack(Gr_Acot21_RNA_normoxia, name ="Normoxia",ylim = c(0,1))),
           type = c("g","histogram"), background.title = "darkblue",
           from = 28941000, to = 28959999)

# 800 x 400

# CpG islands and gene model Hbs
geneModel <- read.csv("Gene_Hbs.txt",sep = "\t")
grtrack_Hb <- GeneRegionTrack(geneModel,name = "Hbs")

Hbs_scale <- data.frame("chr"=3,"start"=55097000,"end"=55106999,"strand"="*")
Gr_hbs_scale <- makeGRangesFromDataFrame(Hbs_scale)

Hbs_CpG <- data.frame("chr"=c(rep(3,7)),"start"=c(55097041,55098477,55098854,55100372,55103466,55103843,55105396),
                          "end"=c(55097426,55098760,55099054,55100675,55103748,55104042,55105699),
                          "strand"="*")
Gr_Hbs_CpG <- makeGRangesFromDataFrame(Hbs_CpG)

Hbs_marks <- data.frame("chr"=c(3,3,3,3,3),"start"=c(55098000,55100000,55102000,55104000,55106000),
                      "end"=c(55098001,55100001,55102001,55104001,55106001),
                      "strand"="*")
Gr_Hbs_marks <- makeGRangesFromDataFrame(Hbs_marks)

data_Hb_bioCAP <- read.csv("hem_probes_BioCAP.txt",sep= "\t")

Gr_Hb_bioCAP_testes <- makeGRangesFromDataFrame(df_Hb_RNA)  
values(Gr_Hb_bioCAP_testes) <- DataFrame(Testes = data_Hb_bioCAP$SRR648831_ordered.bam)

Gr_Hb_bioCAP_liver <- makeGRangesFromDataFrame(df_Hb_RNA)  
values(Gr_Hb_bioCAP_liver) <- DataFrame(Liver = data_Hb_bioCAP$SRR648833_ordered.bam)

Gr_Hb_bioCAP_24hpf <- makeGRangesFromDataFrame(df_Hb_RNA)  
values(Gr_Hb_bioCAP_24hpf) <- DataFrame(x24hpf = data_Hb_bioCAP$SRR648835_ordered.bam)

plotTracks(list(grtrack_Hb,AnnotationTrack(Gr_Hbs_CpG, name ="CpG"),AnnotationTrack(Gr_Hbs_marks,name = "Marks"),
                DataTrack(Gr_Hb_bioCAP_testes, name="Testes",ylim= c(0,100)),
                DataTrack(Gr_Hb_bioCAP_liver, name= "Liver",ylim= c(0,100)),
                DataTrack(Gr_Hb_bioCAP_24hpf, name="24hpf",ylim= c(0,100))),
           type = c("g","histogram"), background.title = "darkblue",
           from = 55097000, to = 55106999) 

# 800 x 150

# CpG islands and gene model Acot21
geneModel <- read.csv("Gene_Acot21.txt",sep = "\t")
grtrack_Acot21 <- GeneRegionTrack(geneModel,name = "Acot21")

Acot21_scale <- data.frame("chr"=20,"start"=28941000,"end"=28959999,"strand"="*")
Gr_Acot21_scale <- makeGRangesFromDataFrame(Acot21_scale)

Acot21_CpG <- data.frame("chr"=c(rep(20,5)),"start"=c(28943620,28944006,28954750,28955945,28960210),
                      "end"=c(28943920,28944287,28954976,28956236,28960678),
                      "strand"="*")
Gr_acot21_CpG <- makeGRangesFromDataFrame(Acot21_CpG)


Acot_21_marks <- data.frame("chr"=c(20,20,20,20),"start"=c(28944000,28948000,28952000,28956000),
                        "end"=c(28944001,28948001,28952001,28956001),
                        "strand"="*")
Gr_Acot21_marks <- makeGRangesFromDataFrame(Acot_21_marks)

data_Acot21_bioCAP <- read.csv("Acot21_probes_BioCAP.txt",sep= "\t")

Gr_Acot21_bioCAP_testes <- makeGRangesFromDataFrame(df_Acot21_RNA)  
values(Gr_Acot21_bioCAP_testes) <- DataFrame(Testes = data_Acot21_bioCAP$SRR648831_ordered.bam)

Gr_Acot21_bioCAP_liver <- makeGRangesFromDataFrame(df_Acot21_RNA)  
values(Gr_Acot21_bioCAP_liver) <- DataFrame(Liver = data_Acot21_bioCAP$SRR648833_ordered.bam)

Gr_Acot21_bioCAP_24hpf <- makeGRangesFromDataFrame(df_Acot21_RNA)  
values(Gr_Acot21_bioCAP_24hpf) <- DataFrame(x24hpf = data_Acot21_bioCAP$SRR648835_ordered.bam)

plotTracks(list(grtrack_Acot21,AnnotationTrack(Gr_acot21_CpG, name ="CpG"),AnnotationTrack(Gr_Acot21_marks,name = "Marks"),
                DataTrack(Gr_Acot21_bioCAP_testes, name="Testes",ylim= c(0,100)),
                DataTrack(Gr_Acot21_bioCAP_liver, name="Liver",ylim= c(0,100)),
                DataTrack(Gr_Acot21_bioCAP_24hpf, name="24hpf",ylim= c(0,100))),
           type = c("g","histogram"), background.title = "darkblue",
           from = 28941000, to = 28959999) 

# Methylation Hb
data_Hb_meth <- read.csv("Met_Hb.txt",sep= "\t")
df_Hb_meth <- data.frame("chr"=data_Hb_meth$Chromosome,"start"=data_Hb_meth$Start,
                            "end"=data_Hb_meth$End,"strand"=rep("*",nrow(data_Hb_meth)))

Gr_Hb_meth_hypoxia <- makeGRangesFromDataFrame(df_Hb_meth)
values(Gr_Hb_meth_hypoxia) <- DataFrame(Hypoxia = data_Hb_meth$Hypoxia_met)

Gr_Hb_meth_normoxia <- makeGRangesFromDataFrame(df_Hb_meth)
values(Gr_Hb_meth_normoxia) <- DataFrame(male = data_Hb_meth$Control_met)

Hb_meth_track_hypoxia <- AnnotationTrack(Gr_Hb_meth_hypoxia, name = "Hypoxia")
Hb_meth_track_normoxia <- AnnotationTrack(Gr_Hb_meth_normoxia, name = "Normoxia")

plotTracks(list(DataTrack(Gr_Hb_meth_hypoxia, name ="Hypoxia",ylim= c(0,100),col = "black"),
                DataTrack(Gr_Hb_meth_normoxia, name ="Normoxia",ylim = c(0,100),col = "black")),
           type = c("p","histogram","g"), background.title = "darkblue",
           from = 55097000, to = 55106999)

# 800 x 400

# Methylation Acot21
data_Acot21_meth <- read.csv("Met_Acot21.txt",sep= "\t")
df_Acot21_meth <- data.frame("chr"=data_Acot21_meth$Chromosome,"start"=data_Acot21_meth$Start,
                         "end"=data_Acot21_meth$End,"strand"=rep("*",nrow(data_Acot21_meth)))

Gr_Acot21_meth_hypoxia <- makeGRangesFromDataFrame(df_Acot21_meth)
values(Gr_Acot21_meth_hypoxia) <- DataFrame(Hypoxia = data_Acot21_meth$Hypoxia_met)

Gr_Acot21_meth_normoxia <- makeGRangesFromDataFrame(df_Acot21_meth)
values(Gr_Acot21_meth_normoxia) <- DataFrame(male = data_Acot21_meth$Control_met)

Acot21_meth_track_hypoxia <- AnnotationTrack(Gr_Acot21_meth_hypoxia, name = "Hypoxia")
Acot21_meth_track_normoxia <- AnnotationTrack(Gr_Acot21_meth_normoxia, name = "Normoxia")

plotTracks(list(DataTrack(Gr_Acot21_meth_hypoxia, name ="Hypoxia",ylim= c(0,100),col = "black"),
                DataTrack(Gr_Acot21_meth_normoxia, name ="Normoxia",ylim = c(0,100),col = "black")),
           type = c("p","histogram","g"), background.title = "darkblue",
           from = 28941000, to = 28959999)


# 800 x 400

# Plot RNA-seq, gene models and CpG predicted regions Hb

plotTracks(list(DataTrack(Gr_Hb_RNA_hypoxia, name ="Hypoxia",ylim= c(0,20)),
                DataTrack(Gr_Hb_RNA_normoxia, name ="Normoxia",ylim = c(0,20)),
                grtrack_Hb,
                AnnotationTrack(Gr_Hbs_CpG, name ="CpG"),
                AnnotationTrack(Gr_Hbs_marks,name = "Marks")),
           type = c("g","histogram"), background.title = "darkblue",
           from = 55097000, to = 55106999)

# 800 x 500 

# Plot CpG BioCap Hb

plotTracks(list(DataTrack(Gr_Hb_bioCAP_testes, name="Testes",ylim= c(0,100))),
           type = c("g","histogram"), background.title = "darkblue",
           from = 55097000, to = 55106999)

# 800 x 100

# Plot methylation Hb

plotTracks(list(DataTrack(Gr_Hb_meth_hypoxia, name ="Hypoxia",ylim= c(0,100),col = "black"),
           DataTrack(Gr_Hb_meth_normoxia, name ="Normoxia",ylim = c(0,100),col = "black")),
           type = c("p","histogram","g"), background.title = "darkblue",
           from = 55097000, to = 55106999)

# 800 x 400

# Plot RNA-seq, gene models and CpG predicted regions Acot21

plotTracks(list(DataTrack(Gr_Acot21_RNA_hypoxia, name ="Hypoxia",ylim= c(0,1)),
                DataTrack(Gr_Acot21_RNA_normoxia, name ="Normoxia",ylim = c(0,1)),
                grtrack_Acot21,
                AnnotationTrack(Gr_acot21_CpG, name ="CpG"),
                AnnotationTrack(Gr_Acot21_marks,name = "Marks")),
           type = c("g","histogram"), background.title = "darkblue",
           from = 28941000, to = 28959999)

# 800 x 500

# Plot CpG BioCap Acot21

plotTracks(list(DataTrack(Gr_Acot21_bioCAP_testes, name="Testes",ylim= c(0,100))),
           type = c("g","histogram"), background.title = "darkblue",
           from = 28941000, to = 28959999)

# 800 x 100

# Plot methylation Acot21

plotTracks(list(DataTrack(Gr_Acot21_meth_hypoxia, name ="Hypoxia",ylim= c(0,100),col = "black"),
                DataTrack(Gr_Acot21_meth_normoxia, name ="Normoxia",ylim = c(0,100),col = "black")),
           type = c("p","histogram","g"), background.title = "darkblue",
           from = 28941000, to = 28959999)

# 800 x 400

# Figure 5 Panel B

# Specific coupling
data = read.csv("Coupling_DEG.txt",sep="\t")

data_clean = data[data$Control_calls>20&data$Hypoxia_calls>20,]
data_clean$expression = ifelse(data_clean$log2FoldChange>0,"Up","Down")

table(data_clean$expression)

data_clean$expression <- factor(data_clean$expression)

colfunc_data <- colorRampPalette(c("#45b6fe", "#ed2939"))
col_data <- colfunc_data(2)[1:2]

plot_control <- ggplot(data_clean, aes(x=expression, y=Control_met,fill=expression)) +
  geom_violin(scale = "area")+
  stat_summary(fun.y=mean, geom="point", size=4, color="black") +
  theme(axis.text = element_text(size=20)) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_blank())+
  scale_fill_manual(values=col_data) + 
  guides(fill=FALSE)

plot_hypoxia <- ggplot(data_clean, aes(x=expression, y=Hypoxia_met,fill=expression)) +
  geom_violin(scale = "area")+
  stat_summary(fun.y=mean, geom="point", size=4, color="black") +
  theme(axis.text = element_text(size=20)) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_blank())+
  scale_fill_manual(values=col_data) + 
  guides(fill=FALSE)


plot_control
plot_hypoxia

# 800 x 600

# Figure 5 Panel C
data = read.csv("Coupling_all_genes.csv",sep = ",",dec=".",stringsAsFactors = FALSE,na.strings = "#N/A")

data_clean = data[data$Control_calls>20&data$Hypoxia_calls>20,]
data_clean$Control_RNA_seq <- as.numeric(data_clean$Control_RNA_seq)
data_clean$Hypoxia_RNA_seq <- as.numeric(data_clean$Hypoxia_RNA_seq)

# Hypoxia
data_hypoxia = data_clean[data_clean$Hypoxia_RNA_seq>0&!is.na(data_clean$Hypoxia_met)&!is.na(data_clean$Hypoxia_RNA_seq),]

quantile_hypoxia <- quantile(data_hypoxia$Hypoxia_RNA_seq, seq(0, 1, 0.25))
data_hypoxia$quantile <- findInterval(data_hypoxia$Hypoxia_RNA_seq, quantile_hypoxia, all.inside = TRUE)

table(data_hypoxia$quantile)

data_hypoxia$quantile <- factor(data_hypoxia$quantile)

colfunc_hypoxia <- colorRampPalette(c("white", "#00048b"))
col_hypoxia <- colfunc_hypoxia(8)[2:5]

plot_hypoxia <- ggplot(data_hypoxia, aes(x=quantile, y=Hypoxia_met,fill=quantile)) +
  geom_violin(scale = "width")+
  stat_summary(fun.y=median, geom="point", size=4, color="black") +
  theme(axis.text = element_text(size=20)) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_blank())+
  scale_fill_manual(values=col_hypoxia) + 
  guides(fill=FALSE)

plot_hypoxia

# Control
data_control = data_clean[data_clean$Control_RNA_seq>0&!is.na(data_clean$Control_met)&!is.na(data_clean$Control_RNA_seq),]

quantile_control <- quantile(data_control$Control_RNA_seq, seq(0, 1, 0.25))
data_control$quantile <- findInterval(data_control$Control_RNA_seq, quantile_control, all.inside = TRUE)

table(data_hypoxia$quantile)
table(data_control$quantile)

data_control$quantile <- factor(data_control$quantile)

colfunc_control <- colorRampPalette(c("white", "#8f0000"))
col_control <- colfunc_control(8)[2:5]

plot_control <- ggplot(data_control, aes(x=quantile, y=Control_met,fill=quantile)) +
  geom_violin(scale = "width")+
  stat_summary(fun.y=median, geom="point", size=4, color="black") +
  theme(axis.text = element_text(size=20)) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_blank())+
  scale_fill_manual(values=col_control) + 
  guides(fill=FALSE)

plot_control
