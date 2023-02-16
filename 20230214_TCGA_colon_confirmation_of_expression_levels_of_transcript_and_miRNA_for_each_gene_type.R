# This script is to classfy plot about correlation between expression level of transcript and miRNA 
# 2023/02/14 made

# activate packages
library(stringr)

# make new directory
setwd("C:/Rdata")
dir.create("20230214_TCGA_colon_confirmation_of_expression_levels_of_transcript_and_miRNA_for_each_gene_type")

# import gencode v36 annotation
# this annotation is located at "//fsw-q02/okamura-lab/Files_related_to_M1_Projects/Hirota/annotations"
setwd("C:/Rdata/annotations")
gen36.anot <-read.table("gencode.v36.annotation.gtf",sep="\t",header = F,stringsAsFactors = F)

# extract only genes and make table about gene name and gene type
gen36.anot <-gen36.anot[gen36.anot[,3]=="gene",]
anot <-str_split(gen36.anot[,9],pattern = "; ",simplify = T)
anot <-anot[,c(3,2)]
anot[,1] <-gsub("gene_name ","",anot[,1])
anot[,2] <-gsub("gene_type ","",anot[,2])
colnames(anot) <-c("gene_name","gene_type")

# import list of transcripts that intersect with miRNAs in gencode v36
# this list is located at "https://github.com/Ryosuke-Hirota/20230110_TCGA_colon_transcriptome_bam_correlation_between_transcript_and_miRNA"
setwd("C:/Rdata")
primir.list <-read.table("TCGA_hg38_transcript_intersect_with_miRNA.txt",sep="\t",header = F,stringsAsFactors = F)

# edit list and remove unnecessary rows and columns
primir.list[,2] <-primir.list[,2]-1
primir.list <-primir.list[primir.list[,2]!=primir.list[,8]&primir.list[,3]!=primir.list[,9],]
primir.list <-primir.list[,c(4,11)]
colnames(primir.list) <-c("miRNA","gene_name")
primir.list <-subset(primir.list,!duplicated(primir.list))

# merge list of transcript that intersect with miRNA and table about gene name and gene type
primir.anot <-merge(primir.list,anot,by=,"gene_name")
primir.anot <-primir.anot[,c(2,1,3)]
colnames(primir.anot)[2] <-"transcript"
primir.anot <-subset(primir.anot,!duplicated(primir.anot))

# import result of correlation analysis between expression level of miRNA and transcript
# this result is lcoated at "https://github.com/Ryosuke-Hirota/20230110_TCGA_colon_transcriptome_bam_correlation_between_transcript_and_miRNA"
setwd("C:/Rdata/20230110_TCGA_colon_transcriptome_bam_correlation_between_transcript_and_miRNA")
cor.result <-read.table("TCGA_colon_transcriptome_summary_of_correlation_between_transcript_and_miRNA.txt",sep="\t",header = T,stringsAsFactors = F)

# merge result of correlation analysis and previously merged table
m.cor.result <-merge(cor.result,primir.anot,by=c("miRNA","transcript"))

# import table about mean of expression leve of miRNA and transcript
# this table is located at "https://github.com/Ryosuke-Hirota/20230206_TCGA_colon_confirmation_of_expression_levels_of_transcript_and_miRNA"
setwd("C:/Rdata/20230206_TCGA_colon_confirmation_of_expression_levels_of_transcript_and_miRNA")
mean.df <-read.table("TCGA_colon_confirmation_of_expression_level_of_transcript_and_miRNA.txt",sep="\t",header = T,stringsAsFactors = F)
mean.df <-mean.df[,c(1,2,9,13)]

# merge table about mean of expression level and previously merged table
m.mean.df <-merge(mean.df,m.cor.result,by=c("miRNA","transcript"))

# set cutoff value
cutoff <-seq(0,200,50)

# make directory for protein coding transcript
setwd("C:/Rdata/20230214_TCGA_colon_confirmation_of_expression_levels_of_transcript_and_miRNA_for_each_gene_type")
dir.create("protein_coding")
setwd("C:/Rdata/20230214_TCGA_colon_confirmation_of_expression_levels_of_transcript_and_miRNA_for_each_gene_type/protein_coding")

# draw plot
for (i in 1:length(cutoff)) {
  
  ### draw volcano plot to show how many combinations have signigicant positive correlation
  # remove NA row and annotate whether combinations have significant positive correlations
  p.vp <-m.mean.df[m.mean.df[,8]=="protein_coding",]
  p.vp <-p.vp[p.vp[,7]>=cutoff[i],]
  p.vp <-p.vp[!is.na(p.vp[,5]),]
  p.vp[,9] <-NA
  p.vp[p.vp[,5]>0&p.vp[,6]<0.05,9] <-"red"
  p.vp[p.vp[,5]<=0|p.vp[,6]>=0.05,9] <-"grey40"

  # draw valcano plot
  pdf(paste0("valcano_plot_about_TCGA_colon_correlation_between_expression_level_of_transcript_and_miRNA_",cutoff[i],".pdf"))
  plot(p.vp[,5],-log10(p.vp[,6]),xlab="correlation coefficient",ylab="-log10 (p.value)",col=p.vp[,9],pch=19,
     main = paste0("r>0 p<0.05 ",nrow(p.vp[p.vp[,9]=="red",]),"/",nrow(p.vp)," (sig.posi / total)"))
  legend("topleft",legend =c("other","r>0, p<0.05"),col=unique(p.vp[,9]),pch=19)
  abline(h=log10(0.05)*-1,v=0,lty=2)
  dev.off()

  ### draw scatter plot to investigate whether combinations without correlation have small mean of miRNA/transcript expression 
  # remove NA row and annotate whether combinations have significant positive correlations
  p.sp <-m.mean.df[m.mean.df[,8]=="protein_coding",]
  p.sp <-p.sp[p.sp[,7]>=cutoff[i],]
  p.sp <-p.sp[!is.na(p.sp[,5]),]
  p.sp[,9] <-NA
  p.sp[p.sp[,5]>0&p.sp[,6]<0.05,9] <-"red"
  p.sp[p.sp[,5]<=0|p.sp[,6]>=0.05,9] <-"grey40"

  # draw scatter plot
  pdf(paste0("plot_about_mean_of_TCGA_colon_expression_level_of_miRNA_and_transcript_",cutoff[i],".pdf"))
  plot(log2(p.sp[,3]),log2(p.sp[,4]),xlab="log2(mean of miRNA expression)",ylab="log2(mean of transcript expression)",pch=19,col=p.sp[,9],
     main = paste0("sig.posi.cor = ",nrow(p.sp[p.sp[,9]=="red",])," other = ",nrow(p.sp)-nrow(p.sp[p.sp[,9]=="red",])))
  legend("topright",legend =c("other","r>0, p<0.05"),col=unique(p.sp[,9]),pch=19)
  abline(h=0,v=0,lty=2)
  dev.off()
}

# make directory for lncRNA transcript
setwd("C:/Rdata/20230214_TCGA_colon_confirmation_of_expression_levels_of_transcript_and_miRNA_for_each_gene_type")
dir.create("lncRNA")
setwd("C:/Rdata/20230214_TCGA_colon_confirmation_of_expression_levels_of_transcript_and_miRNA_for_each_gene_type/lncRNA")

# draw plot
for (i in 1:length(cutoff)) {
  
  ### draw volcano plot to show how many combinations have signigicant positive correlation
  # remove NA row and annotate whether combinations have significant positive correlations
  l.vp <-m.mean.df[m.mean.df[,8]=="lncRNA",]
  l.vp <-l.vp[l.vp[,7]>=cutoff[i],]
  l.vp <-l.vp[!is.na(l.vp[,5]),]
  l.vp[,9] <-NA
  l.vp[l.vp[,5]>0&l.vp[,6]<0.05,9] <-"red"
  l.vp[l.vp[,5]<=0|l.vp[,6]>=0.05,9] <-"grey40"
  
  # draw valcano plot
  pdf(paste0("valcano_plot_about_TCGA_colon_correlation_between_expression_level_of_transcript_and_miRNA_",cutoff[i],".pdf"))
  plot(l.vp[,5],-log10(l.vp[,6]),xlab="correlation coefficient",ylab="-log10 (p.value)",col=l.vp[,9],pch=19,
       main = paste0("r>0 p<0.05 ",nrow(l.vp[l.vp[,9]=="red",]),"/",nrow(l.vp)," (sig.posi / total)"))
  legend("topleft",legend =c("other","r>0, p<0.05"),col=unique(l.vp[,9]),pch=19)
  abline(h=log10(0.05)*-1,v=0,lty=2)
  dev.off()
  
  ### draw scatter plot to investigate whether combinations without correlation have small mean of miRNA/transcript expression 
  # remove NA row and annotate whether combinations have significant positive correlations
  l.sp <-m.mean.df[m.mean.df[,8]=="lncRNA",]
  l.sp <-l.sp[l.sp[,7]>=cutoff[i],]
  l.sp <-l.sp[!is.na(l.sp[,5]),]
  l.sp[,9] <-NA
  l.sp[l.sp[,5]>0&l.sp[,6]<0.05,9] <-"red"
  l.sp[l.sp[,5]<=0|l.sp[,6]>=0.05,9] <-"grey40"
  
  # draw scatter plot
  pdf(paste0("plot_about_mean_of_TCGA_colon_expression_level_of_miRNA_and_transcript_",cutoff[i],".pdf"))
  plot(log2(l.sp[,3]),log2(l.sp[,4]),xlab="log2(mean of miRNA expression)",ylab="log2(mean of transcript expression)",pch=19,col=l.sp[,9],
       main = paste0("sig.posi.cor = ",nrow(l.sp[l.sp[,9]=="red",])," other = ",nrow(l.sp)-nrow(l.sp[l.sp[,9]=="red",])))
  legend("topright",legend =c("other","r>0, p<0.05"),col=unique(l.sp[,9]),pch=19)
  abline(h=0,v=0,lty=2)
  dev.off()
}
