##### Part-1: Genomic Feature Enrichment #####


##### Adapted from Alex's codes
library(limma)
library(reshape)
#library(rstatix)

load('/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/results_all.RData')
load('/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/summary_r2g003_62hits.RData')

load("/data/js95/ALSPAC/ARIES/featureData.Rdata")

price.annot <- read.table("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Annotation/price_annotation_450k.txt", 
                          header=T, sep='\t')

annot <- price.annot[which(price.annot$ID %in% results.jobless.sI$Probe),] #epigenome background
rm(results.Fscore.sI,results.income.sI, results.income_reduction.sI,results.jobless.sI, results.major_F.sI, results.nbhqual.sI)
annot.hits <- annot[which(annot$ID %in% cpg_r2_g003$Probe),] #hits
dim(annot.hits)

featureData.annot <- featureData[which(featureData$TargetID %in% annot$ID),]
featureData.annot$Genomic_location <- strsplit2(featureData.annot$UCSC_REFGENE_GROUP, split=";")[,1] #select the first 
features.hits <- featureData.annot[which(featureData.annot$TargetID %in% cpg_r2_g003$Probe),]
dim(features.hits)

features.summary <- data.frame(hits = summary(as.factor(features.hits$Genomic_location))/
                                 nrow(features.hits),
                               array = summary(as.factor(featureData.annot$Genomic_location))/
                                 nrow(featureData.annot))

features.summary$feature <- rownames(features.summary)
features.summary$feature[which(features.summary$feature=="")] <- "Intergenic"
features.melt <- melt(features.summary)

features.melt$feature<-factor(features.melt$feature,levels=c("TSS1500","TSS200","5'UTR","1stExon","Body","3'UTR","Intergenic"))
features.melt$variable<-factor(features.melt$variable,levels=c("hits","array"),labels = c("CpGs with R²>3% (n=62)", "All tested CpGs (n=412,956)"))



g1<-ggplot(features.melt, aes(y=value*100, x=feature, fill = variable, group=variable))+
  geom_bar(position=position_dodge(), stat='identity', width = 0.5)+
  theme_classic()+
  scale_fill_manual(values = wes_palette("Chevalier1"))+
  ylab("Percent of CpGs")+
  xlab("Genomic feature")+
  theme(
    axis.title.x =element_text(size = 20),
    axis.text.x = element_text(size = 17),
    axis.title.y = element_text(size=20),
    axis.text.y = element_text(size = 17),
    legend.text=element_text(size=18)
  ) +
  labs(fill=" ")

png("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/Enrichment/figure_genomic_enrichment-2021-05-20.png", 1100, 350)
g1
dev.off()


#chi-square test
temp.fea.sum<-features.summary[,1:2]
temp.fea.sum[,1]<-temp.fea.sum[,1]*62
temp.fea.sum[,2]<-temp.fea.sum[,2]*412956
temp.fea.sum[,2]<-temp.fea.sum[,2]-temp.fea.sum[,1]# "non-significant" probes
chisq.test(temp.fea.sum)

row_wise_fisher_test(temp.fea.sum)
row_wise_fisher_test(temp.fea.sum,p.adjust.method = "bonferroni",alternative = "greater")
row_wise_fisher_test(temp.fea.sum,p.adjust.method = "BH",alternative = "greater")


## CGI enrichment 


CGI.summary <- data.frame(hits = summary(as.factor(features.hits$RELATION_TO_UCSC_CPG_ISLAND))/
                            nrow(features.hits),   
                               array = summary(as.factor(featureData.annot$RELATION_TO_UCSC_CPG_ISLAND))/
                                 nrow(featureData.annot))

CGI.summary$feature <- rownames(CGI.summary)
CGI.summary$feature[which(CGI.summary$feature=="")] <- "Open sea"
CGI.melt <- melt(CGI.summary)


CGI.melt$feature<-factor(CGI.melt$feature,levels=c("Open sea","N_Shelf","N_Shore","Island","S_Shore","S_Shelf"),labels = c("Open sea","N.Shelf","N.Shore","Island","S.Shore","S.Shelf"))
CGI.melt$variable<-factor(CGI.melt$variable,levels=c("hits","array"),labels = c("CpGs with R²>3% (n=62)", "All tested CpGs (n=412,956)"))




g1<-ggplot(CGI.melt, aes(y=value*100, x=feature, fill = variable, group=variable))+
  geom_bar(position=position_dodge(), stat='identity', width = 0.5)+
  theme_classic()+
  scale_fill_manual(values = wes_palette("Chevalier1"))+
  ylab("Percent of CpGs")+
  xlab("Relation to CpG Island")+
  theme(
    axis.title.x =element_text(size = 20),
    axis.text.x = element_text(size = 17),
    axis.title.y = element_text(size=20),
    axis.text.y = element_text(size = 17),
    legend.text=element_text(size=18)
  ) +
  labs(fill=" ")

png("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/Enrichment/figure_CGI_enrichment-2021-05-20.png", 1090, 350)
g1
dev.off()


#chi-squre test
temp.CGI.sum<-CGI.summary[,1:2]
temp.CGI.sum[,1]<-temp.CGI.sum[,1]*62
temp.CGI.sum[,2]<-temp.CGI.sum[,2]*412956
temp.CGI.sum[,2]<-temp.CGI.sum[,2]-temp.CGI.sum[,1]# "non-significant" probes
chisq.test(temp.CGI.sum)
row_wise_fisher_test(temp.CGI.sum,p.adjust.method = "BH")
row_wise_fisher_test(temp.CGI.sum,p.adjust.method = "bonferroni")
row_wise_fisher_test(temp.CGI.sum,p.adjust.method = "BH",alternative = "greater")

## enhancer enrichment 

anno<-read.csv("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Annotation/HumanMethylation450_15017482_v1-2.csv",header=T)
featureData.annot.enhancer<-merge(featureData.annot,anno[,c("IlmnID", "Enhancer")],by.x="TargetID",by.y="IlmnID")
table(featureData.annot.enhancer$Enhancer)

features.hits.enhancer <- merge(features.hits,anno[,c("IlmnID", "Enhancer")],by.x="TargetID",by.y="IlmnID")
table(features.hits.enhancer$Enhancer)

# test

16/62-92837/412956
featureData.annot.enhancer
datatable <- matrix(c(16,92837,62-16,412956-92837),nrow=2,ncol=2)
datatable

chisq.test(datatable,correct=FALSE) #p=0.53
fisher.test(datatable,alternative = "greater") #p=0.3

# plot

enhancer.summary <- data.frame(hits = summary(as.factor(features.hits.enhancer$Enhancer))/
                            nrow(features.hits.enhancer),   
                          array = summary(as.factor(featureData.annot.enhancer$Enhancer))/
                            nrow(featureData.annot.enhancer))

enhancer.summary <- enhancer.summary[1,]
enhancer.summary$feature <- " "
enhancer.melt <- melt(enhancer.summary)


enhancer.melt$variable<-factor(enhancer.melt$variable,levels=c("hits","array"),labels = c("CpGs with R²>3% (n=62)", "All tested CpGs (n=412,956)"))




g1<-ggplot(enhancer.melt, aes(y=value*100, x=feature, fill = variable, group=variable))+
  geom_bar(position=position_dodge(), stat='identity', width = 0.5)+
  theme_classic()+
  scale_fill_manual(values = wes_palette("Chevalier1"))+
  ylab("Percent of CpGs")+
  xlab("Enhancer")+
  theme(
    axis.title.x =element_text(size = 20),
    axis.text.x = element_text(size = 17),
    axis.title.y = element_text(size=20),
    axis.text.y = element_text(size = 17),
    legend.text=element_text(size=18)
  ) +
  labs(fill=" ")

png("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/Enrichment/figure_Enhancer_enrichment-2021-05-20.png", 600, 350)
g1
dev.off()




##### Part-2: Biological Pathway Enrichment (no updates on 2021-05-12) #####
##### Gene Set analysis 

library("methylGSA")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)


results.list<-c("results.jobless.sI","results.income_reduction.sI","results.Fscore.sI","results.income.sI","results.major_F.sI","results.nbhqual.sI")

for (i in results.list) {
  df<-get(i)
  df <- df[order(df$`1.p`),]
  p.adver<-df$`1.p`
  names(p.adver)<-rownames(df)
  GSA.KEGG<- methylRRA(cpg.pval = p.adver, method = "GSEA", GS.type="KEGG",
                                minsize = 50, maxsize = 800)
  
  GSA.GO  <- methylRRA(cpg.pval = p.adver , method = "GSEA", GS.type="GO",
                               minsize = 50, maxsize = 800)
  
  adver <-strsplit(i, "[.]")[[1]][2]
  assign(paste0("GSA.KEGG.",adver),GSA.KEGG)
  assign(paste0("GSA.GO.",adver),GSA.GO)
}


GSA.KEGG.jobless[1:20,c("Description","Size","pvalue","padj","enrichmentScore")]
GSA.KEGG.income_reduction[1:20,c("Description","Size","pvalue","padj","enrichmentScore")]
GSA.KEGG.Fscore[1:20,c("Description","Size","pvalue","padj","enrichmentScore")]
GSA.KEGG.income[1:20,c("Description","Size","pvalue","padj","enrichmentScore")]
GSA.KEGG.major_F[1:20,c("Description","Size","pvalue","padj","enrichmentScore")]
GSA.KEGG.nbhqual[1:20,c("Description","Size","pvalue","padj","enrichmentScore")]

GSA.GO.jobless[1:20,c("Description","Size","pvalue","padj","enrichmentScore")]
GSA.GO.income_reduction[1:20,c("Description","Size","pvalue","padj","enrichmentScore")]
GSA.GO.Fscore[1:20,c("Description","Size","pvalue","padj","enrichmentScore")]
GSA.GO.income[1:20,c("Description","Size","pvalue","padj","enrichmentScore")]
GSA.GO.major_F[1:20,c("Description","Size","pvalue","padj","enrichmentScore")]
GSA.GO.nbhqual[1:20,c("Description","Size","pvalue","padj","enrichmentScore")]

save(GSA.KEGG.jobless,GSA.KEGG.income_reduction,GSA.KEGG.Fscore,GSA.KEGG.income,GSA.KEGG.major_F,GSA.KEGG.nbhqual,
     GSA.GO.jobless,GSA.GO.income_reduction,GSA.GO.Fscore,GSA.GO.income,GSA.GO.major_F,GSA.GO.nbhqual,
     file='/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/enrichment_GSA_results.RData')
load('/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/enrichment_GSA_results.RData')

# output the top GO terms for the REVIGO analysis (http://revigo.irb.hr). 
GSA.GO.income[1:100,c("Description","pvalue")]
write.csv(GSA.GO.jobless[1:100,c("Description","pvalue")],file='/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/Enrichment/GSA_jobless_top100_GO.csv', row.names = TRUE)
write.csv(GSA.GO.income[1:100,c("Description","pvalue")],file='/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/Enrichment/GSA_income_top100_GO.csv', row.names = TRUE)
write.csv(GSA.GO.income_reduction[1:100,c("Description","pvalue")],file='/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/Enrichment/GSA_income_reduction_top100_GO.csv', row.names = TRUE)
write.csv(GSA.GO.Fscore[1:100,c("Description","pvalue")],file='/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/Enrichment/GSA_Fscore_top100_GO.csv', row.names = TRUE)
write.csv(GSA.GO.major_F[1:100,c("Description","pvalue")],file='/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/Enrichment/GSA_major_F_top100_GO.csv', row.names = TRUE)
write.csv(GSA.GO.nbhqual[1:100,c("Description","pvalue")],file='/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/Enrichment/GSA_nbhqual_top100_GO.csv', row.names = TRUE)





#make a plot
library(ggpubr)

GO_results<-rbind(data.frame(GO=GSA.GO.jobless$ID, pathway=GSA.GO.jobless$Description,pvalue=GSA.GO.jobless$pvalue,pvalue.adj=GSA.GO.jobless$padj,adversity="Job loss", measure="jobless"),
                  data.frame(GO=GSA.GO.income_reduction$ID, pathway=GSA.GO.income_reduction$Description,pvalue=GSA.GO.income_reduction$pvalue,pvalue.adj=GSA.GO.income_reduction$padj,adversity="Income reduction", measure="income_reduction"),
                  data.frame(GO=GSA.GO.income$ID, pathway=GSA.GO.income$Description,pvalue=GSA.GO.income$pvalue,pvalue.adj=GSA.GO.income$padj,adversity="Low family income", measure="income"),
                  data.frame(GO=GSA.GO.Fscore$ID, pathway=GSA.GO.Fscore$Description,pvalue=GSA.GO.Fscore$pvalue,pvalue.adj=GSA.GO.Fscore$padj,adversity="Financial hardship", measure="Fscore"),
                  data.frame(GO=GSA.GO.major_F$ID, pathway=GSA.GO.major_F$Description,pvalue=GSA.GO.major_F$pvalue,pvalue.adj=GSA.GO.major_F$padj,adversity="Major financial problem", measure="major_F"),
                  data.frame(GO=GSA.GO.nbhqual$ID, pathway=GSA.GO.nbhqual$Description,pvalue=GSA.GO.nbhqual$pvalue,pvalue.adj=GSA.GO.nbhqual$padj,adversity="Neighborhood disadvantage", measure="nbhqual")
)


GO_results.sig<-GO_results[GO_results$pvalue<0.001,]
GO_results.sig$adversity<-as.factor(GO_results.sig$adversity)
GO_results.sig$adversity<-factor(GO_results.sig$adversity,levels=c("Job loss","Income reduction","Low family income","Financial hardship","Major financial problem","Neighborhood disadvantage"))

make_go_plot <- function (adv) {
  GO_results.sig.adv <- GO_results.sig[GO_results.sig$adversity==adv,]
  
  g<-ggplot(GO_results.sig.adv,
             aes(x =pathway, y= -log10(pvalue)))+
    geom_hline(yintercept = -log10(0.05/3167), col='red', linetype=2)+
    geom_bar(stat='identity', width=0.1)+
    geom_point(size=3)+
    scale_x_discrete(" ", limits=rev(GO_results.sig.adv$pathway), labels=rev(GO_results.sig.adv$pathway)) + 
    theme_classic()+
    theme(axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 13, face="bold"),
          strip.text.x = element_text(size = 13, face="bold"),
          plot.background = element_rect(fill = "transparent",colour = NA))+
    ylim(0,6)+
    ylab("-log10(P-value)")+
    facet_wrap(~adversity, nrow=1)+
    coord_flip()
  return(g)
}

g1 <- make_go_plot("Job loss")
g2 <- make_go_plot("Income reduction")
g3 <- make_go_plot("Low family income")
g4 <- make_go_plot("Financial hardship")
g5 <- make_go_plot("Major financial problem")
g6 <- make_go_plot("Neighborhood disadvantage")


png("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/Enrichment/figure_GO_enrichment_methylGSA_pe3_2021-04-04.png", 1550, 700,bg = "transparent")
ggarrange(
  ggarrange(g1,g2,ncol=1, align = "v", heights = c(0.9,0.3)),
  ggarrange(g3,g4,g6,ncol=1,align = "v", heights = c(0.38,1,0.7)),
  g5,
  ncol = 3, widths = c(1.03, 1.03, 1.15))

dev.off()




