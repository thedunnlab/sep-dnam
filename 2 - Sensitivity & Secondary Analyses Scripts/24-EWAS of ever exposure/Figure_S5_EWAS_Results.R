##### Examine EWAS results  #####
library(dplyr)
library(data.table)
library(ggplot2)
library(wesanderson)

##load results
setwd("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/newDNAm/EWAS")

files <- list.files(pattern = "*081820.Rdata")

for(i in files){
  load(i, verbose=T)
}


results.list<-c("results.EWAS.jobless","results.EWAS.income_reduction","results.EWAS.Fscore","results.EWAS.income_bin","results.EWAS.major_F","results.EWAS.nbhqual")

for (i in results.list) {
  df<-get(i)
  df <- df[order(df$p_value),]
  print(head(df))
  assign(i,df)
}




##:::::::Remove sex chr and cross-hybridzed CpGs :::::::##

anno<-read.csv("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Annotation/HumanMethylation450_15017482_v1-2.csv",header=T)
dim(anno) #486428 probes
sex.chr<-anno %>% dplyr::filter(CHR=="X" | CHR=="Y") %>% dplyr::select(IlmnID)
sex.chr<-as.character(unlist(sex.chr$IlmnID)) #11648

# List of hybridized probes
price <- fread("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Annotation/price_annotation_450k.txt") 
dim(price)
hybrid<-price %>% dplyr::filter(XY_Hits=="XY_YES" | Autosomal_Hits=="A_YES")  %>% dplyr::select(ID)
hybrid<-as.character(unlist(hybrid$ID)) #41937

# Remove the probes from the results files (#412956 CpGs left)

for (i in results.list) {
  df<-get(i)
  df <-df[!(df$probe %in% sex.chr | df$probe %in% hybrid),]
  print(dim(df))
  assign(i,df)
}





#### Correction for multiple testing ####
bonf_total <- 0
fdr005_total <- 0
fdr01_total <- 0
fdr015_total <- 0
fdr02_total <- 0

for (i in results.list) {
  df<-get(i)
  fdr.df<-p.adjust(df$`1.p`, method = "BH")
  bonf.df<-p.adjust(df$`1.p`, method = "bonferroni")
  print(paste0(i," - ",sum(bonf.df<0.05)," CpGs passed Bonf.p<0.05"))
  print(paste0(i," - ",sum(fdr.df<0.05)," CpGs passed FDR<0.05"))
  print(paste0(i," - ",sum(fdr.df<0.1)," CpGs passed FDR<0.1"))
  print(paste0(i," - ",sum(fdr.df<0.15)," CpGs passed FDR<0.15"))
  print(paste0(i," - ",sum(fdr.df<0.2)," CpGs passed FDR<0.2"))
  bonf_total <- bonf_total+sum(bonf.df<0.05)
  fdr005_total <- fdr005_total +sum(fdr.df<0.05)
  fdr01_total <- fdr01_total +sum(fdr.df<0.1)
  fdr015_total <- fdr015_total +sum(fdr.df<0.15)
  fdr02_total <- fdr02_total +sum(fdr.df<0.2)
}
print(bonf_total)
print(fdr005_total)
print(fdr01_total)
print(fdr015_total)
print(fdr02_total) # all 0


# Create a summary object
EWAS.summary <- rbind(data.frame(Probe = results.EWAS.Fscore$probe, p.value =results.EWAS.Fscore$p_value, FDR_q= p.adjust(results.EWAS.Fscore$p_value, method = "BH"), 
                               expected = c(1:412956)/412956,measure = "Financial difficulty score",
                               Beta=results.EWAS.Fscore$beta, adversity=results.EWAS.Fscore$adversity),
                    
                    data.frame(Probe = results.EWAS.income_bin$probe, p.value =results.EWAS.income_bin$p_value, FDR_q= p.adjust(results.EWAS.income_bin$p_value, method = "BH"), 
                               expected = c(1:412956)/412956, measure = "Low family income",
                               Beta=results.EWAS.income_bin$beta, adversity=results.EWAS.income_bin$adversity),
                    
                    data.frame(Probe = results.EWAS.income_reduction$probe, p.value =results.EWAS.income_reduction$p_value, FDR_q= p.adjust(results.EWAS.income_reduction$p_value, method = "BH"),
                               expected = c(1:412956)/412956,   measure = "Income reduction",
                               Beta=results.EWAS.income_reduction$beta, adversity=results.EWAS.income_reduction$adversity),
                    
                    data.frame(Probe = results.EWAS.jobless$probe, p.value =results.EWAS.jobless$p_value, FDR_q= p.adjust(results.EWAS.jobless$p_value, method = "BH"),
                               expected = c(1:412956)/412956,  measure = "Job loss",
                               Beta=results.EWAS.jobless$beta, adversity=results.EWAS.jobless$adversity),
                    
                    data.frame(Probe = results.EWAS.major_F$probe, p.value =results.EWAS.major_F$p_value, FDR_q= p.adjust(results.EWAS.major_F$p_value, method = "BH"),
                               expected = c(1:412956)/412956,  measure = "Major financial problem",
                               Beta=results.EWAS.major_F$beta, adversity=results.EWAS.major_F$adversity),
                    
                    data.frame(Probe = results.EWAS.nbhqual$probe, p.value = results.EWAS.nbhqual$p_value, FDR_q= p.adjust(results.EWAS.nbhqual$p_value, method = "BH"),
                               expected = c(1:412956)/412956, measure = "Neighborhood disadvantage",
                               Beta=results.EWAS.nbhqual$beta, adversity=results.EWAS.nbhqual$adversity))




## QQ plot
ggplot(EWAS.summary, aes(y = -log10(p.value), 
                       x = -log10(expected)))+
  geom_point(alpha =0.5, size = 0.5)+
  geom_abline()+
  facet_wrap(~measure)+
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14))+
  ylab("Observed -log10(p)")+
  xlab("Expected -log10(p)")+
  ggtitle("QQ plots for EWAS (excluding sex chr and hybridized CpG)")


## p-value distribution
ggplot(EWAS.summary, aes(p.value))+
  geom_density()+
  facet_wrap(~measure)+
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14))+
  ylab("Density")+
  xlab("P-value")+
  ggtitle("P-value distributions for EWAS (excluding sex chr and hybridized CpG)")


### Save results

save(EWAS.summary,file='/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/newDNAm/EWAS/EWAS_allresults.RData')
#load('/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/newDNAm/EWAS/EWAS_allresults.RData')

# fix naming issue of income
EWAS.summary[EWAS.summary$adversity=="income_bin","adversity"]<-"lowincome"

### Compare with SLCMA results 

# load SLCMA results  (62 R^2 hits)

load("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/Table2/Table2_results_62hits_2021-05-16.Rdata")

### compare effect estimates 

##for SLCMA hits
compare<-left_join(table2.final[,c("CpG","adversity","Beta_SI")], EWAS.summary[,c("Probe","adversity","Beta")], by = c("CpG" = "Probe","adversity"))
compare$Beta_SI<-as.numeric(compare$Beta_SI)

compare$adversity<-factor(compare$adversity,levels=c("jobless","income_reduction","lowincome","Fscore","major_F","nbhqual"),
                          labels = c("Job loss", "Income reduction","Low family income","Financial hardship","Major financial problem","Neighborhood disadvantage"))
compare$adversity

summary(abs(compare$Beta_SI)-abs(compare$Beta)) # effects are greater in magnitude in SLCMA than EWAS for all 62 hits
sum(abs(compare$Beta_SI)>abs(compare$Beta))

g1<-ggplot(compare, aes(x=Beta, y=Beta_SI, shape=adversity, color=adversity)) +
  geom_point(size=6)+
  geom_abline(intercept=0, slope=1, color="red", linetype=2)+
  geom_vline(xintercept = 0, color="grey",linetype="dashed")+
  geom_hline(yintercept = 0, color="grey",linetype="dashed")+
  ylab("SLCMA effect estimates ")+
  xlab("EWAS effect estimates (ever-exposed vs. never-exposed)")+
  xlim(-0.13, 0.10)+
  ylim(-0.13, 0.10)+
  scale_color_manual(values=c("#88CCEE", "#AA4499","#CC6677","#DDCC77","#332288","#44AA99"))+
  scale_shape_manual(values = c(20, 18,8, 15,17,16))+
  labs(color='SEP Measure',shape='SEP Measure')+ theme_bw(base_size = 17)+
  theme(plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        legend.background = element_rect(fill = "transparent",colour = NA),
        legend.key = element_rect(fill = "transparent",colour = NA),
        legend.position="right",legend.title = element_blank(),
        text = element_text(size = 32),
        axis.text.x = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=20))

png("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/EWAS/EWAS_comparebeta_62hits_2021-05-27.png", 900, 500,bg = "transparent")
g1
dev.off()



# subset to the 4 FDR-significant CpGs
load('/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/summary_fdr005_4hits.RData')

table2.final<- table2.final[table2.final$CpG %in% cpg_fdr005$Probe,]
dim(table2.final)


### compare effect estimates 

##for SLCMA hits
compare<-left_join(table2.final[,c("CpG","adversity","Beta_SI")], EWAS.summary[,c("Probe","adversity","Beta")], by = c("CpG" = "Probe","adversity"))
compare$Beta_SI<-as.numeric(compare$Beta_SI)

compare$adversity<-factor(compare$adversity,levels=c("nbhqual"),
                          labels = c("Neighborhood Disadvantage"))
compare$adversity

summary(abs(compare$Beta_SI-compare$Beta)) 


g1<-ggplot(compare, aes(x=Beta, y=Beta_SI, shape=adversity, color=adversity)) +
  geom_point(size=10)+
  geom_abline(intercept=0, slope=1, color="red", linetype=2)+
  geom_vline(xintercept = 0, color="grey",linetype="dashed")+
  geom_hline(yintercept = 0, color="grey",linetype="dashed")+
  ylab("SLCMA effect estimates ")+
  xlab("EWAS effect estimates (ever-exposed vs. never-exposed)")+
  xlim(-0.13, 0.10)+
  ylim(-0.13, 0.10)+
  scale_color_manual(values=c("#44AA99"))+
  scale_shape_manual(values = c(20, 18,8, 15,17,16))+
  labs(color='SEP Measure',shape='SEP Measure')+ theme_bw(base_size = 17)+
  theme(plot.background = element_rect(fill = "transparent",colour = NA),
        legend.background = element_rect(fill = "transparent",colour = NA),
        legend.key = element_rect(fill = "transparent",colour = NA),
        legend.position="right",legend.title = element_blank(),
        text = element_text(size = 32),
        axis.text.x = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=20))

png("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/EWAS/EWAS_comparebeta_2021-05-12.png", 900, 500,bg = "transparent")
g1
dev.off()

