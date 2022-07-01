#####################################
# Prepare Final Results Reporting
# Correct for the p-value threshold error
# 2021-05-07
# Jiaxuan (Jessie) Liu
#####################################



## load packages
library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(dplyr)



## Load results
setwd("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/newDNAm/FWL_Results")

files <- list.files(pattern = "*031720_sI.Rdata")

for(i in files){
  load(i, verbose=T)
}

results.list<-c("results.jobless.sI","results.income_reduction.sI","results.Fscore.sI","results.income.sI","results.major_F.sI","results.nbhqual.sI")

for (i in results.list) {
  df<-get(i)
  df <- df[order(df$`1.p`),]
  print(head(df))
  assign(i,df)
}


##:::::::Remove sex chr and cross-hybridzed CpGs :::::::##

anno<-read.csv("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Annotation/HumanMethylation450_15017482_v1-2.csv",header=T)
dim(anno) #486428 probes
sex.chr<-anno %>% filter(CHR=="X" | CHR=="Y") %>% select(IlmnID)
sex.chr<-as.character(unlist(sex.chr$IlmnID)) #11648

# List of hybridized probes
price <- fread("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Annotation/price_annotation_450k.txt") 
dim(price)
hybrid<-price %>% filter(XY_Hits=="XY_YES" | Autosomal_Hits=="A_YES")  %>% select(ID)
hybrid<-as.character(unlist(hybrid$ID)) #41937


# Remove the probes from the results files (#412956 CpGs left)

for (i in results.list) {
  df<-get(i)
  df <-df[!(rownames(df) %in% sex.chr | rownames(df) %in% hybrid),]
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
print(fdr02_total)
# "results.jobless.sI - 0 CpGs passed Bonf.p<0.05"
# "results.jobless.sI - 0 CpGs passed FDR<0.05"
# "results.jobless.sI - 0 CpGs passed FDR<0.1"
# "results.jobless.sI - 11 CpGs passed FDR<0.2"
# "results.jobless.sI - 24 CpGs passed FDR<0.3"
# "results.income_reduction.sI - 0 CpGs passed Bonf.p<0.05"
# "results.income_reduction.sI - 0 CpGs passed FDR<0.05"
# "results.income_reduction.sI - 0 CpGs passed FDR<0.1"
# "results.income_reduction.sI - 0 CpGs passed FDR<0.2"
# "results.income_reduction.sI - 0 CpGs passed FDR<0.3"
# "results.Fscore.sI - 0 CpGs passed Bonf.p<0.05"
# "results.Fscore.sI - 0 CpGs passed FDR<0.05"
# "results.Fscore.sI - 0 CpGs passed FDR<0.1"
# "results.Fscore.sI - 3 CpGs passed FDR<0.2"
# "results.Fscore.sI - 20 CpGs passed FDR<0.3"
# "results.income.sI - 0 CpGs passed Bonf.p<0.05"
# "results.income.sI - 0 CpGs passed FDR<0.05"
# "results.income.sI - 0 CpGs passed FDR<0.1"
# "results.income.sI - 0 CpGs passed FDR<0.2"
# "results.income.sI - 0 CpGs passed FDR<0.3"
# "results.major_F.sI - 0 CpGs passed Bonf.p<0.05"
# "results.major_F.sI - 0 CpGs passed FDR<0.05"
# "results.major_F.sI - 2 CpGs passed FDR<0.1"
# "results.major_F.sI - 4 CpGs passed FDR<0.2"
# "results.major_F.sI - 4 CpGs passed FDR<0.3"
# "results.nbhqual.sI - 0 CpGs passed Bonf.p<0.05"
# "results.nbhqual.sI - 4 CpGs passed FDR<0.05"
# "results.nbhqual.sI - 5 CpGs passed FDR<0.1"
# "results.nbhqual.sI - 18 CpGs passed FDR<0.2"
# "results.nbhqual.sI - 43 CpGs passed FDR<0.3"


# Create a summary object
sI.summary <- rbind(data.frame(Probe = results.Fscore.sI$Probe, hypothesis = results.Fscore.sI$`1.hypo`, 
                               r2 = results.Fscore.sI$`1.r2`, p.value = results.Fscore.sI$`1.p`, FDR_q= p.adjust(results.Fscore.sI$`1.p`, method = "BH"), 
                               expected = c(1:412956)/412956,measure = "Financial difficulty score"),
                    
                    data.frame(Probe = results.income.sI$Probe, hypothesis = results.income.sI$`1.hypo`, 
                               r2 = results.income.sI$`1.r2`, p.value = results.income.sI$`1.p`, FDR_q= p.adjust(results.income.sI$`1.p`, method = "BH"), 
                               expected = c(1:412956)/412956, measure = "Income"),
                    
                    data.frame(Probe = results.income_reduction.sI$Probe, hypothesis = results.income_reduction.sI$`1.hypo`, 
                               r2 = results.income_reduction.sI$`1.r2`, p.value = results.income_reduction.sI$`1.p`, FDR_q= p.adjust(results.income_reduction.sI$`1.p`, method = "BH"),
                               expected = c(1:412956)/412956,   measure = "Income reduction"),
                    
                    data.frame(Probe = results.jobless.sI$Probe, hypothesis = results.jobless.sI$`1.hypo`, 
                               r2 = results.jobless.sI$`1.r2`, p.value = results.jobless.sI$`1.p`, FDR_q= p.adjust(results.jobless.sI$`1.p`, method = "BH"), 
                               expected = c(1:412956)/412956,  measure = "Job loss"),
                    
                    data.frame(Probe = results.major_F.sI$Probe, hypothesis = results.major_F.sI$`1.hypo`, 
                               r2 = results.major_F.sI$`1.r2`, p.value = results.major_F.sI$`1.p`, FDR_q= p.adjust(results.major_F.sI$`1.p`, method = "BH"), 
                               expected = c(1:412956)/412956,  measure = "Major financial problem"),
                    
                    data.frame(Probe = results.nbhqual.sI$Probe, hypothesis = results.nbhqual.sI$`1.hypo`, 
                               r2 = results.nbhqual.sI$`1.r2`, p.value = results.nbhqual.sI$`1.p`, FDR_q= p.adjust(results.nbhqual.sI$`1.p`, method = "BH"), 
                               expected = c(1:412956)/412956, measure = "Neighborhood disadvantage"))


summary(as.factor(sI.summary$hypothesis))
sI.summary$hypothesis <- gsub("Fscore_newbin_","", sI.summary$hypothesis)
sI.summary$hypothesis <- gsub("income_bin_","", sI.summary$hypothesis)
sI.summary$hypothesis <- gsub("income_reduction_newbin_","", sI.summary$hypothesis)
sI.summary$hypothesis <- gsub("jobless_newbin_","", sI.summary$hypothesis)
sI.summary$hypothesis <- gsub("major_F_newbin_","", sI.summary$hypothesis)
sI.summary$hypothesis <- gsub("nbhqual_newbin_","", sI.summary$hypothesis)

sI.summary$hypothesis <- factor(sI.summary$hypothesis, levels = c("accumulation", "very_early","early","middle",
                                                                  "mobility_D12", "mobility_D23", "mobility_U12", "mobility_U23"))



## QQ plot
ggplot(sI.summary, aes(y = -log10(p.value), 
                       x = -log10(expected)))+
  geom_point(alpha =0.5, size = 0.5)+
  geom_abline()+
  facet_wrap(~measure)+
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14))+
  ylab("Observed -log10(p)")+
  xlab("Expected -log10(p)")+
  ggtitle("QQ plots for selective inference (excluding sex chr and hybridized CpG)")

## p-value distribution
ggplot(sI.summary, aes(p.value))+
  geom_density()+
  facet_wrap(~measure)+
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14))+
  ylab("Density")+
  xlab("P-value")+
  ggtitle("P-value distributions for selective inference (excluding sex chr and hybridized CpG)")


## Barplot of hypothesis by FDR q thresholds 

ggplot(sI.summary[which(sI.summary$FDR_q <0.05),], aes(x= hypothesis, fill= hypothesis))+
  geom_bar(stat='count', position=position_dodge(preserve = 'single'), width =0.5, color = 'grey')+
  xlab("Measure")+
  ylab("# of CpGs under hypothesis")+ 
  facet_wrap(~measure)+
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        axis.text.x = element_text(angle = 45, hjust= 1))+
  #scale_fill_brewer(palette = "RdYlBu")+
  scale_fill_manual(values = c("#fee090","#abd9e9","#74add1","#4575b4"))+
  #scale_fill_got(option = 'Tyrell', discrete=T)+
  ggtitle("FDR q<0.05")+
  geom_text(stat= 'count', aes(label = ..count..), position = position_dodge(width =1), vjust = -0.5)+
  ylim(0,10)


ggplot(sI.summary[which(sI.summary$FDR_q <0.1),], aes(x= hypothesis, fill= hypothesis))+
  geom_bar(stat='count', position=position_dodge(preserve = 'single'), width =0.5, color = 'grey')+
  xlab("Measure")+
  ylab("# of CpGs under hypothesis")+ 
  facet_wrap(~measure)+
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        axis.text.x = element_text(angle = 45, hjust= 1))+
  #scale_fill_brewer(palette = "RdYlBu")+
  scale_fill_manual(values = c("#fee090","#abd9e9","#74add1","#4575b4"))+
  #scale_fill_got(option = 'Tyrell', discrete=T)+
  ggtitle("FDR q<0.1")+
  geom_text(stat= 'count', aes(label = ..count..), position = position_dodge(width =1), vjust = -0.5)+
  ylim(0,10)


ggplot(sI.summary[which(sI.summary$FDR_q <0.2),], aes(x= hypothesis, fill= hypothesis))+
  geom_bar(stat='count', position=position_dodge(preserve = 'single'), width =0.5, color = 'grey')+
  xlab("Measure")+
  ylab("# of CpGs under hypothesis")+ 
  facet_wrap(~measure)+
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        axis.text.x = element_text(angle = 45, hjust= 1))+
  #scale_fill_brewer(palette = "RdYlBu")+
  scale_fill_manual(values = c("#f46d43","#fdae61","#fee090","#abd9e9","#74add1","#4575b4"))+
  #scale_fill_got(option = 'Tyrell', discrete=T)+
  ggtitle("FDR q<0.2")+
  geom_text(stat= 'count', aes(label = ..count..), position = position_dodge(width =1), vjust = -0.5)+
  ylim(0,10)

ggplot(sI.summary[which(sI.summary$FDR_q <0.3),], aes(x= hypothesis, fill= hypothesis))+
  geom_bar(stat='count', position=position_dodge(preserve = 'single'), width =0.5, color = 'grey')+
  xlab("Measure")+
  ylab("# of CpGs under hypothesis")+ 
  facet_wrap(~measure)+
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        axis.text.x = element_text(angle = 45, hjust= 1))+
  scale_fill_brewer(palette = "RdYlBu")+
  #scale_fill_got(option = 'Tyrell', discrete=T)+
  ggtitle("FDR q<0.3")+
  geom_text(stat= 'count', aes(label = ..count..), position = position_dodge(width =1), vjust = -0.5)+
  ylim(0,20)



### Save results

cpg_pe5<-sI.summary[sI.summary$`p.value`<1e-05,]#64
save(cpg_pe5,file='/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/summary_pe5_64hits.RData')


cpg_fdr005<-sI.summary[sI.summary$FDR_q <0.05,]#4
save(cpg_fdr005,file='/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/summary_fdr005_4hits.RData')
print(cpg_fdr005)

cpg_fdr01<-sI.summary[sI.summary$FDR_q <0.1,]#7
save(cpg_fdr01,file='/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/summary_fdr01_7hits.RData')

cpg_fdr015<-sI.summary[sI.summary$FDR_q <0.15,]#14
save(cpg_fdr015,file='/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/summary_fdr015_14hits.RData')


cpg_fdr02<-sI.summary[sI.summary$FDR_q <0.2,]#36
save(cpg_fdr02,file='/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/summary_fdr02_36hits.RData')

cpg_r2_g003<-sI.summary[sI.summary$r2 >0.03,]#62
save(cpg_r2_g003,file='/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/summary_r2g003_62hits.RData')

# sort the 62 hits
cpg.order <- read.csv(file = '/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/Table2/62hits_order.csv')

cpg_r2_g003 <- cpg_r2_g003[match(cpg.order$CpG, cpg_r2_g003$Probe),] # order the CpGs 

write.csv(cpg_r2_g003, "/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/Table2/Table2_62hits_FDR_2021-05-19.csv", row.names = F)


# Hypotheses selected among R^2>0.03 hits
table(cpg_r2_g003$measure)
adv.list <- c("Financial difficulty score","Income","Income reduction","Job loss","Major financial problem","Neighborhood disadvantage")
for (i in adv.list){
  print(i)
  print(table(cpg_r2_g003[cpg_r2_g003$measure==i,"hypothesis"]))
}


### Take a look at the top 10 hits for each SEP

for (adv in adv.list) {
  df <- sI.summary[sI.summary$measure == adv,]
  df <- df [order(df$`p.value`),]
  df.top10 <- df[1:10,]
  print(df.top10)
  
}



# Create a dataset that contains the original SLCMA output for the top 62 hits with R2>0.03.

pull.raw<- function(adver){
  df<-get(adver)
  df <- df[df$'1.r2'>0.03,c("Probe","1.hypo","1.r2","1.p")]
  df$adversity <-strsplit(adver, "[.]")[[1]][2]
  return(df)
}

cpg_r2_g003_raw <- lapply(results.list, pull.raw)
cpg_r2_g003_raw <- data.frame(do.call("rbind", cpg_r2_g003_raw))

save(cpg_r2_g003_raw,file='/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/summary_r2g003_62hits_raw.RData')

# Save all data
save(results.Fscore.sI,results.income.sI, results.income_reduction.sI,results.jobless.sI, results.major_F.sI, results.nbhqual.sI,file='/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/results_all.RData')
