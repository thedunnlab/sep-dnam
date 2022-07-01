#####################################
# Examine mobility-excluded results
# 06-23-2020
# Jiaxuan (Jessie) Liu
#####################################
## load packages
library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(dplyr)

## Load results
load("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/newDNAm/FWL_Results/Ex_mobility/results_Fscore_period_exmob_062320_sI.Rdata")
load("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/newDNAm/FWL_Results/Ex_mobility/results_income_period_exmob_062320_sI.Rdata")
load("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/newDNAm/FWL_Results/Ex_mobility/results_major_F_period_exmob_062320_sI.Rdata")
load("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/newDNAm/FWL_Results/Ex_mobility/results_nbhqual_period_exmob_062320_sI.Rdata")

head(results.Fscore.sI)
results.Fscore.sI <- results.Fscore.sI[order(results.Fscore.sI$`1.p`),]
head(results.income.sI)
results.income.sI <- results.income.sI[order(results.income.sI$`1.p`),]
head(results.major_F.sI)
results.major_F.sI <- results.major_F.sI[order(results.major_F.sI$`1.p`),]
head(results.nbhqual.sI)
results.nbhqual.sI <- results.nbhqual.sI[order(results.nbhqual.sI$`1.p`),]


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



# Remove the probes from the results files
results.Fscore.sI <-results.Fscore.sI[!(rownames(results.Fscore.sI) %in% sex.chr | rownames(results.Fscore.sI) %in% hybrid),]
results.income.sI <- results.income.sI[!(rownames(results.income.sI) %in% sex.chr | rownames(results.income.sI) %in% hybrid),]
results.major_F.sI <- results.major_F.sI[!(rownames(results.major_F.sI) %in% sex.chr | rownames(results.major_F.sI) %in% hybrid),]
results.nbhqual.sI <- results.nbhqual.sI[!(rownames(results.nbhqual.sI) %in% sex.chr | rownames(results.nbhqual.sI) %in% hybrid),]
dim(results.nbhqual.sI) #412956



#Bonferroni
dim(results.Fscore.sI[results.Fscore.sI$`1.p`<0.05/412956,])#1
dim(results.income.sI[results.income.sI$`1.p`<0.05/412956,])#0
dim(results.major_F.sI[results.major_F.sI$`1.p`<0.05/412956,])#0
dim(results.nbhqual.sI[results.nbhqual.sI$`1.p`<0.05/412956,])#1

#FDR
fdr.Fscore<-p.adjust(results.Fscore.sI$`1.p`, method = "BH")
sum(fdr.Fscore<0.1)#10
sum(fdr.Fscore<0.2)#122
sum(fdr.Fscore<0.3)#904
fdr.income<-p.adjust(results.income.sI$`1.p`, method = "BH")
sum(fdr.income<0.1)#0
sum(fdr.income<0.2)#0
sum(fdr.income<0.3)#13
fdr.major_F<-p.adjust(results.major_F.sI$`1.p`, method = "BH")
sum(fdr.major_F<0.1)#0
sum(fdr.major_F<0.2)#0
sum(fdr.major_F<0.3)#0
fdr.nbhqual<-p.adjust(results.nbhqual.sI$`1.p`, method = "BH")
sum(fdr.nbhqual<0.1)#1
sum(fdr.nbhqual<0.2)#39
sum(fdr.nbhqual<0.3)#469


# Make a new summary object using the updated results file
sI.summary <- rbind(data.frame(Probe = results.Fscore.sI$Probe, hypothesis = results.Fscore.sI$`1.hypo`, 
                               r2 = results.Fscore.sI$`1.r2`, p.value = results.Fscore.sI$`1.p`, FDR_q= p.adjust(results.Fscore.sI$`1.p`, method = "BH"), 
                               expected = c(1:412956)/412956,measure = "Financial difficulty score"),
                    
                    data.frame(Probe = results.income.sI$Probe, hypothesis = results.income.sI$`1.hypo`, 
                               r2 = results.income.sI$`1.r2`, p.value = results.income.sI$`1.p`, FDR_q= p.adjust(results.income.sI$`1.p`, method = "BH"), 
                               expected = c(1:412956)/412956, measure = "Income"),
                  
                    data.frame(Probe = results.major_F.sI$Probe, hypothesis = results.major_F.sI$`1.hypo`, 
                               r2 = results.major_F.sI$`1.r2`, p.value = results.major_F.sI$`1.p`, FDR_q= p.adjust(results.major_F.sI$`1.p`, method = "BH"), 
                               expected = c(1:412956)/412956,  measure = "Major financial problem"),
                    
                    data.frame(Probe = results.nbhqual.sI$Probe, hypothesis = results.nbhqual.sI$`1.hypo`, 
                               r2 = results.nbhqual.sI$`1.r2`, p.value = results.nbhqual.sI$`1.p`, FDR_q= p.adjust(results.nbhqual.sI$`1.p`, method = "BH"), 
                               expected = c(1:412956)/412956, measure = "Neighborhood disadvantage"))


summary(as.factor(sI.summary$hypothesis))
sI.summary$hypothesis <- gsub("Fscore_newbin_","", sI.summary$hypothesis)
sI.summary$hypothesis <- gsub("income_bin_","", sI.summary$hypothesis)
sI.summary$hypothesis <- gsub("major_F_newbin_","", sI.summary$hypothesis)
sI.summary$hypothesis <- gsub("nbhqual_newbin_","", sI.summary$hypothesis)

sI.summary$hypothesis <- factor(sI.summary$hypothesis, levels = c("accumulation", "very_early","early","middle"))



#total #probes passing Bonferroni
cpg_bonf<-sI.summary[sI.summary$'p.value'<0.05/412956,]#2
#total #probes passing FDR q<0.1
cpg_fdr01<-sI.summary[sI.summary$FDR_q<0.1,]#11

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

ggplot(sI.summary[which(sI.summary$FDR_q <0.1),], aes(x= hypothesis, fill= hypothesis))+
  geom_bar(stat='count', position=position_dodge(preserve = 'single'), width =0.5, color = 'grey')+
  xlab("Measure")+
  ylab("# of CpGs under hypothesis")+ 
  facet_wrap(~measure)+
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        axis.text.x = element_text(angle = 45, hjust= 1))+
  #scale_fill_brewer(palette = "RdYlBu")+
  scale_fill_manual(values = c("#a50026","#f46d43","#fee090"))+
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
  scale_fill_manual(values = c("#a50026","#f46d43","#fee090"))+
  #scale_fill_got(option = 'Tyrell', discrete=T)+
  ggtitle("FDR q<0.2")+
  geom_text(stat= 'count', aes(label = ..count..), position = position_dodge(width =1), vjust = -0.5)+
  ylim(0,115)

ggplot(sI.summary[which(sI.summary$FDR_q <0.3),], aes(x= hypothesis, fill= hypothesis))+
  geom_bar(stat='count', position=position_dodge(preserve = 'single'), width =0.5, color = 'grey')+
  xlab("Measure")+
  ylab("# of CpGs under hypothesis")+ 
  facet_wrap(~measure)+
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        axis.text.x = element_text(angle = 45, hjust= 1))+
 # scale_fill_brewer(palette = "RdYlBu")+
  scale_fill_manual(values = c("#a50026","#f46d43","#fdae61","#fee090"))+
  ggtitle("FDR q<0.3")+
  geom_text(stat= 'count', aes(label = ..count..), position = position_dodge(width =1), vjust = -0.5)+
  ylim(0,800)

ggplot(sI.summary[which(sI.summary$'p.value' <1/412956),], aes(x= hypothesis, fill= hypothesis))+
  geom_bar(stat='count', position=position_dodge(preserve = 'single'), width =0.5, color = 'grey')+
  xlab("Measure")+
  ylab("# of CpGs under hypothesis")+ 
  facet_wrap(~measure)+
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        axis.text.x = element_text(angle = 45, hjust= 1))+
  #scale_fill_brewer(palette = "RdYlBu")+
  scale_fill_manual(values = c("#a50026","#f46d43","#fee090"))+
  #scale_fill_got(option = 'Tyrell', discrete=T)+
  ggtitle("P value <2.4e-6 (Bonferroni)")+
  geom_text(stat= 'count', aes(label = ..count..), position = position_dodge(width =1), vjust = -0.5)+
  ylim(0,10)


##### Do we get the same hits as when mobility was included?

# load significant results

load('/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/summary_r2g003_62hits.RData')

# Neighborhood quality

nbhql_r2<-cpg_r2_g003[which(cpg_r2_g003$measure=="Neighborhood disadvantage"),] #17
dim(sI.summary[which(sI.summary$Probe %in% nbhql_r2$Probe & sI.summary$measure=="Neighborhood disadvantage" & sI.summary$r2>0.03),])[1] # 6

# Fscore

fscore_r2<-cpg_r2_g003[which(cpg_r2_g003$measure=="Financial difficulty score"),] #9
dim(sI.summary[which(sI.summary$Probe %in% fscore_r2$Probe & sI.summary$measure=="Financial difficulty score"& sI.summary$r2>0.03),])[1] # 9

# Income
income_r2<-cpg_r2_g003[which(cpg_r2_g003$measure=="Income"),] #13
dim(sI.summary[which(sI.summary$Probe %in%income_r2$Probe & sI.summary$measure=="Income"& sI.summary$r2>0.03),])[1] #8

# Major F
majorF_r2<-cpg_r2_g003[which(cpg_r2_g003$measure=="Major financial problem"),] #5
dim(sI.summary[which(sI.summary$Probe %in%majorF_r2$Probe & sI.summary$measure=="Major financial problem"& sI.summary$r2>0.03),])[1]  #1 


#### Export Ex-Mobility results for the main analysis hits

a <- sI.summary[which(sI.summary$Probe %in% nbhql_r2$Probe & sI.summary$measure=="Neighborhood disadvantage"),] 
b <- sI.summary[which(sI.summary$Probe %in% fscore_r2$Probe & sI.summary$measure=="Financial difficulty score"),]
c <- sI.summary[which(sI.summary$Probe %in%income_r2$Probe & sI.summary$measure=="Income"),]
d <- sI.summary[which(sI.summary$Probe %in%majorF_r2$Probe & sI.summary$measure=="Major financial problem"),] 

ex_mob_res<- rbind(a,b,c,d)
ex_mob_res$adversity <- as.character(ex_mob_res$measure)
ex_mob_res[ex_mob_res$measure=="Neighborhood disadvantage","adversity"] <- "nbhqual"
ex_mob_res[ex_mob_res$measure=="Financial difficulty score","adversity"] <- "Fscore"
ex_mob_res[ex_mob_res$measure=="Income","adversity"] <- "lowincome"
ex_mob_res[ex_mob_res$measure=="Major financial problem","adversity"] <- "major_F"

## Calculate beta values
source("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Codes/newDNAm/FWL/LARS-mobility-FWL-CI-20200722.R")

# load raw data
load("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/newDNAm/newDNAm_pheno.Rdata")
colnames(df)[11:13]<-c("lowincome_newbin_very_early","lowincome_newbin_early" ,"lowincome_newbin_middle") #fix the inconsistent variable naming

adv.list <- c("lowincome", "Fscore","major_F", "nbhqual")

lars.fun <- function(adver){ 
  
  outcome.vec <- as.character(ex_mob_res[ex_mob_res$adversity==adver,"Probe"])
     lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood"),
                                          hypos=c("accumulation"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE)
 
  lars.out$adversity<-adver
  return(lars.out)
}


tableS_exMob<- data.frame(do.call("rbind", lapply(adv.list, lars.fun)))

colnames(tableS_exMob)[3]<-"Beta_SI"

save(tableS_exMob, file = "/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/TableS_ex_Mob_2021-05-17.Rdata")
write.csv(tableS_exMob,file='/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/TableS_ex_Mob_2021-05-17.csv', row.names = FALSE)

# Merge with the main table
load("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/Table2/Table2_62hits_annot_2021-05-17.Rdata")
colnames(tableS_exMob)[2:9] <- paste0 (colnames(tableS_exMob[2:9]),"_exMob")
tableS_exMob$Probe <- as.character(tableS_exMob$Probe)
tableS_exMob_combined <- merge(table2.final.annot[,c("CpG","adversity","hypo","Beta_SI","SE","Pvalue","R2")],tableS_exMob,all.y=T,by.y="Probe",by.x="CpG")

write.csv(tableS_exMob_combined,file='/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/TableS_ex_Mob_2021-05-17_combined.csv', row.names = FALSE)

