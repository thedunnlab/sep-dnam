# Sensitivity Analyses 
library(xlsx)
library(tidyverse)
library(dplyr)

source("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Codes/newDNAm/FWL/LARS-mobility-FWL-CI-20200722.R")

# load top hits
load('/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/summary_r2g003_62hits_raw.RData')
dim(cpg_r2_g003_raw)

cpg_r2_g003_raw$'X1.hypo' <- as.character(cpg_r2_g003_raw$'X1.hypo')
table(cpg_r2_g003_raw[cpg_r2_g003_raw$adversity=="income",'X1.hypo'])
cpg_r2_g003_raw[cpg_r2_g003_raw$adversity=="income","adversity"]<-"lowincome" #fix the naming issue with income variable (overlap with income_reduction; if not fix this will get incorrect sample size for income; in the original SLCMA analysis I entered the exact variable names, i.e. "income_bin" and "income_reduction_newbin" so not affected by this)
cpg_r2_g003_raw[cpg_r2_g003_raw$'X1.hypo'=="income_bin_very_early",'X1.hypo']<-"lowincome_newbin_very_early" #fix the naming issue 
cpg_r2_g003_raw[cpg_r2_g003_raw$'X1.hypo'=="income_bin_early",'X1.hypo']<-"lowincome_newbin_early" #fix the naming issue 
cpg_r2_g003_raw[cpg_r2_g003_raw$'X1.hypo'=="income_bin_middle",'X1.hypo']<-"lowincome_newbin_middle" #fix the naming issue 
table(cpg_r2_g003_raw[cpg_r2_g003_raw$adversity=="lowincome",'X1.hypo'])

# read in the order of the CpGs that we want to arrange the output files
cpg.order <- read.csv(file = '/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/Table2/62hits_order.csv')

# load raw data (DNAm + Phenotypes)
load("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/newDNAm/newDNAm_pheno.Rdata")
colnames(df)[11:13]<-c("lowincome_newbin_very_early","lowincome_newbin_early" ,"lowincome_newbin_middle") #fix the inconsistent variable naming

df.hits <- df[,colnames(df) %in% cpg_r2_g003_raw$Probe]
df.hits <- cbind(df[,1:67],df.hits)
df<-df.hits

# load indivdiual time points data 
load("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/SES_period_040419.Rdata")


adv.list <- c("jobless", "income_reduction",   "lowincome", "Fscore","major_F", "nbhqual")

# test the function of primary analysis
lars.fun <- function(adver){ 
  
  outcome.vec <- cpg_r2_g003_raw[cpg_r2_g003_raw$adversity==adver,"Probe"]
  
  if (adver %in% c("jobless", "income_reduction")) {
    lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood"),
                                          hypos=c("accumulation"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE) 
  }
  else {
    lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood"),
                                          hypos=c("accumulation","mobility"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE)
  } 
  lars.out$adversity<-adver
  return(lars.out)
}

table2.ci<- data.frame(do.call("rbind", lapply(adv.list, lars.fun)))

table2.ci <- table2.ci[match(cpg.order$CpG, table2.ci$Probe),] # order the CpGs 


########### Sensitivity Analysis: Adjust for additional covariates ###########


#further adjust for home_owner (created by Miriam )

lars.fun.home <- function(adver){ 
  
  outcome.vec <- cpg_r2_g003_raw[cpg_r2_g003_raw$adversity==adver,"Probe"]
  
  if (adver %in% c("jobless", "income_reduction")) {
    lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood","home_owner"),
                                          hypos=c("accumulation"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE) 
  }
  else {
    lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood","home_owner"),
                                          hypos=c("accumulation","mobility"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE)
  } 
  lars.out$adversity<-adver
  return(lars.out)
}

table(df$home_owner)
df$home_owner <- ifelse(str_trim(df$home_owner)==".",NA,df$home_owner)
df$home_owner<-as.factor(df$home_owner)
sens.home<- data.frame(do.call("rbind", lapply(adv.list, lars.fun.home)))
sens.home <- sens.home[match(cpg.order$CpG, sens.home$Probe),] # order the CpGs 


table2.ci$Variable==sens.home$Variable #Same hypothesis is selected for all probes
1-sum(is.na(df$home_owner))/dim(df)[1]#0.973
write.xlsx(sens.home, file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/RnR/sensitivity_62hits_2022-06-04.xlsx", sheetName="home owner", row.names=TRUE) 



#further adjust for marital status

lars.fun.marital <- function(adver){ 
  
  outcome.vec <- cpg_r2_g003_raw[cpg_r2_g003_raw$adversity==adver,"Probe"]
  
  if (adver %in% c("jobless", "income_reduction")) {
    lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood","marital"),
                                          hypos=c("accumulation"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE) 
  }
  else {
    lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood","marital"),
                                          hypos=c("accumulation","mobility"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE)
  } 
  lars.out$adversity<-adver
  return(lars.out)
}

table(df$marital)
df$marital <- ifelse(str_trim(df$marital)==".",NA,df$marital)
df$marital<-as.factor(df$marital)
sens.marital<- data.frame(do.call("rbind", lapply(adv.list, lars.fun.marital)))
sens.marital<- sens.marital[match(cpg.order$CpG, sens.marital$Probe),] # order the CpGs 

table2.ci$Variable==sens.marital$Variable #Same hypothesis is selected for all probes
1-sum(is.na(df$marital))/dim(df)[1]#0.984
write.xlsx(sens.marital, file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/RnR/sensitivity_62hits_2022-06-04.xlsx", sheetName="marital", append=TRUE, row.names=TRUE) 

#further adjust for maternal education 

lars.fun.edu <- function(adver){ 
  
  outcome.vec <- cpg_r2_g003_raw[cpg_r2_g003_raw$adversity==adver,"Probe"]
  
  if (adver %in% c("jobless", "income_reduction")) {
    lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood","ed_momgest"),
                                          hypos=c("accumulation"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE) 
  }
  else {
    lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood","ed_momgest"),
                                          hypos=c("accumulation","mobility"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE)
  } 
  lars.out$adversity<-adver
  return(lars.out)
}

table(df$ed_momgest)
df$ed_momgest <- ifelse(str_trim(df$ed_momgest)==".",NA,df$ed_momgest)
df$ed_momgest<-as.factor(df$ed_momgest)
sens.edu<- data.frame(do.call("rbind", lapply(adv.list, lars.fun.edu)))
sens.edu<- sens.edu[match(cpg.order$CpG, sens.edu$Probe),] # order the CpGs 


table2.ci$Variable==sens.edu$Variable #Same hypothesis is selected for all probes
1-sum(is.na(df$ed_momgest))/dim(df)[1]#0.979
write.xlsx(sens.edu, file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/RnR/sensitivity_62hits_2022-06-04.xlsx", sheetName="education", append=TRUE, row.names=TRUE) 


#further adjust for Townsend Deprivation Index 

df<-merge(df,SES.period[,c("cidB1471","Townsend_pre32")],by="cidB1471")

lars.fun.town <- function(adver){ 
  
  outcome.vec <- cpg_r2_g003_raw[cpg_r2_g003_raw$adversity==adver,"Probe"]
  
  if (adver %in% c("jobless", "income_reduction")) {
    lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood","Townsend_pre32"),
                                          hypos=c("accumulation"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE) 
  }
  else {
    lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood","Townsend_pre32"),
                                          hypos=c("accumulation","mobility"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE)
  } 
  lars.out$adversity<-adver
  return(lars.out)
}

df$Townsend_pre32<-as.factor(df$Townsend_pre32)

sens.town<- data.frame(do.call("rbind", lapply(adv.list, lars.fun.town)))

sens.town<- sens.town[match(cpg.order$CpG, sens.town$Probe),] # order the CpGs 


table2.ci$Variable==sens.town$Variable #Same hypothesis is selected for all probes
1-sum(is.na(df$Townsend_pre32))/dim(df)[1] #0.897
write.xlsx(sens.town, file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/RnR/sensitivity_62hits_2022-06-04.xlsx", sheetName="Townsend", append=TRUE, row.names=TRUE) 



#further adjust for ever-homeless (variable created by Brigette)


df<-merge(df,SES.period[,c("cidB1471","ever_homeless")],by="cidB1471")

lars.fun.homeless <- function(adver){ 
  
  outcome.vec <- cpg_r2_g003_raw[cpg_r2_g003_raw$adversity==adver,"Probe"]
  
  if (adver %in% c("jobless", "income_reduction")) {
    lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood","ever_homeless"),
                                          hypos=c("accumulation"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE) 
  }
  else {
    lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood","ever_homeless"),
                                          hypos=c("accumulation","mobility"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE)
  } 
  lars.out$adversity<-adver
  return(lars.out)
}

table(df$ever_homeless)
df$ever_homeless<-as.factor(df$ever_homeless)
sens.homeless<- data.frame(do.call("rbind", lapply(adv.list, lars.fun.homeless)))

sens.homeless <- sens.homeless [match(cpg.order$CpG, sens.homeless$Probe),] # order the CpGs 


table2.ci$Variable==sens.homeless$Variable #Same hypothesis is selected for all probes
1-sum(is.na(df$ever_homeless))/dim(df)[1]#0.782

write.xlsx(sens.homeless, file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/RnR/sensitivity_62hits_2022-06-04.xlsx", sheetName="homeless(adj)", append=TRUE, row.names=TRUE) 



#further adjust for all stable SEP indicators 


lars.fun.all <- function(adver){ 
  
  outcome.vec <- cpg_r2_g003_raw[cpg_r2_g003_raw$adversity==adver,"Probe"]
  
  if (adver %in% c("jobless", "income_reduction")) {
    lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood",
                                                                                 "ever_homeless","home_owner","marital","ed_momgest","Townsend_pre32"),
                                          hypos=c("accumulation"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE) 
  }
  else {
    lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood",
                                                                                 "ever_homeless","home_owner","marital","ed_momgest","Townsend_pre32"),
                                          hypos=c("accumulation","mobility"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE)
  } 
  lars.out$adversity<-adver
  return(lars.out)
}


sens.all<- data.frame(do.call("rbind", lapply(adv.list, lars.fun.all)))

sens.all <- sens.all [match(cpg.order$CpG, sens.all$Probe),] # order the CpGs 



table2.ci$Variable==sens.all$Variable #Same hypothesis is selected for all probes
1-sum(is.na(df$ever_homeless) | is.na(df$home_owner) | is.na(df$marital) | is.na(df$ed_momgest) | is.na(df$Townsend_pre32))/dim(df)[1] #0.700

write.xlsx(sens.all, file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/RnR/sensitivity_62hits_2022-06-04.xlsx", sheetName="all inv SEP", append=TRUE, row.names=TRUE) 





#further adjust for epigenetic PCs
lars.fun.ePC <- function(adver){ 
  
  outcome.vec <- cpg_r2_g003_raw[cpg_r2_g003_raw$adversity==adver,"Probe"]
  
  if (adver %in% c("jobless", "income_reduction")) {
    lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood","EWAS.PC.1","EWAS.PC.2",
                                                                                 "EWAS.PC.3","EWAS.PC.4"),
                                          hypos=c("accumulation"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE) 
  }
  else {
    lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood","EWAS.PC.1","EWAS.PC.2",
                                                                                 "EWAS.PC.3","EWAS.PC.4"),
                                          hypos=c("accumulation","mobility"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE)
  } 
  lars.out$adversity<-adver
  return(lars.out)
}

df$'EWAS.PC.1'<-as.numeric(as.character(df$'EWAS.PC.1')) # PC 1 is a factor in the original data
sens.ePC<- data.frame(do.call("rbind", lapply(adv.list, lars.fun.ePC)))
sens.ePC <- sens.ePC [match(cpg.order$CpG, sens.ePC$Probe),] # order the CpGs 


table2.ci$Variable==sens.ePC$Variable #Same hypothesis is selected for all probes
write.xlsx(sens.ePC, file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/RnR/sensitivity_62hits_2022-06-04.xlsx", sheetName="4 ePCs", append=TRUE, row.names=TRUE) 




#further adjust for cord blood DNAm

load("/data/js95/ALSPAC/ARIES/DNAm_2020/Cord/betas_Cord_WIN_20200128.Rdata")
dim(betas.cord.win) #N=905
betas_cord_hits <-betas.cord.win[,colnames(betas.cord.win) %in% cpg_r2_g003_raw$Probe]

samplesheet<-read.table("/data/js95/ALSPAC/ARIES/DNAm_2020/Cord/samplesheet_Cord_20200121.txt",header=T)

cell<-read.table("/data/js95/ALSPAC/ARIES/DNAm_2020/Cord/celltypes_Cord_20200121.txt",header=T)

sample<-merge(samplesheet,betas_cord_hits,by.x="Sample_Name",by.y="row.names")
sample<-merge(cell,sample,by.x="IID",by.y="Sample_Name")

summary(rowSums(sample[,c("Bcell","CD4T","CD8T","Gran","Mono","NK")]))
summary(rowSums(sample[,c("Bcell","CD4T","CD8T","Gran","Mono","NK","nRBC")])) #should include nRBC

# regress the betas on the cell counts first

cord_betas_adj <- do.call("cbind", lapply(cpg_r2_g003_raw$Probe, function(x){
  r <- lm(as.formula(paste0(x, "~ CD8T + CD4T + Bcell + 
                              Mono + NK + Gran + nRBC")), data = sample)$resid
  r <- as.data.frame(r)
  colnames(r) <- paste0("cord_",x)
  return(r)}))  

cord_betas_adj<-cbind(sample[,"ALN"],cord_betas_adj)
colnames(cord_betas_adj)[1] <-"cidB1471"

df<-merge(df,cord_betas_adj,by="cidB1471")


#need to update this
lars.fun.cord <- function(cpg){ 
  
  outcome.vec <- cpg
  adver <-cpg_r2_g003_raw[cpg_r2_g003_raw$Probe==cpg,"adversity"]
  cord_DNAm <- paste0("cord_",cpg)
  
  if (adver %in% c("jobless", "income_reduction")) {
    lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood",cord_DNAm),
                                          hypos=c("accumulation"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE) 
  }
  else {
    lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood",cord_DNAm),
                                          hypos=c("accumulation","mobility"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE)
  } 
  lars.out$adversity<-adver
  return(lars.out)
}


sens.cord<- data.frame(do.call("rbind", lapply(cpg_r2_g003_raw$Probe, lars.fun.cord)))

table2.ci$Probe <- as.character(table2.ci$Probe)
sens.cord$Probe <- as.character(sens.cord$Probe)

sens.cord <- sens.cord [match(cpg.order$CpG, sens.cord$Probe),] # order the CpGs 


table2.ci$Variable==sens.cord$Variable #Same hypothesis is selected for all probes

write.xlsx(sens.cord, file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/RnR/sensitivity_62hits_2022-06-04.xlsx", sheetName="Cord blood", append=TRUE, row.names=TRUE) 


########### Sensitivity Analysis: Adjust for folate intake (R&R) ###########
# read in folate data 
# b143: taking folic acid during this pregnancy (18 weeks gestation) 1=Yes, 2=No
# bc113: taken folic acid in the last 3 months (32 weeks gestation) 1=Yes, 2=No

load("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Folate/folate_data_2022-05-31.Rdata")
dim(folate) #15646 samples
head(folate)

df<-df.hits
df<-merge(df,folate[,c("cidB1471","DEL_P1591","b143","c113")],by="cidB1471")
dim(df) #946

table(df$qlet)

table(df$DEL_P159,useNA = "always") #358/946=37.8% missing
table(df$b143,useNA = "always") #23/946=2.4% missing
table(df$c113,useNA = "always") #37/946=3.9% missing

table(df$b143,df$c113,useNA = "always") 
class(df$b143) #integer
class(df$c113) #integer

df$folate_any <- ifelse ((df$b143==1 | df$c113==1),1,ifelse((df$b143==2 | df$c113==2),0,NA))
table(df$folate_any,useNA = "always") 
prop.table(table(df$folate_any,useNA = "always")) #5.07% missing, 23.2% prevalent

lars.fun.folate <- function(adver){ 
  
  outcome.vec <- cpg_r2_g003_raw[cpg_r2_g003_raw$adversity==adver,"Probe"]
  
  if (adver %in% c("jobless", "income_reduction")) {
    lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood","folate_any"),
                                          hypos=c("accumulation"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE) 
  }
  else {
    lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood","folate_any"),
                                          hypos=c("accumulation","mobility"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE)
  } 
  lars.out$adversity<-adver
  return(lars.out)
}


sens.folate<- data.frame(do.call("rbind", lapply(adv.list, lars.fun.folate)))
sens.folate <- sens.folate [match(cpg.order$CpG, sens.folate$Probe),] # order the CpGs 
sens.folate$Pvalue<-as.numeric(sens.folate$Pvalue)
dim(sens.folate[sens.folate$Pvalue>0.05,]) # all had p<0.05

table2.ci$Variable==sens.folate$Variable #Same hypothesis is selected for all probes
write.xlsx(sens.folate, file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/RnR/sensitivity_62hits_2022-06-04.xlsx", sheetName="folate_any", append=TRUE, row.names=TRUE) 




########### Sensitivity Analysis: Adjust for exposure to other exposure SEP indicators (R&R) ###########
# I went back to the following script to see how I created the ever exposure variables for Table 1:
# /data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Codes/Manuscript/Table1_FigureS1_SEP_descriptives_2021-04-25.R
# I followed Alex's coding of "ever exposure": /ALSPAC- Russell Sage/SES Measurement/Jessie Data/Team paper drafts/Alex Aim 3/Aim3_mutually_adjusted_2022-05-31.Rmd
# Having at least one exposure (regardless of other missing values) -> coded as exposed; no missing values and all unexposed -> coded as unexposed; otherwise missing. 

adv.list <- c("jobless", "income_reduction",   "lowincome", "Fscore","major_F", "nbhqual")

#create ever exposure indicators 
for (exp in adv.list) {
  exp_ever <- paste0(exp,"_ever")
  df$newcol<- rep(NA,dim(df)[1])
  
  df[,"newcol"] <- ifelse(rowSums(df[,c(paste0(exp,"_newbin_very_early"),paste0(exp,"_newbin_early"),paste0(exp,"_newbin_middle"))], na.rm=T)>0, 1,
                             ifelse(rowSums(df[,c(paste0(exp,"_newbin_very_early"),paste0(exp,"_newbin_early"),paste0(exp,"_newbin_middle"))]) ==0, 0, NA))
    
  names(df)[names(df) == "newcol"] <- exp_ever

}

#create indicators of exposure to any other adversity (RESUME HERE)

for (i in 1:6) {
  exp <- adv.list[i]
  exp_ever <- paste0(exp,"_ever_otherSEPs")
  df$newcol<- rep(NA,dim(df)[1])
  
  df[,"newcol"] <- ifelse(rowSums(df[,paste0(adv.list[adv.list != exp],"_ever")], na.rm=T)>0, 1,
                          ifelse(rowSums(df[,paste0(adv.list[adv.list != exp],"_ever")]) ==0, 0, NA))
  
  names(df)[names(df) == "newcol"] <- exp_ever
  
}

table(df$jobless_ever_otherSEPs,useNA = "always")
table(df$income_reduction_ever_otherSEPs,useNA = "always")

summary_everexp <-matrix(ncol=4, nrow=12) 
colnames(summary_everexp)<-c("exp", "prop.missing", "n.avai","prevalence") 

for (i in 1:6) {
  summary_everexp[i,"exp"] <- adv.list[i]
  exp_ever <- paste0(adv.list[i],"_ever")
  summary_everexp[i,"prop.missing"] <- round(prop.table(table(df[,exp_ever] ,useNA = "always"))[3],4)
  summary_everexp[i,"n.avai"] <- sum(!is.na(df[,exp_ever]))
  summary_everexp[i,"prevalence"] <- round(prop.table(table(df[,exp_ever]))[2],4)
  
} # the numbers are not exact the same as in Table 1 (prevalence of ever exposure is slightly higher here because we used less stringent coding)


for (i in 1:6) {
  summary_everexp[i+6,"exp"] <- paste0(adv.list[i],"_ever_otherSEPs")
  exp_ever <- paste0(adv.list[i],"_ever_otherSEPs")
  summary_everexp[i+6,"prop.missing"] <- round(prop.table(table(df[,exp_ever] ,useNA = "always"))[3],4)
  summary_everexp[i+6,"n.avai"] <- sum(!is.na(df[,exp_ever]))
  summary_everexp[i+6,"prevalence"] <- round(prop.table(table(df[,exp_ever]))[2],4)
  
} 

write.xlsx(summary_everexp, file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/RnR/ever_exp_summary_2022-06-07.xlsx") 


# mutually adjusted analysis

lars.fun.mutadj <- function(adver){ 
  
  outcome.vec <- cpg_r2_g003_raw[cpg_r2_g003_raw$adversity==adver,"Probe"]
  addcovr <- paste0(adver,"_ever_otherSEPs")
  
  if (adver %in% c("jobless", "income_reduction")) {
    lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood",addcovr),
                                          hypos=c("accumulation"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE) 
  }
  else {
    lars.out<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood",addcovr),
                                          hypos=c("accumulation","mobility"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE)
  } 
  lars.out$adversity<-adver
  return(lars.out)
}


sens.mutadj<- data.frame(do.call("rbind", lapply(adv.list, lars.fun.mutadj)))
sens.mutadj <- sens.mutadj [match(cpg.order$CpG, sens.mutadj$Probe),] # order the CpGs 


table2.ci$Variable==sens.mutadj$Variable #Same hypothesis is selected for all probes

sens.mutadj$Pvalue<-as.numeric(sens.mutadj$Pvalue)
dim(sens.mutadj[sens.mutadj$Pvalue>0.05,]) # all had p<0.05


write.xlsx(sens.mutadj , file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/RnR/sensitivity_62hits_2022-06-04.xlsx", sheetName="mutual adj-any exp", append=TRUE, row.names=TRUE) 
