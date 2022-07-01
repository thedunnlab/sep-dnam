## Sensitivty analysis after controlling for genotypes 
## for top CpGs linked to mQTLs 
## Adapted from Yiwen's script:/Dropbox (Partners HealthCare)/Results Yiwen/BiologicalPsychiatry-2019/scripts/rev-2/SuppTables/TableS5_Sensitivity_mQTL_20181120.R


## packages 

## load Bonferroni significant CpGs 
load('/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/summary_r2g003_62hits_raw.RData')
dim(cpg_r2_g003_raw)


cpg_r2_g003_raw$'X1.hypo' <- as.character(cpg_r2_g003_raw$'X1.hypo')
table(cpg_r2_g003_raw[cpg_r2_g003_raw$adversity=="income",'X1.hypo'])
cpg_r2_g003_raw[cpg_r2_g003_raw$adversity=="income","adversity"]<-"lowincome" #fix the naming issue with income variable (overlap with income_reduction; if not fix this will get incorrect sample size for income; in the original SLCMA analysis I entered the exact variable names, i.e. "income_bin" and "income_reduction_newbin" so not affected by this)
cpg_r2_g003_raw[cpg_r2_g003_raw$'X1.hypo'=="income_bin_very_early",'X1.hypo']<-"lowincome_newbin_very_early" #fix the naming issue 
cpg_r2_g003_raw[cpg_r2_g003_raw$'X1.hypo'=="income_bin_early",'X1.hypo']<-"lowincome_newbin_early" #fix the naming issue 
cpg_r2_g003_raw[cpg_r2_g003_raw$'X1.hypo'=="income_bin_middle",'X1.hypo']<-"lowincome_newbin_middle" #fix the naming issue 
table(cpg_r2_g003_raw[cpg_r2_g003_raw$adversity=="lowincome",'X1.hypo'])


F7.mqtl <- read_delim(file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/F7.ALL.M.tab",delim = "\t")

# filter to 62 CpGs:
F7.mqtl.subset <- F7.mqtl %>% filter(gene %in% cpg_r2_g003_raw$Probe)
rm(F7.mqtl)
# how many associations?
nrow(F7.mqtl.subset) # 2523
# how many CpGs?
length(unique(F7.mqtl.subset$gene)) # 24 
# how many SNPs?
length(unique(F7.mqtl.subset$SNP)) # 2523
## length(unique(F7.mqtl.subset$SNP)) SNPs associated with DNAm at 8 CpGs (presumably both cis and trans)



# export SNP list to lookup in PLINK (not doing GRASP for this table yet)
F7.mqtl.subset$SNP <- as.character(F7.mqtl.subset$SNP)
F7.mqtl.subset$gene <- as.character(F7.mqtl.subset$gene)

F7.mqtl.subset <- F7.mqtl.subset[order(F7.mqtl.subset$SNP), ]

save(F7.mqtl.subset, file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/R2_CpG_mQTLs_2021-05-17.Rdata")
# plink:
write(F7.mqtl.subset$SNP, sep="", file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/R2CpGmQTLSNPsforPLINK_2021-05-17.txt")

## IDs for plink
## load pheno data
load("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/period_phenotypes_050719.Rdata")
pheno.period$ID <- paste(pheno.period$cidB1471,pheno.period$qlet,sep="")
write(pheno.period$ID, sep="", file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/ARIES.IDs_forPLINK_2021-05-17.txt")


## now cut the dosage in plink
## see /data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Codes/Manuscript/mQTL_dosage_2021-05-17.sh


#### START HERE for adjusted analysis #### 
library(dplyr)
library(tibble)


# for each CpG, need to extract mQTL SNP names, extract SNP dosages,
# and then run the regression for that particular adversity, controlling
# for the dosages.



# subset list to CpGs affected by mQTLs:
CpG.mQTL <- cpg_r2_g003_raw[cpg_r2_g003_raw$Probe %in% unique(F7.mqtl.subset$gene), ] 



## analysis controlling for genotypes
## loading dosage data downloaded from cluster:
chr <- c('01',"02", "03", "05", "06","07", "09", "10","11","12","13","16","17","20","21") # only available at these chromosomes

load.dat <- function(x){
  mfile <- paste0("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/mQTLdosages/ARIES.F7_dosages_2021-05-17-chr", x, ".raw")
  dat <- read.delim(file=mfile, 
                    sep="", header=T, stringsAsFactors = F)
  return(dat)
}
mylist <- lapply(chr, load.dat)
rm(load.dat)

# how many subjects with genotype data?
sapply(mylist, nrow) # 893 
dosages <- Reduce(merge, mylist)
dim(dosages) # 893 x 2529 (2523+6)
dosages <- dosages[, -c(3,4,6)] # rm uncessary vars
rm(mylist, chr)
head(dosages)
save(dosages, file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/ARIES.F7_mQTL_dosages_2021-05-17.Rdata")

# clean dosages:
dosages <- dosages[, -c(2:3)] # remove extraneous columns
## do they all start with rs? 
sum(!grepl("^rs", colnames(dosages)[-1])) # yup




# how many SNPs per CpG
count <- table(F7.mqtl.subset$gene)
#cg00863271 cg03305840 cg03764134 cg04913057 cg05060427 cg05465572 cg06957310 cg08638097 cg08641963 cg09087363 cg10498926 cg10639395 cg11256802 cg11308211 cg13294509 cg14062745 cg14613617 cg19498110 cg19543471 cg20102336 cg23685969 cg26112574 cg26114043 
#1          2          1          1          2          1          3          2       1001         33          2        553          2          3          1          1        205          1          1          3          1        310        215 
#cg27025925 
#178 
names(count[count >= 5]) #7

# also, issue with missing data.
# how many SNPs have missing data?
check.miss <- apply(dosages[, -1], 2, 
                    function(x) sum(is.na(x)))
table(check.miss)
# what happens if we exclude SNPs with missing rates > 1%?
miss.cutoff <- nrow(dosages) * 0.01
length(check.miss[check.miss > miss.cutoff])
#[1] 1066


## restrict missingness < 1% for those PCs
dosages2 <- dosages[, c("FID", names(check.miss[check.miss < miss.cutoff]))]
rownames(dosages2) <- dosages$FID
keep.SNPs <- gsub("_.*", "", names(dosages2))
F7.mqtl.subset.nomiss <- F7.mqtl.subset[F7.mqtl.subset$SNP %in% keep.SNPs, ]
length(unique(F7.mqtl.subset.nomiss$gene)) # 10 CpGs

# hopefully we still have all 2 for those with # SNP > 4:
sum(CpG.mQTL$Probe %in% names(count[count >= 5]))#7 yup!
length(intersect(CpG.mQTL[ CpG.mQTL$Probe %in% names(count[count > 4]), "Probe"], 
                 unique(F7.mqtl.subset.nomiss$gene)))
#[1] 7

# how many SNPs per CpG?
F7.mqtl.subset.nomiss <- F7.mqtl.subset.nomiss[F7.mqtl.subset.nomiss$gene %in% 
                                                 CpG.mQTL[CpG.mQTL$Probe %in% names(count[count > 4]), "Probe"], ]
count2<-table(F7.mqtl.subset.nomiss$gene)

#cg08641963 cg09087363 cg10639395 cg14613617 cg26112574 cg26114043 cg27025925 
#700         31        283        137        175         31         97 







## SLCMA analysis

source("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Codes/newDNAm/FWL/LARS-mobility-FWL-CI-20200722.R")

# load raw data (DNAm + Phenotypes)
load("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/newDNAm/newDNAm_pheno.Rdata")
colnames(df)[11:13]<-c("lowincome_newbin_very_early","lowincome_newbin_early" ,"lowincome_newbin_middle") #fix the inconsistent variable naming

df.mQTL <- df[,colnames(df) %in% CpG.mQTL$Probe]
df.mQTL <- cbind(df[,1:67],df.mQTL)
df<-df.mQTL

df$ID <- paste(df$cidB1471,df$qlet,sep="")




probe.list <- unique(F7.mqtl.subset$gene)


lars.fun.mQTL <- function(probe){ 
  outcome.vec <- probe
  adver=cpg_r2_g003_raw[cpg_r2_g003_raw$Probe==probe,"adversity"]
  SNP=F7.mqtl.subset[F7.mqtl.subset$gene==probe,"SNP"]$SNP
  SNP=unlist(lapply(SNP, grep,x=colnames(dosages), value=T))
  
  
  # Direct adjustment for CpGs with less than 4 SNPs
  if (length(SNP) <= 4) { 
    df.dosage<-merge(df,dosages,by.x="ID",by.y="FID") # Use the unfiltered SNP data (otherwise many CpGs would have no SNP)
    covar=c( "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
          "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood",SNP)
  
  
    if (adver %in% c("jobless", "income_reduction")) {
      lars.out<-select.LARS.complete.addmob(df.dosage, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=covar,
                                          hypos=c("accumulation"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE) 
   }
    else {
      lars.out<-select.LARS.complete.addmob(df.dosage, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=covar,
                                          hypos=c("accumulation","mobility"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE)
    } 
    lars.out$adversity<-adver
    return(lars.out)
  }
  
  
  # For those CpGs with more than 4 SNPs, adjust for the top 4 PCs
  if (length(SNP) > 4) { 
    df.dosage <- merge(df,dosages2,by.x="ID",by.y="FID")
    SNP_nomiss <- SNP[SNP %in% colnames(dosages2)]
    pcs <- prcomp(as.formula(paste0("~ ", 
                                    paste(SNP_nomiss, 
                                          collapse=" + "))), data=df.dosage[,SNP_nomiss], center = TRUE,scale. = TRUE,na.action = na.omit)
    df.dosage_pcs <-merge(df.dosage,pcs$x[,1:4],by="row.names")
    covar=c( "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
             "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood","PC1","PC2","PC3","PC4")
    
    if (adver %in% c("jobless", "income_reduction")) {
      lars.out<-select.LARS.complete.addmob(df.dosage_pcs, outcome.vec=outcome.vec, 
                                            adver=paste0(adver,"_newbin"),covars=covar,
                                            hypos=c("accumulation"), recency.vec =NULL, 
                                            inf.method = "sI",exposures="default",FWL=TRUE) 
    }
    else {
      lars.out<-select.LARS.complete.addmob(df.dosage_pcs, outcome.vec=outcome.vec, 
                                            adver=paste0(adver,"_newbin"),covars=covar,
                                            hypos=c("accumulation","mobility"), recency.vec =NULL, 
                                            inf.method = "sI",exposures="default",FWL=TRUE)
    } 
    
    
    lars.out$adversity<-adver
    return(lars.out)
    
    }
  
}

#probe="cg08641963"
#lars.fun.mQTL("cg00863271")
sens.mQTL<- data.frame(do.call("rbind", lapply(probe.list, lars.fun.mQTL)))

# read in the order of the CpGs that we want to arrange the output files
cpg.order <- read.csv(file = '/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/Table2/62hits_order.csv')
sens.mQTL.all <- merge(cpg.order,sens.mQTL[,-9],by.x="CpG",by.y="Probe",all.x=T)
cpg_r2_g003_raw<- cpg_r2_g003_raw[match(cpg.order$CpG, cpg_r2_g003_raw$Probe),] # order the CpGs 
sens.mQTL.all <- sens.mQTL.all[match(cpg.order$CpG, sens.mQTL.all$CpG),] # order the CpGs 

sens.mQTL.all$Variable<- as.character(sens.mQT.all$Variable)
cpg_r2_g003_raw$'X1.hypo'<- as.character(cpg_r2_g003_raw$'X1.hypo')
sum(cpg_r2_g003_raw$'X1.hypo'==sens.mQTL.all$Variable,na.rm = T) #24.Same hypothesis is selected for all probes

# calculate sample size

mQTL_samplesize <- function(probe){ 
  outcome.vec <- probe
  adver=cpg_r2_g003_raw[cpg_r2_g003_raw$Probe==probe,"adversity"]
  SNP=F7.mqtl.subset[F7.mqtl.subset$gene==probe,"SNP"]$SNP
  SNP=unlist(lapply(SNP, grep,x=colnames(dosages), value=T))
  
  if (length(SNP) <= 4) { 
    
    df.dosage<-merge(df,dosages,by.x="ID",by.y="FID")
    covar=c( "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
           "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood",SNP)
    df.analysis <- df.dosage[,c(covar,  "cidB1471", grep(paste0("^",adver), colnames(df.dosage), value=T))]
    df.analysis <- df.analysis[complete.cases(df.analysis),] #remove subject with missing data
    n <- nrow(df.analysis)
  }
  
  if (length(SNP) > 4) { 
    df.dosage <- merge(df,dosages2,by.x="ID",by.y="FID")
    SNP_nomiss <- SNP[SNP %in% colnames(dosages2)]
    pcs <- prcomp(as.formula(paste0("~ ", 
                                    paste(SNP_nomiss, 
                                          collapse=" + "))), data=df.dosage[,SNP_nomiss], center = TRUE,scale. = TRUE,na.action = na.omit)
    df.dosage_pcs <-merge(df.dosage,pcs$x[,1:4],by="row.names")
    covar=c( "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
             "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood","PC1","PC2","PC3","PC4")
    df.analysis <- df.dosage_pcs[,c(covar,  "cidB1471", grep(paste0("^",adver), colnames(df.dosage_pcs), value=T))]
    df.analysis <- df.analysis[complete.cases(df.analysis),] #remove subject with missing data
    n <- nrow(df.analysis)

  }
  
  return(c(probe,n))
}


n.mQTL<- data.frame(do.call("rbind", lapply(probe.list, mQTL_samplesize)))
colnames(n.mQTL) <-c("probe","samplesize")

# number of mQTLs
num.mQTL <- data.frame(count)
colnames(num.mQTL) <- c("probe","n_mQTL")
num.mQTL$probe <- as.character(num.mQTL$probe)
count2<-table(F7.mqtl.subset.nomiss$gene)
num.mQTL.ge5 <- data.frame(count2)
colnames(num.mQTL.ge5) <- c("probe","n_mQTL")
num.mQTL.ge5$probe <- as.character(num.mQTL.ge5$probe)


for (i in 1:dim(num.mQTL)[1]) {
  probe=num.mQTL[i,"probe"] 
  if (probe %in% num.mQTL.ge5$probe) {
    num.mQTL[i,"n_mQTL"] <- num.mQTL.ge5[num.mQTL.ge5$probe==probe,"n_mQTL"] # Update the sample size for those CpGs with more than 5 mQTLs
  }
}
# save results

sens.mQTL.all <- merge(sens.mQTL.all,n.mQTL,by.x="CpG",by.y="probe",all.x=T)
sens.mQTL.all <- merge(sens.mQTL.all,num.mQTL,by.x="CpG",by.y="probe",all.x=T)
sens.mQTL.all <- sens.mQTL.all[match(cpg.order$CpG, sens.mQTL.all$CpG),] # order the CpGs 


write.xlsx(sens.mQTL.all, file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/Sensitivity/sensitivity_mQTL_2021-05-17.xlsx")

save(sens.mQTL.all, file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/mQTL_results_2021-05-17.Rdata")
