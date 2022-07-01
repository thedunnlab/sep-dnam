# Table 2

#### packages ---- 
library(dplyr)
library(tibble)


#### load top hits ####
load('/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/summary_r2g003_62hits_raw.RData')
cpg_r2_g003_raw$'X1.hypo' <- as.character(cpg_r2_g003_raw$'X1.hypo')
table(cpg_r2_g003_raw[cpg_r2_g003_raw$adversity=="income",'X1.hypo'])
cpg_r2_g003_raw[cpg_r2_g003_raw$adversity=="income","adversity"]<-"lowincome" #fix the naming issue with income variable (overlap with income_reduction; if not fix this will get incorrect sample size for income; in the original SLCMA analysis I entered the exact variable names, i.e. "income_bin" and "income_reduction_newbin" so no need to worry)
cpg_r2_g003_raw[cpg_r2_g003_raw$'X1.hypo'=="income_bin_very_early",'X1.hypo']<-"lowincome_newbin_very_early" #fix the naming issue 
cpg_r2_g003_raw[cpg_r2_g003_raw$'X1.hypo'=="income_bin_early",'X1.hypo']<-"lowincome_newbin_early" #fix the naming issue 
cpg_r2_g003_raw[cpg_r2_g003_raw$'X1.hypo'=="income_bin_middle",'X1.hypo']<-"lowincome_newbin_middle" #fix the naming issue 
table(cpg_r2_g003_raw[cpg_r2_g003_raw$adversity=="lowincome",'X1.hypo'])


## load DNAm & phenotype data
load("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/newDNAm/newDNAm_pheno.Rdata")

## get the DNAm for hits
df.hits <- df[,colnames(df) %in% cpg_r2_g003_raw$Probe]
df.hits <- cbind(df[,1:67],df.hits)
colnames(df.hits)[11:13]<-c("lowincome_newbin_very_early","lowincome_newbin_early" ,"lowincome_newbin_middle") #fix the inconsistent variable naming
rm(df)

# Estimate beta coefficients

#adv.list <- c("Job loss", "Income reduction", "Income", "Financial difficulty score", "Major financial problem", "Neighborhood disadvantage")
adv.list <- c("jobless", "income_reduction",   "lowincome", "Fscore","major_F", "nbhqual")

#### The function of generating table-2 ####
run.stage2 <- function(adver){
  
  ### step 1: construct the dataset --- 
  
  results <- cpg_r2_g003_raw %>% filter(adversity == adver)
  results <- results %>% 
    dplyr::select(Probe, p.val = `X1.p`, hypo = `X1.hypo`, r2 = `X1.r2`)
  list.cpg <- results$Probe  
  
  ## regress the betas on the cell counts first
  betas <- df.hits[,c(list.cpg, "cidB1471", "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran")]
  betas[,colnames(betas) %in% list.cpg] <- do.call("cbind", lapply(list.cpg, function(x){
    r <- lm(as.formula(paste0(x, "~ CD8T + CD4T + Bcell + 
                              Mono + NK + Gran")), data = betas)$resid
    r <- as.data.frame(r)
    colnames(r) <- x
    return(r)}))  
  
  covars <- c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
              "age_DNAm7","wholeblood") 
  
  df <- df.hits[,c(covars,  "cidB1471", grep(paste0("^",adver), colnames(df.hits), value=T))]
  df <- df[complete.cases(df),] #remove subject with missing data 
  
  # create accumulation and mobility variables
  adv <- df[,grep(paste0("^",adver), colnames(df))]
  adv <- sapply(adv, function(x) as.numeric(as.character(x)))
  df$accumulation <- rowSums(adv)
  x1=paste0(adver,"_newbin_very_early")
  x2=paste0(adver,"_newbin_early")
  x3=paste0(adver,"_newbin_middle")
  df$mobility_D12 <- (1-df[,x2])*df[,x1]
  df$mobility_U12 <- (1-df[,x1])*df[,x2]
  df$mobility_D23 <- (1-df[,x3])*df[,x2]
  df$mobility_U23 <- (1-df[,x2])*df[,x3]
  
  # residuals, used in regression
  df.adj <- left_join(df, betas, by = "cidB1471")
  # raw betas, used in calculating group DNAm differences
  betas <- df.hits[,c(list.cpg, "cidB1471")]
  df.raw <- left_join(df, betas, by = "cidB1471")
  df.raw[, grep(paste0("^",adver), colnames(df.raw))] <- sapply(df.raw[, grep(paste0("^",adver), colnames(df.raw))], function(x) as.numeric(as.character(x))) 
  
  
  ### step 2 run stage 2 analysis ---- 
  stage2.est <- do.call("rbind", lapply(list.cpg, function(x){
    var <- as.character(results[results$Probe == x,"hypo"])
    form <- as.formula(paste0(x, "~ WHITE + Female + mom_birthage + ppregnum + birth.weight + sustained.smoke + age_DNAm7 + wholeblood +", var))
    fit <- lm(form, data = df.adj)
    b <- summary(fit)$coef[var, "Estimate"]
    s <- summary(fit)$coef[var, "Std. Error"]
    p <- results[results$Probe == x,"p.val"]
    lcl <- b + qnorm((0.025 - p/4)/(1-p/2))*s
    ucl <- b + qnorm((0.975 - p/4)/(1-p/2))*s
    
    # group differences
    DNAm.expos <- mean(df.raw %>% filter(!!as.name(var) > 0) %>% dplyr::select(!!x) %>% .[,1] )
    DNAm.unexpos <- mean(df.raw %>% filter(!!as.name(var) == 0) %>% dplyr::select(!!x) %>% .[,1] )
    
    r <- c(x, b, s, lcl, ucl, DNAm.expos, DNAm.unexpos)
  }))
  
  
  stage2.est <- as.data.frame(stage2.est)
  colnames(stage2.est) <- c("CpG", "beta", "SE", "CI.lo", "CI.up", "DNAm.expos", "DNAm.unexpos")
  
  ### step 3 add other statistics --- 
  table1 <- left_join(results, stage2.est, by = c("Probe" = "CpG")) %>% 
    mutate(adversity = adver, 
           samplesize = nrow(df)) %>% 
    dplyr::select(CpG = Probe, 
                  adversity, 
                  samplesize, 
                  hypo, 
                  DNAm.unexpos,
                  DNAm.expos, 
                  r2, 
                  p.val, 
                  beta, 
                  SE,
                  CI.lo,
                  CI.up)
  
  return(table1)
}


#### END OF FUNCTION ####

table2 <- data.frame(do.call("rbind", lapply(adv.list, run.stage2)))



#### Calculate CI using the updated pipeline  ####
source("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Codes/newDNAm/FWL/LARS-mobility-FWL-CI-20200722.R")

adv.list <- c("jobless", "income_reduction",   "lowincome", "Fscore","major_F", "nbhqual")

lars.fun <- function(adver){ 
  
  outcome.vec <- cpg_r2_g003_raw[cpg_r2_g003_raw$adversity==adver,"Probe"]
  
  if (adver %in% c("jobless", "income_reduction")) {
    lars.out<-select.LARS.complete.addmob(df.hits, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood"),
                                          hypos=c("accumulation"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE) 
  }
  else {
    lars.out<-select.LARS.complete.addmob(df.hits, outcome.vec=outcome.vec, 
                                          adver=paste0(adver,"_newbin"),covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                                                 "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood"),
                                          hypos=c("accumulation","mobility"), recency.vec =NULL, 
                                          inf.method = "sI",exposures="default",FWL=TRUE)
  } 
  lars.out$adversity<-adver
  return(lars.out)
}


table2.ci<- data.frame(do.call("rbind", lapply(adv.list, lars.fun)))
colnames(table2.ci)[3]<-"Beta_SI"

table2.final<-merge(table2[,c("CpG","adversity","samplesize","hypo","DNAm.unexpos","DNAm.expos")],
                    table2.ci[,c("Probe","Beta_SI","SE","Pvalue","CI.lo","CI.up","R2")],
                    by.x="CpG",by.y="Probe")

table2.final<-table2.final[order(table2.final$adversity),]

save(table2.final, file = "/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/Table2/Table2_results_62hits_2021-05-16.Rdata")




#### run annotation ####

## packages
library(dplyr)
#BiocManager::install("FDb.InfiniumMethylation.hg19")
library(FDb.InfiniumMethylation.hg19)
library(tibble)

hm450 <- get450k()
probes <- hm450[table2.final$CpG]
x <- getNearestGene(probes)
summary(x$distance)
x <- x %>% rownames_to_column("CpG")

## output: 
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.0      0.0      0.0   7865.1    489.8 115163.0 

table2.final.annot <- left_join(table2.final, x[,c("CpG", "nearestGeneSymbol", "distance")], by = "CpG")

# Other annotation

load("/data/js95/ALSPAC/ARIES/featureData.Rdata")

fData <- featureData[row.names(featureData) %in% table2.final.annot$CpG, ]
rm(featureData)
fData <- fData %>% rownames_to_column("CpG")
table2.final.annot<-merge(table2.final.annot,fData[,c("CpG", "CHR37", "COORDINATE_37")],by="CpG")

anno<-read.csv("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Annotation/HumanMethylation450_15017482_v1-2.csv",header=T)
table2.final.annot<-merge(table2.final.annot,anno[,c("IlmnID","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island", "Enhancer","Regulatory_Feature_Name","DHS")],by.x="CpG",by.y="IlmnID")

## clean up the format 
## it seems FDb.InfiniumMethylation.hg19 used the h19 build, so we should use CHR37 and COORDINATE_37 to be consistent
table2.final.annot <- table2.final.annot %>% 
  dplyr::rename(Chr = CHR37, Coordinate = COORDINATE_37, NearestGene = nearestGeneSymbol) ## lookin' good 

table2.final.annot$Chr <- unlist(table2.final.annot$Chr)
table2.final.annot$Coordinate <- unlist(table2.final.annot$Coordinate)
lapply(table2.final.annot, class)

## save 
save(table2.final.annot, file = "/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/Table2/Table2_62hits_annot_2021-05-17.Rdata")
write.csv(table2.final.annot, "/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/Table2/Table2_62hits_annot_2021-05-17.csv", row.names = F)





