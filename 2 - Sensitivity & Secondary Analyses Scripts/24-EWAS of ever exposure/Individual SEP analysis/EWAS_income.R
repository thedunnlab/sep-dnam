##### EWAS of ever-exposure  #####
library(dplyr,lib.loc = "/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/R_package")

### load data
load("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/newDNAm/newDNAm_pheno.Rdata") #dataset name: df

adver<-"income_bin"
outcome.vec <- colnames(df)[68:ncol(df)]
covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
         "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood")

### Run EWAS
# Restrict analysis to people with complete data, to be consistent with SLCMA analysis 
DNAm <- df[, c("ID", outcome.vec)]
df <- df[, !grepl("^cg", colnames(df)) & !grepl("^ch", colnames(df))]
exposures <- grep(paste0("^", adver), colnames(df), value = T)


# subset to those with complete data on this exposure and all covariates
df <- df %>% filter(rowSums(is.na(.[,c(exposures,covars)])) == 0)
DNAm <- DNAm[DNAm$ID %in% df$ID, ]

# define sample size:
n <- nrow(df) 

# create ever-exposure
exp <- do.call("cbind", lapply(df[, exposures], function(x) as.numeric(as.character(x))))
df <- df %>% mutate(ever = ifelse(rowSums(exp) > 0, 1, 0))

# what are the categorical covariates
cat.covars <- do.call("c", lapply(1:length(covars), function(x) {
  cov1 <- covars[x]
  if (is.factor(df[, cov1])) {
    return(cov1)
  } else {
    return(NULL)
  }
}))
# dummy code the categorical covars
cat.covariates <- model.matrix(as.formula(paste0("~", paste(cat.covars, collapse = "+"))), df)[,-1]
predictors <- as.matrix(cbind(df$ever,cat.covariates, as.matrix(df[,setdiff(covars, cat.covars)])))
colnames(predictors)[1]<-"ever"

# Linear regression
results<-do.call("rbind", lapply(outcome.vec, function(x) {
  lm.res<-lm(DNAm[,x] ~ predictors)
  result<-summary(lm.res)$coefficients[2,]
  return(result)
}))
results<-data.frame(results)
colnames(results)<-c("beta","SE","t_value","p_value")
results$probe<-outcome.vec
results$adversity<-adver
results$samplesize<-n
results$ever_preval<-prop.table(table(df$ever))[2]


# output results 
results.EWAS.income_bin<-results
save(results.EWAS.income_bin,file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/newDNAm/EWAS/results_EWAS_income_bin_081820.Rdata")

