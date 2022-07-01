

source("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Codes/newDNAm/FWL/LARS-mobility-FWL-20200316.R")

load("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/newDNAm/newDNAm_pheno.Rdata")

outcome.vec <- colnames(df)[68:ncol(df)]

results.nbhqual.sI<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                            adver="nbhqual_newbin",covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                            "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood"),
                            hypos=c("accumulation","mobility"), recency.vec =NULL, 
                            inf.method = "sI",exposures="default",FWL=TRUE)

save(results.nbhqual.sI,file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/newDNAm/FWL_Results/results_nbhqual_period_031720_sI.Rdata")



