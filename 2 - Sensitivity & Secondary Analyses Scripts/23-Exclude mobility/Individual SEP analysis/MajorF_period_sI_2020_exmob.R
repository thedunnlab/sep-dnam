
source("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Codes/newDNAm/FWL/LARS-mobility-FWL-20200316.R")

load("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/newDNAm/newDNAm_pheno.Rdata")

outcome.vec <- colnames(df)[68:ncol(df)]

results.major_F.sI<-select.LARS.complete.addmob(df, outcome.vec=outcome.vec, 
                            adver="major_F_newbin",covars=c("WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                                            "CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran","age_DNAm7","wholeblood"),
                            hypos=c("accumulation"), recency.vec =NULL, 
                            inf.method = "sI",exposures="default",FWL=TRUE)

save(results.major_F.sI,file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/newDNAm/FWL_Results/Ex_mobility/results_major_F_period_exmob_062320_sI.Rdata")




