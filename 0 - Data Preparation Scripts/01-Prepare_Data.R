############################
#Create analysis dataset
############################

load("/data/js95/ALSPAC/ARIES/DNAm_2020/F7/betas_F7_WIN_20200128.Rdata")
load("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/period_phenotypes_050719.Rdata")

samplesheet<-read.table("/data/js95/ALSPAC/ARIES/DNAm_2020/F7/samplesheet_F7_20200121.txt",header=T)
samplesheet$wholeblood<-ifelse(samplesheet$sample_type=="wholeblood",1,0)

cell<-read.table("/data/js95/ALSPAC/ARIES/DNAm_2020/F7/celltypes_F7_20200121.txt",header=T)

sample<-merge(samplesheet,betas.f7.win,by.x="Sample_Name",by.y="row.names")
sample<-merge(cell,sample,by="Sample_Name")

df<-merge(pheno.period[,-c(43:48)],sample,by.x="cidB1471",by.y = "ALN")

save(df,file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/newDNAm/newDNAm_pheno.Rdata")

#check overlaps with old DNAm data
load("/data/js95/ALSPAC/ARIES/F7/tostBetas_F7.rv2_winsor.Rdata")
sample$cids=paste(sample$ALN,"A",sep="")
overlap<-merge(sample,tbetas.winsor, by.x="cids",by.y="row.names") #950 overlaps

#test the pipeline in the small subset
load("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/newDNAm/newDNAm_pheno.Rdata")
df.sb=df[,1:80]
save(df.sb,file="/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/newDNAm/newDNAm_pheno_sub.Rdata")

sum(is.na(df$jobless_newbin_early)==F & is.na(df$jobless_newbin_middle)==F 
    & is.na( df$jobless_newbin_very_early)==F )

sum(is.na(df$income_bin_early)==F & is.na(df$income_bin_middle)==F 
    & is.na( df$income_bin_very_early)==F  )

sum(is.na(df$major_F_newbin_early)==F & is.na(df$major_F_newbin_middle)==F 
    & is.na( df$major_F_newbin_very_early)==F  )


