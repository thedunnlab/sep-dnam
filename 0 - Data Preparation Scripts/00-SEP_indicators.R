#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ARIES SES measures
# Notes
# - This is a cleaned version of coding script for SEP measures.
# - I removed redundant codes or exploratory analyses that are no longer included in the manuscript.
# - The previous version (reviewed by Brooke):
#   /Users/jiaxliu/Dropbox (Partners HealthCare)/ALSPAC- Russell Sage/SES Measurement/Jessie Data/Codes/SEP_clean_2021-05-10.R


# Created by: Jiaxuan (Jessie) Liu
# Originally created on: Dec 11 2018
# Last checked on: July 2, 2021
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(stringr)
library(dplyr)
library(xlsx)
library(magrittr)
############################################################
# Load the very first dataset  (from ARIES_Table1_Jessie.R)
############################################################

#load data, 971 subjects with DNAm at age 7
load("/Users/jiaxliu/Dropbox (Partners HealthCare)/ALSPAC- Russell Sage/Data/ARIES_SES_20181207.Rdata")
dim(ARIES_SES)
colnames(ARIES_SES)
head(ARIES_SES)#971 149

# What's the missing rate of Fscore at this stage? (2021-05-09)
dim(na.omit(ARIES_SES[,c("FScore_8m","FScore_21m","FScore_33m","FScore_61m","FScore_85m")])) # 773

##What's the missing rate of neighborhood quality at this stage? (2021-05-09)
table(ARIES_SES$nbhqual_85)# continuous
table(ARIES_SES$nbhqual_85_bin)# binary
dim(na.omit(ARIES_SES[,c("nbhqual_21" ,"nbhqual_33","nbhqual_61" ,"nbhqual_85")])) # 769


#take a peek at the data
#jobless
table(ARIES_SES$jbls_8wks,exclude=F)
table(ARIES_SES$jbls_8mos,exclude=F)
table(ARIES_SES$jbls_21mos,exclude=F)
table(ARIES_SES$jbls_33mos,exclude=F) 
table(ARIES_SES$jbls_47mos,exclude=F) 
table(ARIES_SES$jbls_61mos,exclude=F) 
table(ARIES_SES$jbls_73mos,exclude=F) 
 
#income reduction
table(ARIES_SES$inrdx_8wks,exclude=F)
table(ARIES_SES$inrdx_8mos,exclude=F)
table(ARIES_SES$inrdx_21mos,exclude=F)
table(ARIES_SES$inrdx_33mos,exclude=F) 
table(ARIES_SES$inrdx_47mos,exclude=F) 
table(ARIES_SES$inrdx_61mos,exclude=F) 
table(ARIES_SES$inrdx_73mos,exclude=F) 
#Neighborhood 
table(ARIES_SES$nbhqual_21,exclude=F) #0 to 16
table(ARIES_SES$nbhqual_33,exclude=F) 
table(ARIES_SES$nbhqual_61,exclude=F) 
table(ARIES_SES$nbhqual_85,exclude=F) 

# Townsend deprivation index (included in sensitivty analysis)
table(ARIES_SES$gTownsendq5)#1 to 5, and missing is coded as "."
table(ARIES_SES$cTownsendq5)#1 to 5, and missing is coded as "."

#recodes the Townsend variables:

ARIES_SES$Townsend_pre32<-ifelse(str_trim(ARIES_SES$cTownsendq5)==".",NA,str_trim(ARIES_SES$cTownsendq5))
table(ARIES_SES$Townsend_pre32,ARIES_SES$cTownsendq5)
ARIES_SES$Townsend_8wks<-ifelse(str_trim(ARIES_SES$eTownsendq5)==".",NA,str_trim(ARIES_SES$eTownsendq5))
table(ARIES_SES$Townsend_8wks,ARIES_SES$eTownsendq5)
ARIES_SES$Townsend_8mos<-ifelse(str_trim(ARIES_SES$fTownsendq5)==".",NA,str_trim(ARIES_SES$fTownsendq5))
table(ARIES_SES$Townsend_8mos,ARIES_SES$fTownsendq5)
ARIES_SES$Townsend_21mos<-ifelse(str_trim(ARIES_SES$gTownsendq5)==".",NA,str_trim(ARIES_SES$gTownsendq5))
table(ARIES_SES$Townsend_21mos,ARIES_SES$gTownsendq5)
ARIES_SES$Townsend_33mos<-ifelse(str_trim(ARIES_SES$hTownsendq5)==".",NA,str_trim(ARIES_SES$hTownsendq5))
table(ARIES_SES$Townsend_33mos,ARIES_SES$hTownsendq5)
ARIES_SES$Townsend_47mos<-ifelse(str_trim(ARIES_SES$jTownsendq5)==".",NA,str_trim(ARIES_SES$jTownsendq5))
table(ARIES_SES$Townsend_47mos,ARIES_SES$jTownsendq5)
ARIES_SES$Townsend_61mos<-ifelse(str_trim(ARIES_SES$kTownsendq5)==".",NA,str_trim(ARIES_SES$kTownsendq5))
table(ARIES_SES$Townsend_61mos,ARIES_SES$kTownsendq5)
ARIES_SES$Townsend_73mos<-ifelse(str_trim(ARIES_SES$lTownsendq5)==".",NA,str_trim(ARIES_SES$lTownsendq5))
table(ARIES_SES$Townsend_73mos,ARIES_SES$lTownsendq5)
ARIES_SES$Townsend_85mos<-ifelse(str_trim(ARIES_SES$mTownsendq5)==".",NA,str_trim(ARIES_SES$mTownsendq5))
table(ARIES_SES$Townsend_85mos,ARIES_SES$mTownsendq5)



# recode missingness for maternal edu and home ownership
table(ARIES_SES$ed_momgest)
ARIES_SES$ed_momgest<-ifelse(str_trim(ARIES_SES$ed_momgest)==".",NA,str_trim(ARIES_SES$ed_momgest))
table(ARIES_SES$ed_momgest)

table(ARIES_SES$home_owner)
ARIES_SES$home_owner<-ifelse(str_trim(ARIES_SES$home_owner)==".",NA,str_trim(ARIES_SES$home_owner))
table(ARIES_SES$home_owner)

save(ARIES_SES,file="/Users/jiaxliu/Dropbox (Partners HealthCare)/ALSPAC- Russell Sage/Data/ARIES_SES_combined_Jessie-2021-07-02.Rdata") # I added today's date "-2021-07-02" to avoid overwriting the original file

load("/Users/jiaxliu/Dropbox (Partners HealthCare)/ALSPAC- Russell Sage/Data/ARIES_SES_combined_Jessie-2021-07-02.Rdata") # I added today's date "-2021-07-02" to avoid overwriting the original file



######################################################
### Load additional SEP data (from additional_SES.R)
######################################################

### load new data ###
load("/Users/jiaxliu/Dropbox (Partners HealthCare)/ALSPAC- Russell Sage/Data/SES_additional_vars_20190301.Rdata")
dim(df.SES.addvars)#971 113
colnames(df.SES.addvars)

### rename variable ###
# get the variable name
rename_list=read.csv(file="/Users/jiaxliu/Dropbox (Partners HealthCare)/ALSPAC- Russell Sage/SES Measurement/Jessie Data/Data/rename_list_add_var.csv", header=T, stringsAsFactors=F)
# rename
df.SES.addvars.rename<-df.SES.addvars %>% rename_at(vars(rename_list$Variable_Name), ~ rename_list$rename) 
colnames(df.SES.addvars.rename)
# keep only the new variables
df.SES.addvars.rename<-df.SES.addvars.rename[,colnames(df.SES.addvars.rename) %in% c("cidB1471","qlet",rename_list$rename)]
dim(df.SES.addvars.rename)#971 56


### Take a peek at the data
## major financial problem
table(df.SES.addvars.rename$major_F_73mos,exclude=F)
table(df.SES.addvars.rename$major_F_8mos,exclude=F)
table(df.SES.addvars.rename$major_F_8mos,exclude=F)
table(df.SES.addvars.rename$major_F_33mos,exclude=F)

## household weekly income
table(df.SES.addvars.rename$income_33mos,exclude =F)
table(df.SES.addvars.rename$income_47mos,exclude =F)
table(df.SES.addvars.rename$income_85mos,exclude =F)
#recode missing for 85mos
df.SES.addvars.rename$income_85mos<-ifelse(df.SES.addvars.rename$income_85mos==9,NA,df.SES.addvars.rename$income_85mos)
table(df.SES.addvars.rename$income_85mos, exclude=F)


save(df.SES.addvars.rename,file="/Users/jiaxliu/Dropbox (Partners HealthCare)/ALSPAC- Russell Sage/Data/ARIES_SES_add_Jessie_2021-07-02.Rdata")


######################################################
### Recode binary indicators (from recode_binary_SES.R)
######################################################

### load the two datasets prepared above ###

load("/Users/jiaxliu/Dropbox (Partners HealthCare)/ALSPAC- Russell Sage/Data/ARIES_SES_add_Jessie_2021-07-02.Rdata")
load("/Users/jiaxliu/Dropbox (Partners HealthCare)/ALSPAC- Russell Sage/Data/ARIES_SES_combined_Jessie-2021-07-02.Rdata")

dim(df.SES.addvars.rename) #971 56
dim(ARIES_SES)# 971 158

## What's the missing rate of Fscore at this stage? (2021-05-09)
dim(na.omit(ARIES_SES[,c("FScore_8m","FScore_21m","FScore_33m","FScore_61m","FScore_85m")])) # 773

## What's the missing rate of neighborhood quality at this stage? (2021-05-09)
table(ARIES_SES$nbhqual_85)# continuous
table(ARIES_SES$nbhqual_85_bin)# binary
dim(na.omit(ARIES_SES[,c("nbhqual_21" ,"nbhqual_33","nbhqual_61" ,"nbhqual_85")])) # 769

## merge two datasets
SES.merge=merge(df.SES.addvars.rename,ARIES_SES,by="cidB1471")

dim(SES.merge) #971 213


### create new binary indicators

attach(SES.merge)

# job loss
jbls_newbin_8wks=ifelse(jbls_8wks<4,1,0)
table(jbls_newbin_8wks,jbls_8wks,exclude = FALSE)
table(jbls_newbin_8wks,jbls_bin_8wks,exclude = FALSE)
jbls_newbin_8mos=ifelse(jbls_8mos<4,1,0)
jbls_newbin_21mos=ifelse(jbls_21mos<4,1,0)
jbls_newbin_33mos=ifelse(jbls_33mos<4,1,0)
jbls_newbin_47mos=ifelse(jbls_47mos<4,1,0)
jbls_newbin_61mos=ifelse(jbls_61mos<4,1,0)
jbls_newbin_73mos=ifelse(jbls_73mos<4,1,0)
jobless_newbin=c("jbls_newbin_8wks","jbls_newbin_8mos","jbls_newbin_21mos","jbls_newbin_33mos","jbls_newbin_47mos","jbls_newbin_61mos","jbls_newbin_73mos")

SES.merge.newbin=data.frame(SES.merge,jbls_newbin_8wks,jbls_newbin_8mos,jbls_newbin_21mos,jbls_newbin_33mos,jbls_newbin_47mos,jbls_newbin_61mos,jbls_newbin_73mos)

# income reduction
inrdx_newbin_8wks=ifelse(inrdx_8wks<4,1,0)
table(inrdx_newbin_8wks,inrdx_8wks,exclude = FALSE)
table(inrdx_newbin_8wks,inrdx_bin_8wks,exclude = FALSE)
inrdx_newbin_8mos=ifelse(inrdx_8mos<4,1,0)
inrdx_newbin_21mos=ifelse(inrdx_21mos<4,1,0)
inrdx_newbin_33mos=ifelse(inrdx_33mos<4,1,0)
inrdx_newbin_47mos=ifelse(inrdx_47mos<4,1,0)
inrdx_newbin_61mos=ifelse(inrdx_61mos<4,1,0)
inrdx_newbin_73mos=ifelse(inrdx_73mos<4,1,0)
income_reduction_newbin=c("inrdx_newbin_8wks","inrdx_newbin_8mos","inrdx_newbin_21mos","inrdx_newbin_33mos","inrdx_newbin_47mos","inrdx_newbin_61mos","inrdx_newbin_73mos")

SES.merge.newbin=data.frame(SES.merge.newbin,inrdx_newbin_8wks,inrdx_newbin_8mos,inrdx_newbin_21mos,inrdx_newbin_33mos,inrdx_newbin_47mos,inrdx_newbin_61mos,inrdx_newbin_73mos)

# major finanical problem
major_F_newbin_8mos=ifelse(major_F_8mos<4,1,0)
table(major_F_newbin_8mos,major_F_8mos,exclude = FALSE)
table(major_F_newbin_8mos,major_F_bin_8mos,exclude = FALSE)
major_F_newbin_33mos=ifelse(major_F_33mos<4,1,0)
major_F_newbin_61mos=ifelse(major_F_61mos<4,1,0)
major_F_newbin_73mos=ifelse(major_F_73mos<4,1,0)
major_F_newbin=c("major_F_newbin_8mos","major_F_newbin_33mos","major_F_newbin_61mos","major_F_newbin_73mos")

SES.merge.newbin=data.frame(SES.merge.newbin,major_F_newbin_8mos,major_F_newbin_33mos,major_F_newbin_61mos,major_F_newbin_73mos)

# family income
income_bin_33mos=ifelse(income_33mos<=2,1,0)
table(income_bin_33mos,income_33mos,exclude = F)
income_bin_47mos=ifelse(income_47mos<=2,1,0)
income_bin_85mos=ifelse(income_85mos<=2,1,0)
income_bin=c("income_bin_33mos","income_bin_47mos","income_bin_85mos")

SES.merge.newbin=data.frame(SES.merge.newbin,income_bin_33mos,income_bin_47mos,income_bin_85mos)

# neighborhood quality
quantile(SES.merge.newbin$nbhqual_21,0.95,na.rm=T)#9
quantile(SES.merge.newbin$nbhqual_21,0.93,na.rm=T)#8
quantile(SES.merge.newbin$nbhqual_21,0.9,na.rm=T)#7
nbhqual_newbin_21mos=ifelse(SES.merge.newbin$nbhqual_21>=8,1,0)
nbhqual_newbin_33mos=ifelse(SES.merge.newbin$nbhqual_33>=8,1,0)
nbhqual_newbin_61mos=ifelse(SES.merge.newbin$nbhqual_61>=8,1,0)
nbhqual_newbin_85mos=ifelse(SES.merge.newbin$nbhqual_85>=8,1,0)
table(nbhqual_newbin_85mos,SES.merge.newbin$nbhqual_85_bin)#all correct
SES.merge.newbin=data.frame(SES.merge.newbin,nbhqual_newbin_21mos,nbhqual_newbin_33mos,nbhqual_newbin_61mos,nbhqual_newbin_85mos)
nbhqual_newbin=c("nbhqual_newbin_21mos","nbhqual_newbin_33mos","nbhqual_newbin_61mos","nbhqual_newbin_85mos")

## What's the missing rate of neighborhood quality at this stage? (2021-05-09)
dim(na.omit(SES.merge.newbin[,nbhqual_newbin])) # 769! 



### create the new Fscore exposure variable from imputation data
load("/Users/jiaxliu/Dropbox (Partners HealthCare)/ALSPAC- Russell Sage/Results/Yiwen/RussellSage-SES/data/Fscore_imputed51_mean_for_subsidies_20190328.Rdata")
dim(df.mean) #971 21
imp.items.list=c(colnames(df.mean)[2:11])

##What's the missing rate of Fscore at this stage? (2021-05-09)
dim (na.omit(df.mean)) # 753


table(df.mean$rent_21mos,SES.merge.newbin$rent_21mos)# for people without subsidy, the imputed value is the same as the observed value
table(df.mean$rent_33mos,SES.merge.newbin$rent_33mos)
table(df.mean$heat_21mos,SES.merge.newbin$heat_21mos)

colnames(df.mean)[2:11]=paste0(colnames(df.mean)[2:11],"_imp_mean")
df.mean.imp=df.mean[,1:11]
df.mean.imp %<>% mutate_all(funs(as.numeric(.)))
table(df.mean.imp$heat_21mos_imp_mean,SES.merge.newbin$heat_21mos)

SES.merge.newbin.imp=merge(SES.merge.newbin,df.mean.imp,by="cidB1471")
table(SES.merge.newbin.imp$rent_21mos_imp_mean,SES.merge.newbin.imp$rent_21mos,exclude=F) # Those with NA were orginally missing


for (i in 1:10) { #fix NA values (they should not be imputed)
  SES.merge.newbin.imp[is.na(SES.merge.newbin.imp[,imp.items.list[i]])==T,paste0(imp.items.list[i],"_imp_mean")]=NA
}

detach(SES.merge)
attach(SES.merge.newbin.imp)

#creat item-specific exposure vairables
timepoints <- c("_8mos","_21mos", "_33mos", "_61mos","_85mos")
items <- c("rent","heat","child","food","cloth")

SES.merge.newbin.imp.mean=SES.merge.newbin.imp

for (i in 1:5) {
  for (j in 1:5) {
    var.name=paste0(items[i],timepoints[j])
    if (var.name %in% imp.items.list) {
      var.name=paste0(var.name,"_imp_mean")
    }
    x=ifelse(SES.merge.newbin.imp[,var.name]<=2,1,0)
    SES.merge.newbin.imp.mean=data.frame(SES.merge.newbin.imp.mean,x)
    colnames(SES.merge.newbin.imp.mean)[dim(SES.merge.newbin.imp.mean)[2]]=paste0(items[i],"_newbin",timepoints[j])
  }
}

table(SES.merge.newbin.imp.mean$rent_61mos_imp_mean,SES.merge.newbin.imp.mean$rent_newbin_61mos,exclude=F)


#create Fscore variables 

# Code Fscore 
for (i in 1:5) {
  var.name=paste0("Fscore_imp",timepoints[i])
  data=SES.merge.newbin.imp.mean[,paste0(items,timepoints[i])]
  #a=rowSums(SES.merge.newbin.imp.mean[,paste0(items,"_newbin",timepoints[i])],na.rm=F)#people with at least one NA wwill be assigned missing
  b=rowSums(SES.merge.newbin.imp.mean[,paste0(items,"_newbin",timepoints[i])],na.rm=T) #number of (reported) exposed items, ignoring NA values
  c=rowSums(1-SES.merge.newbin.imp.mean[,paste0(items,"_newbin",timepoints[i])],na.rm=T) #number of (reported) unexposed items, ignoring NA values
  x=ifelse(b >= 3, 1, ifelse(c>=3, 0, NA)) # coded as Fscore=1 if 3 or more items were 1 (regardless of missing values), and coded as Fscore=0 if 3 or more items were 0.
  SES.merge.newbin.imp.mean=data.frame(SES.merge.newbin.imp.mean,x)
  colnames(SES.merge.newbin.imp.mean)[dim(SES.merge.newbin.imp.mean)[2]]=paste0("Fscore_newbin",timepoints[i])
}



# save data
save(SES.merge.newbin.imp.mean,file="/Users/jiaxliu/Dropbox (Partners HealthCare)/ALSPAC- Russell Sage/SES Measurement/Jessie Data/Data/Recoded_binary_imputed_mean-2021-07-02.Rdata")


#######################################################################
### Combine time points within the same period (from period_collapsed.R)
#######################################################################

### load data ###
load("/Users/jiaxliu/Dropbox (Partners HealthCare)/ALSPAC- Russell Sage/SES Measurement/Jessie Data/Data/Recoded_binary_imputed_mean-2021-07-02.Rdata")


## create period-measures: very early (before 3 years), early (3-5 years), middle childhood (6-7 years)

# "period" specific global measures 
jobless_newbin_very_early=c("jbls_newbin_8wks","jbls_newbin_8mos","jbls_newbin_21mos","jbls_newbin_33mos")
jobless_newbin_early=c("jbls_newbin_47mos","jbls_newbin_61mos")
jobless_newbin_middle=c("jbls_newbin_73mos")

income_reduction_newbin_very_early=c("inrdx_newbin_8wks","inrdx_newbin_8mos","inrdx_newbin_21mos","inrdx_newbin_33mos")
income_reduction_newbin_early=c("inrdx_newbin_47mos","inrdx_newbin_61mos")
income_reduction_newbin_middle=c("inrdx_newbin_73mos")

major_F_newbin_very_early=c("major_F_newbin_8mos","major_F_newbin_33mos")
major_F_newbin_early=c("major_F_newbin_61mos")
major_F_newbin_middle=c("major_F_newbin_73mos")


income_bin_very_early=c("income_bin_33mos")
income_bin_early=c("income_bin_47mos")
income_bin_middle=c("income_bin_85mos")

nbhqual_newbin_very_early=c("nbhqual_newbin_21mos","nbhqual_newbin_33mos")
nbhqual_newbin_early=c("nbhqual_newbin_61mos")
nbhqual_newbin_middle=c("nbhqual_newbin_85mos")


Fscore_newbin_very_early=c("Fscore_newbin_8mos","Fscore_newbin_21mos","Fscore_newbin_33mos")
Fscore_newbin_early=c("Fscore_newbin_61mos")
Fscore_newbin_middle=c("Fscore_newbin_85mos")

period=c("_very_early","_early","_middle")

create.period=function(varname,data) {
  for (i in 1:3) {
    if (length(get(paste0(varname,period[i])))>1) {
      a=rowSums(data[,get(paste0(varname,period[i]))],na.rm=F) # sum the # of exposed time points; people with any missing value will be assgined as NA- 
      b=rowSums(data[,get(paste0(varname,period[i]))],na.rm=T) # sum the # of exposed time points if not missing
      x=ifelse(b >0, 1, ifelse(a==0, 0, NA)) # people with one or more missing time points will get NA, unless reported at least one exposed time point
      SES.period=data.frame(SES.period,x)
    }
    else {
      SES.period=data.frame(SES.period,data[,get(paste0(varname,period[i]))])
    }
    colnames(SES.period)[dim(SES.period)[2]]=paste0(varname,period[i])
  }
  return(SES.period)
}

SES.period=SES.merge.newbin.imp.mean
SES.period=create.period(varname="jobless_newbin",data=SES.merge.newbin.imp.mean)

table(SES.period$jobless_newbin_middle,SES.period$jbls_newbin_73mos)
table(SES.period$jobless_newbin_early,SES.period$jbls_newbin_47mos)

SES.period=create.period(varname="income_reduction_newbin",data=SES.merge.newbin.imp.mean)

SES.period=create.period(varname="major_F_newbin",data=SES.merge.newbin.imp.mean)
SES.period=create.period(varname="income_bin",data=SES.merge.newbin.imp.mean)

SES.period=create.period(varname="nbhqual_newbin",data=SES.merge.newbin.imp.mean)
SES.period=create.period(varname="Fscore_newbin",data=SES.merge.newbin.imp.mean)


# I didn't resave this on 2021-07-02
#save(SES.period,file="/Users/JiaxuanLIU/Dropbox (Partners HealthCare)/ALSPAC- Russell Sage/SES Measurement/Jessie Data/Data/SES_period_040119.Rdata")

############ Compare the data with those I previously prepared (the version used for analysis)
#### By far, the period-collapsed exposure variables are all ready.
## rename the dataset name (since the old version has the same name)
SES.period.2021 <- SES.period
dim(SES.period.2021) #971 296


# load the dataset previouly prepared
load ("/Users/Jiaxliu/Dropbox (Partners HealthCare)/ALSPAC- Russell Sage/SES Measurement/Jessie Data/Data/SES_period_040119.Rdata")
dim(SES.period) #971 320 (I deleted many redundant codes today for those variables that were not used in analysis, that's why the old version has many more variables )

# check if the period-collapsed exposure variables are the same
var_list <- c("jobless_newbin_very_early", "jobless_newbin_early","jobless_newbin_middle",
              "income_reduction_newbin_very_early", "income_reduction_newbin_early","income_reduction_newbin_middle",
              "major_F_newbin_very_early", "major_F_newbin_early","major_F_newbin_middle",
              "income_bin_very_early", "income_bin_early","income_bin_middle",
              "nbhqual_newbin_very_early", "nbhqual_newbin_early","nbhqual_newbin_middle",
              "Fscore_newbin_very_early", "Fscore_newbin_early","Fscore_newbin_middle",
              )

all (SES.period.2021[,var_list] == SES.period[,var_list], na.rm = T) #Yep@ They are identifcal.

# I didn't save today's version since the data is the same with the existing one; (this script is a cleaned version of my old scripts; just to check I coded everything correctly ) 
# The codes below used the old version, and I didn't update it.

#######################################################################
### Prepare for analysis (from prepare_phenotype.R)
#######################################################################


load("/Users/jiaxliu/Dropbox (Partners HealthCare)/ALSPAC- Russell Sage/SES Measurement/Jessie Data/Data/SES_period_040419.Rdata")

pheno.period=SES.period[,c(1,303:320)] # keep only the combined period SEP exposures

load("/Users/jiaxliu/Dropbox (Partners HealthCare)/ALSPAC- Russell Sage/Data/ARIES_SES_20190506.Rdata")
pheno.period=merge(pheno.period, ARIES_SES[,-(13:149)],by="cidB1471") # get covariates from the original dataset
pheno.period$ID<-paste0(pheno.period$cidB1471,pheno.period$qlet)

pheno.period$SES_parent=as.factor(pheno.period$SES_parent)

table(pheno.period$SES_parent,exclude=F)
table(pheno.period$WHITE,exclude=F)
table(pheno.period$Female,exclude=F)
table(pheno.period$sustained.smoke,exclude=F)
table(pheno.period$mom_birthage,exclude=F)
table(pheno.period$ppregnum,exclude=F)


pheno.period$WHITE=as.factor(ifelse(str_trim(pheno.period$WHITE)==".",NA, pheno.period$WHITE)) # fix the coding for missingness

pheno.period$Female=as.factor(pheno.period$Female)
pheno.period$sustained.smoke=as.factor(pheno.period$sustained.smoke)
pheno.period$mom_birthage=as.factor(pheno.period$mom_birthage)
pheno.period$ppregnum=as.factor(pheno.period$ppregnum)


# I didn't resave this on 2021-07-02
#save(pheno.period,file="/Users/JiaxuanLIU/Dropbox (Partners HealthCare)/ALSPAC- Russell Sage/SES Measurement/Jessie Data/Data/period_phenotypes_050719.Rdata")

# Check twin ID
load("/Users/jiaxliu/Dropbox (Partners HealthCare)/ALSPAC- Russell Sage/SES Measurement/Jessie Data/Data/ALSPAC-twinsIDs204.Rdata")
all(pheno.period$cidB1471 %in% twins.IDs==FALSE) # This data contains no twins


