# Check for replications in Janine's review paper EWAS studies
library(xlsx)
library(dplyr)


# read in our top 62 CpGs (R2 hits)
top62 <- read.csv(file="/Users/jiaxliu/Dropbox (Partners HealthCare)/ALSPAC- Russell Sage/SES Measurement/Jessie Data/Results/Manuscript/Supplementary/Table2_62hits_annot_2021-05-17.csv")
top_probes <- as.character(top62$CpG)

# read in Janine's list
review_list <- read.csv(file="/Users/jiaxliu/Dropbox (Partners HealthCare)/ALSPAC- Russell Sage/Paper Draft - Review Article/Janine's Data/Summary stats/Jessie's paper/StackedSummaryStats_RussSage_841CpGs_2021-05-18.csv")
table(review_list$SEP)

table(review_list$Study,review_list$SEP)
# dplyr::count_(review_list, vars = c('Study','SEP','Tissue','SEP_Age',"DNAm_AgePer"))
#review_list <- review_list[order(review_list$CpG),]
review_list[review_list$Study=="Suderman" & review_list$SEP=="Composite",]
sum(review_list$FDR<0.05) # 3
min(review_list$FDR)
min(review_list$pval)
review_list[review_list$FDR<0.05,]

# exclude "Neighborhood" and "Other", which correspond to neighborhood quality and Fscore in the BP paper (overlap with our measures, avoid direct comparision considerent the difference in data and methods)
review_list <- review_list[!(review_list$SEP %in% c("Neighborhood","Other")),]
table(review_list$Study)


# Merge
replication_merge <- merge(top62[,c("CpG","Beta_SI","adversity")],review_list[,c("CpG","db","pval","SEP","SEP_Age","Tissue","DNAm_AgePer","Study")],by="CpG")
table(replication_merge$CpG)
length(unique(replication_merge$CpG)) # 62 hits

replication_merge$CpG <- as.character(replication_merge$CpG)
replication_merge <- na.omit(replication_merge) 
table(replication_merge$CpG)


# Check direction of effects and p values

replication_merge$same_sign <- ifelse(replication_merge$Beta_SI*replication_merge$db>0,1,0)
replication_merge$p005 <-ifelse(replication_merge$pval<0.05,1,0)

sum(replication_merge$same_sign)/length(replication_merge$same_sign) #0.5104603
sum(replication_merge$p005)/length(replication_merge$p005) #0.08786611


# CpG-level results:

cpg.order <- read.csv(file = '/Users/jiaxliu/Dropbox (Partners HealthCare)/ALSPAC- Russell Sage/SES Measurement/Jessie Data/Results/Manuscript/Supplementary/62hits_order.csv') # read in the order of the CpGs that we want to arrange the output files

# write a function to calcualte the summary statistics
check_by_CpG <- function(df) {
  res <- data.frame(CpG=character(),prop_same_sign=double(),prop_p005=double(),n_compare=integer(),stringsAsFactors=FALSE)
  i=1
  for (probe in top_probes) {
    temp <- df[df$CpG==probe,]
    res[i,"CpG"] <- probe
    res[i,"n_compare"] <-dim(temp)[1]
    res[i,"prop_same_sign"] <- sum(temp$same_sign)/length(temp$same_sign) # proportion of EWAS having an estiamte in the same direction
    res[i,"prop_p005"] <-sum(temp$p005)/length(temp$p005) # proportion of EWAS having an estiamte with p<0.05
    i=i+1
  }
  res<- res[match(cpg.order$CpG, res$CpG),] # order the CpGs as desired
  
  # add some overall statistics at the end (average across all the CpGs)
  res[i,"CpG"] <- "Average"
  res[i,"prop_same_sign"] <- mean(res[1:62,"prop_same_sign"])
  res[i,"prop_p005"] <-mean(res[1:62,"prop_p005"])
  res[i+1,"CpG"] <- "Fraction"
  res[i+1,"prop_same_sign"] <- sum(ifelse(res[1:62,"prop_same_sign"]>0.5,1,0))/62 # fraction of CpGs showing consistent findings with majority (>50%) of studies 
  res[i+1,"prop_p005"] <-sum(ifelse(res[1:62,"prop_p005"]>0.05,1,0))/62 # fraction of CpGs showing more than 20% significant results
  

  
  return(res)
}

# Run the analsysi for different criteria of comparison.

all_byCpG <- check_by_CpG(replication_merge)
all_byCpG


write.xlsx(all_byCpG, file="/Users/jiaxliu/Dropbox (Partners HealthCare)/ALSPAC- Russell Sage/SES Measurement/Jessie Data/Results/Manuscript/Supplementary/Replication/summary_2021-06-07.xlsx", sheetName="All Studies", row.names=TRUE)




