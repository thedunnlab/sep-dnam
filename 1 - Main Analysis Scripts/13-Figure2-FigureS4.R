# Figure 2
library(dplyr)
library(plyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(scales)

# Use colorblind accessible color scales: friendly : https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible


# FDR <0.05
load('/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/summary_fdr005_4hits.RData')

f2<-cpg_fdr005
f2$facet.flag <- factor(ifelse(f2$hypothesis %in% c("very_early", "early","middle"), 
                               "Sensitive period","Mobility"), 
                        levels=c("Sensitive period","Mobility"))

f2$hypothesis<-as.factor(f2$hypothesis)
f2$hypothesis<-recode_factor(f2$hypothesis,very_early="Very early", early="Early", middle="Middle", mobility_D12="Early improvement",
                             mobility_U12="Early worsening", mobility_D23="Later improvement", mobility_U23="Later worsening")




f2 <- f2[order(f2$measure, decreasing=F), ]

f2 <- f2 %>% 
  group_by(measure, hypothesis, facet.flag) %>% 
  dplyr::summarise(obs = n())

f2$measure <- factor(f2$measure)


adv.labels <- c( "Neighborhood disadvantage")

g2 <- ggplot(f2, aes(x=hypothesis, group=measure)) + 
  geom_bar(aes(y=obs, fill=measure), stat="identity") +
  #scale_fill_manual("measure", values=greys, labels=adv.labels) + 
  #scale_fill_manual("measure", values=cbbPalette, labels=adv.labels) + 
  scale_fill_manual("measure", values=c("#44AA99"), labels=adv.labels) + 
  # scale_fill_discrete("measure", 
  #                    labels=adv.labels, 
  #                     guide = guide_legend(reverse=FALSE)) + 
  facet_grid(.~ facet.flag, scales="free_x", space="free") + 
  xlab("Theoretical Model") + 
  ylab("Number of CpGs with FDR<0.05") +
  theme_bw() + theme(text = element_text(size = 32),
                     axis.text.x = element_text(size=26),
                     axis.title.x = element_blank(),
                     axis.text.y = element_text(size=28),
                     axis.title.y = element_text(size=25),
                     legend.text = element_text(size=28),
                     legend.title = element_blank(),
                     legend.position="top",
                     strip.text = element_text(size=28),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     plot.background = element_rect(fill = "transparent",colour = NA),
                     legend.background = element_rect(fill = "transparent",colour = NA)) +
  ggtitle(" ") + 
  guides(fill=guide_legend(
    keywidth=0.3,
    keyheight=0.5,
    default.unit="inch",
    title="SEP Measure",nrow=2,byrow=TRUE)
  )



# CpGs with R2>0.03 (62 hits)


load('/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/summary_r2g003_62hits.RData')


# Panel A-Instability Measures

f1<-cpg_r2_g003[cpg_r2_g003$measure %in% c("Job loss","Income reduction"),]

f1$facet.flag <- factor(ifelse(f1$hypothesis %in% c("very_early", "early","middle"), 
                               "Sensitive period",ifelse(f1$hypothesis == "accumulation","Accumulation","Mobility")), 
                        levels=c("Sensitive period","Mobility","Accumulation"))

f1$hypothesis<-as.factor(f1$hypothesis)
f1$hypothesis<-recode_factor(f1$hypothesis, very_early="Very early", early="Early", middle="Middle", mobility_D12="Early improvement",
                             mobility_U12="Early worsening", mobility_D23="Later improvement", mobility_U23="Later worsening",accumulation="Accumulation")




f1 <- f1[order(f1$measure, decreasing=F), ]

f1 <- f1 %>% 
  group_by(measure, hypothesis, facet.flag) %>% 
  dplyr::summarise(obs = n())

f1$measure <- factor(f1$measure, levels = c("Job loss", "Income reduction"))


adv.labels <- c("Job loss", 
                "Income reduction")




g1 <- ggplot(f1, aes(x=hypothesis, group=measure)) + 
  geom_bar(aes(y=obs, fill=measure), stat="identity", width = 0.6) +
  #scale_fill_manual("measure", values=greys, labels=adv.labels) + 
  #scale_fill_manual("measure", values=cbbPalette, labels=adv.labels) + 
  scale_fill_manual("measure", values=c("#88CCEE", "#AA4499"), labels=adv.labels) + 
  #scale_fill_discrete("measure", 
  #                    labels=adv.labels, 
  #                    guide = guide_legend(reverse=FALSE)) + 
  facet_grid(.~ facet.flag, scales="free_x", space="free") + 
  xlab("Theoretical Model") + 
  ylab(expression(Number~of~CpGs~with~R^2~">3%")) +
  theme_bw() + theme(text = element_text( size = 32),
                     axis.text.x = element_text(size=28),
                     axis.title.x = element_blank(),
                     axis.text.y = element_text(size=28),
                     axis.title.y = element_text(size=28),
                     legend.text = element_text(size=28),
                     legend.title = element_blank(),
                     legend.position="top",
                     strip.text = element_text(size=28),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     plot.background = element_rect(fill = "transparent",colour = NA),
                     legend.background = element_rect(fill = "transparent",colour = NA)) +
  scale_y_continuous(breaks = pretty_breaks())+
  ggtitle(" ") + 
  guides(fill=guide_legend(
    keywidth=0.3,
    keyheight=0.5,
    default.unit="inch",
    title="SEP Measure")
  )

png("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/figure2/figure2_r2_panelA_2021-05-20.png", 1100, 700,bg = "transparent")
g1
dev.off()

# Panel B-Adversity Measures

f2<-cpg_r2_g003[!(cpg_r2_g003$measure %in% c("Job loss","Income reduction")),]

f2$facet.flag <- factor(ifelse(f2$hypothesis %in% c("very_early", "early","middle"), 
                               "Sensitive period",ifelse(f2$hypothesis == "accumulation","Accumulation","Mobility")), 
                        levels=c("Sensitive period","Mobility","Accumulation"))

f2$hypothesis<-as.factor(f2$hypothesis)
f2$hypothesis<-recode_factor(f2$hypothesis,very_early="Very early", early="Early", middle="Middle", mobility_D12="Early improvement",
                             mobility_U12="Early worsening", mobility_D23="Later improvement", mobility_U23="Later worsening",accumulation="Accumulation")





f2 <- f2[order(f2$measure, decreasing=F), ]

f2 <- f2 %>% 
  group_by(measure, hypothesis, facet.flag) %>% 
  dplyr::summarise(obs = n())

f2$measure <- factor(f2$measure, levels = c("Income","Financial difficulty score",  "Major financial problem", "Neighborhood disadvantage"))


adv.labels <- c("Low family income", "Finanical hardship", 
                "Major financial problem", 
                "Neighborhood disadvantage")


g2 <- ggplot(f2, aes(x=hypothesis, group=measure)) + 
  geom_bar(aes(y=obs, fill=measure), stat="identity", width = 0.6) +
  #scale_fill_manual("measure", values=greys, labels=adv.labels) + 
  #scale_fill_manual("measure", values=cbbPalette, labels=adv.labels) + 
  scale_fill_manual("measure", values=c("#CC6677","#DDCC77","#332288","#44AA99"), labels=adv.labels) + 
  # scale_fill_discrete("measure", 
  #                    labels=adv.labels, 
  #                     guide = guide_legend(reverse=FALSE)) + 
  facet_grid(.~ facet.flag, scales="free_x", space="free") + 
  xlab("Theoretical Model") + 
  ylab(expression(Number~of~CpGs~with~R^2~">3%")) +
  theme_bw() + theme(text = element_text( size = 32),
                     axis.text.x = element_text(size=28),
                     axis.title.x = element_blank(),
                     axis.text.y = element_text(size=28),
                     axis.title.y = element_text(size=28),
                     legend.text = element_text(size=28),
                     legend.title = element_blank(),
                     legend.position="top",
                     strip.text = element_text(size=28),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     plot.background = element_rect(fill = "transparent",colour = NA),
                     legend.background = element_rect(fill = "transparent",colour = NA)) +
  scale_y_continuous(breaks = pretty_breaks())+
  ggtitle(" ") + 
  guides(fill=guide_legend(
    keywidth=0.3,
    keyheight=0.5,
    default.unit="inch",
    title="SEP Measure")
  )

png("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/figure2/figure2_R2_panelB_2021-05-19.png", 1800, 700,bg = "transparent")
g2
dev.off()




#### Genome-wide ##
load('/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Data/Manuscript/results_all.RData')

sI.summary <- rbind(data.frame(Probe = results.Fscore.sI$Probe, hypothesis = results.Fscore.sI$`1.hypo`, 
                               r2 = results.Fscore.sI$`1.r2`, p.value = results.Fscore.sI$`1.p`, FDR_q= p.adjust(results.Fscore.sI$`1.p`, method = "BH"), 
                               expected = c(1:412956)/412956,measure = "Financial difficulty score"),
                    
                    data.frame(Probe = results.income.sI$Probe, hypothesis = results.income.sI$`1.hypo`, 
                               r2 = results.income.sI$`1.r2`, p.value = results.income.sI$`1.p`, FDR_q= p.adjust(results.income.sI$`1.p`, method = "BH"), 
                               expected = c(1:412956)/412956, measure = "Income"),
                    
                    data.frame(Probe = results.income_reduction.sI$Probe, hypothesis = results.income_reduction.sI$`1.hypo`, 
                               r2 = results.income_reduction.sI$`1.r2`, p.value = results.income_reduction.sI$`1.p`, FDR_q= p.adjust(results.income_reduction.sI$`1.p`, method = "BH"),
                               expected = c(1:412956)/412956,   measure = "Income reduction"),
                    
                    data.frame(Probe = results.jobless.sI$Probe, hypothesis = results.jobless.sI$`1.hypo`, 
                               r2 = results.jobless.sI$`1.r2`, p.value = results.jobless.sI$`1.p`, FDR_q= p.adjust(results.jobless.sI$`1.p`, method = "BH"), 
                               expected = c(1:412956)/412956,  measure = "Job loss"),
                    
                    data.frame(Probe = results.major_F.sI$Probe, hypothesis = results.major_F.sI$`1.hypo`, 
                               r2 = results.major_F.sI$`1.r2`, p.value = results.major_F.sI$`1.p`, FDR_q= p.adjust(results.major_F.sI$`1.p`, method = "BH"), 
                               expected = c(1:412956)/412956,  measure = "Major financial problem"),
                    
                    data.frame(Probe = results.nbhqual.sI$Probe, hypothesis = results.nbhqual.sI$`1.hypo`, 
                               r2 = results.nbhqual.sI$`1.r2`, p.value = results.nbhqual.sI$`1.p`, FDR_q= p.adjust(results.nbhqual.sI$`1.p`, method = "BH"), 
                               expected = c(1:412956)/412956, measure = "Neighborhood disadvantage"))


summary(as.factor(sI.summary$hypothesis))
sI.summary$hypothesis <- gsub("Fscore_newbin_","", sI.summary$hypothesis)
sI.summary$hypothesis <- gsub("income_bin_","", sI.summary$hypothesis)
sI.summary$hypothesis <- gsub("income_reduction_newbin_","", sI.summary$hypothesis)
sI.summary$hypothesis <- gsub("jobless_newbin_","", sI.summary$hypothesis)
sI.summary$hypothesis <- gsub("major_F_newbin_","", sI.summary$hypothesis)
sI.summary$hypothesis <- gsub("nbhqual_newbin_","", sI.summary$hypothesis)

sI.summary$hypothesis <- factor(sI.summary$hypothesis, levels = c("accumulation", "very_early","early","middle",
                                                                  "mobility_D12", "mobility_D23", "mobility_U12", "mobility_U23"))


# Panel A-Instability Measures

f1<-sI.summary[sI.summary$measure %in% c("Job loss","Income reduction"),]

f1$facet.flag <- factor(ifelse(f1$hypothesis %in% c("very_early", "early","middle"), 
                               "Sensitive period",ifelse(f1$hypothesis == "accumulation","Accumulation","Mobility")), 
                        levels=c("Sensitive period","Mobility","Accumulation"))

f1$hypothesis<-as.factor(f1$hypothesis)
f1$hypothesis<-recode_factor(f1$hypothesis, very_early="Very early", early="Early", middle="Middle", mobility_D12="Early improvement",
                             mobility_U12="Early worsening", mobility_D23="Later improvement", mobility_U23="Later worsening",accumulation="Accumulation")







f1 <- f1[order(f1$measure, decreasing=F), ]

f1 <- f1 %>% 
  group_by(measure, hypothesis, facet.flag) %>% 
  dplyr::summarise(obs = n())

f1$measure <- factor(f1$measure, levels = c("Job loss", "Income reduction"))


adv.labels <- c("Job loss", 
                "Income reduction")



g1 <- ggplot(f1, aes(x=hypothesis, group=measure)) + 
  geom_bar(aes(y=obs, fill=measure), stat="identity") +
  #scale_fill_manual("measure", values=greys, labels=adv.labels) + 
  #scale_fill_manual("measure", values=cbbPalette, labels=adv.labels) + 
  scale_fill_manual("measure", values=c("#88CCEE", "#AA4499"), labels=adv.labels) + 
  #scale_fill_discrete("measure", 
  #                    labels=adv.labels, 
  #                    guide = guide_legend(reverse=FALSE)) + 
  facet_grid(.~ facet.flag, scales="free_x", space="free") + 
  xlab("Theoretical Model") + 
  ylab("Number of CpGs across epigenome") +
  scale_y_continuous(labels = scales::comma)+
  theme_bw() + theme(text = element_text( size = 32),
                     axis.text.x = element_text(size=26),
                     axis.title.x = element_blank(),
                     axis.text.y = element_text(size=28),
                     axis.title.y = element_text(size=25),
                     legend.text = element_text(size=28),
                     legend.title = element_blank(),
                     legend.position="top",
                     strip.text = element_text(size=28),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     plot.background = element_rect(fill = "transparent",colour = NA),
                     legend.background = element_rect(fill = "transparent",colour = NA)) +
  ggtitle(" ") + 
  guides(fill=guide_legend(
    keywidth=0.3,
    keyheight=0.5,
    default.unit="inch",
    title="SEP Measure")
  )

png("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/figure2/figure2_genome_panelA_2021-04-25.png", 1300, 700,bg = "transparent")
g1
dev.off()

# Panel B-Adversity Measures

f2<-sI.summary[!(sI.summary$measure %in% c("Job loss","Income reduction")),]

f2$facet.flag <- factor(ifelse(f2$hypothesis %in% c("very_early", "early","middle"), 
                               "Sensitive period",ifelse(f2$hypothesis == "accumulation","Accumulation","Mobility")), 
                        levels=c("Sensitive period","Mobility","Accumulation"))

f2$hypothesis<-as.factor(f2$hypothesis)
f2$hypothesis<-recode_factor(f2$hypothesis,very_early="Very early", early="Early", middle="Middle", mobility_D12="Early improvement",
                             mobility_U12="Early worsening", mobility_D23="Later improvement", mobility_U23="Later worsening",accumulation="Accumulation")





f2 <- f2[order(f2$measure, decreasing=F), ]

f2 <- f2 %>% 
  group_by(measure, hypothesis, facet.flag) %>% 
  dplyr::summarise(obs = n())

f2$measure <- factor(f2$measure, levels = c("Income","Financial difficulty score",  "Major financial problem", "Neighborhood disadvantage"))


adv.labels <- c("Low family income", "Finanical hardship", 
                "Major financial problem", 
                "Neighborhood disadvantage")


g2 <- ggplot(f2, aes(x=hypothesis, group=measure)) + 
  geom_bar(aes(y=obs, fill=measure), stat="identity") +
  #scale_fill_manual("measure", values=greys, labels=adv.labels) + 
  #scale_fill_manual("measure", values=cbbPalette, labels=adv.labels) + 
  scale_fill_manual("measure", values=c("#CC6677","#DDCC77","#332288","#44AA99"), labels=adv.labels) + 
  # scale_fill_discrete("measure", 
  #                    labels=adv.labels, 
  #                     guide = guide_legend(reverse=FALSE)) + 
  facet_grid(.~ facet.flag, scales="free_x", space="free") + 
  xlab("Theoretical Model") + 
  ylab("Number of CpGs across epigenome") +
  scale_y_continuous(labels = scales::comma)+
  theme_bw() + theme(text = element_text( size = 32),
                     axis.text.x = element_text(size=26),
                     axis.title.x = element_blank(),
                     axis.text.y = element_text(size=28),
                     axis.title.y = element_text(size=25),
                     legend.text = element_text(size=28),
                     legend.title = element_blank(),
                     legend.position="top",
                     strip.text = element_text(size=28),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     plot.background = element_rect(fill = "transparent",colour = NA),
                     legend.background = element_rect(fill = "transparent",colour = NA)) +
  ggtitle(" ") + 
  guides(fill=guide_legend(
    keywidth=0.3,
    keyheight=0.5,
    default.unit="inch",
    title="SEP Measure",nrow=2,byrow=TRUE)
  )

png("/data/js95/ALSPAC/ARIES/Russell Sage/Jessie/Results/Manuscript/figure2/figure2_genome_panelB_2021-04-25.png", 2000, 810,bg = "transparent")
g2
dev.off()


