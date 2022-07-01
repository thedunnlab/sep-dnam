Documentation for the paper "Socioeconomic changes predict genome-wide DNA methylation in childhood".

* Author: 	Jiaxuan (Jessie) Liu
* Last updated: Jun 30, 2022

In this folder, you can find the following information regarding our manuscript.


## 0 - Data Preparation Scripts
This folder includes scripts to create SEP indicators and prepare data for SLCMA analysis. 

* `00-SEP_indicators.R`  Creates the SEP indicators.
* `01-Prepare_Data.R`   Combine SEP data with DNAm data, and prepare the datasets for SLCMA analysis.


## 1 - Main Analysis Scripts
This folder includes scripts for the main SLCMA analysis. 

* `10-SLCMA analyses` This folder includes SLCMA codes.
	* `LARS-mobility-FWL-20200316.R` The SLCMA pipeline codes.
	* `Individual SEP analysis` Codes for individual SEP analysis.
* `11-Examine_Primary_Results_Figure2_FigureS4.R` Examine the main SLCMA analysis results and explore the top hits.
* `12-TableS2-Summary_62hits.R` Create Table S2 which includes summary statistics, hypothesis selected, and gene annotations for the top 62 CpGs with R^2>3%.
* `13-Figure2-FigureS4.R` Plot CpGs by hypothesis selected (Figure 2 and Figure S4).


## 2 - Sensitivity & Secondary Analyses Scripts
This folder includes scripts for the sensitivity and secondary analyses.

* `20-TableS3_Replication.R` Checking potential replication in prior published EWAS of SEP (Table S3).
* `21-TableS5_Additional_Covar_Adj.R` Sensitivity analysis adjusting for additional covariates (Table S5)
* `22-TableS6_mQTL_Analysis.R` Sensitivity analysis adjusting for genetic variation for mQTLs (Table S6)
* `23-Exclude mobility` This folder includes scripts for the sensitivity analysis of excluding mobility hypothesis in SLCMA
	* `Individual SEP analysis` Scripts for individual SEP analysis.
	* `TableS7_Exclude_Mobility.R` Examine analysis results and create Table S7.
* `24-EWAS of ever exposure` This folder includes scripts for the EWAS of ever exposure to SEP indicators.
	* `Individual SEP analysis` Scripts for individual SEP analysis.
	* `Figure_S5_EWAS_Results.R` Examine EWAS results and create Figure S5
* `25-FigureS6_S7_S8_Enrichment.R` Secondary enrichment analysis. Creates Figure S6, S7, and S8.
