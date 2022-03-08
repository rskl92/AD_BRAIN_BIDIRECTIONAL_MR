R
install.packages("devtools")
devtools::install_github("MRCIEU/TwoSampleMR")
devtools::install_github("MRCIEU/MRInstruments")
devtools::install_github("MRCIEU/xlsx")
install.packages("xlsx")
install.packages("psych")
install.packages("plyr")
install.packages("ggplot2")
install.packages("microbenchmark")
install.packages("data.table")

## This calls the packages you just installed so that you can use them ##
library(devtools)
library(TwoSampleMR)
library(MRInstruments)
library(xlsx)
library(psych)
library(ggplot2)
library(plyr)
library(dplyr)
library(microbenchmark)
library(data.table)

##Set working directory
setwd("~/PRS_AD_BRAIN_STRUCTURES/prs_alzheimers_brain_structures/Data/Grasby_cerebral_cortex/Full_sumstats/ENIGMA3_withGlobal/")
temp = list.files()
##Create a list of files and subset genome-wide significant snps, only select regions that have SNPs
for (i in 1:length(temp)) assign(temp[i], fread(temp[i],quote=""))
datalist_thickness <- list(ENIGMA3_mixed_se_wo_Mean_Full_Thickness_20190429.txt,ENIGMA3_mixed_se_wTHICK_Mean_bankssts_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_caudalanteriorcingulate_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_caudalmiddlefrontal_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_cuneus_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_entorhinal_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_frontalpole_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_fusiform_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_inferiorparietal_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_inferiortemporal_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_insula_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_isthmuscingulate_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_lateraloccipital_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_lateralorbitofrontal_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_lingual_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_medialorbitofrontal_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_middletemporal_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_paracentral_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_parahippocampal_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_parsopercularis_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_parsorbitalis_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_parstriangularis_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_pericalcarine_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_postcentral_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_posteriorcingulate_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_precentral_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_precuneus_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_rostralanteriorcingulate_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_rostralmiddlefrontal_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_superiorfrontal_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_superiorparietal_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_superiortemporal_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_supramarginal_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_temporalpole_thickavg_20200522.txt,ENIGMA3_mixed_se_wTHICK_Mean_transversetemporal_thickavg_20200522.txt)
names(datalist_thickness) <- c("ENIGMA3_mixed_se_wo_Mean_Full_Thickness","ENIGMA3_mixed_se_wTHICK_Mean_bankssts_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_caudalanteriorcingulate_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_caudalmiddlefrontal_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_cuneus_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_entorhinal_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_frontalpole_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_fusiform_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_inferiorparietal_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_inferiortemporal_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_insula_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_isthmuscingulate_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_lateraloccipital_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_lateralorbitofrontal_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_lingual_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_medialorbitofrontal_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_middletemporal_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_paracentral_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_parahippocampal_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_parsopercularis_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_parsorbitalis_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_parstriangularis_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_pericalcarine_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_postcentral_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_posteriorcingulate_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_precentral_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_precuneus_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_rostralanteriorcingulate_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_rostralmiddlefrontal_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_superiorfrontal_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_superiorparietal_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_superiortemporal_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_supramarginal_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_temporalpole_thickavg","ENIGMA3_mixed_se_wTHICK_Mean_transversetemporal_thickavg")
datasets_cerebral_thickness_snps <- lapply(datalist_thickness, function(x) x[x$P<=5e-08,])
datasets_cerebral_thickness_snps <- datasets_cerebral_thickness_snps[sapply(datasets_cerebral_thickness_snps,function(x) nrow(x)[1])>0]

##Bind datasets
my_list_thickness <- Map(cbind, datasets_cerebral_thickness_snps, brainStructure= names(datasets_cerebral_thickness_snps))
my_list_thickness_full <- Map(cbind, datalist_thickness, brainStructure= names(datasets_cerebral_thickness_snps))

##Bind all list to make one dataset for all regions
outcome_data_thickness <- do.call(rbind,my_list_thickness)

##Standardise coefficients by Grasby, using the formula by Chauhan et al (2015)
outcome_data_thickness$Zscore <- outcome_data_thickness$BETA1/outcome_data_thickness$SE
outcome_data_thickness$FreqA2 <- (1-outcome_data_thickness$FREQ1)
outcome_data_thickness$FreqA1A2 <- (outcome_data_thickness$FREQ1*outcome_data_thickness$FreqA2)
outcome_data_thickness$FreqA1A2_X2 <- 2*(outcome_data_thickness$FREQ1*outcome_data_thickness$FreqA2)
outcome_data_thickness$SE_conv <- sqrt(1/(outcome_data_thickness$N*outcome_data_thickness$FreqA1A2_X2))
outcome_data_thickness$beta_conv <- (outcome_data_thickness$SE_conv*outcome_data_thickness$Zscore)
outcome_data_thickness$A1 <- toupper(outcome_data_thickness$A1)
outcome_data_thickness$A2 <- toupper(outcome_data_thickness$A2)


##SAVE EFFECTS ON CORTICAL SURFACE AREA FOR ALZHEIMER'S SNPS
setwd("~/PRS_AD_BRAIN_STRUCTURES/prs_alzheimers_brain_structures/Data/Grasby_cerebral_cortex/Cortical_AD_snp_effects")
by(outcome_data_thickness, outcome_data_thickness$brainStructure, FUN=function(i) write.csv(i, paste0(i$brainStructure[1],".csv"),quote=F,row.names=F))



##read in exposure (cortical regions)
setwd("~/PRS_AD_BRAIN_STRUCTURES/prs_alzheimers_brain_structures/Data/Grasby_cerebral_cortex/Cortical_AD_snp_effects/Thickness")
filelist<- list.files()
l<- list()
for(i in filelist)
{
  l[[i]] <- read_exposure_data(filename=i,sep=",",snp_col="SNP",phenotype_col="brainStructure",beta_col="beta_conv",se_col="SE_conv", eaf_col="FREQ1", effect_allele_col = "A1", other_allele_col="A2",pval_col="P", samplesize_col="N",clump=T) 
  l[[i]]$Outcome <- i
}
exposure_dat <- bind_rows(l)

##Read in outcome (Alzheimer's disease)
outcome_dat <- read_outcome_data("~/PRS_AD_BRAIN_STRUCTURES/prs_alzheimers_brain_structures/Data/DiscoverySample/Kunkle_etal_2019_IGAP_Summary_statistics.with_allelefreqs.txt",exposure_dat$SNP,snp_col = "MarkerName", beta_col = "Beta", se_col = "SE", effect_allele_col = "Effect_allele", other_allele_col = "Non_Effect_allele", pval_col = "Pvalue",eaf="Effect_allele_freq")


#harmonising exposure and outcome datasets (making sure alleles aligned, making sure no pallindromic SNPs etc)
dat <- harmonise_data(exposure_dat, outcome_dat) 



#Running the MR
mr_results <- mr(dat)

#Remove MR simple mode estimates
mr_results <- mr_results[!mr_results$method=="simple mode",]

#Set results directory
setwd("~/PRS_AD_BRAIN_STRUCTURES/prs_alzheimers_brain_structures/Results/CORTICAL_GWAS_GRASBY/CORTICAL_to_AD/WITH_OUTCOME_EAFS")
#Putting MR results into a HTML report
mr_report(dat, output_path =paste("dat$brainStructure",""), author="Analyst", study = paste("brain structures.txt","-Alzheimer's disease",sep=""))
results<-cbind.data.frame(mr_results$outcome,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval)

#Write results to a file
write.csv(mr_results,"mr_results_cortical_to_ad_cortical_withglobadj_witheafs.txt",quote=F,row.names = F)


#Prepare columns for Steiger filtering
dat$exposure <- dat$Outcome
dat$ncase.outcome=21982
dat$ncontrol.outcome=41944
dat$samplesize.outcome=63926	
dat$units.outcome <- "log odds"
dat$units.exposure <- "SD"
dat$prevalence.outcome=0.07
dat$eaf.outcome=dat$eaf.outcome
dat2<- steiger_filtering(dat)

directionality_test <- directionality_test(dat)

#pleiotropy test
pleiotropy <- mr_pleiotropy_test(dat)
pleiotropy <- pleiotropy[!duplicated(pleiotropy),]
write.csv(pleiotropy,"pleiotropy_mr.csv",row.names=F,quote=F)

#heterogeneity test
heterogeneity <- mr_heterogeneity(dat)
heterogeneity <- heterogeneity[heterogeneity$method=="Inverse variance weighted",]
write.csv(heterogeneity,"heterogeneity_mr.csv",row.names=F,quote=F)

#Directionality test
directionality <- directionality_test(dat2)
write.table(directionality,"directionality_Test_cortex_to_Ad.txt",quote=FALSE, row.names=FALSE)


# Calculate F statistics
dat2$F   = dat2$BXG^2/dat2$seBetaXG^2


#Write dataset
write.table(dat2,"steiger_filtering_cortical_to_ad_eafs.csv", row.names=F,quote=F)
