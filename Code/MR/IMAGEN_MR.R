#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR")
#devtools::install_github("MRCIEU/MRInstruments")
#devtools::install_github("MRCIEU/xlsx")
#install.packages("xlsx")
#install.packages("psych")
#install.packages("plyr")
#install.packages("ggplot2")
#install.packages("microbenchmark")

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
##############################################################################################################################
########################################## MR  ################################################################
##############################################################################################################################

## OR, instead of the above you can use an MR-BASE exposure file. The ID in quotation marks is the ID of the exposure you want in the MR-BASE catalogue. You can find these on the MR-BASE website 

# get all the available outcomes from MR-BASE
ao <- available_outcomes(access_token=NULL)


#Setting your working directory where everything will be called from and also saved
rm(list=ls())

exposure <- read_exposure_data("~/BRAIN_ANALYSES/IMAGEN/DATA/ALZHEIMERS_GWAS/alzheimer_alleles_assignments_to_use_in_mr.csv",sep=",",snp_col = "SNP", beta_col = "BETA", se_col = "SE", eaf_col = "EAF", effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",ncase_col = "ncase.exposre",ncontrol_col = "ncontrol.exposure")
#naming exposure 
exposure$exposure<- "Alzheimer's disease"
#seeing how many rows are int he datset
head(exposure)
dim(exposure)

setwd("~/BRAIN_ANALYSES/IMAGEN/DATA/BRAIN_GWAS/Summary_stats")
filelist<- list.files()
l<- list()
for(i in filelist)
{
  l[[i]] <- read_outcome_data(filename=i, phenotype_col = "Outcome",beta_col ="BETA",se_col="SE", eaf_col="EAF", effect_allele_col = "A1", other_allele_col = "A2",samplesize_col = "N",pval_col = "P")
  l[[i]]$Outcome <- i
}
outcome_dat <- bind_rows(l)
outcome_dat$Outcome<- gsub(pattern = "\\.txt$", "", outcome_dat$Outcome)

#harmonising exposure and outcome datasets (making sure alleles aligned, making sure no pallindromic SNPs etc)
dat <- harmonise_data(exposure, outcome_dat) 
dat$outcome<- dat$Outcome


#Running the MR
mr_results <- mr(dat)

#Remove MR analysis of simple mode
mr_results <- mr_results[!mr_results$method=="simple mode",]

#displayng MR results 
mr_results


setwd("~/BRAIN_ANALYSES/IMAGEN/RESULTS/MR/")

#Putting MR results into a HTML report
results<-cbind.data.frame(mr_results$outcome,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval)
write.csv(mr_results,"mr_results_cortical_imagen_cohort.txt",quote=F,row.names = F)



pleiotropy <- mr_pleiotropy_test(dat)
heterogeneity <- mr_heterogeneity(dat)
write.csv(heterogeneity,"heterogeneity_imagen_cortical.csv",row.names=F,quote=F)
write.csv(pleiotropy,"pleiotropy_imagen_cortical.csv",row.names=F,quote=F)


#Putting MR results into a HTML report
mr_report(dat, output_path = "cortical_brain.txt", author="Analyst", study = paste("cortical brain structures.txt","-Alzheimer's disease",sep=""))



#columns required for steiger filtering
dat$samplesize.exposure=82771
dat$ncase.exposure=30344
dat$ncontrol.exposure=52427
dat$ncase.exposure[dat$SNP=="rs7412"]=21982
dat$ncontrol.exposure[dat$SNP=="rs7412"]=41944
dat$samplesize.exposure[dat$SNP=="rs7412"]=63926
dat$ncase.exposure[dat$SNP=="rs429358"]=21982
dat$ncontrol.exposure[dat$SNP=="rs429358"]=41944
dat$samplesize.exposure[dat$SNP=="rs429358"]=63926
dat$ncase.exposure[dat$SNP=="rs62039712"]=24078
dat$ncontrol.exposure[dat$SNP=="rs62039712"]=45820
dat$samplesize.exposure[dat$SNP=="rs62039712"]=69898
dat$ncase.exposure[dat$SNP=="rs114812713"]=59163
dat$ncontrol.exposure[dat$SNP=="rs114812713"]=35274
dat$samplesize.exposure[dat$SNP=="rs114812713"]=94437
dat$ncase.exposure[dat$SNP=="rs2830500"]=59163
dat$ncontrol.exposure[dat$SNP=="rs2830500"]=35274
dat$samplesize.exposure[dat$SNP=="rs2830500"]=94437
dat$ncase.exposure[dat$SNP=="rs7920721"]=59163
dat$ncontrol.exposure[dat$SNP=="rs7920721"]=35274
dat$samplesize.exposure[dat$SNP=="rs7920721"]=94437
dat$units.exposure <- "log odds"
dat$units.outcome <- "SD"
dat$prevalence.exposure=0.07
dat$prevalence.outcome=1

##Perform steiger filtering
dat<- steiger_filtering(dat)


#Save results
write.csv(dat, "full_dat_steiger_apoe_imagen.csv", row.names = FALSE,quote=FALSE)
