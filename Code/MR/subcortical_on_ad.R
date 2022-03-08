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
library(knitr)

setwd("~/PRS_AD_BRAIN_STRUCTURES/prs_alzheimers_brain_structures/Data/Subcortical_GWAS/Full_sumstats/")
temp = list.files( )
for (i in 1:length(temp)) assign(temp[i], fread(temp[i]))
icv <- read.table("CHARGE_ENIGMA_ICV_METAANALYSIS_201311141.txt",he=T)
icv <- icv[,c("RSNUMBERS","N","Allele1","Allele2","Zscore","Freq1","P.value")]
hippocampal <-read.table("CHARGE_ENIGMA_HV_METAANALYSIS_201311141.txt",he=T)
hippocampal<- hippocampal[,c("RSNUMBERS","N","Allele1","Allele2","Zscore","Freq1","P.value")]
colnames(hippocampal) <- c("rsid","N","A1","A2","Zscore","FreqA1","P")
colnames(icv) <- c("rsid","N","A1","A2","Zscore","FreqA1","P")
accumbens <- accumbens_eur_z_ldsc_restricted_NG05SEP19[,c("rsid","N","A1","A2","Zscore","FreqA1","P")]
amygdala <- amygdala_eur_z_ldsc_restricted_NG05SEP19[,c("rsid","N","A1","A2","Zscore","FreqA1","P")]
brainstem <- brainstem_eur_z_ldsc_restricted_NG05SEP19[,c("rsid","N","A1","A2","Zscore","FreqA1","P")]
caudate <- caudate_eur_z_ldsc_restricted_NG05SEP19[,c("rsid","N","A1","A2","Zscore","FreqA1","P")]
pallidum <- pallidum_eur_z_ldsc_restricted_NG05SEP19[,c("rsid","N","A1","A2","Zscore","FreqA1","P")]
putamen <- putamen_eur_z_ldsc_restricted_NG05SEP19[,c("rsid","N","A1","A2","Zscore","FreqA1","P")]
thalamus <- thalamus_eur_z_ldsc_restricted_NG05SEP19[,c("rsid","N","A1","A2","Zscore","FreqA1","P")]

temp2 <- list(accumbens,amygdala,brainstem,caudate,pallidum,putamen,thalamus,icv,hippocampal)
names(temp2) <- c("accumbens","amygdala","brainstem","caudate","pallidum","putamen","thalamus","icv","hippocampal")
alzheimers_snps <- c("rs4844610","rs6733839","rs10933431","rs9271058","rs9473117","rs12539172","rs10808026","rs73223431","rs9331896","rs3740688","rs7933202","rs3851179","rs11218343","rs17125924","rs12881735","rs3752246","rs429358","rs7412","rs6024870","rs7920721","rs593742","rs7185636","rs138190086","rs2830500","rs114812713","rs62039712")

##no proxies identified for HLA in satizabal gwas and proxies for WOX rs62039712
#keep alzheimer's snps##
ad_subcortical_snps <- lapply(temp2, function(x) x[x$rsid %in% alzheimers_snps,])
ad_subcortical_snps<- Map(cbind, ad_subcortical_snps, brainStructure= names(ad_subcortical_snps))
ad_subcortical_snps_dataframe <- do.call(rbind,ad_subcortical_snps)
##no proxies identified for HLA in satizabal gwas and proxies for WOX rs62039712

##SAVE EFFECTS ON SUBCORTICAL VOLUMES FOR ALZHEIMER'S SNPS
ad_subcortical_snps_dataframe$FreqA2 <- (1-ad_subcortical_snps_dataframe$FreqA1)
ad_subcortical_snps_dataframe$FreqA1A2 <- (ad_subcortical_snps_dataframe$FreqA1*ad_subcortical_snps_dataframe$FreqA2)

ad_subcortical_snps_dataframe$SE_conv <- sqrt(1/(ad_subcortical_snps_dataframe$N)*2*(ad_subcortical_snps_dataframe$FreqA1A2))
ad_subcortical_snps_dataframe$beta_conv <- (ad_subcortical_snps_dataframe$SE_conv*ad_subcortical_snps_dataframe$Zscore)
ad_subcortical_snps_dataframe$A1 <- toupper(ad_subcortical_snps_dataframe$A1)
ad_subcortical_snps_dataframe$A2 <- toupper(ad_subcortical_snps_dataframe$A2)

##Set GWAS directory
setwd("~/PRS_AD_BRAIN_STRUCTURES/prs_alzheimers_brain_structures/Data/Subcortical_GWAS/AD_subcortical_snps")
write.csv(ad_subcortical_snps_dataframe,"ad_subcortical_snps.csv",quote=F,row.names=F)


##read in exposure
exposure <- read_exposure_data("~/PRS_AD_BRAIN_STRUCTURES/prs_alzheimers_brain_structures/Data/DiscoverySample/alzheimer_alleles_assignments_to_use_satizabal.csv", sep =",", snp_col = "SNP", beta_col = "BETA", se_col = "SE", eaf_col = "EAF", effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P")


setwd("~/PRS_AD_BRAIN_STRUCTURES/prs_alzheimers_brain_structures/Data/Subcortical_GWAS/AD_subcortical_snps")
outcome_dat <- read_outcome_data("ad_subcortical_snps.csv",exposure$SNP,sep=",",snp_col="rsid",phenotype_col="brainStructure",beta_col="beta_conv",se_col="SE_conv", eaf_col="FreqA1", effect_allele_col = "A1", other_allele_col="A2",pval_col="P", samplesize_col="N") 



#harmonising exposure and outcome datasets (making sure alleles aligned, making sure no pallindromic SNPs etc)
dat <- harmonise_data(exposure, outcome_dat) 


#Running the MR
mr_results <- mr(dat)


mr_results <- mr_results[!mr_results$method=="simple mode",]

setwd("~/PRS_AD_BRAIN_STRUCTURES/prs_alzheimers_brain_structures/Results/Subcortical_GWAS/AD_to_Subcortical")
#Putting MR results into a HTML report
mr_report(dat, output_path =paste("Alzheimer's disease","dat$brainStructure"), author="Analyst", study = paste("Alzheimer's disease","dat$brainStructure",sep=""))
results<-cbind.data.frame(mr_results$outcome,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval)

write.csv(mr_results,"mr_results_ad_to_subcortical_structures.txt",quote=F,row.names = F)



#prepare columns for steiger filtering
dat<- harmonise_data(exposure,outcome_dat)
dat <- dat[dat$mr_keep=="TRUE",]
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
dat$ncase.exposure[dat$SNP=="rs138190086"]=59163
dat$ncontrol.exposure[dat$SNP=="rs138190086"]=35274
dat$samplesize.exposure[dat$SNP=="rs138190086"]=94437
dat$units.exposure <- "log odds"
dat$units.outcome <- "SD"
dat$prevalence.exposure=0.07


#steiger filtering
dat2<- steiger_filtering(dat)
write.table(dat2,"steiger_filtering_ad_to_subcortical.csv", row.names=F,quote=F)

#pleiotropy test
pleiotropy <- mr_pleiotropy_test(dat)
pleiotropy <- pleiotropy[!duplicated(pleiotropy),]
write.csv(pleiotropy,"pleiotropy_mr.csv",row.names=F,quote=F)

#heterogeneity test
heterogeneity <- mr_heterogeneity(dat)
heterogeneity <- heterogeneity[heterogeneity$method=="Inverse variance weighted",]
write.csv(heterogeneity,"heterogeneity_mr.csv",row.names=FALSE,quote=FALSE)

# Calculate F statistics
dat$F   = BXG^2/seBetaXG^2


#Save results
write.csv(dat, file="full_dat.csv", row.names = FALSE)


