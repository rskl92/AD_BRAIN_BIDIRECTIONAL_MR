library(meta)
setwd("~/BRAIN_ANALYSES/Metaanalysis")
meta_analysis_data <- read_excel("children_dataset_for_metanalysing.xlsx")
meta_analysis_data<- meta_analysis_data[meta_analysis_data$Method=="Inverse variance weighted",]
metaanalysis <- metagen(TE,seTE,data=meta_analysis_data,studlab = Study,comb.random=TRUE, method.tau="DL",sm="SMD")
metanalysis <- update.meta(metaanalysis,byvar=meta_analysis_data$Outcome,comb.random=TRUE,comb.fixed=FALSE)
setwd("~/BRAIN_ANALYSES/Metaanalysis/Childhood_metanalysis")
sink("results_metaanalysis_childhood240222_ivw.txt")
print(metanalysis)
sink()



setwd("~/BRAIN_ANALYSES/Metaanalysis")
meta_analysis_data <- read_excel("children_dataset_for_metanalysing.xlsx")
meta_analysis_data<- meta_analysis_data[meta_analysis_data$Method=="MR Egger",]
metaanalysis <- metagen(TE,seTE,data=meta_analysis_data,studlab = Study,comb.random=TRUE, method.tau="DL",sm="SMD")
metanalysis <- update.meta(metaanalysis,byvar=meta_analysis_data$Outcome,comb.random=TRUE,comb.fixed=FALSE)
setwd("~/BRAIN_ANALYSES/Metaanalysis/Childhood_metanalysis")
sink("results_metaanalysis_childhood240222_mregger.txt")
print(metanalysis)
sink()


setwd("~/BRAIN_ANALYSES/Metaanalysis")
meta_analysis_data <- read_excel("children_dataset_for_metanalysing.xlsx")
meta_analysis_data<- meta_analysis_data[meta_analysis_data$Method=="Weighted median",]
metaanalysis <- metagen(TE,seTE,data=meta_analysis_data,studlab = Study,comb.random=TRUE, method.tau="DL",sm="SMD")
metanalysis <- update.meta(metaanalysis,byvar=meta_analysis_data$Outcome,comb.random=TRUE,comb.fixed=FALSE)
setwd("~/BRAIN_ANALYSES/Metaanalysis/Childhood_metanalysis")
sink("results_metaanalysis_childhood240222_wmedian.txt")
print(metanalysis)
sink()

setwd("~/BRAIN_ANALYSES/Metaanalysis")
meta_analysis_data <- read_excel("children_dataset_for_metanalysing.xlsx")
meta_analysis_data<- meta_analysis_data[meta_analysis_data$Method=="Weighted mode",]
metaanalysis <- metagen(TE,seTE,data=meta_analysis_data,studlab = Study,comb.random=TRUE, method.tau="DL",sm="SMD")
metanalysis <- update.meta(metaanalysis,byvar=meta_analysis_data$Outcome,comb.random=TRUE,comb.fixed=FALSE)
setwd("~/BRAIN_ANALYSES/Metaanalysis/Childhood_metanalysis")
sink("results_metaanalysis_childhood240222_wmode.txt")
print(metanalysis)
sink()


setwd("~/BRAIN_ANALYSES/Metaanalysis")
meta_analysis_data <- read_excel("children_dataset_for_metanalysing_noapoe.xlsx")
metaanalysis <- metagen(TE,seTE,data=meta_analysis_data,studlab = Study,comb.random=TRUE,method.tau = "DL",sm="SMD")
metanalysis <- update.meta(metaanalysis,byvar=meta_analysis_data$Outcome,comb.random=TRUE,comb.fixed=FALSE)
setwd("~/BRAIN_ANALYSES/Metaanalysis/Childhood_metanalysis")
sink("results_metaanalysis_updated_noapoe_childhood28022.txt")
print(metanalysis)
sink()
