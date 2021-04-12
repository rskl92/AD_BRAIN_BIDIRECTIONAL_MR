# AD_BRAIN_BIDIRECTIONAL_MR

## Analysis of the polygenic risk score for Alzheimer's disease on brain structures in different cohorts of different ages

This repository accompanies the paper:
Korologou-Linden RS, et al. The bidirectional causal effects of brain structures across the life course on risk of Alzheimer’s disease: A cross-cohort comparison and Mendelian randomization meta-analysis

# Environment details
I use the following language versions: R-3.3.1-ATLAS, plink xx/

The code uses some environment variables that need to be set in your linux environment.

I set the results directory and project data directory temporarily with:

export RES_DIR="${HOME}/PRS_AD_BRAIN_STRUCTURES"
export PROJECT_DATA="${HOME}PRS_AD_BRAIN_STRUCTURES/Data"
export UKB_PHENO="/path/to/UKBapplication/"

I set the IEU shared ALSPAC and UK Biobank data directories by adding the following to my ~/.bash_profile file

export UKB_DATA="/path/to/ukb/data"
export ALSPAC_DATA="/path/to/ALSPAC/data"



