setwd("~/Working/GitHub/twin-data-models")
rm(list = ls())

# Load libraries
library(umx)
library(tidyverse)

# Set up session 01 data
famData.ses01 <- readRDS("QTAB_familywise_ses01.RDS")
mzData.ses01 <- famData.ses01 %>% filter(zyg == 1 | zyg == 2)
dzData.ses01 <- famData.ses01 %>% filter(zyg == 3 | zyg == 4 | zyg == 5)

#### Variance-based ####
# See Verhulst et al., 2019. https://doi.org/10.1007/s10519-018-9942-y
# May produce negative variance component estimates
# Include age & sex in the model as covariates (note: umx will drop rows with missing covariates)
m1vc <- umxACEv(selDVs = "SMFQ_score", sep = "_0", mzData = mzData.ses01, dzData = dzData.ses01, selCovs = c("ses01_age_months", "sex_female"))
m1vc$top$A_std$result # % variance attributed to A
m1vc$top$C_std$result # % variance attributed to C
m1vc$top$E_std$result # % variance attributed to E
umxSummaryACEv(m1vc, std = T)

# Submodels
m1vc.AE <- umxModify(m1vc, update = "C_r1c1", name = "AE", comparison = TRUE)
m1vc.CE <- umxModify(m1vc, update = "A_r1c1", name = "CE", comparison = TRUE)
m1vc.E <- umxModify(m1vc.CE, update = "C_r1c1", name = "E", comparison = TRUE)
umxCompare(base = m1vc, comparison = c(m1vc.AE, m1vc.CE, m1vc.E))

# Covariates (age, sex) regressed out prior to twin modelling
famData.ses01.resid = umx_residualize(c("SMFQ_score"), cov = c("ses01_age_months", "sex_female"), data = famData.ses01, suffixes = c("_01", "_02")) # regress out age & sex
mzData.ses01.resid <- famData.ses01.resid %>% filter(zyg == 1 | zyg == 2)
dzData.ses01.resid <- famData.ses01.resid %>% filter(zyg == 3 | zyg == 4 | zyg == 5)
m2vc <- umxACEv(selDVs = "SMFQ_score", sep = "_0", mzData = mzData.ses01.resid, dzData = dzData.ses01.resid)
m2vc$top$A_std$result # % variance attributed to A
m2vc$top$C_std$result # % variance attributed to C
m2vc$top$E_std$result # % variance attributed to E
umxSummaryACEv(m2vc, std = T)

# Submodels
m2vc.AE <- umxModify(m2vc, update = "C_r1c1", name = "AE", comparison = TRUE)
m2vc.CE <- umxModify(m2vc, update = "A_r1c1", name = "CE", comparison = TRUE)
m2vc.E <- umxModify(m2vc.CE, update = "C_r1c1", name = "E", comparison = TRUE)
umxCompare(base = m2vc, comparison = c(m2vc.AE, m2vc.CE, m2vc.E))

# Bivariate
m3vc <- umxACEv(selDVs = c("SMFQ_score", "SCAS_score"), sep = "_0", mzData = mzData.ses01, dzData = dzData.ses01, selCovs = c("ses01_age_months", "sex_female"))
umxSummaryACEv(m3vc, std = T, showRg = T)
m3vc$top$A_std$result[1,1] # % variance attributed to A - phenotype 1
m3vc$top$A_std$result[2,2] # % variance attributed to A - phenotype 2
m3vc$top$C_std$result[1,1] # % variance attributed to C - phenotype 1
m3vc$top$C_std$result[2,2] # % variance attributed to C - phenotype 2
m3vc$top$E_std$result[1,1] # % variance attributed to E - phenotype 1
m3vc$top$E_std$result[2,2] # % variance attributed to E - phenotype 2

m3vc.AE <- umxModify(m3vc, update = c("C_r1c1", "C_r2c1", "C_r2c2"), name = "AE", comparison = TRUE)
umxSummaryACEv(m3vc.AE, std = T, showRg = T)

#### Path-based ####
# Include age & sex in the model as covariates (note: umx will drop rows with missing covariates)
m1path <- umxACE(selDVs = "SMFQ_score", sep = "_0", mzData = mzData.ses01, dzData = dzData.ses01, selCovs = c("ses01_age_months", "sex_female"))
m1path$top$A_std$result # % variance attributed to A
m1path$top$C_std$result # % variance attributed to C
m1path$top$E_std$result # % variance attributed to E
umxSummaryACE(m1path, std = T)

## Covariates (age, sex) regressed out prior to twin modelling
famData.ses01.resid = umx_residualize(c("SMFQ_score"), cov = c("ses01_age_months", "sex_female"), data = famData.ses01, suffixes = c("_01", "_02")) # regress out age & sex
mzData.ses01.resid <- famData.ses01.resid %>% filter(zyg == 1 | zyg == 2)
dzData.ses01.resid <- famData.ses01.resid %>% filter(zyg == 3 | zyg == 4 | zyg == 5)
m2path <- umxACE(selDVs = "SMFQ_score", sep = "_0", mzData = mzData.ses01.resid, dzData = dzData.ses01.resid)
m2path$top$A_std$result # % variance attributed to A
m2path$top$C_std$result # % variance attributed to C
m2path$top$E_std$result # % variance attributed to E
umxSummaryACE(m2path, std = T)

# Bivariate
m3path <- umxACE(selDVs = c("SMFQ_score", "SCAS_score"), sep = "_0", mzData = mzData.ses01, dzData = dzData.ses01, selCovs = c("ses01_age_months", "sex_female"))
umxSummaryACE(m3path, std = T, showRg = T)
m3path$top$A_std$result[1,1] # % variance attributed to A - phenotype 1
m3path$top$A_std$result[2,2] # % variance attributed to A - phenotype 2
m3path$top$C_std$result[1,1] # % variance attributed to C - phenotype 1
m3path$top$C_std$result[2,2] # % variance attributed to C - phenotype 2
m3path$top$E_std$result[1,1] # % variance attributed to E - phenotype 1
m3path$top$E_std$result[2,2] # % variance attributed to E - phenotype 2

m3path.AE <- umxModify(m3path, update = c("c_r1c1", "c_r2c1", "c_r2c2"), name = "AE", comparison = TRUE)
umxSummaryACE(m3path.AE, std = T, showRg = T)

m3path.AE$top$A_std$result[1,1] # % variance attributed to A - phenotype 1
m3path.AE$top$A_std$result[2,2] # % variance attributed to A - phenotype 2
m3path.AE$top$E_std$result[1,1] # % variance attributed to E - phenotype 1
m3path.AE$top$E_std$result[2,2] # % variance attributed to E - phenotype 1

#### Compare Variance & Path estimates ####
m1vc$top$A_std$result # 0.5078013
m1path$top$A_std$result # 0.3834879

m1vc$top$C_std$result # -0.1129604. Negative C - sampling error, or an ADE model is a better fit? rMZ is more than twice rDZ.
m1path$top$C_std$result # 0

m1vc$top$E_std$result # 0.6051591
m1path$top$E_std$result # 0.6165121

m1vc.ADE <- umxACEv(selDVs = "SMFQ_score", sep = "_0", mzData = mzData.ses01, dzData = dzData.ses01, selCovs = c("ses01_age_months", "sex_female"), dzCr = .25)

m1vc$top$A_std$result # 0.5078013
m1vc.ADE$top$A_std$result # 0.5078013
m1vc.ADE$top$C_std$result # 0.5078013
m1vc.ADE$top$E_std$result # 0.5078013

umxSummarizeTwinData(data = famData.ses01, selVars = "SMFQ_score", sep = "_0", zyg = "zyg", age = "ses01_age_months", MZ = c(1,2), DZ = c(3,4,5))
