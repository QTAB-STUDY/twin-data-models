# Lachlan T Strike
# Based on scripts provided by Hermine Maes (https://hermine-maes.squarespace.com)
# Saturated model used to check assumptions of twin data. Is setup for QTAB data (i.e. 5 zyg groups, with opposite sex twin pairs coded females as 01 and males as 02).
# Take rMZ/DZ & covariate estimates from model with equated means/variances, but separate covariances for MZ and DZ (i.e. CModel2)
# Sex_female - 1 = yes (female), 0 = no (male)
# Expects phenotype data to be Z-score, with "Z" appended to the end of the variable name (if not, comment out 'Get raw means & sd for table' section)
# ###############################

rm(list = ls())
setwd("/qtab/twin-data-models")

#### Load libraries ####
library(OpenMx)
library(stringr)
library(dplyr)
source("miFunctions.R")
mxOption(NULL, "Default optimizer", "CSOLNP") # NPSOL SLSQP CSOLNP

#### Models ####
Saturated_covariate <- function(phenotype, twin.data) {
  covariate <- c("sex_female", "ses01_age_months")
  nc <- length(covariate)

  # OpenMx does not tolerate missing values for definition variables.
  # Recode any missing definition variables as -999
  # BUT!!! Make sure there are not any cases of missing definition variables
  # with a phenotype present
  for (x in covariate) {
    twin01.missing <- twin.data[, paste0(phenotype, "_01")][is.na(twin.data[, paste0(x, "_01")])]
    twin02.missing <- twin.data[, paste0(phenotype, "_02")][is.na(twin.data[, paste0(x, "_02")])]
    stopifnot(is.na(twin01.missing))
    stopifnot(is.na(twin02.missing))
    twin.data[, paste0(x, "_01")][is.na(twin.data[, paste0(x, "_01")])] <- -999
    twin.data[, paste0(x, "_02")][is.na(twin.data[, paste0(x, "_02")])] <- -999
  }

  # Select variables
  selVars <- c(paste0(phenotype, "_01"), paste0(phenotype, "_02"))
  covVars <- c(paste0(covariate, "_01"), paste0(covariate, "_02"))
  useVars <- c(selVars, covVars)

  # Select data for analysis
  mzfData <- subset(twin.data, zyg == 1, useVars)
  mzmData <- subset(twin.data, zyg == 2, useVars)
  dzfData <- subset(twin.data, zyg == 3, useVars)
  dzmData <- subset(twin.data, zyg == 4, useVars)
  dzosData <- subset(twin.data, zyg == 5, useVars)

  # Set Starting Values
  nt <- 2
  nv <- length(phenotype)
  ntv <- nv * nt
  svMe <- mean(unlist(twin.data[, selVars]), na.rm = T)
  svVa <- var(unlist(twin.data[, selVars]), na.rm = T)
  svBe <- 0.01

  # Get raw means & sd for table (needed if using Z-score data)
  selVars.2 <- str_sub(string = phenotype, start = 1, end = -2)
  selVars.3 <- c(paste0(selVars.2, "_01"), paste0(selVars.2, "_02"))
  mean.val <- mean(unlist(twin.data[, c(selVars.3[1], selVars.3[2])]), na.rm = T)
  sd.val <- sd(unlist(twin.data[, c(selVars.3[1], selVars.3[2])]), na.rm = T)

  # ------------------------------------------------------------------------------
  # PREPARE MODEL
  # Saturated Model
  # Create Matrices for Covariates and Linear Regression Coefficients
  defMZF <- mxMatrix(type = "Full", nrow = nc, ncol = ntv, free = F, labels = c(paste0("data.", covariate, "_01"), paste0("data.", covariate, "_02")), name = "defMZF")
  defMZM <- mxMatrix(type = "Full", nrow = nc, ncol = ntv, free = F, labels = c(paste0("data.", covariate, "_01"), paste0("data.", covariate, "_02")), name = "defMZM")
  defDZF <- mxMatrix(type = "Full", nrow = nc, ncol = ntv, free = F, labels = c(paste0("data.", covariate, "_01"), paste0("data.", covariate, "_02")), name = "defDZF")
  defDZM <- mxMatrix(type = "Full", nrow = nc, ncol = ntv, free = F, labels = c(paste0("data.", covariate, "_01"), paste0("data.", covariate, "_02")), name = "defDZM")
  defDZOS <- mxMatrix(type = "Full", nrow = nc, ncol = ntv, free = F, labels = c(paste0("data.", covariate, "_01"), paste0("data.", covariate, "_02")), name = "defDZOS")

  betaMZF <- mxMatrix(type = "Full", nrow = 1, ncol = nc, free = TRUE, values = svBe, labels = c(paste0("beta", covariate)), name = "betaMZF")
  betaMZM <- mxMatrix(type = "Full", nrow = 1, ncol = nc, free = TRUE, values = svBe, labels = c(paste0("beta", covariate)), name = "betaMZM")
  betaDZF <- mxMatrix(type = "Full", nrow = 1, ncol = nc, free = TRUE, values = svBe, labels = c(paste0("beta", covariate)), name = "betaDZF")
  betaDZM <- mxMatrix(type = "Full", nrow = 1, ncol = nc, free = TRUE, values = svBe, labels = c(paste0("beta", covariate)), name = "betaDZM")
  betaDZOS <- mxMatrix(type = "Full", nrow = 1, ncol = nc, free = TRUE, values = svBe, labels = c(paste0("beta", covariate)), name = "betaDZOS")

  # Algebra for expected Mean Matrices in MZ & DZ twins
  meanMZF <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, values = svMe, labels = c("mMZF1", "mMZF2"), name = "meanMZF")
  meanMZM <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, values = svMe, labels = c("mMZM1", "mMZM2"), name = "meanMZM")
  meanDZF <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, values = svMe, labels = c("mDZF1", "mDZF2"), name = "meanDZF")
  meanDZM <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, values = svMe, labels = c("mDZM1", "mDZM2"), name = "meanDZM")
  meanDZOS <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, values = svMe, labels = c("mDZOS1", "mDZOS2"), name = "meanDZOS")

  expMeanMZF <- mxAlgebra(expression = meanMZF + betaMZF %*% defMZF, name = "expMeanMZF")
  expMeanMZM <- mxAlgebra(expression = meanMZM + betaMZM %*% defMZM, name = "expMeanMZM")
  expMeanDZF <- mxAlgebra(expression = meanDZF + betaDZF %*% defDZF, name = "expMeanDZF")
  expMeanDZM <- mxAlgebra(expression = meanDZM + betaDZM %*% defDZM, name = "expMeanDZM")
  expMeanDZOS <- mxAlgebra(expression = meanDZOS + betaDZOS %*% defDZOS, name = "expMeanDZOS")

  # Twin Correlations
  expCorMZF <- mxAlgebra(cov2cor(expCovMZF), name = "expCorMZF")
  expCorMZM <- mxAlgebra(cov2cor(expCovMZM), name = "expCorMZM")
  expCorDZF <- mxAlgebra(cov2cor(expCovDZF), name = "expCorDZF")
  expCorDZM <- mxAlgebra(cov2cor(expCovDZM), name = "expCorDZM")
  expCorDZOS <- mxAlgebra(cov2cor(expCovDZOS), name = "expCorDZOS")

  # Create Algebra for expected Variance/Covariance Matrices
  expCovMZF <- mxMatrix(type = "Symm", nrow = ntv, free = T, values = valDiag(svVa, ntv), labels = c("vMZF1", "cMZF21", "vMZF2"), name = "expCovMZF")
  expCovMZM <- mxMatrix(type = "Symm", nrow = ntv, free = T, values = valDiag(svVa, ntv), labels = c("vMZM1", "cMZM21", "vMZM2"), name = "expCovMZM")
  expCovDZF <- mxMatrix(type = "Symm", nrow = ntv, free = T, values = valDiag(svVa, ntv), labels = c("vDZF1", "cDZF21", "vDZF2"), name = "expCovDZF")
  expCovDZM <- mxMatrix(type = "Symm", nrow = ntv, free = T, values = valDiag(svVa, ntv), labels = c("vDZM1", "cDZM21", "vDZM2"), name = "expCovDZM")
  expCovDZOS <- mxMatrix(type = "Symm", nrow = ntv, free = T, values = valDiag(svVa, ntv), labels = c("vDZOS1", "cDZOS21", "vDZOS2"), name = "expCovDZOS")

  # Data objects for Multiple Groups
  dataMZF <- mxData(observed = mzfData, type = "raw")
  dataMZM <- mxData(observed = mzmData, type = "raw")
  dataDZF <- mxData(observed = dzfData, type = "raw")
  dataDZM <- mxData(observed = dzmData, type = "raw")
  dataDZOS <- mxData(observed = dzosData, type = "raw")

  # Objective objects for Multiple Groups
  funML <- mxFitFunctionML()
  objMZF <- mxExpectationNormal(covariance = "expCovMZF", means = "expMeanMZF", dimnames = selVars)
  objMZM <- mxExpectationNormal(covariance = "expCovMZM", means = "expMeanMZM", dimnames = selVars)
  objDZF <- mxExpectationNormal(covariance = "expCovDZF", means = "expMeanDZF", dimnames = selVars)
  objDZM <- mxExpectationNormal(covariance = "expCovDZM", means = "expMeanDZM", dimnames = selVars)
  objDZOS <- mxExpectationNormal(covariance = "expCovDZOS", means = "expMeanDZOS", dimnames = selVars)

  # Combine Groups
  modelMZF <- mxModel("MZF", meanMZF, betaMZF, defMZF, expMeanMZF, expCorMZF, expCovMZF, dataMZF, objMZF, funML)
  modelMZM <- mxModel("MZM", meanMZM, betaMZM, defMZM, expMeanMZM, expCorMZM, expCovMZM, dataMZM, objMZM, funML)
  modelDZF <- mxModel("DZF", meanDZF, betaDZF, defDZF, expMeanDZF, expCorDZF, expCovDZF, dataDZF, objDZF, funML)
  modelDZM <- mxModel("DZM", meanDZM, betaDZM, defDZM, expMeanDZM, expCorDZM, expCovDZM, dataDZM, objDZM, funML)
  modelDZOS <- mxModel("DZOS", meanDZOS, betaDZOS, defDZOS, expMeanDZOS, expCorDZOS, expCovDZOS, dataDZOS, objDZOS, funML)
  minus2ll <- mxAlgebra(MZF.objective + MZM.objective + DZF.objective + DZM.objective + DZOS.objective, name = "minus2sumloglikelihood")
  obj <- mxFitFunctionAlgebra("minus2sumloglikelihood")
  ciCor <- mxCI(c("MZF.expCorMZF[2,1]", "DZF.expCorDZF[2,1]"))
  twinSatModel <- mxModel("twinSat", modelMZF, modelMZM, modelDZF, modelDZM, modelDZOS, minus2ll, obj, ciCor)

  # Run Saturated Model
  twinSatFit <- mxTryHard(twinSatModel, intervals = T, extraTries = 100)
  twinSatSumm <- summary(twinSatFit)

  #### Sub-models ####
  # MODELS TO TEST HETEROGENIETY OF MEANS
  # Constrain expected Means to be equal across twin/sibling order
  MModel1 <- twinSatModel
  MModel1 <- omxSetParameters(MModel1, label = c("mMZF1", "mMZF2"), newlabels = "mMZF")
  MModel1 <- omxSetParameters(MModel1, label = c("mMZM1", "mMZM2"), newlabels = "mMZM")
  MModel1 <- omxSetParameters(MModel1, label = c("mDZF1", "mDZF2"), newlabels = "mDZF")
  MModel1 <- omxSetParameters(MModel1, label = c("mDZM1", "mDZM2"), newlabels = "mDZM")
  MModel1Fit <- mxTryHard(MModel1, intervals = F, extraTries = 50)

  # Constrain expected Means to be equal across same sex groups
  MModel2 <- MModel1
  MModel2 <- omxSetParameters(MModel2, label = c("mMZF", "mDZF"), newlabels = "mFemale")
  MModel2 <- omxSetParameters(MModel2, label = c("mMZM", "mDZM"), newlabels = "mMale")
  MModel2Fit <- mxTryHard(MModel2, intervals = F, extraTries = 50)

  # Constrain expected Means to be equal across opposite sex groups. DZOS indid 1 = female, indid 2 = male
  MModel3 <- MModel2
  MModel3 <- omxSetParameters(MModel3, label = c("mDZOS1", "mFemale"), newlabels = "mF")
  MModel3 <- omxSetParameters(MModel3, label = c("mDZOS2", "mMale"), newlabels = "mM")
  MModel3Fit <- mxTryHard(MModel3, intervals = F, extraTries = 50)

  # Constrain expected Means to be equal across sex
  MModel4 <- MModel3
  MModel4 <- omxSetParameters(MModel4, label = c("mF", "mM"), newlabels = "m")
  MModel4Fit <- mxTryHard(MModel4, intervals = F, extraTries = 50)

  #################################################################################################
  # MODELS TO TEST HETEROGENIETY OF VARIANCES
  # Constrain expected Variances to be equal across twin/sibling order
  VModel1 <- MModel4
  VModel1 <- omxSetParameters(VModel1, label = c("vMZF1", "vMZF2"), newlabels = "vMZF")
  VModel1 <- omxSetParameters(VModel1, label = c("vMZM1", "vMZM2"), newlabels = "vMZM")
  VModel1 <- omxSetParameters(VModel1, label = c("vDZF1", "vDZF2"), newlabels = "vDZF")
  VModel1 <- omxSetParameters(VModel1, label = c("vDZM1", "vDZM2"), newlabels = "vDZM")
  VModel1Fit <- mxTryHard(VModel1, intervals = F, extraTries = 50)

  # Constrain expected Variances to be equal across same sex groups
  VModel2 <- VModel1
  VModel2 <- omxSetParameters(VModel2, label = c("vMZF", "vDZF"), newlabels = "vFemale")
  VModel2 <- omxSetParameters(VModel2, label = c("vMZM", "vDZM"), newlabels = "vMale")
  VModel2Fit <- mxTryHard(VModel2, intervals = F, extraTries = 50)

  # Constrain expected Variances to be equal across opposite sex groups
  VModel3 <- VModel2
  VModel3 <- omxSetParameters(VModel3, label = c("vDZOS1", "vFemale"), newlabels = "vF")
  VModel3 <- omxSetParameters(VModel3, label = c("vDZOS2", "vMale"), newlabels = "vM")
  VModel3Fit <- mxTryHard(VModel3, intervals = F, extraTries = 50)

  # Constrain expected Variances to be equal across sex
  VModel4 <- VModel3
  VModel4 <- omxSetParameters(VModel4, label = c("vF", "vM"), newlabels = "v")
  VModel4Fit <- mxTryHard(VModel4, intervals = F, extraTries = 50)

  # MODELS TO TEST HETEROGENIETY OF COVARIANCES
  # Constrain expected covs to be equal within zygosity
  # H1c equated the covariances of MZ twins, and equated the covariance of same-sex DZ twins. Thus, the comparision between H1c and H0c tested for the presence of scalar sex limitation.
  CModel1 <- VModel4
  CModel1 <- omxSetParameters(CModel1, label = c("cMZF21", "cMZM21"), values = 0, newlabels = "cMZ")
  CModel1 <- omxSetParameters(CModel1, label = c("cDZF21", "cDZM21"), values = 0, newlabels = "cDZ21")
  CModel1Fit <- mxTryHard(CModel1, intervals = F, extraTries = 50)

  # Constrain expected covs to be equal across same and & opposite sex groups
  # H2c equated the covariances of MZ twins, and the covariances of all DZ twins. The comparision H2c therefore tested for the presence of non scalar sex limitation.
  CModel2 <- CModel1
  CModel2 <- omxSetParameters(CModel2, label = c("cDZOS21", "cDZ21"), values = 0, newlabels = "cDZ")
  CModel2Fit <- mxTryHard(CModel2, intervals = T, extraTries = 50)
  CModel2Summ <- summary(CModel2Fit)

  # Constrain expected covs to be equal between MZ and DZ groups
  # H3c equated all covariances The comparison tested whether variance on a trait was influenced by genetic factors.
  CModel3 <- CModel2
  CModel3 <- omxSetParameters(CModel3, label = c("cMZ", "cDZ"), values = 0, newlabels = "c")
  CModel3Fit <- mxTryHard(CModel3, intervals = F, extraTries = 50)

  # Drop cov to zero
  # H4c set all covariances to zero and tested for effects of familial aggregation
  CModel4 <- CModel3
  CModel4 <- omxSetParameters(CModel4, label = "c", free = FALSE, values = 0)
  CModel4Fit <- mxTryHard(CModel4, intervals = F, extraTries = 50)

  # Model fitting
  H1m <- mxCompare(twinSatFit, MModel1Fit)
  H2m <- mxCompare(MModel1Fit, MModel2Fit)
  H3m <- mxCompare(MModel2Fit, MModel3Fit)
  H4m <- mxCompare(MModel3Fit, MModel4Fit)

  H1v <- mxCompare(MModel4Fit, VModel1Fit)
  H2v <- mxCompare(VModel1Fit, VModel2Fit)
  H3v <- mxCompare(VModel2Fit, VModel3Fit)
  H4v <- mxCompare(VModel3Fit, VModel4Fit)

  H1c <- mxCompare(VModel4Fit, CModel1Fit)
  H2c <- mxCompare(CModel1Fit, CModel2Fit)
  H3c <- mxCompare(CModel2Fit, CModel3Fit)
  H4c <- mxCompare(CModel3Fit, CModel4Fit)


  #########################################################################################
  # Test Covariates
  # Drop sex_female
  noSex <- CModel2
  noSex <- omxSetParameters(model = noSex, labels = "betasex_female", free = FALSE, values = 0, name = "noSexModel")
  noSexFit <- mxTryHard(noSex, extraTries = 100, intervals = F)
  noSexSumm <- summary(noSexFit)

  # Drop Age
  noAge <- CModel2
  noAge <- omxSetParameters(model = noAge, labels = "betases01_age_months", free = FALSE, values = 0, name = "noAgeModel")
  noAgeFit <- mxTryHard(noAge, extraTries = 100, intervals = F)
  noAgeSumm <- summary(noAgeFit)

  # Model fitting
  noSexLRT <- mxCompare(CModel2Fit, noSexFit)
  noAgeLRT <- mxCompare(CModel2Fit, noAgeFit)
  
  Results <- data.frame(
    Variable = phenotype,
    ObservedStatistics = twinSatSumm$observedStatistics,
    Mean = mean.val,
    SD = sd.val,
    rMZF = VModel4Fit$MZF$expCorMZF$result[2, 1],
    rMZM = VModel4Fit$MZM$expCorMZM$result[2, 1],
    rDZF = VModel4Fit$DZF$expCorDZF$result[2, 1],
    rDZM = VModel4Fit$DZM$expCorDZM$result[2, 1],
    rDZOS = VModel4Fit$DZOS$expCorDZOS$result[2, 1],
    rMZ = CModel2Fit$MZF$expCorMZF$result[2, 1],
    rDZ = CModel2Fit$DZF$expCorDZF$result[2, 1],
    rMZ_95CI = paste0(sprintf("%.2f", round(CModel2Summ$CI$estimate[1], 2)), " (", sprintf("%.2f", round(CModel2Summ$CI$lbound[1], 2)), ", ", sprintf("%.2f", round(CModel2Summ$CI$ubound[1], 2)), ")"),
    rDZ_95CI = paste0(sprintf("%.2f", round(CModel2Summ$CI$estimate[2], 2)), " (", sprintf("%.2f", round(CModel2Summ$CI$lbound[2], 2)), ", ", sprintf("%.2f", round(CModel2Summ$CI$ubound[2], 2)), ")"),
    sexEstimate = CModel2Summ$parameters$Estimate[2],
    ageEstimate = CModel2Summ$parameters$Estimate[3],
    NoSex_pval = noSexLRT$p[2],
    NoAge_pval = noAgeLRT$p[2],
    AIC_Sat = twinSatSumm$AIC.Mx,
    CODE_Sat = twinSatFit$output$status$code,
    H1m_pval = H1m$p[2],
    H2m_pval = H2m$p[2],
    H3m_pval = H3m$p[2],
    H4m_pval = H4m$p[2],
    H1v_pval = H1v$p[2],
    H2v_pval = H2v$p[2],
    H3v_pval = H3v$p[2],
    H4v_pval = H4v$p[2],
    H1c_pval = H1c$p[2],
    H2c_pval = H2c$p[2],
    H3c_pval = H3c$p[2],
    H4c_pval = H4c$p[2],
    stringsAsFactors = FALSE
  )
  return(Results)
}

#### Run models ####
twin.data <- readRDS("QTAB_familywise.RDS")

# Run for single phenotype
Saturated_covariate(phenotype = "ProcSpeed_rawZ", twin.data = twin.data)

# Run for list of phenotypes
variable_list <- c("ProcSpeed_rawZ", "CJOLOZ", "TotalComposite_agecorr_stanZ")
results.sat <- lapply(variable_list, Saturated_covariate, twin.data = twin.data) %>% bind_rows()

# Assumption testing (sig after FDR adjustment)
results.sat$H1m_pval_fdr <- p.adjust(p = results.sat$H1m_pval, method = "fdr")
results.sat$H2m_pval_fdr <- p.adjust(p = results.sat$H2m_pval, method = "fdr")
results.sat$H3m_pval_fdr <- p.adjust(p = results.sat$H3m_pval, method = "fdr")
results.sat$H4m_pval_fdr <- p.adjust(p = results.sat$H4m_pval, method = "fdr")

results.sat$H1v_pval_fdr <- p.adjust(p = results.sat$H1v_pval, method = "fdr")
results.sat$H2v_pval_fdr <- p.adjust(p = results.sat$H2v_pval, method = "fdr")
results.sat$H3v_pval_fdr <- p.adjust(p = results.sat$H3v_pval, method = "fdr")
results.sat$H4v_pval_fdr <- p.adjust(p = results.sat$H4v_pval, method = "fdr")

results.sat$H1c_pval_fdr <- p.adjust(p = results.sat$H1c_pval, method = "fdr")
results.sat$H2c_pval_fdr <- p.adjust(p = results.sat$H2c_pval, method = "fdr")

select(results.sat %>% filter(H1m_pval_fdr < 0.05), c(Variable, H1m_pval_fdr, H1m_pval))
select(results.sat %>% filter(H2m_pval_fdr < 0.05), c(Variable, H2m_pval_fdr, H2m_pval))
select(results.sat %>% filter(H3m_pval_fdr < 0.05), c(Variable, H3m_pval_fdr, H3m_pval))
select(results.sat %>% filter(H4m_pval_fdr < 0.05), c(Variable, H4m_pval_fdr, H4m_pval))

select(results.sat %>% filter(H1v_pval_fdr < 0.05), c(Variable, H1v_pval_fdr, H1v_pval))
select(results.sat %>% filter(H2v_pval_fdr < 0.05), c(Variable, H2v_pval_fdr, H2v_pval))
select(results.sat %>% filter(H3v_pval_fdr < 0.05), c(Variable, H3v_pval_fdr, H3v_pval))
select(results.sat %>% filter(H4v_pval_fdr < 0.05), c(Variable, H4v_pval_fdr, H4v_pval))

select(results.sat %>% filter(H1c_pval_fdr < 0.05), c(Variable, H1c_pval_fdr, H1c_pval))
select(results.sat %>% filter(H2c_pval_fdr < 0.05), c(Variable, H2c_pval_fdr, H2c_pval))

# Covariate effects (sig after FDR adjustment)
results.sat$NoSex_pval_fdr <- p.adjust(p = results.sat$NoSex_pval, method = "fdr")
results.sat$NoAge_pval_fdr <- p.adjust(p = results.sat$NoAge_pval, method = "fdr")

select(results.sat %>% filter(NoSex_pval_fdr < 0.05), c(Variable, NoSex_pval_fdr, sexEstimate))
select(results.sat %>% filter(NoAge_pval_fdr < 0.05), c(Variable, NoAge_pval_fdr, ageEstimate))

results.sat$Mean_SD <- paste0(sprintf("%.2f", round(results.sat$Mean, 2)), " (", sprintf("%.2f", round(results.sat$SD, 2)), ")")
write.csv(results.sat, "03_sat_covariate_OpenMX_results.csv", row.names = F)
