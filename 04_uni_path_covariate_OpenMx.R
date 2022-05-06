# Lachlan T Strike
# Based on scripts provided by Hermine Maes (https://hermine-maes.squarespace.com)
# Univariate model specifying path coefficients and covariates

rm(list = ls())
setwd("~/GitHub/twin-data-models")

#### Load libraries ####
library(OpenMx)
library(dplyr)
source("miFunctions.R")
mxOption(NULL, "Default optimizer", "SLSQP") # NPSOL SLSQP CSOLNP

#### Models ####
Univariate_path_covariate <- function(phenotype, twin.data) {
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
  mzData <- subset(twin.data, zyg < 3, useVars)
  dzData <- subset(twin.data, zyg > 2, useVars)

  # Set Starting Values
  nt <- 2
  nv <- length(phenotype)
  ntv <- nv * nt
  svPa <- sqrt(var(c(twin.data[, paste0(phenotype, "_01")], twin.data[, paste0(phenotype, "_02")]), na.rm = T) / 3)
  svMe <- mean(unlist(twin.data[, selVars]), na.rm = T)

  # ------------------------------------------------------------------------------
  # PREPARE GENETIC MODEL
  # ------------------------------------------------------------------------------
  # Cholesky Decomposition ACE Model
  # ------------------------------------------------------------------------------
  pathA <- mxMatrix(
    type = "Lower", nrow = nv, ncol = nv, free = TRUE,
    values = svPa, labels = labLower("a", nv), name = "a"
  )
  pathC <- mxMatrix(
    type = "Lower", nrow = nv, ncol = nv, free = TRUE,
    values = svPa, labels = labLower("c", nv), name = "c"
  )
  pathE <- mxMatrix(
    type = "Lower", nrow = nv, ncol = nv, free = TRUE,
    values = svPa, labels = labLower("e", nv), name = "e"
  )

  StPathA <- mxAlgebra(iSD %*% a, name = "sta")
  StPathC <- mxAlgebra(iSD %*% c, name = "stc")
  StPathE <- mxAlgebra(iSD %*% e, name = "ste")

  # Matrices generated to hold A, C, and E computed variance components
  covA <- mxAlgebra(expression = a %*% t(a), name = "A")
  covC <- mxAlgebra(expression = c %*% t(c), name = "C")
  covE <- mxAlgebra(expression = e %*% t(e), name = "E")

  # Algebra to compute total variances and standard deviations
  # (diagonal only)
  covP <- mxAlgebra(expression = A + C + E, name = "V") # Total within-twin covariance for phenotypes
  matI <- mxMatrix(type = "Iden", nrow = nv, ncol = nv, name = "I")
  invSD <- mxAlgebra(expression = solve(sqrt(I * V)), name = "iSD")

  stVA <- mxAlgebra(A / V, name = "StVA")
  stVC <- mxAlgebra(C / V, name = "StVC")
  stVE <- mxAlgebra(E / V, name = "StVE")

  # Covariate matrix for expected Mean
  meanG <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, values = svMe, label = "mean", name = "InterceptMean")
  defVars <- mxMatrix(type = "Full", nrow = nc, ncol = ntv, free = F, labels = c(paste0("data.", covariate, "_01"), paste0("data.", covariate, "_02")), name = "DefVars")
  beta <- mxMatrix(type = "Full", nrow = 1, ncol = nc, free = TRUE, values = 0.1, labels = c(paste0("beta", covariate)), name = "Beta")
  expMean <- mxAlgebra(expression = InterceptMean + Beta %*% DefVars, name = "expMean")

  # Algebra for expected Mean and Variance/Covariance Matrices
  # in MZ & DZ twins
  covMZ <- mxAlgebra(expression = rbind(
    cbind(V, A + C),
    cbind(A + C, V)
  ), name = "expCovMZ")
  covDZ <- mxAlgebra(expression = rbind(
    cbind(V, 0.5 %x% A + C),
    cbind(0.5 %x% A + C, V)
  ), name = "expCovDZ")

  # Data objects for Multiple Groups
  dataMZ <- mxData(observed = mzData, type = "raw")
  dataDZ <- mxData(observed = dzData, type = "raw")

  # Objective objects for Multiple Groups
  objMZ <- mxExpectationNormal(covariance = "expCovMZ", means = "expMean", dimnames = selVars)
  objDZ <- mxExpectationNormal(covariance = "expCovDZ", means = "expMean", dimnames = selVars)
  fitFunction <- mxFitFunctionML()

  # Confidence intervals for standardised variance components
  CIs <- mxCI(c("StVA", "StVC", "StVE"))

  # Combine Groups
  pars <- list(
    pathA, pathC, pathE, covA, covC, covE, covP, matI,
    invSD, StPathA, StPathC, StPathE, stVA, stVC, stVE
  )
  modelMZ <- mxModel(pars, beta, defVars, meanG, expMean, covMZ, dataMZ, fitFunction, objMZ, name = "MZ")
  modelDZ <- mxModel(pars, beta, defVars, meanG, expMean, covDZ, dataDZ, fitFunction, objDZ, name = "DZ")
  minus2ll <- mxAlgebra(expression = MZ.objective + DZ.objective, name = "m2LL")
  obj <- mxFitFunctionAlgebra("m2LL")
  CholAceModel <- mxModel("CholACE", pars, modelMZ, modelDZ, minus2ll, obj, CIs)

  # ------------------------------------------------------------------------------
  CholAceFit <- mxTryHard(CholAceModel, intervals = T, extraTries = 50)
  CholAceSum <- summary(CholAceFit)

  # Run AE model
  CholAeModel <- CholAceModel
  CholAeModel <- mxModel(CholAceModel, name = "CholAE")
  CholAeModel <- omxSetParameters(CholAeModel,
    labels = c("c11"),
    free = FALSE, values = 0
  )
  CholAeFit <- mxTryHard(CholAeModel, intervals = T, extraTries = 50)
  CholAeSum <- summary(CholAeFit)

  # Run CE model
  CholCeModel <- CholAceModel
  CholCeModel <- mxModel(CholCeModel, name = "CholCE")
  CholCeModel <- omxSetParameters(CholCeModel,
    labels = c("a11"),
    free = FALSE, values = 0
  )
  CholCeFit <- mxTryHard(CholCeModel, intervals = T, extraTries = 50)
  CholCeSumm <- summary(CholCeFit)

  # Run E model
  CholEModel <- CholAeModel
  CholEModel <- mxModel(CholEModel, name = "CholE")
  CholEModel <- omxSetParameters(CholEModel,
    labels = c("a11"),
    free = FALSE, values = 0
  )
  CholEFit <- mxTryHard(CholEModel, intervals = T, extraTries = 50)
  CholESumm <- summary(CholEFit)

  ACENested <- list(CholAeFit, CholCeFit, CholEFit)
  Results_modelfit <- mxCompare(CholAceFit, ACENested)

  AW <- omxAkaikeWeights(models = list(CholAceFit, CholAeFit, CholCeFit, CholEFit))

  #### Create results table ####
  Results <- NULL
  Results$Variable <- phenotype
  Results$ObservedStatistics <- CholAceSum$observedStatistics
  Results$ACE_AIC <- CholAceSum$AIC.Mx
  Results$CE_AIC <- CholCeSumm$AIC.Mx
  Results$AE_AIC <- CholAeSum$AIC.Mx
  Results$E_AIC <- CholESumm$AIC.Mx

  Results$ACE_noC_pval <- Results_modelfit$p[2]
  Results$ACE_noA_pval <- Results_modelfit$p[3]
  Results$ACE_noAC_pval <- Results_modelfit$p[4]

  Results$ACE_model_conf_set <- AW[which(AW$model == "CholACE"), "inConfidenceSet"]
  Results$AE_model_conf_set <- AW[which(AW$model == "CholAE"), "inConfidenceSet"]
  Results$CE_model_conf_set <- AW[which(AW$model == "CholCE"), "inConfidenceSet"]
  Results$E_model_conf_set <- AW[which(AW$model == "CholE"), "inConfidenceSet"]

  Results$ACE_A_95CI <- paste0(sprintf("%.2f", round(CholAceSum$CI$estimate[1], 2)), " (", sprintf("%.2f", round(CholAceSum$CI$lbound[1], 2)), ", ", sprintf("%.2f", round(CholAceSum$CI$ubound[1], 2)), ")")
  Results$ACE_C_95CI <- paste0(sprintf("%.2f", round(CholAceSum$CI$estimate[2], 2)), " (", sprintf("%.2f", round(CholAceSum$CI$lbound[2], 2)), ", ", sprintf("%.2f", round(CholAceSum$CI$ubound[2], 2)), ")")
  Results$ACE_E_95CI <- paste0(sprintf("%.2f", round(CholAceSum$CI$estimate[3], 2)), " (", sprintf("%.2f", round(CholAceSum$CI$lbound[3], 2)), ", ", sprintf("%.2f", round(CholAceSum$CI$ubound[3], 2)), ")")
  Results$AE_A_95CI <- paste0(sprintf("%.2f", round(CholAeSum$CI$estimate[1], 2)), " (", sprintf("%.2f", round(CholAeSum$CI$lbound[1], 2)), ", ", sprintf("%.2f", round(CholAeSum$CI$ubound[1], 2)), ")")
  Results$AE_E_95CI <- paste0(sprintf("%.2f", round(CholAeSum$CI$estimate[3], 2)), " (", sprintf("%.2f", round(CholAeSum$CI$lbound[3], 2)), ", ", sprintf("%.2f", round(CholAeSum$CI$ubound[3], 2)), ")")
  Results$CE_C_95CI <- paste0(sprintf("%.2f", round(CholCeSumm$CI$estimate[2], 2)), " (", sprintf("%.2f", round(CholCeSumm$CI$lbound[2], 2)), ", ", sprintf("%.2f", round(CholCeSumm$CI$ubound[2], 2)), ")")
  Results$CE_E_95CI <- paste0(sprintf("%.2f", round(CholCeSumm$CI$estimate[3], 2)), " (", sprintf("%.2f", round(CholCeSumm$CI$lbound[3], 2)), ", ", sprintf("%.2f", round(CholCeSumm$CI$ubound[3], 2)), ")")
  Results$E_E_95CI <- paste0(sprintf("%.2f", round(CholESumm$CI$estimate[3], 2)), " (", sprintf("%.2f", round(CholESumm$CI$lbound[3], 2)), ", ", sprintf("%.2f", round(CholESumm$CI$ubound[3], 2)), ")")

  Results$CODE_ACE <- CholAceFit$output$status$code
  Results$CODE_AE <- CholAeFit$output$status$code
  Results$CODE_CE <- CholCeFit$output$status$code
  Results$CODE_E <- CholEFit$output$status$code

  Results$ACE_A <- CholAceFit$StVA@result[1, 1]
  Results$ACE_C <- CholAceFit$StVC@result[1, 1]
  Results$ACE_E <- CholAceFit$StVE@result[1, 1]

  Results$AE_A <- CholAeFit$StVA@result[1, 1]
  Results$AE_E <- CholAeFit$StVE@result[1, 1]

  Results$CE_C <- CholCeFit$StVC@result[1, 1]
  Results$CE_E <- CholCeFit$StVE@result[1, 1]

  Results$E_E <- CholEFit$StVE@result[1, 1]
  return(Results)
}

#### Run models ####
mydata <- readRDS("QTAB_familywise.RDS")

# Run for single phenotype
Univariate_path_covariate(phenotype = "ProcSpeed_rawZ", twin.data = mydata)

# Run for list of phenotypes
variable_list <- c("ProcSpeed_rawZ", "CJOLOZ", "TotalComposite_agecorr_stanZ")
results <- lapply(variable_list, Univariate_path_covariate, twin.data = mydata) %>% bind_rows()
results <- as.data.frame(results)

write.csv(results, "04_uni_path_covariate_OpenMx_results.csv", row.names = F)
