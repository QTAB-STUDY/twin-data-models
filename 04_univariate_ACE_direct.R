# Lachlan T Strike
# Based on scripts provided by Hermine Maes (https://hermine-maes.squarespace.com)
# Direct estimates of variance components

#### Load libraries ####
library(OpenMx)
library(dplyr)
source("miFunctions.R")
#mxOption(NULL, "Default optimizer", "CSOLNP") # NPSOL SLSQP CSOLNP

Univariate_ACE <- function(phenotype, twin.data, covariate) {
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
  # Create Matrices for Covariates and linear Regression Coefficients
  defL <- mxMatrix(type = "Full", nrow = nc, ncol = ntv, free = F, labels = c(paste0("data.", covariate, "_01"), paste0("data.", covariate, "_02")), name = "defL")
  pathBl <- mxMatrix(type = "Full", nrow = 1, ncol = nc, free = TRUE, values = 0.1, labels = c(paste0("beta", covariate)), name = "bl")
  
  # Create Algebra for expected Mean Matrices
  meanG <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, values = svMe, label = "mean", name = "meanG")
  expMean <- mxAlgebra(expression = meanG + bl %*% defL, name = "expMeanG")
  
  # Create Matrices for Variance Components
  covA <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VA11", name="VA" )
  covC <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VC11", name="VC" )
  covE <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VE11", name="VE" )
  
  # Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
  covP <- mxAlgebra( expression= VA+VC+VE, name="V" )
  covMZ <- mxAlgebra( expression= VA+VC, name="cMZ" )
  covDZ <- mxAlgebra( expression= 0.5%x%VA+VC, name="cDZ" )
  expCovMZ <- mxAlgebra( expression= rbind( cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ" )
  expCovDZ <- mxAlgebra( expression= rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )
  
  # Create Data Objects for Multiple Groups
  dataMZ <- mxData( observed=mzData, type="raw" )
  dataDZ <- mxData( observed=dzData, type="raw" )
  
  # Create Expectation Objects for Multiple Groups
  expMZ <- mxExpectationNormal( covariance="expCovMZ", means="expMeanG", dimnames=selVars )
  expDZ <- mxExpectationNormal( covariance="expCovDZ", means="expMeanG", dimnames=selVars )
  funML <- mxFitFunctionML()
  
  # Create Model Objects for Multiple Groups
  pars <- list( pathBl, meanG, covA, covC, covE, covP )
  defs <- list( defL )
  modelMZ <- mxModel( pars, defs, expMean, covMZ, expCovMZ, dataMZ, expMZ, funML, name="MZ" )
  modelDZ <- mxModel( pars, defs, expMean, covDZ, expCovDZ, dataDZ, expDZ, funML, name="DZ" )
  multi <- mxFitFunctionMultigroup( c("MZ","DZ") )
  
  # Create Algebra for Variance Components
  rowUS <- rep('US',nv)
  colUS <- rep(c('VA','VC','VE','SA','SC','SE'),each=nv)
  estUS <- mxAlgebra( expression=cbind(VA,VC,VE,VA/V,VC/V,VE/V), name="US", dimnames=list(rowUS,colUS) )
  
  # Create Confidence Interval Objects
  ciACE <- mxCI( "US[1,1:6]" )
  
  # Build Model with Confidence Intervals
  modelACE <- mxModel( "oneACEvca", pars, modelMZ, modelDZ, multi, estUS, ciACE )
  
  # ------------------------------------------------------------------------------
  fitACE <- mxTryHard(modelACE, intervals = T, extraTries = 50)
  sumACE <- summary(fitACE)
  
  # Run AE model
  modelAE <- mxModel( modelACE, name="oneAEvca" )
  modelAE <- omxSetParameters( modelAE, labels="VC11", free=FALSE, values=0 )
  fitAE <- mxTryHard( modelAE, intervals=T, extraTries = 50)
  sumAE <- summary(fitAE)
  
  # Run CE model
  modelCE <- mxModel( modelACE, name="oneCEvca" )
  modelCE <- omxSetParameters( modelCE, labels="VA11", free=FALSE, values=0 )
  fitCE <- mxTryHard( modelCE, intervals=T, extraTries = 50)
  sumCE <- summary(fitCE)
  
  # Run E model
  modelE <- mxModel( modelAE, name="oneEvca" )
  modelE <- omxSetParameters( modelE, labels="VA11", free=FALSE, values=0)
  fitE <- mxTryHard( modelE, intervals=T, extraTries = 50)
  sumE <- summary(fitE)
  
  # Model Comparison
  AW <- omxAkaikeWeights(models = list(fitACE, fitAE, fitCE, fitE))
  model.comparison <- mxCompare( fitACE, nested <- list(fitAE, fitCE, fitE))
  
  #### Create results table ####
  NameCols <- c("Variable", "ObservedStatistics", "ACE_A_95CI", "ACE_C_95CI", "ACE_E_95CI", "AE_A_95CI", "AE_E_95CI", "CE_C_95CI", "CE_E_95CI", "AIC_BestFit",
                "dropC", "dropA", "dropAC", "ACE_A", "ACE_C", "ACE_E", "AE_A", "AE_E", "CE_C", "CE_E")
  
  Results <- data.frame(matrix(nrow = 1, ncol = length(NameCols)))
  colnames(Results) = NameCols
  Results[, "Variable"] <- phenotype
  Results[, "ObservedStatistics"] <- sumACE$observedStatistics
  Results[, "ACE_A_95CI"] <- paste0(sprintf("%.2f", round(sumACE$CI$estimate[4], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[4], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[4], 2)), ")")
  Results[, "ACE_C_95CI"] <- paste0(sprintf("%.2f", round(sumACE$CI$estimate[5], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[5], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[5], 2)), ")")
  Results[, "ACE_E_95CI"] <- paste0(sprintf("%.2f", round(sumACE$CI$estimate[6], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[6], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[6], 2)), ")")
  Results[, "AE_A_95CI"] <- paste0(sprintf("%.2f", round(sumAE$CI$estimate[4], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[4], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[4], 2)), ")")
  Results[, "AE_E_95CI"] <- paste0(sprintf("%.2f", round(sumAE$CI$estimate[6], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[6], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[6], 2)), ")")
  Results[, "CE_C_95CI"] <- paste0(sprintf("%.2f", round(sumCE$CI$estimate[5], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[5], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[5], 2)), ")")
  Results[, "CE_E_95CI"] <- paste0(sprintf("%.2f", round(sumCE$CI$estimate[6], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[6], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[6], 2)), ")")
  Results[, "ACE_AIC"] <- sumACE$AIC.Mx
  Results[, "CE_AIC"] <- sumCE$AIC.Mx
  Results[, "AE_AIC"] <- sumAE$AIC.Mx
  Results[, "E_AIC"] <- sumE$AIC.Mx
  Results[, "AIC_BestFit"] <- AW$model[1]
  Results[, "dropC"] <- model.comparison$p[2]
  Results[, "dropA"] <- model.comparison$p[3]
  Results[, "dropAC"] <- model.comparison$p[4]
  Results[, "ACE_A"] <- fitACE$US$result[4]
  Results[, "ACE_C"] <- fitACE$US$result[5]
  Results[, "ACE_E"] <- fitACE$US$result[6]
  Results[, "AE_A"] <- fitAE$US$result[4]
  Results[, "AE_E"] <- fitAE$US$result[6]
  Results[, "CE_C"] <- fitCE$US$result[5]
  Results[, "CE_E"] <- fitCE$US$result[6]
  return(Results)
}

#### Run ####
setwd("C:/GitHub/twin-data-models")
my.data <- readRDS("QTAB_familywise_ses01.RDS")

# Single phenotype
result.single <- as_tibble(Univariate_ACE(phenotype = "ProcSpeed_rawZ", twin.data = my.data, covariate = c("sex_female", "ses01_age_months")))
result.single %>% select(Variable, AIC_BestFit, ACE_A_95CI, ACE_C_95CI, ACE_E_95CI, AE_A_95CI, AE_E_95CI, ACE_A, ACE_C, ACE_E, AE_A, AE_E, CE_C, CE_E)

# Multiple phenotypes
phenotype.array <-  c("ProcSpeed_rawZ", "CJOLOZ", "FAS_60Z", "DSfZ", "SCAS_scoreZ",
                      "SMFQ_scoreZ", "SPHERE_anxdep_scoreZ", "SPHERE_fat_scoreZ",
                      "CRSQ_rum_scoreZ", "CASQ_comp_scoreZ", "MSPSS_scoreZ", "DLSS_scoreZ",
                      "GBSvictim_scoreZ", "PDSS_scoreZ", "pPRMQ_scoreZ", "pSCAS_scoreZ",
                      "pSMFQ_scoreZ", "pSDQ_total_scoreZ", "pEATQ_act_scoreZ", "kJwithDFZ")
result.multiple <- as_tibble(lapply(phenotype.array, Univariate_ACE, twin.data = my.data, covariate = c("sex_female", "ses01_age_months")) %>% bind_rows())
result.multiple %>% select(Variable, AIC_BestFit, ACE_A_95CI, ACE_C_95CI, ACE_E_95CI, AE_A_95CI, AE_E_95CI, ACE_A, ACE_C, ACE_E, AE_A, AE_E, CE_C, CE_E)

write.csv(result.multiple, "04_univariate_direct_variance_ACE_estimates.csv", row.names = F)
