setwd("~/GitHub/twin-data-models")
rm(list = ls())

library(tidyverse)

#### Results Table ####
source("04a_uni_path_covariate_OpenMx.R")
qtab.ses01 <- readRDS("QTAB_familywise_ses01.RDS")

# Single phenotype
result.single <- as_tibble(Univariate_path_covariate(phenotype = "ProcSpeed_rawZ", twin.data = qtab.ses01, covariate = c("ses01_age_months", "sex_female")))
result.single %>% select(Variable, ACE_A_95CI, ACE_C_95CI, ACE_E_95CI, ACE_AIC, AE_AIC, CE_AIC, E_AIC)

# Phenotype array
phenotype.array <- c("ProcSpeed_rawZ", "CJOLOZ", "FAS_60Z", "DSfZ", "SCAS_scoreZ",
                   "SMFQ_scoreZ", "SPHERE_anxdep_scoreZ", "SPHERE_fat_scoreZ",
                   "CRSQ_rum_scoreZ", "CASQ_comp_scoreZ", "MSPSS_scoreZ", "DLSS_scoreZ",
                   "GBSvictim_scoreZ", "PDSS_scoreZ", "pPRMQ_scoreZ", "pSCAS_scoreZ",
                   "pSMFQ_scoreZ", "pSDQ_total_scoreZ", "pEATQ_act_scoreZ", "kJwithDFZ")
result.multiple <- lapply(phenotype.array, Univariate_path_covariate, twin.data = qtab.ses01, covariate = c("ses01_age_months", "sex_female")) %>% bind_rows()
result.multiple <- as_tibble(result.multiple)
write.csv(result.multiple, "qtab_ses01_ace_estimates.csv", row.names = F)

result.multiple %>% select(Variable, ACE_A_95CI, ACE_C_95CI, ACE_E_95CI, ACE_AIC, AE_AIC, CE_AIC, E_AIC)

#### Results Plot ####
library(viridis)
library(ggplot2)
my_theme <- function() {
  theme_minimal(base_size = 8) +
    theme(
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      text=element_text(family="sans"),
      legend.position="bottom"
    )
}

results <- read.csv("qtab_ses01_ace_estimates.csv")
nvars <- length(results$Variable)
results_table <- as.data.frame(c(results$ACE_A, results$ACE_C, results$ACE_E))
colnames(results_table) <- "VC"
results_table$variance <- c(rep("A", nvars), rep("C", nvars), rep("E", nvars))
results_table$Variable <- factor(results$Variable) # specify levels here for ordering
results_table$variance <- factor(results_table$variance, levels = c("E", "C", "A"))

gACE <- ggplot() +
  geom_bar(aes(y = VC, x = Variable, fill=variance), data = results_table, stat="identity") + 
  scale_fill_manual(values=viridis(n = 3, option = "D", direction = 1), name = "", breaks=c("A","C","E")) +
  labs(x = "", y = "Variance Explained")  +
  ggtitle("ses-01 Variance Components") +
  coord_flip() + 
  my_theme()

pdf(file = "ses01_ace.pdf", paper = "a4")
gACE
dev.off()

