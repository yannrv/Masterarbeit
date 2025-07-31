#### GAM ANALYSIS OF ALPHA DIVERSITY (16S & 18S) VS ENVIRONMENTAL PARAMETERS ####


# Load packages
library(dplyr)
library(tidyverse)
library(mgcv)
library(MASS)
library(caret)
library(ggplot2)

# ===========
# 16S Dataset
# ===========

# Load data 16S
sample.data_16S <- read.csv("~/Metadata_16S.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
sample.data_16S$Evenness_16S <- sample.data_16S$Shannon_16S / log(sample.data_16S$Richness_16S)
filtered_data_16S <- sample.data_16S %>% dplyr::select(Sea_Surface_Temperature, Salinity, Orbital_Velocity_Mean, Bathymetry, Sheltered_Index, Richness_16S, Evenness_16S)

# Define variable
env_vars <- c("Sea_Surface_Temperature", "Salinity", "Orbital_Velocity_Mean", "Bathymetry", "Sheltered_Index")
response_richness <- "Richness_16S"
response_evenness <- "Evenness_16S"

# Correlation test 
predictor_cor_matrix <- cor(filtered_data_16S[, env_vars], method = "spearman")
high_cor_env <- findCorrelation(predictor_cor_matrix, cutoff = 0.50)
cleaned_env_vars <- if(length(high_cor_env) > 0) env_vars[-high_cor_env] else env_vars

# Standardization 
pp_env <- preProcess(filtered_data_16S[, cleaned_env_vars], method = c("center", "scale"))
scaled_env <- predict(pp_env, filtered_data_16S[, cleaned_env_vars])
final_data_env <- cbind(scaled_env, filtered_data_16S[, c(response_richness, response_evenness)])

# Optimal K
calculate_optimal_k <- function(n, min_k = 4, max_k = 18, ratio = 0.7) {
  k <- round(n * ratio)
  k <- pmin(pmax(k, min_k), max_k)
  k <- pmin(k, n - 1)
  return(k)
}

# Richness 16S
unique_counts <- sapply(cleaned_env_vars, function(v) length(unique(final_data_env[[v]])))
unique_counts
k_vals <- sapply(unique_counts, calculate_optimal_k)
gam_terms <- mapply(function(v, k) paste0("s(", v, ", k=", k, ")"), cleaned_env_vars, k_vals)

# Define family
formula_richness <- as.formula(paste(response_richness, "~", paste(gam_terms, collapse = " + ")))
overdispersion <- var(final_data_env[[response_richness]]) / mean(final_data_env[[response_richness]])
family_used <- if (overdispersion > 2) quasipoisson("log") else poisson("log")

# Full model Richness 16S
gam_richness <- gam(formula_richness, data = final_data_env, family = family_used, method = "REML")
summary(gam_richness)

# Parsimonious model Richness 16S
formula_env_pars <- as.formula("Richness_16S ~ s(Sea_Surface_Temperature, k=12) + s(Salinity, k=12) + s(Orbital_Velocity_Mean, k=14) + s(Bathymetry, k=12)")
gam_env_pars <- gam(formula_env_pars, data = final_data_env, family = quasipoisson("log"), method = "REML")
print(summary(gam_env_pars))

# Save GAM diagnostics for Richness 16S
pdf("gam_diagnostics_Env_richness_16S.pdf", width = 12, height = 10)
par(mfrow = c(2, 2))  
gam.check(gam_env_pars, rep = 0)
par(mfrow = c(1, 1))  
dev.off()

# Full model for evenness 16S
formula_evenness <- as.formula("Evenness_16S ~ s(Sea_Surface_Temperature, k=18) + s(Salinity, k=18) + s(Orbital_Velocity_Mean, k=18) + s(Bathymetry, k=17) + s(Sheltered_Index, k = 4)")
gam_evenness <- gam(formula_evenness, data = final_data_env, family = betar(link = "logit"), method = "REML")
print(summary(gam_evenness))

# Parsimonious model for evenness 16S
formula_evenness_pars <- as.formula("Evenness_16S ~ s(Salinity, k=12)")
gam_evenness_pars <- gam(formula_evenness_pars, data = final_data_env, family = betar(link = "logit"), method = "REML")
print(summary(gam_evenness_pars))

# Save GAM diagnostics for Evenness 16S
pdf("gam_diagnostics_Env_evenness_16S.pdf", width = 12, height = 10)
par(mfrow = c(2, 2))  
gam.check(gam_env_pars, rep = 0)
par(mfrow = c(1, 1))  
dev.off()

# Remove List
rm(list = ls())

# ===========
# 18S Dataset
# ===========

# Load data 18S
sample.data_18S <- read.csv("~/Metadata_18S.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
sample.data_18S$Evenness_18S <- sample.data_18S$Shannon_18S / log(sample.data_18S$Richness_18S)
filtered_data_18S <- sample.data_18S %>% 
  dplyr::select(Sea_Surface_Temperature, Salinity, Orbital_Velocity_Mean, Bathymetry, Sheltered_Index, Richness_18S, Evenness_18S) %>% 
  na.omit()

# Define variable
env_vars <- c("Sea_Surface_Temperature", "Salinity", "Orbital_Velocity_Mean", "Bathymetry", "Sheltered_Index")
response_richness <- "Richness_18S"
response_evenness <- "Evenness_18S"

# Correlation test
predictor_cor_matrix <- cor(filtered_data_18S[, env_vars], method = "spearman")
high_cor_env <- findCorrelation(predictor_cor_matrix, cutoff = 0.50)
cleaned_env_vars <- if(length(high_cor_env) > 0) env_vars[-high_cor_env] else env_vars
cleaned_env_vars

# Standardization
pp_env <- preProcess(filtered_data_18S[, cleaned_env_vars], method = c("center", "scale"))
scaled_env <- predict(pp_env, filtered_data_18S[, cleaned_env_vars])
final_data_env <- cbind(scaled_env, filtered_data_18S[, c(response_richness, response_evenness)])

# Optimal K
calculate_optimal_k <- function(n, min_k = 4, max_k = 18, ratio = 0.7) {
  k <- round(n * ratio)
  k <- pmin(pmax(k, min_k), max_k)
  k <- pmin(k, n - 1)
  return(k)
}

# Richness 18S
unique_counts <- sapply(cleaned_env_vars, function(v) length(unique(final_data_env[[v]])))
k_vals <- sapply(unique_counts, calculate_optimal_k)
gam_terms <- mapply(function(v, k) paste0("s(", v, ", k=", k, ")"), cleaned_env_vars, k_vals)

# define family 
formula_richness <- as.formula(paste(response_richness, "~", paste(gam_terms, collapse = " + ")))
overdispersion <- var(final_data_env[[response_richness]]) / mean(final_data_env[[response_richness]])
family_used <- if (overdispersion > 2) quasipoisson("log") else poisson("log")

# Full model Richness 18S
gam_richness <- gam(formula_richness, data = final_data_env, family = family_used, method = "REML")
print(summary(gam_richness))

# Parsimonious model Richness 18S
formula_env_pars <- as.formula("Richness_18S ~ s(Salinity, k=12) + s(Bathymetry, k = 12) + s(Sheltered_Index, k = 4)")
gam_env_pars <- gam(formula_env_pars, data = final_data_env, family = quasipoisson("log"), method = "REML")
print(summary(gam_env_pars))

# Save GAM diagnostics for Richness 18S
pdf("gam_diagnostics_Env_richness_18S.pdf", width = 12, height = 10)
par(mfrow = c(2, 2))  
gam.check(gam_env_pars, rep = 0)
par(mfrow = c(1, 1))  
dev.off()

# Full model for evenness 18S
formula_evenness <- as.formula("Evenness_18S ~ s(Salinity, k=18) + s(Orbital_Velocity_Mean, k=18) + s(Bathymetry, k=17) + s(Sheltered_Index, k = 4)")
gam_evenness <- gam(formula_evenness, data = final_data_env, family = betar(link = "logit"), method = "REML")
print(summary(gam_evenness))

# Parsimonious model for evenness 18S
formula_evenness_pars <- as.formula("Evenness_18S ~ s(Orbital_Velocity_Mean, k=14)+ s(Bathymetry, k=14) + s(Sheltered_Index, k = 4)")
gam_evenness_pars <- gam(formula_evenness_pars, data = final_data_env, family = betar(link = "logit"), method = "REML")
print(summary(gam_evenness_pars))

# Save GAM diagnostics for Evenness 18S
pdf("gam_diagnostics_Env_evenness_18S.pdf", width = 12, height = 10)
par(mfrow = c(2, 2))  
gam.check(gam_evenness_pars, rep = 0)
par(mfrow = c(1, 1))  
dev.off()
