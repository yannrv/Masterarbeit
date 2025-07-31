#### GAM ANALYSIS OF ALPHA DIVERSITY (16S & 18S) VS MORPHOLOGICAL TRAITS ####


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
sample.data_16S <- read.csv("Metadata_16S.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
sample.data_16S$Evenness_16S <- sample.data_16S$Shannon_16S / log(sample.data_16S$Richness_16S)
filtered_data_16S <- sample.data_16S %>% dplyr::select(Shoot_Length, Leaf_Width, No_Leaves, Richness_16S, Evenness_16S) %>% na.omit()

# Define variable
mor_var <- c("Shoot_Length", "Leaf_Width", "No_Leaves")
response_richness <- "Richness_16S"
response_evenness <- "Evenness_16S"

# Standardization 
pp <- preProcess(filtered_data_16S[, mor_var], method = c("center", "scale"))
scaled_morpho <- predict(pp, filtered_data_16S[, mor_var])
final_data_16S <- cbind(scaled_morpho, filtered_data_16S[, c(response_richness, response_evenness)])

# Optimal k
calculate_optimal_k <- function(n_unique, min_k = 4, max_k = 15, ratio = 0.6) {
  optimal_k <- round(n_unique * ratio)
  optimal_k <- pmax(optimal_k, min_k)
  optimal_k <- pmin(optimal_k, max_k)
  optimal_k <- pmin(optimal_k, n_unique - 1)
  return(optimal_k)
}

# Richness 16S
unique_counts <- sapply(mor_var, function(v) length(unique(final_data_16S[[v]])))
k_vals <- sapply(unique_counts, calculate_optimal_k)
gam_terms <- mapply(function(var, k) paste0("s(", var, ", k=", k, ")"), names(k_vals), k_vals, USE.NAMES = FALSE)
formula_full <- as.formula(paste(response_richness, "~", paste(gam_terms, collapse = " + ")))

# Define family
overdispersion <- var(final_data_16S[[response_richness]]) / mean(final_data_16S[[response_richness]])
family_used <- if (overdispersion > 2) quasipoisson("log") else poisson("log")
cat("\nSelected family:", if (overdispersion > 2) "quasipoisson" else "poisson", "\n")

# Full Richness 16S
gam_full <- gam(formula_full, data = final_data_16S, family = family_used, method = "REML")
print(summary(gam_full))

# Parsimonious model Richness 16S
formula_pars <- as.formula("Richness_16S ~ s(Shoot_Length, k=15) + s(Leaf_Width, k=9)")
gam_pars <- gam(formula_pars, data = final_data_16S, family = quasipoisson("log"), method = "REML")
print(summary(gam_pars))

# Save GAM diagnostics for Richness 16S
pdf("gam_diagnostics_Morpho_richness_16S.pdf", width = 12, height = 10)
par(mfrow = c(2, 2))  
gam.check(gam_pars, rep = 0)
par(mfrow = c(1, 1))  
dev.off()

# Evenness 16S
formula_evenness_full <- as.formula("Evenness_16S ~ s(Shoot_Length, k=15) + s(Leaf_Width, k=9) + s(No_Leaves, k=7)")
gam_evenness_full <- gam(formula_evenness_full, data = final_data_16S, family = betar(link = "logit"), method = "REML")
print(summary(gam_evenness_full))

# Parsimonious Evenness 16S
formula_evenness <- as.formula("Evenness_16S ~ s(Shoot_Length, k=15)")
gam_evenness <- gam(formula_evenness, data = final_data_16S, family = betar(link = "logit"), method = "REML")
print(summary(gam_evenness))

# Save GAM diagnostics for Evenness 16S
pdf("gam_diagnostics_Morpho_Evenness_16S.pdf", width = 12, height = 10)
par(mfrow = c(2, 2))  
gam.check(gam_evenness, rep = 0)
par(mfrow = c(1, 1))  
dev.off()


# ===========
# 18S Dataset
# ===========

# Load data 18S
sample.data_18S <- read.csv("Metadata_18S.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
sample.data_18S$Evenness_18S <- sample.data_18S$Shannon_18S / log(sample.data_18S$Richness_18S)
filtered_data_18S <- sample.data_18S %>% dplyr::select(Shoot_Length, Leaf_Width, No_Leaves, Richness_18S, Evenness_18S) %>% na.omit()

# Define variable
mor_var_18S <- c("Shoot_Length", "Leaf_Width", "No_Leaves")
response_richness_18S <- "Richness_18S"
response_evenness_18S <- "Evenness_18S"

# Standardization 
pp_18S <- preProcess(filtered_data_18S[, mor_var_18S], method = c("center", "scale"))
scaled_morpho_18S <- predict(pp_18S, filtered_data_18S[, mor_var_18S])
final_data_18S <- cbind(scaled_morpho_18S, filtered_data_18S[, c(response_richness_18S, response_evenness_18S)])

# Richness 18S
unique_counts_18S <- sapply(mor_var_18S, function(v) length(unique(final_data_18S[[v]])))
k_vals_18S <- sapply(unique_counts_18S, calculate_optimal_k)
gam_terms_18S <- mapply(function(var, k) paste0("s(", var, ", k=", k, ")"), names(k_vals_18S), k_vals_18S, USE.NAMES = FALSE)
formula_full_18S <- as.formula(paste(response_richness_18S, "~", paste(gam_terms_18S, collapse = " + ")))

# Define family
overdispersion_18S <- var(final_data_18S[[response_richness_18S]]) / mean(final_data_18S[[response_richness_18S]])
family_used_18S <- if (overdispersion_18S > 2) quasipoisson("log") else poisson("log")
cat("\nSelected family:", if (overdispersion_18S > 2) "quasipoisson" else "poisson", "\n")

# Full Richness 18S
gam_full_18S <- gam(formula_full_18S, data = final_data_18S, family = family_used_18S, method = "REML")
print(summary(gam_full_18S))

# Parsimonious Richness 18S
formula_pars_18S <- as.formula("Richness_18S ~ s(Shoot_Length, k=15)")
gam_pars_18S <- gam(formula_pars_18S, data = final_data_18S, family = quasipoisson("log"), method = "REML")
print(summary(gam_pars_18S))

# Save GAM diagnostics for Richness 18S
pdf("gam_diagnostics_Morpho_richness_18S.pdf", width = 12, height = 10)
par(mfrow = c(2, 2))  
gam.check(gam_pars_18S, rep = 0)
par(mfrow = c(1, 1))
dev.off()


# Evenness 18S
formula_evenness_full_18S <- as.formula("Evenness_18S ~ s(Shoot_Length, k=15) + s(Leaf_Width, k=9) + s(No_Leaves, k=7)")
gam_evenness_full_18S <- gam(formula_evenness_full_18S, data = final_data_18S, family = betar(link = "logit"), method = "REML")
print(summary(gam_evenness_full_18S))

# Parsimonious Evenness 18S
formula_evenness_pars_18S <- as.formula("Evenness_18S ~ s(Shoot_Length, k=15)")
gam_evenness_pars_18S <- gam(formula_evenness_pars_18S, data = final_data_18S, family = betar(link = "logit"), method = "REML")
print(summary(gam_evenness_pars_18S))


#Save GAM diagnostics for Evenness 18S
pdf("gam_diagnostics_Morpho_evenness_18S.pdf", width = 12, height = 10)
par(mfrow = c(2, 2))
gam.check(gam_evenness_pars_18S, rep = 0)
par(mfrow = c(1, 1))
dev.off()
