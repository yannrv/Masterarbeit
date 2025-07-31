#### Samples Variables Effects on Microbial Communities (16S and 18S) ####

# Load required packages
library(vegan)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(grid)
library(dplyr)

# ===========
# 16S Dataset
# ===========

# Load ASV and Metadata
ASV.table_16S <- read.csv("~/ASV_tab_16S_mix.csv", row.names = 1, stringsAsFactors = FALSE)
sample.data_16S <- read.csv("~/Metadata_16S.csv", row.names = 1, stringsAsFactors = FALSE)

# Select sampling variables
Sampling_var <- sample.data_16S %>%
  select(Sample_ID, Sample_type, Dislodgment_method, Dataset) %>%
  na.omit()
rownames(Sampling_var) <- Sampling_var$Sample_ID

# Filter to shared samples
common_samples <- intersect(rownames(ASV.table_16S), rownames(Sampling_var))
ASV.table.filtered <- ASV.table_16S[common_samples, ]
Sampling_var.filtered <- Sampling_var[common_samples, ]

# Hellinger transformation and Bray-Curtis distance
ASV.table.nor <- decostand(ASV.table.filtered, method = "hellinger")
dist_bray <- vegdist(ASV.table.nor, method = "bray")

# NMDS ordination
set.seed(123)
nmds_result <- metaMDS(dist_bray, k = 3, trymax = 100)
nmds_coords <- as.data.frame(scores(nmds_result))

# Add metadata
nmds_coords$Sample_type <- Sampling_var.filtered$Sample_type
nmds_coords$Dislodgment_method <- Sampling_var.filtered$Dislodgment_method
nmds_coords$Dataset <- Sampling_var.filtered$Dataset

# Define colors and shapes
sample_type_colors <- c("old_leaf" = "red", "whole_plant" = "blue", "young_leaf" = "darkgreen")
dislodgment_colors <- c("swabbing" = "#33A02C", "leaf_piece" = "#FF7F00")
dataset_shapes <- c("Germination_Experiment" = 16, "Hiddensee_Study" = 17, "Salinity_Gradient" = 15, "Seagrass_Vibrio" = 3, "Seastore_Project" = 8)

# NMDS Plots by sample type
nmds_plot_16S_sample <- ggplot(nmds_coords, aes(NMDS1, NMDS2, color = Sample_type, shape = Dataset)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Sample_type), level = 0.95, linetype = 2) +
  scale_color_manual(values = sample_type_colors) +
  scale_shape_manual(values = dataset_shapes) +
  theme_classic() +
  labs(title = "16S NMDS by Sample Type", subtitle = paste("Stress:", round(nmds_result$stress, 4))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# NMDS Plots by dislodgment method
nmds_plot_16S_dislodgment <- ggplot(nmds_coords, aes(NMDS1, NMDS2, color = Dislodgment_method, shape = Dataset)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Dislodgment_method), level = 0.95, linetype = 2) +
  scale_color_manual(values = dislodgment_colors) +
  scale_shape_manual(values = dataset_shapes) +
  theme_classic() +
  labs(title = "16S NMDS by Dislodgment Method", subtitle = paste("Stress:", round(nmds_result$stress, 4))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Save plots
ggsave("NMDS_16S_Sample_Type.png", nmds_plot_16S_sample, width = 8, height = 6, dpi = 300)
ggsave("NMDS_16S_Dislodgment_Method.png", nmds_plot_16S_dislodgment, width = 8, height = 6, dpi = 300)

# PERMANOVA test
adonis_sample <- adonis2(dist_bray ~ Sample_type, data = Sampling_var.filtered, permutations = 999)
adonis_dislodgment <- adonis2(dist_bray ~ Dislodgment_method, data = Sampling_var.filtered, permutations = 999)

#Print results
print(adonis_sample)
print(adonis_dislodgment)

# Betadisper and Tukey tests
bd_sample <- betadisper(dist_bray, Sampling_var.filtered$Sample_type)
bd_dislodgment <- betadisper(dist_bray, Sampling_var.filtered$Dislodgment_method)

#Print results
print(permutest(bd_sample, pairwise = TRUE, permutations = 999))
print(permutest(bd_dislodgment, pairwise = TRUE, permutations = 999))
print(TukeyHSD(bd_sample))
print(TukeyHSD(bd_dislodgment))

#remove Environment
rm(list = ls())

# ===========
# 18S Dataset
# ===========

# Load ASV and Metadata
ASV.table_18S <- read.csv("~/ASV_tab_18S_mix.csv", row.names = 1, stringsAsFactors = FALSE)
sample.data_18S <- read.csv("~/Metadata_18S.csv", row.names = 1, stringsAsFactors = FALSE)

# Select sampling variables
Sampling_var <- sample.data_18S %>%
  select(Sample_ID, Sample_type, Dislodgment_method, Dataset) %>%
  na.omit()
rownames(Sampling_var) <- Sampling_var$Sample_ID

# Filter to shared samples
common_samples <- intersect(rownames(ASV.table_18S), rownames(Sampling_var))
ASV.table.filtered <- ASV.table_18S[common_samples, ]
Sampling_var.filtered <- Sampling_var[common_samples, ]

# Hellinger transformation and Bray-Curtis distance
ASV.table.nor <- decostand(ASV.table.filtered, method = "hellinger")
dist_bray <- vegdist(ASV.table.nor, method = "bray")

# NMDS ordination
set.seed(123)
nmds_result <- metaMDS(dist_bray, k = 3, trymax = 100)
nmds_coords <- as.data.frame(scores(nmds_result))

# Add metadata
nmds_coords$Sample_type <- Sampling_var.filtered$Sample_type
nmds_coords$Dislodgment_method <- Sampling_var.filtered$Dislodgment_method
nmds_coords$Dataset <- Sampling_var.filtered$Dataset

# Define colors and shapes
sample_type_colors <- c("old_leaf" = "red", "whole_plant" = "blue", "young_leaf" = "darkgreen")
dislodgment_colors <- c("swabbing" = "#33A02C", "leaf_piece" = "#FF7F00")
dataset_shapes <- c("Germination_Experiment" = 16, "Hiddensee_Study" = 17, "Salinity_Gradient" = 15, "Seagrass_Vibrio" = 3, "Seastore_Project" = 8)

# NMDS Plots by sample type
nmds_plot_18S_sample <- ggplot(nmds_coords, aes(NMDS1, NMDS2, color = Sample_type, shape = Dataset)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Sample_type), level = 0.95, linetype = 2) +
  scale_color_manual(values = sample_type_colors) +
  scale_shape_manual(values = dataset_shapes) +
  theme_classic() +
  labs(title = "18S NMDS by Sample Type", subtitle = paste("Stress:", round(nmds_result$stress, 4))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# NMDS Plots by dislodgment method
nmds_plot_18S_dislodgment <- ggplot(nmds_coords, aes(NMDS1, NMDS2, color = Dislodgment_method, shape = Dataset)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Dislodgment_method), level = 0.95, linetype = 2) +
  scale_color_manual(values = dislodgment_colors) +
  scale_shape_manual(values = dataset_shapes) +
  theme_classic() +
  labs(title = "18S NMDS by Dislodgment Method", subtitle = paste("Stress:", round(nmds_result$stress, 4))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Save plot
ggsave("NMDS_18S_Sample_Type.png", nmds_plot_18S_sample, width = 8, height = 6, dpi = 300)
ggsave("NMDS_18S_Dislodgment_Method.png", nmds_plot_18S_dislodgment, width = 8, height = 6, dpi = 300)

# PERMANOVA test
adonis_sample <- adonis2(dist_bray ~ Sample_type, data = Sampling_var.filtered, permutations = 999)
adonis_dislodgment <- adonis2(dist_bray ~ Dislodgment_method, data = Sampling_var.filtered, permutations = 999)

# Print results
print(adonis_sample)
print(adonis_dislodgment)

# Betadisper and Tukey tests
bd_sample <- betadisper(dist_bray, Sampling_var.filtered$Sample_type)
bd_dislodgment <- betadisper(dist_bray, Sampling_var.filtered$Dislodgment_method)

#Print results
print(permutest(bd_sample, pairwise = TRUE, permutations = 999))
print(permutest(bd_dislodgment, pairwise = TRUE, permutations = 999))
print(TukeyHSD(bd_sample))
print(TukeyHSD(bd_dislodgment))