#### EFFECTS OF SAMPLING VARIABLES ON ALPHA DIVERSITY ####

# Load required libraries
library(readr) 
library(vegan)
library(ggplot2)
library(dplyr)
library(FSA)
library(ggsignif)  
library(patchwork)
library(cowplot)

# Load sample data
Metadata_16S <- read.csv("~/Metadata_16S.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
Metadata_18S <- read.csv("~/Metadata_18S.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Calculate Evenness
Metadata_16S$Evenness_16S <- Metadata_16S$Shannon_16S / log(Metadata_16S$Richness_16S)
Metadata_18S$Evenness_18S <- Metadata_18S$Shannon_18S / log(Metadata_18S$Richness_18S)

# ========================
# Analysis by Sample type
# ========================

## Richness 16S

# Perform Kruskal-Wallis test
kruskal_test_16S <- kruskal.test(Richness_16S ~ Sample_type, data = Metadata_16S)
print(kruskal_test_16S)

# Perform post hoc Dunn test with Benjamini-Hochberg correction
dunn_test_16S <- dunnTest(Richness_16S ~ Sample_type, data = Metadata_16S, method = "bh")
print(dunn_test_16S)

# Richness plot by sample type
richness_sample_type_plot_16S <- ggplot(data = Metadata_16S, 
                                        aes(x = Sample_type, 
                                            y = Richness_16S, 
                                            fill = Sample_type)) +
  geom_violin(trim = FALSE, color = "black", width = 0.7, linewidth = 0.8) +
  geom_boxplot(width = 0.1, color = "white", outlier.shape = NA) +
  geom_jitter(aes(shape = Dataset), 
              color = "black", 
              position = position_jitter(0.2), 
              size = 2, 
              alpha = 0.6) +
  stat_summary(fun.data = function(x) data.frame(y = -100, 
                                                 label = paste("n =", length(x))),
               geom = "text", 
               color = "black", 
               size = 4, 
               vjust = 1.2) +
  scale_fill_manual(values = c("old_leaf" = "#FDD166", 
                               "whole_plant" = "#118AB2", 
                               "young_leaf" = "#87986A")) +
  scale_shape_manual(values = c(16, 17, 15, 3, 8)) +
  labs(title = "16S Richness by Sample Type", 
       x = "Sample Type", 
       y = "16S Richness (rarefied)") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.3)
  ) +
  scale_y_continuous(limits = c(-100, 1100 * 1.2), breaks = seq(0, 1000 * 1.2, by = 200)) +
  # Add significance annotations
  geom_signif(comparisons = list(c("old_leaf", "whole_plant"),
                                 c("old_leaf", "young_leaf"),
                                 c("whole_plant", "young_leaf")),
              map_signif_level = TRUE,
              textsize = 5, 
              y_position = c(1100, 1170, 1250),
              tip_length = 0.01)

print(richness_sample_type_plot_16S)

## Richness 18S

# Perform Kruskal-Wallis test
kruskal_test_18S <- kruskal.test(Richness_18S ~ Sample_type, data = Metadata_18S)
print(kruskal_test_18S)

# Perform post hoc Dunn test
dunn_test_18S <- dunnTest(Richness_18S ~ Sample_type, data = Metadata_18S, method = "bh")
print(dunn_test_18S)

# Create 18S richness plot by sample type
richness_sample_type_plot_18S <- ggplot(data = Metadata_18S, 
                                        aes(x = Sample_type, 
                                            y = Richness_18S, 
                                            fill = Sample_type)) +
  geom_violin(trim = FALSE, color = "black", width = 0.7, linewidth = 0.8) +
  geom_jitter(aes(shape = Dataset), 
              color = "black", 
              position = position_jitter(0.2), 
              size = 2, 
              alpha = 0.6) +
  stat_summary(fun.data = function(x) data.frame(y = -100, 
                                                 label = paste("n =", length(x))),
               geom = "text", 
               color = "black", 
               size = 4, 
               vjust = 1.2) +
  scale_fill_manual(values = c("old_leaf" = "#FDD166", 
                               "whole_plant" = "#118AB2", 
                               "young_leaf" = "#87986A")) +
  scale_shape_manual(values = c(16, 17, 15, 3, 8)) +
  labs(title = "18S Richness by Sample Type", 
       x = "Sample Type", 
       y = "18S Richness (rarefied)") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.3)
  ) +
  scale_y_continuous(limits = c(-100, 400 * 1.2), breaks = seq(0, 400 * 1.2, by = 100)) +
  # Add significance annotations
  geom_signif(comparisons = list(c("old_leaf", "whole_plant"),
                                 c("old_leaf", "young_leaf"),
                                 c("whole_plant", "young_leaf")),
              map_signif_level = TRUE,
              textsize = 5, 
              y_position = c(350, 400, 450),
              tip_length = 0.01)

print(richness_sample_type_plot_18S)

# Combine richness plots for sample type
richness_sample_type_plot_16S_no_legend <- richness_sample_type_plot_16S + theme(legend.position = "none")
richness_sample_type_plot_18S_no_legend <- richness_sample_type_plot_18S + theme(legend.position = "none")

shared_legend <- get_legend(richness_sample_type_plot_16S)

final_richness_plot <- plot_grid(
  plot_grid(
    richness_sample_type_plot_16S_no_legend, 
    richness_sample_type_plot_18S_no_legend, 
    ncol = 2, 
    labels = c("A", "B"),
    label_size = 14
  ),
  shared_legend,
  ncol = 2, 
  rel_widths = c(3, 1)
)

print(final_richness_plot)
ggsave("Richness_by_Sample_type.png", plot = final_richness_plot, width = 15, height = 7, bg = "white")

## Evenness 16S

# Perform Kruskal-Wallis test
kruskal_test_evenness_16S <- kruskal.test(Evenness_16S ~ Sample_type, data = Metadata_16S)
print(kruskal_test_evenness_16S)

# Perform post hoc Dunn test
dunn_test_evenness_16S <- dunnTest(Evenness_16S ~ Sample_type, data = Metadata_16S, method = "bh")
print(dunn_test_evenness_16S)

# Create 16S evenness plot by sample type
evenness_sample_type_plot_16S <- ggplot(data = Metadata_16S, 
                                        aes(x = Sample_type, 
                                            y = Evenness_16S, 
                                            fill = Sample_type)) +
  geom_violin(trim = FALSE, color = "black", width = 0.7, linewidth = 0.8) +
  geom_boxplot(width = 0.1, color = "white", outlier.shape = NA) +
  geom_jitter(aes(shape = Dataset), 
              color = "black", 
              position = position_jitter(0.2), 
              size = 2, 
              alpha = 0.6) +
  stat_summary(fun.data = function(x) data.frame(y = 0, 
                                                 label = paste("n =", length(x))),
               geom = "text", 
               color = "black", 
               size = 4, 
               vjust = 1.2) +
  scale_fill_manual(values = c("old_leaf" = "#FDD166", 
                               "whole_plant" = "#118AB2", 
                               "young_leaf" = "#87986A")) +
  scale_shape_manual(values = c(16, 17, 15, 3, 8)) +
  labs(title = "Evenness by Sample Type (16S)", 
       x = "Sample Type", 
       y = "16S Evenness") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.3)
  ) +
  scale_y_continuous(limits = c(-0.1, 1 * 1.3), breaks = seq(0, 1 * 1.3, by = 0.5)) +
  geom_signif(comparisons = list(c("old_leaf", "whole_plant"),
                                 c("old_leaf", "young_leaf"),
                                 c("whole_plant", "young_leaf")),
              map_signif_level = TRUE,
              textsize = 5, 
              y_position = c(1, 1.1, 1.2),
              tip_length = 0.01)

print(evenness_sample_type_plot_16S)

## Evenness 18S

# Perform Kruskal-Wallis test
kruskal_test_evenness_18S <- kruskal.test(Evenness_18S ~ Sample_type, data = Metadata_18S)
print(kruskal_test_evenness_18S)

# Perform post hoc Dunn test
dunn_test_evenness_18S <- dunnTest(Evenness_18S ~ Sample_type, data = Metadata_18S, method = "bh")
print(dunn_test_evenness_18S)

# Prepare custom annotations for 18S evenness
dunn_results_18S <- dunn_test_evenness_18S$res %>%
  mutate(
    group1 = sub(" - .*", "", Comparison),
    group2 = sub(".*- ", "", Comparison)
  )

annotations_18S <- ifelse(
  dunn_results_18S$P.adj > 0.05, 
  "NS",  # Non-significant
  sprintf("p = %.3f", dunn_results_18S$P.adj)
)

comparisons_18S <- list(
  c("old_leaf", "whole_plant"),
  c("old_leaf", "young_leaf"),
  c("whole_plant", "young_leaf")
)

# Evenness 18S plot by sample type
evenness_sample_type_plot_18S <- ggplot(data = Metadata_18S, 
                                        aes(x = Sample_type, 
                                            y = Evenness_18S, 
                                            fill = Sample_type)) +
  geom_violin(trim = FALSE, color = "black", width = 0.7, linewidth = 0.8) +
  geom_boxplot(width = 0.1, color = "white", outlier.shape = NA) +
  stat_summary(fun.data = function(x) data.frame(y = -0.3, 
                                                 label = paste("n =", length(x))),
               geom = "text", 
               color = "black", 
               size = 4, 
               vjust = 1.2) +
  geom_jitter(aes(shape = Dataset), 
              color = "black", 
              position = position_jitter(0.2), 
              size = 2, 
              alpha = 0.6) +
  scale_fill_manual(values = c("old_leaf" = "#FDD166", 
                               "whole_plant" = "#118AB2", 
                               "young_leaf" = "#87986A")) +
  scale_shape_manual(values = c(16, 17, 15, 3, 8)) +
  labs(title = "Evenness by Sample Type (18S)", 
       x = "Sample Type", 
       y = "18S Evenness") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.3)
  ) +
  scale_y_continuous(limits = c(-0.5, 1.15 * 1.3), breaks = seq(0, 1.5 * 1.3, by = 0.5)) +
  geom_signif(
    comparisons = comparisons_18S,
    annotations = annotations_18S,
    textsize = 4,
    y_position = c(1.15, 1.25, 1.35),
    tip_length = 0.01
  )

print(evenness_sample_type_plot_18S)

# Combine evenness plots for sample type
evenness_sample_type_plot_16S_no_legend <- evenness_sample_type_plot_16S + theme(legend.position = "none")
evenness_sample_type_plot_18S_no_legend <- evenness_sample_type_plot_18S + theme(legend.position = "none")

shared_legend_evenness <- get_legend(evenness_sample_type_plot_16S)

final_evenness_plot <- plot_grid(
  plot_grid(
    evenness_sample_type_plot_16S_no_legend, 
    evenness_sample_type_plot_18S_no_legend, 
    ncol = 2, 
    labels = c("A", "B"),
    label_size = 14
  ),
  shared_legend_evenness,
  ncol = 2, 
  rel_widths = c(3, 1)
)

print(final_evenness_plot)
ggsave("Evenness_by_Sample_type.png", plot = final_evenness_plot, width = 15, height = 7, bg = "white")

# ==============================
# Analysis by Dislodgment method
# ==============================

## Richness Analysis

# Richness 16S by dislodgment method
wilcox_test_16S_richness <- wilcox.test(Richness_16S ~ Dislodgment_method, 
                                        p.adjust.method = "BH", 
                                        data = Metadata_16S)
print("16S Richness by Dislodgment Method - Wilcoxon test:")
print(wilcox_test_16S_richness)

# Richness 18S by dislodgment method
wilcox_test_18S_richness <- wilcox.test(Richness_18S ~ Dislodgment_method, 
                                        p.adjust.method = "BH", 
                                        data = Metadata_18S)
print("18S Richness by Dislodgment Method - Wilcoxon test:")
print(wilcox_test_18S_richness)

# Create richness 16S plot by dislodgment method
richness_dislodgment_plot_16S <- ggplot(Metadata_16S, 
                                        aes(x = Dislodgment_method, 
                                            y = Richness_16S, 
                                            fill = Dislodgment_method)) +
  geom_violin(trim = FALSE, color = "black", width = 0.7) +
  geom_boxplot(width = 0.1, color = "white", outlier.shape = NA, linewidth = 0.5) +
  geom_jitter(aes(shape = Dataset), color = "black", position = position_jitter(0.2), size = 2, alpha = 0.6) +
  stat_summary(fun.data = function(x) data.frame(y = -100, label = paste("n =", length(x))),
               geom = "text", color = "black", size = 4, vjust = 1) +
  scale_fill_manual(values = c("swabbing" = "#ffab41", "leaf_piece" = "#06d6a0")) +
  scale_shape_manual(values = c(16, 17, 15, 3, 8)) +
  labs(title = "16S Richness by Dislodgment Method", 
       x = "Dislodgment Method", 
       y = "16S Richness (rarefied)") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.3)
  ) +
  scale_y_continuous(limits = c(-100, max(Metadata_16S$Richness_16S, na.rm = TRUE) * 1.4), 
                     breaks = seq(0, max(Metadata_16S$Richness_16S, na.rm = TRUE) * 1.4, by = 200)) +
  geom_signif(comparisons = list(c("swabbing", "leaf_piece")),
              annotations = "***",  
              y_position = max(Metadata_16S$Richness_16S, na.rm = TRUE) * 1.3,  
              tip_length = 0.01, 
              textsize = 4)

print(richness_dislodgment_plot_16S)

# Create richness 18S plot by dislodgment method
richness_dislodgment_plot_18S <- ggplot(Metadata_18S, 
                                        aes(x = Dislodgment_method, 
                                            y = Richness_18S, 
                                            fill = Dislodgment_method)) +
  geom_violin(trim = FALSE, color = "black", width = 0.7) +
  geom_boxplot(width = 0.1, color = "white", outlier.shape = NA, linewidth = 0.5) +
  geom_jitter(aes(shape = Dataset), color = "black", position = position_jitter(0.2), size = 2, alpha = 0.6) +
  stat_summary(fun.data = function(x) data.frame(y = -100, label = paste("n =", length(x))),
               geom = "text", color = "black", size = 4, vjust = 1) +
  scale_fill_manual(values = c("swabbing" = "#ffab41", "leaf_piece" = "#06d6a0")) +
  scale_shape_manual(values = c(16, 17, 15, 3, 8)) +
  labs(title = "18S Richness by Dislodgment Method", 
       x = "Dislodgment Method", 
       y = "18S Richness (rarefied)") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.3)
  ) +
  scale_y_continuous(limits = c(-100, max(Metadata_18S$Richness_18S, na.rm = TRUE) * 1.4), 
                     breaks = seq(0, max(Metadata_18S$Richness_18S, na.rm = TRUE) * 1.4, by = 100)) +
  geom_signif(comparisons = list(c("swabbing", "leaf_piece")),
              annotations = "***",  
              y_position = max(Metadata_18S$Richness_18S, na.rm = TRUE) * 1.3,  
              tip_length = 0.01, 
              textsize = 4)

print(richness_dislodgment_plot_18S)

# Combine richness plots for dislodgment method
combined_richness_dislodgment_plot <- (richness_dislodgment_plot_16S + richness_dislodgment_plot_18S) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

print(combined_richness_dislodgment_plot)
ggsave("Richness_by_Dislodgment_method.png", plot = combined_richness_dislodgment_plot,
       width = 15, height = 7, dpi = 300, units = "in", bg = "white")

## Evenness Analysis

# Evenness 16S by dislodgment method
wilcox_test_16S_evenness <- wilcox.test(Evenness_16S ~ Dislodgment_method, 
                                        p.adjust.method = "BH", 
                                        data = Metadata_16S)
print(wilcox_test_16S_evenness)

# 18S evenness by dislodgment method
wilcox_test_18S_evenness <- wilcox.test(Evenness_18S ~ Dislodgment_method, 
                                        p.adjust.method = "BH", 
                                        data = Metadata_18S)
print(wilcox_test_18S_evenness)

# Create evenness 16S plot by dislodgment method
evenness_dislodgment_plot_16S <- ggplot(Metadata_16S, 
                                        aes(x = Dislodgment_method, 
                                            y = Evenness_16S, 
                                            fill = Dislodgment_method)) +
  geom_violin(trim = FALSE, color = "black", width = 0.7) +
  geom_boxplot(width = 0.1, color = "white", outlier.shape = NA, linewidth = 0.5) +
  geom_jitter(aes(shape = Dataset), color = "black", position = position_jitter(0.2), size = 2, alpha = 0.6) +
  stat_summary(fun.data = function(x) data.frame(y = 0, label = paste("n =", length(x))),
               geom = "text", color = "black", size = 4, vjust = 1) +
  scale_fill_manual(values = c("swabbing" = "#ffab41", "leaf_piece" = "#06d6a0")) +
  scale_shape_manual(values = c(16, 17, 15, 3, 8)) +
  labs(title = "16S Evenness by Dislodgment Method", 
       x = "Dislodgment Method", 
       y = "16S Evenness") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.3)
  ) +
  scale_y_continuous(limits = c(-0.1, max(Metadata_16S$Evenness_16S, na.rm = TRUE) * 1.4), 
                     breaks = seq(0, max(Metadata_16S$Evenness_16S, na.rm = TRUE) * 1.4, by = 0.4)) +
  geom_signif(comparisons = list(c("swabbing", "leaf_piece")),
              annotations = "***",
              y_position = max(Metadata_16S$Evenness_16S, na.rm = TRUE) * 1.3,
              tip_length = 0.01,
              textsize = 4)

print(evenness_dislodgment_plot_16S)

# Create evenness 18S plot by dislodgment method
evenness_dislodgment_plot_18S <- ggplot(Metadata_18S, 
                                        aes(x = Dislodgment_method, 
                                            y = Evenness_18S, 
                                            fill = Dislodgment_method)) +
  geom_violin(trim = FALSE, color = "black", width = 0.7) +
  geom_boxplot(width = 0.1, color = "white", outlier.shape = NA, linewidth = 0.5) +
  geom_jitter(aes(shape = Dataset), color = "black", position = position_jitter(0.2), size = 2, alpha = 0.6) +
  stat_summary(fun.data = function(x) data.frame(y = -0.25, label = paste("n =", length(x))),
               geom = "text", color = "black", size = 4, vjust = 1) +
  scale_fill_manual(values = c("swabbing" = "#ffab41", "leaf_piece" = "#06d6a0")) +
  scale_shape_manual(values = c(16, 17, 15, 3, 8)) +
  labs(title = "18S Evenness by Dislodgment Method", 
       x = "Dislodgment Method", 
       y = "18S Evenness") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.3)
  ) +
  scale_y_continuous(limits = c(-0.1, max(Metadata_18S$Evenness_18S, na.rm = TRUE) * 1.4), 
                     breaks = seq(0, max(Metadata_18S$Evenness_18S, na.rm = TRUE) * 1.4, by = 0.4)) +
  geom_signif(comparisons = list(c("swabbing", "leaf_piece")),
              annotations = "***",
              y_position = max(Metadata_18S$Evenness_18S, na.rm = TRUE) * 1.3,
              tip_length = 0.01,
              textsize = 4)

print(evenness_dislodgment_plot_18S)


# Combine evenness plots for dislodgment method
combined_evenness_dislodgment_plot <- (evenness_dislodgment_plot_16S + evenness_dislodgment_plot_18S) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

print(combined_evenness_dislodgment_plot)
ggsave("Evenness_by_Dislodgment_method.png", plot = combined_evenness_dislodgment_plot,
       width = 15, height = 7, dpi = 300, units = "in", bg = "white")


