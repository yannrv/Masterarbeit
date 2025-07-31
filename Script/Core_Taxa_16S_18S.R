#### CORE TAXA ANALYSIS: 16S AND 18S rRNA GENE SEQUENCING DATA ####

# Load Required Libraries
library(ggplot2)
library(ggupset)
library(tidyverse)
library(readr)
library(ggrepel)
library(knitr)
library(phyloseq)

# ===========
# 16S Dataset
# ===========

# Load 16S Data
ASV.table_16S <- read.csv("~/Documents/cvs file/ASV_tab_16S_mix.csv", header = TRUE, row.names = 1)
Phyloseq_16S_mix_clean <- readRDS("~/Documents/phyloseq/Phyloseq_16S_mix_clean.rds")

# Compute Core ASVs
freq_ASVs <- apply(ASV.table_16S > 0, 2, sum)
abund_ASVs <- colSums(ASV.table_16S)
total_samples <- nrow(ASV.table_16S)
threshold_50 <- total_samples * 0.5
total_ASVs <- length(freq_ASVs)

# Dataframe for core ASVs
plot_data_16S <- data.frame(
  ASV = names(freq_ASVs),
  Abundance = abund_ASVs,
  Frequency = freq_ASVs
) %>%
  mutate(is_core = Frequency >= threshold_50,
         label = ifelse(is_core, ASV, NA))

# Summary
n_core <- sum(plot_data_16S$is_core)
percent_core <- round(n_core / total_ASVs * 100, 1)

# Abundance vs Frequency plot
core_plot_16S <- ggplot(plot_data_16S, aes(x = Abundance, y = Frequency, color = is_core)) +
  geom_point(aes(size = is_core), alpha = 0.7) +
  geom_text_repel(aes(label = label), size = 3) +
  geom_hline(yintercept = threshold_50, color = "red", linetype = "dashed") +
  scale_x_log10() +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#2166AC")) +
  labs(title = "Core ASVs (16S)", subtitle = paste(n_core, "core ASVs (", percent_core, "%)"), 
       x = "Total Abundance (log10)", y = "Frequency (Number of Samples)") +
  theme_minimal()
ggsave("Core_16S_plot.png", core_plot_16S, width = 10, height = 6)

# Extract core ASVs by Dataset
datasets <- unique(sample_data(Phyloseq_16S_mix_clean)$Dataset)
core_asvs <- plot_data_16S$ASV[plot_data_16S$is_core]
asv_by_dataset_16S <- list()
for (dataset in datasets) {
  ps_subset <- subset_samples(Phyloseq_16S_mix_clean, Dataset == dataset)
  asv_table <- otu_table(ps_subset)
  if (!is.null(dim(asv_table))) {
    asv_counts <- apply(asv_table > 0, 2, sum)
    threshold_ds <- nrow(sample_data(ps_subset)) * 0.5
    core_in_ds <- names(asv_counts[asv_counts >= threshold_ds])
    asv_by_dataset_16S[[dataset]] <- intersect(core_in_ds, core_asvs)
  }
}
# UpSet input from core ASVs
presence_df_16S <- tibble(ASV = unique(unlist(asv_by_dataset_16S))) %>%
  mutate(Datasets = map(ASV, function(asv) {
    keep(names(asv_by_dataset_16S), function(dataset) {
      asv %in% asv_by_dataset_16S[[dataset]]
    })
  }))

# UpSet plot
upset_plot_16S <- ggplot(presence_df_16S, aes(x = Datasets)) +
  geom_bar(fill = "#4A90E2") +
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5) +
  scale_x_upset(n_intersections = 30) +
  theme_minimal() +
  labs(title = "ASV Intersections Across Datasets (16S)", y = "Number of ASVs")
ggsave("Core_16S_upset.png", upset_plot_16S, width = 10, height = 6)

# Taxonomy Table
tax_core_16S <- tax_table(Phyloseq_16S_mix_clean)[core_asvs, ]
df_tax_16S <- data.frame(
  ASV = core_asvs,
  Phylum = tax_core_16S[, "phylum"],
  Class = tax_core_16S[, "class"],
  Order = tax_core_16S[, "order"],
  Family = tax_core_16S[, "family"],
  Genus = tax_core_16S[, "genus"]
)
df_tax_16S <- df_tax_16S[order(df_tax_16S$Phylum, df_tax_16S$Class), ]
print(knitr::kable(df_tax_16S, format = "markdown"))

# Remove List
rm(list = ls())

# ============
# 18S Dataset
# ============

# Load 18S ASV abundance table and phyloseq object
ASV_18S <- read_csv("~/Documents/cvs file/ASV_tab_18S_mix.csv")
Phyloseq_18S_mix <- readRDS("~/Documents/phyloseq/Phyloseq_SS_18S_mix_clean.rds")
ASV_matrix <- as.matrix(ASV_18S[,-1])  

# Compute Core ASVs
freq_ASVs <- apply(ASV_matrix > 0, 2, sum)
abund_ASVs <- colSums(ASV_matrix)
total_samples <- nrow(ASV_matrix)
threshold_50 <- total_samples * 0.5
total_ASVs <- length(freq_ASVs)

# Dataframe for core ASVs
plot_data <- data.frame(
  ASV = names(freq_ASVs),
  Abundance = abund_ASVs,
  Frequency = freq_ASVs
) %>%
  mutate(
    is_core = Frequency >= threshold_50,
    label = ifelse(is_core, ASV, NA)
  )

# Summary
n_core <- sum(plot_data$is_core)
percent_core <- round(n_core / total_ASVs * 100, 1)

# Abundance vs Frequency plot for 18S
core_plot_18S <- ggplot(plot_data, aes(x = Abundance, y = Frequency, color = is_core)) +
  geom_point(aes(size = is_core), alpha = 0.7) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = Inf, box.padding = 0.5) +
  geom_hline(yintercept = threshold_50, color = "#E41A1C", linetype = "dashed", size = 0.8) +
  scale_x_log10(labels = scales::comma) +
  scale_y_continuous(limits = c(0, max(freq_ASVs) * 1.1)) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#2166AC"), labels = c("Other ASVs", "Core Taxa")) +
  scale_size_manual(values = c("FALSE" = 2, "TRUE" = 3), guide = "none") +
  annotate("text", x = min(abund_ASVs), y = threshold_50 * 1.05,
           label = paste0("50% threshold\n(", round(threshold_50), " samples)"),
           hjust = 0, size = 3.5) +
  labs(
    x = "Abundance (log scale)",
    y = "Frequency (number of samples)",
    title = "ASV Frequency vs Abundance Distribution (18S)",
    subtitle = paste0(n_core, " core ASVs (", percent_core, "% of total) present in >50% of samples"),
    color = "ASV Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "top"
  )
core_plot_18S
ggsave("Core_18S.png", plot = core_plot_18S, width = 20, height = 10, dpi = 300, bg = "white")

# Extract core ASVs by Dataset
core_asvs <- names(freq_ASVs[freq_ASVs >= threshold_50])
datasets <- c("Seastore_Project", "Hiddensee_Study", "Salinity_Gradient", "Seagrass_Vibrio", "Germination_Experiment")
asv_by_dataset <- list()
for (dataset in datasets) {
  ps_subset <- subset_samples(Phyloseq_18S_mix, Dataset == dataset)
  asv_table <- otu_table(ps_subset)
  if (!is.null(dim(asv_table))) {
    asv_counts <- apply(asv_table > 0, 2, sum)
    dataset_samples <- nrow(sample_data(ps_subset))
    dataset_threshold <- dataset_samples * 0.5
    core_in_dataset <- names(asv_counts[asv_counts >= dataset_threshold])
    asv_by_dataset[[dataset]] <- intersect(core_in_dataset, core_asvs)
  } else {
    asv_by_dataset[[dataset]] <- character(0)
  }
}

# UpSet input from core ASVs
all_core_asvs <- unique(unlist(asv_by_dataset))
presence_df <- tibble(ASV = all_core_asvs) %>%
  mutate(Datasets = map(ASV, function(asv) {
    names(asv_by_dataset)[map_lgl(asv_by_dataset, ~ asv %in% .)]
  }))

# UpSet plot 
upset_plot_18S <- ggplot(presence_df, aes(x = Datasets)) +
  geom_bar(fill = "#4A90E2", alpha = 0.8) +
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5, size = 5, fontface = "bold") +
  scale_x_upset(n_intersections = 30) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    panel.grid.major = element_line(colour = "grey90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  labs(title = "ASV Intersections Across Datasets (18S)", y = "Number of ASVs")
upset_plot_18S
ggsave("Core_upset_plot_18S.png", plot = upset_plot_18S, width = 15, height = 10, dpi = 300, bg = "white")

# Taxonomy table
tax_table <- tax_table(Phyloseq_18S_mix)
core_tax <- tax_table[all_core_asvs, ]

df_tax <- data.frame(
  ASV = all_core_asvs,
  Phylum = core_tax[, "rank_III"],
  Class = core_tax[, "rank_IV"],
  Order = core_tax[, "rank_IX"],
  Family = core_tax[, "rank_X"],
  Genus = core_tax[, "rank_XI"]
)
df_tax <- df_tax[order(df_tax$Phylum, df_tax$Class, df_tax$Order, df_tax$Family, df_tax$Genus), ]

knitr::kable(df_tax,
             col.names = c("ASV ID", "Phylum", "Class", "Order", "Family", "Genus"),
             format = "markdown",
             align = c("l", "l", "l", "l", "l", "l"))
