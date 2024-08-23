## Microbiome Analysis Script for 16S Data

# Import required libraries
library(vegan)
library(ggplot2)
library(gridExtra)
library(dplyr)

# Set seed for reproducibility
set.seed(100)

# Load the necessary data tables (ASV table, taxonomy table, and sample metadata)
ASV.table <- read.csv("~/ASV_tab_16S.csv", header = TRUE, row.names = 1, sep = ",", stringsAsFactors = FALSE)
sample.data <- read.csv("~/Metadata.csv", header = TRUE, row.names = 1, sep = ",", stringsAsFactors = FALSE)
taxa.table <- read.csv("~/Tax_tab_16S.csv", header = TRUE, row.names = 1, sep = ",", stringsAsFactors = FALSE)

# Data Summary and Initial Checks
dim(ASV.table)
min(rowSums(ASV.table))

# Rarefy data 
ASV.table.rarefy <- rrarefy(ASV.table, min(rowSums(ASV.table)))


# Alpha Diversity Calculation and Visualization
# By Site and Area


#Calculate species richness and Shannon diversity
sample.data$richness <- specnumber(ASV.table.rarefy)  
sample.data$shannon <- diversity(ASV.table.rarefy, index = "shannon")

# Create a data frame with the selected metrics
alpha_metrics <- sample.data %>%
  dplyr::select(Site, Area, richness, shannon)


# Summarize alpha diversity by Site and Area
alpha_summary <- alpha_metrics %>%
  group_by(Site, Area) %>%
  summarise(Richness = mean(richness, na.rm = TRUE),
            Shannon = mean(shannon, na.rm = TRUE),
            .groups = 'drop')  # This removes grouping after summarization

# Save the summary to CSV file
write.csv(alpha_summary, "~/Alpha_summary_by_Site_and_Area.csv", row.names = FALSE)


# Plot for richness by site
richness_plot <- ggplot(sample.data, aes(x = Site, y = richness)) +
  stat_boxplot(geom = "errorbar", width = 0.6, color = "black") + 
  geom_boxplot(fill = "#66C2A5", color = "black") +
  labs(title = "16S Richness across Sites", x = "Site", y = "16S Richness (rarefied)") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_line(linewidth = 0.2, linetype = "dashed", color = "gray"),
    plot.margin = margin(10, 10, 10, 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.3)  
  ) +
  scale_y_continuous(limits = c(0, max(sample.data$richness, na.rm = TRUE) * 1.2), 
                     breaks = seq(0, max(sample.data$richness, na.rm = TRUE) * 1.2, by = 2))

# Plot for Shannon index by site
shannon_plot <- ggplot(sample.data, aes(x = Site, y = shannon)) +
  stat_boxplot(geom = "errorbar", width = 0.6, color = "black") +  
  geom_boxplot(fill = "#66C2A5", color = "black") +
  labs(title = "Shannon Diversity across Sites", x = "Site", y = "Shannon Diversity") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_line(linewidth = 0.2, linetype = "dashed", color = "gray"),
    plot.margin = margin(10, 10, 10, 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.3)
  ) +
  scale_y_continuous(limits = c(0, max(sample.data$shannon, na.rm = TRUE) * 1.2), 
                     breaks = seq(0, max(sample.data$shannon, na.rm = TRUE) * 1.2, by = 0.5))

# Arrange the plots side by side and save
final_plot <- grid.arrange(richness_plot, shannon_plot, ncol = 2)
ggsave(filename = "~/Documents/plot/alpha_diversity_by_Site.png", plot = final_plot,
       width = 12, height = 6, dpi = 300, units = "in")


# Plot for richness by area
richness_area_plot <- ggplot(sample.data, aes(x = Area, y = richness)) +
  stat_boxplot(geom = "errorbar", width = 0.3, color = "black") +  
  geom_boxplot(fill = "#8DA0CB", color = "black") +
  labs(title = "16S Richness across Areas", x = "Area", y = "16S Richness (rarefied)") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_line(linewidth = 0.2, linetype = "dashed", color = "gray"),
    plot.margin = margin(10, 10, 10, 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.3)
  ) +
  scale_y_continuous(limits = c(0, max(sample.data$richness, na.rm = TRUE) * 1.2), 
                     breaks = seq(0, max(sample.data$richness, na.rm = TRUE) * 1.2, by = 2))

# Plot for Shannon index by area
shannon_area_plot <- ggplot(sample.data, aes(x = Area, y = shannon)) +
  stat_boxplot(geom = "errorbar", width = 0.3, color = "black") +  
  geom_boxplot(fill = "#8DA0CB", color = "black") +
  labs(title = "Shannon Diversity across Areas", x = "Area", y = "Shannon Diversity") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_line(linewidth = 0.2, linetype = "dashed", color = "gray"),
    plot.margin = margin(10, 10, 10, 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.3)  
  ) +
  scale_y_continuous(limits = c(0, max(sample.data$shannon, na.rm = TRUE) * 1.2), 
                     breaks = seq(0, max(sample.data$shannon, na.rm = TRUE) * 1.2, by = 0.5))

# Arrange the plots side by side and save
final_plot_area <- grid.arrange(richness_area_plot, shannon_area_plot, ncol = 2)
ggsave(filename = "~/Documents/plot/alpha_diversity_plot_by_area.png", plot = final_plot_area,
       width = 12, height = 6, dpi = 300, units = "in")


##check the significance of difference between Areas
# Visual check for normality using histograms and Q-Q plots
par(mfrow = c(2, 2))

# Histogram and Q-Q plot for Shannon Diversity
hist(sample.data$shannon, main = "Histogram of Shannon Diversity", xlab = "Shannon Diversity", col = "lightblue")
qqnorm(sample.data$shannon, main = "Q-Q Plot of Shannon Diversity")
qqline(sample.data$shannon, col = "red")

# Histogram and Q-Q plot for Species Richness
hist(sample.data$richness, main = "Histogram of Species Richness", xlab = "Species Richness", col = "lightgreen")
qqnorm(sample.data$richness, main = "Q-Q Plot of Species Richness")
qqline(sample.data$richness, col = "red")

# Shapiro-Wilk test for normality
shapiro_shannon <- shapiro.test(sample.data$shannon)
shapiro_richness <- shapiro.test(sample.data$richness)

# Print the results of the Shapiro-Wilk tests
print(shapiro_shannon)
print(shapiro_richness) #the p-value is less than 0.05, the data is not normally distributed, appropriate to use non-parametric Kruskal-Wallis test.


# Kruskal-Wallis test for Shannon Diversity across Areas
kruskal_shannon <- kruskal.test(shannon ~ Area, data = sample.data)
print(kruskal_shannon) 

# Kruskal-Wallis test for Species Richness across Areas
kruskal_richness <- kruskal.test(richness ~ Area, data = sample.data)
print(kruskal_richness)



# By Sample_type


# Group by Sample_Type and calculate the mean diversity indices
mean_diversity_by_sample_type <- sample.data %>%
  group_by(Sample_type) %>%
  summarise(Richness = mean(richness, na.rm = TRUE))

# Display the results
print(mean_diversity_by_sample_type)

# save the results to CSV file
write.csv(mean_diversity_by_sample_type, "~/Documents/Diversity_by_sample_type.csv", row.names = FALSE)


# Boxplot for Species Richness by Sample Type
richness_sample_type_plot <- ggplot(sample.data, aes(x = Sample_type, y = richness, fill = Sample_type)) +
  stat_boxplot(geom = "errorbar", width = 0.1, color = "black") +  
  geom_boxplot(color = "black", outlier.shape = NA, width = 0.4) +  
  geom_jitter(position = position_jitter(width = 0.15), size = 1, color = "gray40") +  
  scale_fill_manual(values = c("old_leaf" = "#FDD835", "whole plant" = "#80CBC4", "young_leaf" = "#66BB6A")) +  
  labs(title = "16S Richness by Sample Type", x = "Sample Type", y = "16S Richness (rarefied)") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 12),  
    axis.text.y = element_text(size = 12),
    legend.position = "right",  
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major.y = element_line(linewidth = 0.3, linetype = "dotted", color = "gray"),
    panel.grid.minor.y = element_blank(),  
    plot.margin = margin(10, 10, 10, 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.3)  
  ) +
  scale_y_continuous(limits = c(0, max(sample.data$richness, na.rm = TRUE) * 1.1),  
                     breaks = seq(0, max(sample.data$richness, na.rm = TRUE) * 1.1, by = 2))

# Display the plot and save
print(richness_sample_type_plot)
ggsave(filename = "~/Documents/plot/richness_by_sample_type.png", plot = richness_sample_type_plot,
       width = 8, height = 6, dpi = 300, units = "in")


# Kruskal-Wallis test for Species Richness across Sample Types
kruskal_richness <- kruskal.test(richness ~ Sample_Type, data = sample.data)
print(kruskal_richness)



# By Dislodgment_method



# Group by Dislodgment_method and calculate the mean richness
mean_diversity_by_dislodgment <- sample.data %>%
  group_by(Dislodgment_method) %>%
  summarise(Richness = mean(richness, na.rm = TRUE))

# Display the results
print(mean_diversity_by_dislodgment)

# save the results to CSV
write.csv(mean_diversity_by_dislodgment, "~/Documents/mean_diversity_by_dislodgment.csv", row.names = FALSE)


# Create the plot 
richness_violin_boxplot <- ggplot(sample.data, aes(x = Dislodgment_method, y = richness, fill = Dislodgment_method)) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.8, linewidth = 1.3) +  
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA, alpha = 0.7,linewidth = 0.5) +  
  labs(title = "16S Richness by Dislodgment Method", x = "Dislodgment Method", y = "16S Richness (rarefied)") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major.y = element_line(linewidth = 0.3, linetype = "dotted", color = "gray"),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.3)  
  ) +
  scale_y_continuous(limits = c(0, max(sample.data$richness, na.rm = TRUE) * 1.1), 
                     breaks = seq(0, max(sample.data$richness, na.rm = TRUE) * 1.1, by = 2)) +
  scale_fill_manual(values = c("swabbing" = "#377EB8", "whole_leaf" = "#4DAF4A"))

# Display the plot and save
print(richness_violin_boxplot)
ggsave(filename = "~/Documents/plot/richness_violin_boxplot.png", plot = richness_violin_boxplot,
       width = 10, height = 6, dpi = 300, units = "in")

# Perform Kruskal-Wallis test for Species Richness across Dislodgment Methods
kruskal_test_result <- kruskal.test(richness ~ Dislodgment_method, data = sample.data)

# Display the results
print(kruskal_test_result)



# By Season


# Group by Season and calculate the mean richness
mean_richness_by_season <- sample.data %>%
  group_by(Season) %>%
  summarise(Richness = mean(richness, na.rm = TRUE))

# Display the results
print(mean_richness_by_season)


# Create the boxplot with jitter for richness by Season
richness_by_season_plot <- ggplot(sample.data, aes(x = Season, y = richness, fill = Season))+
  stat_boxplot(geom = "errorbar", width = 0.1, color = "black") + 
  geom_boxplot(color = "black", outlier.shape = NA, width = 0.4) +  
  geom_jitter(position = position_jitter(width = 0.15), size = 1.5, color = "gray40", alpha = 0.6) +  
  scale_fill_manual(values = c("Summer" = "#FDD835", "Autumn" = "#FFA726", "Spring" = "#66BB6A")) +  
  labs(title = "16S Richness by Season", x = "Season", y = "16S Richness (rarefied)") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major.y = element_line(linewidth = 0.3, linetype = "dotted", color = "gray"),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.3)
  ) +
  scale_y_continuous(limits = c(0, max(sample.data$richness, na.rm = TRUE) * 1.1), 
                     breaks = seq(0, max(sample.data$richness, na.rm = TRUE) * 1.1, by = 2))  

# Display the plot and save
print(richness_by_season_plot)
ggsave(filename = "~/Documents/plot/richness_by_season.png", plot = richness_by_season_plot,
       width = 10, height = 6, dpi = 300, units = "in")


# Perform Kruskal-Wallis test for Species Richness across Seasons
kruskal_test_richness <- kruskal.test(richness ~ Season, data = sample.data)


# Display the results
print(kruskal_test_richness)




#Across Julian day


#Group by Julian day and calculate the mean richness
mean_diversity_by_julian <- sample.data %>%
  group_by(Julian_day) %>%
  summarise(
    Richness = mean(richness, na.rm = TRUE))

# Check if the summarization worked correctly
print(mean_diversity_by_julian)

# Plot mean richness across Julian days 
richness_julian_plot <- ggplot(mean_diversity_by_julian, aes(x = Julian_day, y = mean_richness)) +
  geom_line(color = "darkgreen", linewidth = 1.2) +  
  geom_point(size = 3, color = "darkgreen") + 
  labs(title = "16S Richness across Julian Days", x = "Julian Day", y = "16S Richness (rarefied)") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none",  
    panel.grid.major.y = element_line(linewidth = 0.6, linetype = "dotted", color = "gray"),  
    panel.grid.major.x = element_line(linewidth = 0.6, linetype = "dotted", color = "gray"),  
    panel.grid.minor = element_blank(),  
    plot.margin = margin(10, 10, 10, 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)  
  ) +
  scale_y_continuous(limits = c(0, max(mean_diversity_by_julian$mean_richness) * 1.2),  
                     breaks = seq(0, max(mean_diversity_by_julian$mean_richness) * 1.2, by = 2)) +  
  scale_x_continuous(breaks = seq(min(mean_diversity_by_julian$Julian_day), max(mean_diversity_by_julian$Julian_day), by = 10)) +  
  annotate("text", x = 250, y = 17.5, label = "16S Richness", color = "darkgreen", size = 4, fontface = "bold")  

# Display the plot and save
print(richness_julian_plot)
ggsave(filename = "~/Documents/plot/Richness_julian_days.png", plot = richness_julian_plot,
       width = 10, height = 6, dpi = 300, units = "in")