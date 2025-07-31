#### CORRELATION ANALYSIS BETWEEN ALPHA DIVERSITY, MORPHOLOGICAL TRIATS AND ENVRIONMENTAL VARIABLES ####

# Load required libraries
library(GGally)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)

# Load sample metadata
sample.data_16S <- read.csv("~/Metadata_16S.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
sample.data_18S <- read.csv("~/Metadata_18S.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Calculate Pielou's Evenness 
sample.data_16S$Evenness_16S <- sample.data_16S$Shannon_16S / log(sample.data_16S$Richness_16S)
sample.data_18S$Evenness_18S <- sample.data_18S$Shannon_18S / log(sample.data_18S$Richness_18S)

# Function to compute correlation and adjusted p-value matrices 
compute_correlation_matrices <- function(data, method = "spearman") {
  cor_matrix <- cor(data, method = method, use = "pairwise.complete.obs")
  n <- ncol(data)
  p_matrix <- matrix(NA, n, n)
  colnames(p_matrix) <- rownames(p_matrix) <- colnames(data)
  diag(p_matrix) <- 0

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      test_result <- cor.test(data[, i], data[, j], method = method, exact = FALSE)
      p_matrix[i, j] <- p_matrix[j, i] <- test_result$p.value
    }
  }

  p_adjusted <- matrix(p.adjust(p_matrix, method = "BH"), n, n)
  colnames(p_adjusted) <- rownames(p_adjusted) <- colnames(data)
  return(list(cor_matrix = cor_matrix, p_matrix = p_matrix, p_adjusted = p_adjusted))
}

# Custom function for displaying only significant correlation 
my_custom_cor <- function(data, mapping, cor_results, ...) {
  x_name <- as.character(mapping$x)[2]
  y_name <- as.character(mapping$y)[2]
  corr <- cor_results$cor_matrix[x_name, y_name]
  adj_pvalue <- cor_results$p_adjusted[x_name, y_name]
  if (adj_pvalue < 0.05) {
    stars <- if (adj_pvalue < 0.001) "***" else if (adj_pvalue < 0.01) "**" else "*"
    ggplot() + 
      geom_text(aes(x = 0.5, y = 0.5), label = sprintf("%.2f\n%s", corr, stars), size = 6) +
      theme_void()
  } else {
    ggplot() + theme_void()
  }
}

# Custom function for only significant scatterplots with LOESS 
my_custom_scatter <- function(data, mapping, cor_results, ...) {
  x_name <- as.character(mapping$x)[2]
  y_name <- as.character(mapping$y)[2]
  adj_pvalue <- cor_results$p_adjusted[x_name, y_name]
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point(alpha = 0.5, size = 0.8, color = "darkblue")
  if (adj_pvalue < 0.05) {
    p <- p + geom_smooth(method = "loess", color = "red", se = TRUE, size = 0.7)
  }
  p + theme_minimal()
}

# Create correlation matrix plots 
create_correlation_matrix <- function(data, cor_results, title) {
  local_custom_cor <- function(data, mapping, ...) {
    my_custom_cor(data, mapping, cor_results, ...)
  }
  local_custom_scatter <- function(data, mapping, ...) {
    my_custom_scatter(data, mapping, cor_results, ...)
  }
  ggpairs(
    data,
    diag = list(continuous = wrap("densityDiag", alpha = 0.5, fill = "lightblue", color = "darkblue")),
    upper = list(continuous = local_custom_cor),
    lower = list(continuous = local_custom_scatter),
    title = title
  ) + theme_bw(base_size = 10)
}

# ====================================================
# Morphological traits vs Alpha diversity (16S & 18S)
# ====================================================

filtered_data_16S <- sample.data_16S %>% select(Richness_16S, Evenness_16S, Shoot_Length, Leaf_Width, No_Leaves) %>% na.omit()
filtered_data_18S <- sample.data_18S %>% select(Richness_18S, Evenness_18S, Shoot_Length, Leaf_Width, No_Leaves) %>% na.omit()

cor_results_morph_16S <- compute_correlation_matrices(filtered_data_16S)
cor_results_morph_18S <- compute_correlation_matrices(filtered_data_18S)

plot_morph_16S <- create_correlation_matrix(filtered_data_16S, cor_results_morph_16S, "Correlation Matrix of Morphological Traits and Alpha Diversity Index (16S)")
plot_morph_18S <- create_correlation_matrix(filtered_data_18S, cor_results_morph_18S, "Correlation Matrix of Morphological Traits and Alpha Diversity Index (18S)")

# Convert in gtable object 
plot_morph_16S_gtable <- ggmatrix_gtable(plot_morph_16S)
plot_morph_18S_gtable <- ggmatrix_gtable(plot_morph_18S)

# Combine plots
combined_morph_plot <- grid.arrange(
  plot_morph_16S_gtable,
  plot_morph_18S_gtable,
  ncol = 1,
  heights = c(1, 1),
  bottom = textGrob(
    "",
    gp = gpar(fontsize = 8),
    just = "left",
    x = 0.02,
    y = 0.5
  )
)

# print and save 
print(combined_morph_plot)

ggsave(filename = "Correlation_Matrix_Morpho_Alpha.png", plot = combined_morph_plot, width = 20, height = 10, dpi = 300, units = "in", bg = "white")

# =======================================================
# Environmental variables vs Alpha diversity (16S & 18S)
# =======================================================
 
env_data_16S <- sample.data_16S %>% select(Richness_16S, Evenness_16S, Julian_day, Latitude, Longitude, Sheltered_Index, Sea_Surface_Temperature, Salinity, Orbital_Velocity_Mean, Orbital_Velocity_Max, Bathymetry)
env_data_18S <- sample.data_18S %>% select(Richness_18S, Evenness_18S, Julian_day, Latitude, Longitude, Sheltered_Index, Sea_Surface_Temperature, Salinity, Orbital_Velocity_Mean, Orbital_Velocity_Max, Bathymetry)

cor_results_env_16S <- compute_correlation_matrices(env_data_16S)
cor_results_env_18S <- compute_correlation_matrices(env_data_18S)

plot_env_16S <- create_correlation_matrix(env_data_16S, cor_results_env_16S, "Correlation Matrix of Environmental Variables and Alpha Diversity Index (16S)")
plot_env_18S <- create_correlation_matrix(env_data_18S, cor_results_env_18S, "Correlation Matrix of Environmental Variables and Alpha Diversity Index (18S)")

ggsave("~/Correlation_Matrix_Env_16S.png", plot = plot_env_16S, width = 20, height = 10, dpi = 300)
ggsave("~/Correlation_Matrix_Env_18S.png", plot = plot_env_18S, width = 20, height = 10, dpi = 300)

# ================================================
# Environmental variables vs Morphological traits
# ================================================

morpho_env_data <- sample.data_16S %>% select(Shoot_Length, Leaf_Width, No_Leaves, Julian_day, Latitude, Longitude, Sheltered_Index, Sea_Surface_Temperature, Salinity, Orbital_Velocity_Mean, Bathymetry)
cor_results_morpho_env <- compute_correlation_matrices(morpho_env_data)
plot_morpho_env <- create_correlation_matrix(morpho_env_data, cor_results_morpho_env, "Correlation Matrix of Morphological Traits and Environmental Factors")

ggsave("~/Correlation_Matrix_Morpho_Envi.png", plot = plot_morpho_env, width = 20, height = 10, dpi = 300)
