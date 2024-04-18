library(tidyverse)
library(patchwork)  # for 2x2 plot layout
library(Rtsne)

my_theme <- theme_minimal() +
  theme(legend.position = c(0.95, 0.95), # Position the legend inside the plot at the top-right corner
        legend.justification = c("right", "top"), # Anchor point for legend position
        legend.box.just = "right",
        legend.background = element_blank(),  # Remove legend background
        legend.box.background = element_rect(colour = "grey", fill = "white"),  # Add rectangle around the legend without fill
        legend.box.margin = margin(5, 5, 5, 5), # Adjust legend box margin to move closer to corner
        legend.margin = margin(0, 0, 0, 0))

dataset <- read_csv("data/Diffusion-all_knots.zip")
features <- dataset %>%
  select(-ID, -Sequence, -Label, -Tool)

# Factor the 'Tool' column with the specified order
df_tsne$Tool <- factor(df_tsne$Tool, levels = desired_order)

# plot A ------------------------------------------------------------------

# Assuming 'features' is a dataframe or matrix with your data
# Convert features to a matrix if it's not already
features_matrix <- as.matrix(features)

# Perform t-SNE
set.seed(123) # For reproducibility
tsne_result <- Rtsne(features_matrix, dims = 2, perplexity = 30, verbose = TRUE)

# Create a dataframe from the t-SNE result
df_tsne <- as.data.frame(tsne_result$Y)
colnames(df_tsne) <- c("comp-1", "comp-2")

# Assuming 'dataset' is your data frame containing the 'Tool' column
# Map the replacement as per your requirement
desired_order <- c("Evodiff", "RFdiffusion + MPNN", "Real")
df_tsne$Tool <- factor(map_chr(dataset$Tool, ~ str_replace(.x, "RFdiffusion", "RFdiffusion + MPNN")), levels = desired_order)

# Plot using ggplot
plotA <- ggplot(df_tsne, aes(x = `comp-1`, y = `comp-2`, color = Tool)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("#5f57db", "#db5f57", "#57db5f"), name = NULL) +
  labs(title = "T-SNE projection of ProtBert-BFD embeddings") +
  my_theme

plotA

# plot B ------------------------------------------------------------------

plot_tsne <- function(features, labels, title = "T-SNE projection of embeddings") {
  set.seed(123)  # For reproducibility
  tsne_result <- Rtsne(features, dims = 2, perplexity = 30, verbose = TRUE)
  
  df_tsne <- as.data.frame(tsne_result$Y)
  colnames(df_tsne) <- c("comp-1", "comp-2")
  
  # Assuming labels are binary, convert to factor with levels 'knotted' and 'unknotted'
  df_tsne$Label <- factor(ifelse(labels == 1, 'knotted', 'unknotted'), levels = c('knotted', 'unknotted'))
  
  p <- ggplot(df_tsne, aes(x = `comp-1`, y = `comp-2`, color = Label)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c('knotted' = 'blue', 'unknotted' = '#ff7f00'), name = NULL) +
    labs(title = title) +
    my_theme
  
  return(p)
}

# Subset the dataset for 'Evodiff'
e <- dataset %>% filter(Tool == 'Evodiff')
features <- e %>% select(-ID, -Sequence, -Label, -Tool) %>% as.matrix()
labels <- e$Label

plotB <- plot_tsne(features, labels, 'T-SNE projection of Evodiff ProtBert-BFD embeddings')

# plot C ------------------------------------------------------------------

# Subset the dataset for 'Evodiff'
e <- dataset %>% filter(Tool == 'RFdiffusion')
features <- e %>% select(-ID, -Sequence, -Label, -Tool) %>% as.matrix()
labels <- e$Label

plotC <- plot_tsne(features, labels, 'T-SNE projection of RFdiffusion + MPNN ProtBert-BFD embeddings')

# plot D ------------------------------------------------------------------

# Subset the dataset for 'Evodiff'
e <- dataset %>% filter(Tool == 'Real')
features <- e %>% select(-ID, -Sequence, -Label, -Tool) %>% as.matrix()
labels <- e$Label

plotD <- plot_tsne(features, labels, 'T-SNE projection of real protein ProtBert-BFD embeddings')

# Combined plot -----------------------------------------------------------

combined_plot <- plotA + plotB + plotC + plotD +
  plot_layout(ncol = 2, nrow = 2) +  # Arrange in a 2x2 grid
  plot_annotation(tag_levels = 'A')  # Add labels (A-D) to each plot

ggsave("figures/figure4.pdf", combined_plot, width = 12, height = 12)
