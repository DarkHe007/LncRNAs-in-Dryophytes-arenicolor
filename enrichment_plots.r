# Create top 25 GO terms data frame by filtering and sorting
top_GO <- significant_GO_df %>%
  arrange(p.adjust) %>%  # Sort by adjusted p-value (most significant first)
  dplyr::distinct(Description, .keep_all = T) %>%  # Remove duplicate GO term descriptions
  head(25)  # Keep only the top 25 most significant terms

# Create bubble plot for top GO term enrichments
ggplot(top_GO, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = -log10(p.adjust))) +
  geom_point(aes(size = Count)) +  # Point size represents number of genes in term
  scale_color_gradient(low = "blue", high = "red") +  # Color gradient: blue (less sig) to red (more sig)
  facet_grid(rows = vars(Module), scales = "free_y", space = "free_y") +  # Separate panels by module
  labs(title = "Top 25 GO Term Enrichments",
       x = "Gene Ratio",  # Proportion of genes in term vs module
       y = "GO Term",
       color = "-log10(adj.PValue)",  # Transform p-value for better visualization
       size = "Gene Count") + 
  theme(
    axis.text = element_text(size = 12, colour = "black", family = "serif"),
    strip.text.y = element_text(angle = 0, face = "bold", hjust = 0, family = "serif"),  # Module labels
    strip.background = element_rect(fill = "darkgray"),  # Background color for module labels
    strip.text = element_text(color = "black", face = "bold"),  # Module label text styling
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold", family = "serif"),  # Title styling
    axis.title = element_text(size = 15, face = "bold", family = "serif"),  # Axis title styling
    legend.title = element_text(size = 15, face = "bold", family = "serif"),  # Legend title styling
    legend.text = element_text(size = 12, family = "serif")  # Legend text styling
  ) 

# Save high-resolution GO enrichment plot
ggsave("Top30GO.png", width = 3840, height = 2160, units = "px", dpi = 600)

# Create top 25 KEGG pathways data frame
top_KEGG <- significant_KEGG_df %>%
  arrange(p.adjust) %>%  # Sort by adjusted p-value
  dplyr::distinct(Description, .keep_all = T) %>%  # Remove duplicate pathway descriptions
  head(25)  # Keep top 25 pathways

# Create bubble plot for top KEGG pathway enrichments
ggplot(top_KEGG, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = -log10(p.adjust))) +
  geom_point(aes(size = Count)) +  # Size = number of genes in pathway
  scale_color_gradient(low = "blue", high = "red") +  # Color by significance
  facet_grid(rows = vars(Module), scales = "free_y", space = "free_y") +  # Separate by module
  labs(title = "Top 25 KEGG Term Enrichments",
       x = "Gene Ratio",
       y = "KEGG Term",
       color = "-log10(adj.PValue)",
       size = "Gene Count") + 
  theme(
    axis.text = element_text(size = 12, colour = "black", family = "serif"),
    strip.text.y = element_text(angle = 0, face = "bold", hjust = 0, family = "serif"),
    strip.background = element_rect(fill = "darkgray"),
    strip.text = element_text(color = "black", face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold", family = "serif"),
    axis.title = element_text(size = 15, face = "bold", family = "serif"),
    legend.title = element_text(size = 15, face = "bold", family = "serif"),
    legend.text = element_text(size = 12, family = "serif")
  ) 

# Examine gene IDs in top KEGG results
head(top_KEGG$geneID)

# Load tidyr for data manipulation
library(tidyr)

# Separate concatenated gene IDs into multiple columns for easier analysis
df_separated <- top_KEGG %>%
  separate(geneID, into = paste0("V", 1:100), sep = "/", fill = "right", convert = TRUE)
# This splits the geneID column (format: "gene1/gene2/gene3") into multiple columns
# fill = "right" fills missing values with NA, convert = TRUE attempts to convert to appropriate data types

