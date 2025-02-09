# 
##############  INSTRUCTIONS ##################
#The plot size should be W = 545, H = 350 
## DO THIS FOR GREENHOUSE DATA: 
##install.packages("FSA", dependencies = TRUE)  
#library(FSA)  # Load package
# Load package for Dunn's test
#dunnTest(Shannon ~ TreatmentGroup, data = alpha_div, method = "bonferroni")

library(phyloseq)
library(pheatmap)
library(biomformat)
library(ggplot2)
library(vegan)
library(dplyr)
library(tidyverse)
library(tibble)
library(pheatmap)
####### Import metadata #################
metadata_file <- "../data/metadata_its.txt"
metadata <- read.table(metadata_file, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE, quote = "")
colnames(metadata) <- c("SampleID", "Replication", "Species", "Diet", "Treatment", "Industry")

print(ncol(metadata))
print(colnames(metadata))

# Ensure row names are set correctly
rownames(metadata) <- metadata$SampleID
metadata$SampleID <- NULL  # Remove SampleID column after setting row names

# Trim whitespace in the Species column
metadata$Species <- trimws(metadata$Species)
colnames(metadata)

# Import phyloseq data
biom_file <- "../data/taxa.biom"
physeq <- import_biom(biom_file)
sample_metadata <- sample_data(metadata)
physeq_object <- merge_phyloseq(physeq, sample_metadata)

######## Rarefaction ##########
physeq_object <- rarefy_even_depth(physeq_object, rngseed=42)
colnames(tax_table(physeq_object)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

########  Relative abundance   ###########  
# Convert counts to relative abundance percentages
physeq_rel <- transform_sample_counts(physeq_object, function(x) x / sum(x) * 100)

# Extract abundance data as dataframe
data_df <- psmelt(physeq_rel)

# Merge metadata with taxonomic data
data_df <- merge(data_df, metadata, by.x = "Sample", by.y = "row.names")

head(tax_table(physeq_rel))

# Generate Stacked Bar Plot (Similar to Second Image)
ggplot(data_df, aes(x = sample_Species, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill") +  # Ensures y-axis is scaled to 100%
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +  # Converts 0.25 to 25%
  theme_minimal() +
  labs(x = "Phylum", y = "Relative Abundance (%)", fill = "Phylum") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
       panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
  # Add border
################################################################


#############################################################
#########Alpha Diversity (Shannon)###############
alpha_div <- estimate_richness(physeq_object, measures = c("Shannon", "Simpson", "Observed","Chao1" ))
#before merging set make sure SampleID is set correctly
alpha_div$SampleID <- rownames(alpha_div)
alpha_div2 <- merge(alpha_div, metadata, by.x = "SampleID", by.y = "row.names")
colnames(alpha_div2)


# Statistical Analysis
kruskal.test(Shannon ~ Species, data = alpha_div2)
kruskal_result2 <- kruskal.test(Shannon ~ Species, data = alpha_div2)

p_value <- signif(kruskal_result2$p.value, 3)


str(alpha_div2)


ggplot(data = alpha_div2, aes(x = Species, y = Shannon, fill = Species)) +
  stat_boxplot(geom = 'errorbar') +  # Error bars
  geom_boxplot(outlier.shape = NA) +  # Removes outlier points
  annotate("text", 
           x = 3.2, y = max(alpha_div2$Shannon) * 1.05,  # Position above the highest point
           label = paste("p =", p_value), 
           size = 4.5, fontface = "bold", hjust = 0) +  # Adds p-value annotation
  labs(title = "",
       x = "Insect Frass Type",
       y = "Shannon Index",
       fill = "Insect Frass Source") +  
  scale_fill_viridis_d() +
  theme(
    legend.position = "right",  # Moves legend to the right
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )
###########  Alpha Diversity (Chao1) ####################################
#alpha_div <- estimate_richness(physeq_object, measures = c("Chao1"))
alpha_div3 <- merge(alpha_div, metadata, by.x = "SampleID", by.y = "row.names")
alpha_div$SampleID <- rownames(alpha_div)
colnames(alpha_div3)
kruskal.test(Chao1 ~ Species, data = alpha_div3)
print(kruskal.test)

kruskal_result <- kruskal.test(Chao1 ~ Species, data = alpha_div3)

p_value <- signif(kruskal_result$p.value, 3) 

ggplot(data = alpha_div3, aes(x = Species, y = Chao1, fill = Species)) +
  stat_boxplot(geom = 'errorbar') +  # Error bars
  geom_boxplot(outlier.shape = NA) +  # Removes outlier points
  annotate("text", 
           x = 3.2, y = max(alpha_div3$Chao1) * 1.05,  # Position above the highest point
           label = paste("p =", p_value), 
           size = 4.5, fontface = "bold", hjust = 0) +  # Adds p-value annotation
  labs(title = "",
       x = "Insect Frass Type",
       y = "Chao1 Index",
       fill = "Insect Frass Source") +  
  scale_fill_viridis_d() +
  theme(
    legend.position = "right",  # Moves legend to the right
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

############ Alpha Diversity (Simpson) ##############
alpha_div <- estimate_richness(physeq_object, measures = c("Shannon", "Simpson", "Observed","Chao1" ))
#before merging set make sure SampleID is set correctly
alpha_div$SampleID <- rownames(alpha_div)
alpha_div4 <- merge(alpha_div, metadata, by.x = "SampleID", by.y = "row.names")
colnames(alpha_div4)
# Statistical Analysis
kruskal.test(Simpson ~ Species, data = alpha_div4)
kruskal_result4 <- kruskal.test(Simpson ~ Species, data = alpha_div4)

p_value <- signif(kruskal_result4$p.value, 3) 
#install.packages("FSA", dependencies = TRUE)  # DO THIS FOR GREENHOUSE DATA
#library(FSA)  # Load package
# Load package for Dunn's test
#dunnTest(Shannon ~ TreatmentGroup, data = alpha_div, method = "bonferroni")
str(alpha_div4)


ggplot(data = alpha_div4, aes(x = Species, y = Simpson, fill = Species)) +
  stat_boxplot(geom = 'errorbar') +  # Error bars
  geom_boxplot(outlier.shape = NA) +  # Removes outlier points
  annotate("text", 
           x = 3.2, y = max(alpha_div4$Simpson) * 1.05,  # Position above the highest point
           label = paste("p =", p_value), 
           size = 4.5, fontface = "bold", hjust = 0) +  # Adds p-value annotation
  labs(title = "",
       x = "Insect Frass Type",
       y = "Simpson Index",
       fill = "Insect Frass Source") +  
  scale_fill_viridis_d() +
  theme(
    legend.position = "right",  # Moves legend to the right
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )
############# Beta diversity (PCoA) #################
# Load necessary libraries
library(phyloseq)
library(ggplot2)

# Perform PCoA ordination using Bray-Curtis distance
ordination <- ordinate(physeq_object, method = "PCoA", distance = "bray")

# Convert ordination result to a dataframe
ordination_df <- as.data.frame(ordination$vectors)  # Extract ordination vectors
ordination_df$SampleID <- rownames(ordination_df)   # Create SampleID column

# Merge Sample Metadata (to include species info)
ordination_df <- merge(ordination_df, as.data.frame(sample_data(physeq_object)), 
                       by.x = "SampleID", by.y = "row.names")

# Verify that the "Species" column exists
colnames(ordination_df)


# Convert OTU table to distance matrix
bray_dist <- phyloseq::distance(physeq_object, method = "bray")

# Extract metadata
metadata_df <- as(sample_data(physeq_object), "data.frame")

sum(is.na(metadata_df$Species))  # Count NAs in Species column
nrow(metadata_df)  # Check number of rows in metadata
nrow(as.matrix(bray_dist))  # Check number of rows in the distance matrix
metadata_df <- metadata_df[rownames(metadata_df) %in% rownames(as.matrix(bray_dist)), ]
metadata_df <- metadata_df[!is.na(metadata_df$Species), ]



# Perform PERMANOVA (test for significant differences between groups)
permanova_result <- adonis2(bray_dist ~ Species, data = metadata_df, permutations = 999)
print(permanova_result)
permanova_result_1 <- adonis2(bray_dist ~ Treatment, data = metadata_df, permutations = 999)
print(permanova_result_1)
permanova_result_2 <- adonis2(bray_dist ~ Industry, data = metadata_df, permutations = 999)
print(permanova_result_2)
permanova_result_3 <- adonis2(bray_dist ~ Diet, data = metadata_df, permutations = 999)
print(permanova_result_3)
# Extract R² and p-value


# Define color and shape palettes
color_palette <- c("black", "red", "blue", "green", "purple", "orange", "brown") # Adjust for sample groups
shape_palette <- c(15, 16, 17, 18, 19, 21, 22, 23, 24, 25) # Unique symbols for species

# Verify the correct column name for sample groups
colnames(ordination_df)  # Check available columns

# Replace "Sample_Group" with the correct column name 
#### PCoA : Treatment ##########
r2_value_1 <- round(permanova_result_1$R2[1], 2)
p_value_1 <- signif(permanova_result_1$Pr[1], 3)

ggplot(ordination_df, aes(x = Axis.1, y = Axis.2, color = Treatment, shape = Treatment)) +
  geom_point(size = 3, alpha = 0.8) +  # Set point size and transparency
  scale_color_manual(values = color_palette) +  # Apply color palette
  scale_shape_manual(values = shape_palette) +  # Apply unique shapes for species
  labs(x = paste("Axis 1 [", round(ordination$values$Relative_eig[1] * 100, 1), "%]", sep = ""),
       y = paste("Axis 2 [", round(ordination$values$Relative_eig[2] * 100, 1), "%]", sep = ""),
       title = "Fungal Communities",
       subtitle = paste("PERMANOVA: R² =", r2_value_1,"; p =", p_value_1)) +  # Add R2 and p-value
  theme_minimal() +
  theme(legend.position = "right", 
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))  # Add border

#### PCoA : Industry ##########
r2_value_2 <- round(permanova_result_2$R2[1], 2)
p_value_2 <- signif(permanova_result_2$Pr[1], 3) 
ggplot(ordination_df, aes(x = Axis.1, y = Axis.2, color = Industry, shape = Industry)) +
  geom_point(size = 3, alpha = 0.8) +  # Set point size and transparency
  scale_color_manual(values = color_palette) +  # Apply color palette
  scale_shape_manual(values = shape_palette) +  # Apply unique shapes for industry
  labs(x = paste("Axis 1 [", round(ordination$values$Relative_eig[1] * 100, 1), "%]", sep = ""),
       y = paste("Axis 2 [", round(ordination$values$Relative_eig[2] * 100, 1), "%]", sep = ""),
       title = "Fungal Communities",
       subtitle = paste("PERMANOVA: R² =", r2_value_2, "; p =", p_value_2)) +  # Add R2 and p-value
  theme_minimal() +
  theme(legend.position = "right", 
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))  # Add border



############## PCoA : Species ##################
r2_value <- round(permanova_result$R2[1], 2)
p_value <- signif(permanova_result$Pr[1], 3)

ggplot(ordination_df, aes(x = Axis.1, y = Axis.2, color = Species, shape = Species)) +
  geom_point(size = 3, alpha = 0.8) +  # Set point size and transparency
  scale_color_manual(values = color_palette) +  # Apply color palette
  scale_shape_manual(values = shape_palette) +  # Apply unique shapes for species
  labs(x = paste("Axis 1 [", round(ordination$values$Relative_eig[1] * 100, 1), "%]", sep = ""),
       y = paste("Axis 2 [", round(ordination$values$Relative_eig[2] * 100, 1), "%]", sep = ""),
       title = "Fungal Communities",
       subtitle = paste("PERMANOVA: R² =", r2_value,"; p =", p_value)) +  # Add R2 and p-value
  theme_minimal()  +
  theme(legend.position = "right", 
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))  # Add border


############## PCoA : Diet ##################
# Extract R² and p-value
r2_value_3 <- round(permanova_result_3$R2[1], 2)
p_value_3 <- signif(permanova_result_3$Pr[1], 3)

# 
ggplot(ordination_df, aes(x = Axis.1, y = Axis.2, color = Diet, shape = Diet)) +
  geom_point(size = 3, alpha = 0.8) +  # Set point size and transparency
  scale_color_manual(values = color_palette) +  # Apply color palette
  scale_shape_manual(values = shape_palette) +  # Apply unique shapes for species
  labs(x = paste("Axis 1 [", round(ordination$values$Relative_eig[1] * 100, 1), "%]", sep = ""),
       y = paste("Axis 2 [", round(ordination$values$Relative_eig[2] * 100, 1), "%]", sep = ""),
       title = "PCoA Ordination of Microbial Communities",
       subtitle = paste("PERMANOVA: R² =", r2_value_3, "; p =", p_value_3)) +  # Add R2 and p-value
  theme_minimal()  +
  theme(legend.position = "right", 
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))  # Add border


  
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Heat Map @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


heatmap_data <- data_df %>%
  group_by(sample_Species, Phylum) %>%
  summarise(Relative_Abundance = sum(Abundance, na.rm = TRUE)) %>%
  ungroup()

# Convert into a matrix format for heatmap
heatmap_matrix <- heatmap_data %>%
  pivot_wider(names_from = sample_Species, values_from = Relative_Abundance, values_fill = list(Relative_Abundance = 0)) %>%
  column_to_rownames(var = "Phylum")

# Perform hierarchical clustering
row_dist <- dist(heatmap_matrix)  # Distance matrix for rows (Phylum)
row_clust <- hclust(row_dist)  # Clustering rows

col_dist <- dist(t(heatmap_matrix))  # Distance matrix for columns (Species)
col_clust <- hclust(col_dist)  # Clustering columns
#heatmap coloring
heat_colors <- colorRampPalette(c("red", "black", "green"))(100)  

# Define custom breaks ensuring black is at -0.5
breaks_list <- c(seq(-1.5, -0.51, length.out = 49), -0.5, seq(-0.49, 1.5, length.out = 50))

# Plot heatmap with adjusted color scale
pheatmap(
  mat = as.matrix(heatmap_matrix),  # Convert to matrix
  color = heat_colors,  # Apply custom color palette
  breaks = breaks_list,  # Set custom breaks to enforce black at -0.5
  cluster_rows = TRUE,  # Apply row clustering
  cluster_cols = TRUE,  # Apply column clustering
  show_colnames = TRUE,
  show_rownames = TRUE,
  border_color = "black",  # Add grid lines
  annotation_names_row = TRUE,
  annotation_names_col = TRUE,
  scale = "row"  # Standardize values for better visualization
)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@




#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@




