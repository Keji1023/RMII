#Package Load ----------
# Loading the necessary package for the analysis
install.packages(c("tidyverse", "patchwork", "viridis", "svglite", "vegan", "ggpubr", "randomForest", "caret", "pROC"))
library(tidyverse)
library(svglite)
library(vegan)
library(ggthemes)
library(patchwork)
library(viridis)
library(ggplot2)
library(ggpubr)
library(randomForest)
library(caret)
library(pROC)

# Data File Load -------
# Loading the files for the aggregated taxonomic abundance data and sample descriptions
microbiome_df <- read.csv("https://raw.githubusercontent.com/Keji1023/RMII/refs/heads/main/Microbiome_data.csv")
metadata_df <- read.csv("https://raw.githubusercontent.com/Keji1023/RMII/refs/heads/main/Metadata_Microbiome.csv")

# Data Transformation --------
# Inspecting column headings
colnames(microbiome_df)
colnames(metadata_df)

# Filter to only include known class-level taxa; inconsistent taxa classification levels
# in data. Select most occurring taxa level: class
class_level_taxa <- c(
  "Bacilli", "Actinobacteria", "Clostridia", "Verrucomicrobiae", "Bacteroidia",
  "Coriobacteriia", "Gammaproteobacteria", "Alphaproteobacteria", "Desulfitobacteriia",
  "Cyanobacteriia", "Negativicutes", "Gemmatimonadetes", "Acidobacteriae", 
  "Fusobacteriia", "Anaerolineae", "Limnochordia", "Elusimicrobia"
)

microbe_class_counts <- microbiome_df %>%
  select(all_of(class_level_taxa)) %>%
  as.data.frame()

rm(microbiome_df)

# Make sure all values are numeric
microbe_class_counts[] <- lapply(microbe_class_counts, as.numeric)

# Check for sparsity
colSums(microbe_class_counts == 0)  # How many samples have zero for each class

# Filter out columns with all zeros; there were none
microbe_class_counts <- microbe_class_counts[, colSums(microbe_class_counts) > 0]

#Normalize counts to relative abundance (TSS)
microbe_rel_abund <- microbe_class_counts / rowSums(microbe_class_counts)

# Keep only relevant metadata
metadata_df <- metadata_df %>%
  select(gender, study, location)

# Combining dataframes & splitting by intestinal location
combined_df_list <- bind_cols(metadata_df, microbe_class_counts) %>%
  group_split(location) %>%
  setNames(unique(bind_cols(metadata_df, microbe_class_counts)$location))
# Split
combined_df_list$L_Intestine  
combined_df_list$Sm_Intestine

# Relative Abundance Combined
combined_rel_df_list <- bind_cols(metadata_df, microbe_rel_abund) %>%
  group_split(location) %>%
  setNames(unique(bind_cols(metadata_df, microbe_rel_abund)$location))
# Split
combined_rel_df_list$L_Intestine  
combined_rel_df_list$Sm_Intestine

# Raw Counts & Relative Abundance Subsetting
#combined_df_list$L_Intestine[, 4:ncol(combined_df_list$L_Intestine)]
#combined_df_list$Sm_Intestine[, 4:ncol(combined_df_list$Sm_Intestine)]

#combined_rel_df_list$L_Intestine[, 4:ncol(combined_rel_df_list$L_Intestine)]
#combined_rel_df_list$Sm_Intestine[, 4:ncol(combined_rel_df_list$Sm_Intestine)]

# Alpha Diversity Analysis--------
# Combining Diversity Metrics with Metadata
Ldiversity_df <- data.frame(
  Sample         = rownames(combined_df_list$L_Intestine[, 4:ncol(combined_df_list$L_Intestine)]), 
  Shannon        = diversity(combined_df_list$L_Intestine[, 4:ncol(combined_df_list$L_Intestine)], index = "shannon"),     # accounts for both richness and evenness
  Simpson        = diversity(combined_df_list$L_Intestine[, 4:ncol(combined_df_list$L_Intestine)], index = "shannon"),     # emphasizes abundant taxa
  Richness       = specnumber(combined_df_list$L_Intestine[, 4:ncol(combined_df_list$L_Intestine)]), # simply the number of taxa present per sample
  Gender         = combined_df_list$L_Intestine$gender,
  Study          = combined_df_list$L_Intestine$study
)

Smdiversity_df <- data.frame(
  Sample         = rownames(combined_df_list$Sm_Intestine[, 4:ncol(combined_df_list$Sm_Intestine)]), 
  Shannon        = diversity(combined_df_list$Sm_Intestine[, 4:ncol(combined_df_list$Sm_Intestine)], index = "shannon"),     # accounts for both richness and evenness
  Simpson        = diversity(combined_df_list$Sm_Intestine[, 4:ncol(combined_df_list$Sm_Intestine)], index = "shannon"),     # emphasizes abundant taxa
  Richness       = specnumber(combined_df_list$Sm_Intestine[, 4:ncol(combined_df_list$Sm_Intestine)]), # simply the number of taxa present per sample
  Gender         = combined_df_list$Sm_Intestine$gender,
  Study          = combined_df_list$Sm_Intestine$study
)

# Check the combined data
head(Ldiversity_df)
head(Smdiversity_df)

# Shapiro-Wilk Test (test for normality) - reject H0, data is non-normal
# While alpha diversity metrics often violate assumptions of normality 
# and homoscedasticity, I am using made up data so I had to confirm
shapiro.test(Ldiversity_df$Shannon)
shapiro.test(Ldiversity_df$Simpson)
shapiro.test(Ldiversity_df$Richness)

shapiro.test(Smdiversity_df$Shannon)
shapiro.test(Smdiversity_df$Simpson)
shapiro.test(Smdiversity_df$Richness)

# Kruskal-Wallis test (non-parametric test for multi-group statistically 
# significant difference) - non significant differences
# Wilcoxn test for two groups
kruskal.test(Shannon ~ Study, data = Ldiversity_df)
kruskal.test(Simpson ~ Study, data = Ldiversity_df)
kruskal.test(Richness ~ Study, data = Ldiversity_df)

kruskal.test(Shannon ~ Study, data = Smdiversity_df)
kruskal.test(Simpson ~ Study, data = Smdiversity_df)
kruskal.test(Richness ~ Study, data = Smdiversity_df)

# Visualization 
# Boxplot of Shannon diversity by Study (exposure/treatment group)
ggplot(Ldiversity_df, aes(x = Study, y = Shannon)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  labs(title = "Alpha Diversity (Shannon Index) by Treatment Group",
       x = "Treatment Group", 
       y = "Shannon Diversity Index") +
  theme_bw()

# Simpson or Richness: can skip
ggplot(Smdiversity_df, aes(x = Gender, y = Simpson)) +
  geom_boxplot(fill = "lightgreen", color = "darkgreen") +
  labs(title = "Simpson by Sex",
       x = "Sex", 
       y = "Number of Taxa (Richness)") +
  theme_bw()

# Shannon comparison across treatment and sex
max_y <- max(Ldiversity_df$Shannon, Smdiversity_df$Shannon, na.rm = TRUE)

plot2 <- ggplot(Ldiversity_df, aes(x = Study, y = Shannon, fill = Study)) +
  geom_boxplot() +
  facet_wrap(~Gender, labeller = as_labeller(c("F" = "Female", "M" = "Male"))) +
  labs(title = "Large Intestine Alpha Diversity (Class Level)",
       x = NULL,
       y = "Shannon Diversity",
       fill = "Group") +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 14, hjust = 0.5), 
        plot.caption = element_text(size = 7.4, hjust = 0),
        axis.title.y = element_text(size = 10)) + 
  scale_fill_viridis(discrete = TRUE, option = "B") 
  coord_cartesian(ylim = c(0, max_y))  # enforce same y-scale

plot1 <- ggplot(Smdiversity_df, aes(x = Study, y = Shannon, fill = Study)) +
  geom_boxplot() +
  facet_wrap(~Gender, labeller = as_labeller(c("F" = "Female", "M" = "Male"))) +
  labs(title = "Small Intestine Alpha Diversity (Class Level)",
       x = NULL,
       y = "Shannon Diversity",
       fill = "Group",
       caption = "Shannon diversity index (class level) across treatment groups and stratified by sex. Groups include: 'Con_Acute' (acute control), 'Con_SC' (subchronic control), 'Dose_Acute' \n(acute PAH exposure), and 'Dose_SC' (subchronic PAH exposure)."
  ) +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 14, hjust = 0.5), 
    plot.caption = element_text(size = 7.4, hjust = 0),
    axis.title.y = element_text(size = 10)) + 
  scale_fill_viridis(discrete = TRUE, option = "B") +
  guides(fill = "none") + # Remove legend from this plot
  coord_cartesian(ylim = c(0, max_y))  # match y-scale

# Combine side-by-side
plot1 + plot2
ggsave(filename = "alpha.png", width = 10, height = 6, units = "in", dpi = 300)
rm(max_y)

# Beta Diversity and Ordination--------
# Compute Bray–Curtis dissimilarity matrix
Lbray_dist <- vegdist(combined_rel_df_list$L_Intestine[, 4:ncol(combined_rel_df_list$L_Intestine)], method = "bray")
Smbray_dist <- vegdist(combined_rel_df_list$Sm_Intestine[, 4:ncol(combined_rel_df_list$Sm_Intestine)], method = "bray")
# PCoA
Lpcoa_result <- cmdscale(Lbray_dist, eig = TRUE, k = 2)  # k = number of dimensions to retain
Smpcoa_result <- cmdscale(Smbray_dist, eig = TRUE, k = 2)  # k = number of dimensions to retain
#Get coordinates and merge with metadata
Lpcoa_points <- as.data.frame(Lpcoa_result$points)
colnames(Lpcoa_points) <- c("PCoA1", "PCoA2")
Lpcoa_points$study <- combined_rel_df_list$L_Intestine$study
Lpcoa_points$sex <- combined_rel_df_list$L_Intestine$gender

Smpcoa_points <- as.data.frame(Smpcoa_result$points)
colnames(Smpcoa_points) <- c("PCoA1", "PCoA2")
Smpcoa_points$study <- combined_rel_df_list$Sm_Intestine$study
Smpcoa_points$sex <- combined_rel_df_list$Sm_Intestine$gender

# Get percent variance explained
Lvariance_explained <- round(100 * Lpcoa_result$eig / sum(Lpcoa_result$eig), 2)
Smvariance_explained <- round(100 * Smpcoa_result$eig / sum(Smpcoa_result$eig), 2)
# Adding PERMANOVA to plot
permanova <- adonis2(Lbray_dist ~ study + gender, data = combined_rel_df_list$L_Intestine[, 1:2]) # additive
r2 <- round(permanova$R2[1], 3)
pval <- signif(permanova$`Pr(>F)`[1], 3)

# Convert sex to factor to ensure proper interpretation
Lpcoa_points$sex <- as.factor(Lpcoa_points$sex)
Smpcoa_points$sex <- as.factor(Smpcoa_points$sex)

# PCoA plot: can remove gender grouping
plot2 <- ggplot(Lpcoa_points, aes(x = PCoA1, y = PCoA2, color = study, shape = sex)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = study), type = "t", linetype = "dashed", linewidth = 1) +
  labs(
    title = "PCoA of Large Intestine Microbial Community",
    subtitle = paste0("PERMANOVA: R² = ", r2, ", p = ", pval),
    x = paste0("PCoA1 (", Lvariance_explained[1], "%)"),
    y = paste0("PCoA2 (", Lvariance_explained[2], "%)"),
    shape = "Sex",
    color = "Group"
  ) +
  theme_bw() + 
  scale_color_viridis(discrete = TRUE, option = "A") 

# Adding PERMANOVA to plot
permanova <- adonis2(Smbray_dist ~ study + gender, data = combined_rel_df_list$Sm_Intestine[, 1:2]) # additive
# Extract key values
r2 <- round(permanova$R2[1], 3)
pval <- signif(permanova$`Pr(>F)`[1], 3)

plot1 <- ggplot(Smpcoa_points, aes(x = PCoA1, y = PCoA2, color = study, shape = sex)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = study), type = "t", linetype = "dashed", linewidth = 1) +
  labs(
    title = "PCoA of Small Intestine Microbial Community",
    subtitle = paste0("PERMANOVA: R² = ", r2, ", p = ", pval),
    x = paste0("PCoA1 (", Smvariance_explained[1], "%)"),
    y = paste0("PCoA2 (", Smvariance_explained[2], "%)"),
    shape = "Sex",
    color = "Group"
  ) +
  theme_bw() +
  scale_color_viridis(discrete = TRUE, option = "A") +
  guides(color = "none", shape = "none") # Remove legend from this plot

# Combine side-by-side
plot1 + plot2
ggsave(filename = "beta.png", width = 10, height = 6, units = "in", dpi = 300)

# Test whether groupings explain community differences
adonis2(Lbray_dist ~ study, data =  combined_rel_df_list$L_Intestine[, 1:2])
adonis2(Lbray_dist ~ study * gender, data = combined_rel_df_list$L_Intestine[, 1:2]) # interaction
adonis2(Smbray_dist ~ study, data =  combined_rel_df_list$Sm_Intestine[, 1:2])
adonis2(Smbray_dist ~ study * gender, data = combined_rel_df_list$Sm_Intestine[, 1:2]) # interaction

# Differential Abundance Analysis--------
# Group-wise comparison to identify differentially abundant classes 
# Using Kruskal Wallis test
combined_rel_df_list$L_Intestine[, 4:ncol(combined_rel_df_list$L_Intestine)]
combined_rel_df_list$Sm_Intestine[, 4:ncol(combined_rel_df_list$Sm_Intestine)]

Lkruskal_results <- map_dfr(names(combined_rel_df_list$L_Intestine[, 4:ncol(combined_rel_df_list$L_Intestine)]), function(class_name) { # iterates over each class name, applying the function defined within
  test <- kruskal.test(
    combined_rel_df_list$L_Intestine[[class_name]] ~ combined_rel_df_list$L_Intestine$study # compare between study groups
  )
  data.frame( # Combine all iterations into a single dataframe
    Class = class_name,
    p.value = test$p.value
  )
}) %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%  # Benjamini–Hochberg (FDR) correction
  arrange(p.adj) #help identify the most significant easily

Smkruskal_results <- map_dfr(names(combined_rel_df_list$Sm_Intestine[, 4:ncol(combined_rel_df_list$L_Intestine)]), function(class_name) { # iterates over each class name, applying the function defined within
  test <- kruskal.test(
    combined_rel_df_list$Sm_Intestine[[class_name]] ~ combined_rel_df_list$Sm_Intestine$study # compare between study groups
  )
  data.frame( # Combine all iterations into a single dataframe
    Class = class_name,
    p.value = test$p.value
  )
}) %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%  # Benjamini–Hochberg (FDR) correction
  arrange(p.adj) #help identify the most significant easily

rm(permanova)
rm(pval)
rm(r2)

# Visualize top differential abundant classes
# Not used since there are no differential abundance classes that pass the FDR correction
# Pick top 6 from full-group analysis
top_classes <- Lkruskal_results %>%
  filter(p.adj < 0.05) %>%
  slice_head(n = 6) %>%
  pull(Class)

# Make for plotting
Lplot_df <- combined_rel_df_list$L_Intestine %>%
  select(study, gender, all_of(Lkruskal_results$Class[1:6])) %>%
  pivot_longer(-c(study, gender), names_to = "Class", values_to = "RelAbund")

Smplot_df <- combined_rel_df_list$Sm_Intestine %>%
  select(study, gender, all_of(Smkruskal_results$Class[1:6])) %>%
  pivot_longer(-c(study, gender), names_to = "Class", values_to = "RelAbund")

# Boxplots
ggplot(Lplot_df, aes(x = study, y = RelAbund, fill = study)) +
  geom_boxplot() +
  facet_wrap(~Class, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 16, hjust = 0.5)) +
  labs(title = "Top Differentially Abundant Microbial Classes in the Large Intestine",
       fill= "Group",
       x = NULL) +
  scale_fill_viridis(discrete = TRUE, option = "B") 
ggsave(filename = "TDL.svg", width = 10, height = 6, units = "in", dpi = 300)

ggplot(Smplot_df, aes(x = study, y = RelAbund, fill = study)) +
  geom_boxplot() +
  facet_wrap(~Class, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 16, hjust = 0.5)) +
  labs(title = "Top Differentially Abundant Microbial Classes in the Small Intestine",
       fill= "Group",
       x = NULL) +
  scale_fill_viridis(discrete = TRUE, option = "B") 
ggsave(filename = "TDS.svg", width = 10, height = 6, units = "in", dpi = 300)

# Diff Abund Side plots----------
# Using ANCOM-BC (bias corrected, designed to work with raw counts).
# Result: 8 taxa at the Class level passed the internal filtering criteria (non-zero abundance, prevalence, etc.), 
# and were used to model the sample-specific bias. However, ANCOM-BC2 models and corrects for compositional bias using many features,
# ideally >50. Having only 8 reduces statistical power and may impact reliability.

# Visualizations-------------------
# Relative Abundance Density Plots
# Reshape for plotting
Lrel_abund_long <- combined_rel_df_list$L_Intestine %>%
  pivot_longer(cols = all_of(class_level_taxa), names_to = "Class", values_to = "RelAbundance")
Smrel_abund_long <- combined_rel_df_list$Sm_Intestine %>%
  pivot_longer(cols = all_of(class_level_taxa), names_to = "Class", values_to = "RelAbundance")

# Trying replot by grouping others
# Get top 10 classes
Ltop12 <- Lrel_abund_long %>%
  group_by(Class) %>%
  summarise(mean_abundance = mean(RelAbundance)) %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 11) %>%
  pull(Class)

Smtop12 <- Smrel_abund_long %>%
  group_by(Class) %>%
  summarise(mean_abundance = mean(RelAbundance)) %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 11) %>%
  pull(Class)

# Recode all others as "Other"
Lrel_abund_long <- Lrel_abund_long %>%
  mutate(Class = ifelse(Class %in% Ltop12, Class, "Other"))
Smrel_abund_long <- Smrel_abund_long %>%
  mutate(Class = ifelse(Class %in% Smtop12, Class, "Other"))

# Plot (fewer categories, so no palette issue)
ggplot(Lrel_abund_long, aes(x = study, y = RelAbundance, fill = Class)) +
  geom_bar(stat = "identity", position = "fill") +
  #facet_wrap(~gender) +
  theme_minimal() +
  labs(title = "Class-Level Relative Abundances in the Large Intestine", y = "Relative Abundance", x = "Group") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Paired")
ggsave(filename = "RAL.svg", width = 8, height = 6, units = "in", dpi = 300)

ggplot(Smrel_abund_long, aes(x = study, y = RelAbundance, fill = Class)) +
  geom_bar(stat = "identity", position = "fill") +
  #facet_wrap(~gender) +
  theme_minimal() +
  labs(title = "Class-Level Relative Abundances in the Small Intestine", y = "Relative Abundance", x = "Group") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Paired")
ggsave(filename = "RAS.svg", width = 8, height = 6, units = "in", dpi = 300)

# HeatMap of Top Microbial Classes (REMOVE)
# Prepare matrix for heatmap
Lheatmap_matrix <- as.data.frame( combined_rel_df_list$L_Intestine[, 4:ncol(combined_rel_df_list$L_Intestine)][Ltop12])
rownames(Lheatmap_matrix) <- rownames(as.data.frame(combined_rel_df_list$L_Intestine[, 1:2]))
# Ensure the annotation row is ordered
annotation_data <- combined_rel_df_list$L_Intestine[, 1:2] %>%
  arrange(Study, Gender)

# Reorder heatmap matrix row names to match the new annotation order
Lheatmap_matrix <- Lheatmap_matrix[rownames(annotation_data), ]

# Plot
pheatmap(Lheatmap_matrix,
         scale = "row",
         annotation_row = annotation_data,
         show_rownames = FALSE,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100))

# Machine Learning for Microbe Classification------------
# Preparing input data and classification labels
Lml_df <- cbind(combined_rel_df_list$L_Intestine[, 4:ncol(combined_rel_df_list$L_Intestine)], study = as.factor(combined_rel_df_list$L_Intestine$study))
Smml_df <- cbind(combined_rel_df_list$Sm_Intestine[, 4:ncol(combined_rel_df_list$Sm_Intestine)], study = as.factor(combined_rel_df_list$Sm_Intestine$study))

# 2. Split into train/test
set.seed(123) # Setting a random seed to ensure reproducibility
Ltrain_idx <- createDataPartition(Lml_df$study, p = 0.8, list = FALSE) # Split the dataset into training (80%) and testing (20%), indices returned as vectors
set.seed(123) # Setting a random seed to ensure reproducibility
Smtrain_idx <- createDataPartition(Smml_df$study, p = 0.8, list = FALSE) # Split the dataset into training (80%) and testing (20%), indices returned as vectors
#train_data <- ml_df[train_idx, ] # Training data
#test_data  <- ml_df[-train_idx, ] # Test data

# Train Random Forest model with training data set
Lrf_model <- randomForest(study ~ ., data = Lml_df[Ltrain_idx, ], importance = TRUE, ntree = 500)
Smrf_model <- randomForest(study ~ ., data = Smml_df[Smtrain_idx, ], importance = TRUE, ntree = 500)

# Predict on test set
Lpreds <- predict(Lrf_model, newdata = Lml_df[-Ltrain_idx, ])
Smpreds <- predict(Smrf_model, newdata = Smml_df[-Smtrain_idx, ])

# Evaluate performance: check the accuracy of the model by comparing predictions with actual labels
Lconf_matrix <- confusionMatrix(Lpreds, Lml_df[-Ltrain_idx, ]$study)
print(Lconf_matrix)

Smconf_matrix <- confusionMatrix(Smpreds, Smml_df[-Ltrain_idx, ]$study)
print(Smconf_matrix)

# Convert to tidy data frame
# truelabel: actual classification, predictedlabel: predicted classification, 
#count: number of occurences for each combo
Lconf_df <- as.data.frame(as.table(Lconf_matrix$table))
colnames(Lconf_df) <- c("TrueLabel", "PredictedLabel", "Count") 

Smconf_df <- as.data.frame(as.table(Smconf_matrix$table))
colnames(Smconf_df) <- c("TrueLabel", "PredictedLabel", "Count") 

ggplot(Lconf_df, aes(x = PredictedLabel, y = TrueLabel, fill = Count)) +
  geom_tile(color = "white") +                    # Creates a heatmap-style grid
  geom_text(aes(label = Count), color = "black", size = 5) +  # Adds counts to tiles
  scale_fill_gradient(low = "white", high = "steelblue") +  # Color scale
  labs(
    title = "Confusion Matrix: Study Group Classification in the Large Intestine",
    x = "Predicted Label",
    y = "True Label"
  ) +
  theme_minimal()
ggsave(filename = "LCM.svg", width = 8, height = 6, units = "in", dpi = 300)

ggplot(Smconf_df, aes(x = PredictedLabel, y = TrueLabel, fill = Count)) +
  geom_tile(color = "white") +                    # Creates a heatmap-style grid
  geom_text(aes(label = Count), color = "black", size = 5) +  # Adds counts to tiles
  scale_fill_gradient(low = "white", high = "steelblue") +  # Color scale
  labs(
    title = "Confusion Matrix: Study Group Classification in the Small Intestine",
    x = "Predicted Label",
    y = "True Label"
  ) +
  theme_minimal()
ggsave(filename = "SmCM.svg", width = 8, height = 6, units = "in", dpi = 300)

# Generates a feature importance plot
png("Lvar_importance.png", width = 1000, height = 800)
varImpPlot(Lrf_model,
           type = 1,                        # Mean decrease in accuracy
           main = "Class Importance for Group Classification in the Large Intestine",
           pch = 16,                        # Use filled circles
           cex = 1.25,                      # Increase point size (default is 1)
           lwd = 1.5                         # Increase line width (default is 1)
)
dev.off()

png("Smvar_importance.png", width = 1000, height = 800)
varImpPlot(Smrf_model,
           type = 1,                        # Mean decrease in accuracy
           main = "Class Importance for Group Classification in the Small Intestine",
           pch = 16,                        # Use filled circles
           cex = 1.25,                      # Increase point size (default is 1)
           lwd = 1.5                         # Increase line width (default is 1)
)
dev.off()
