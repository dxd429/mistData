setwd("./Data")
library(ggplot2)
library(pROC)

# General parameters
n_iterations <- 20
n_gene <- 10000
common_fpr <- seq(0, 1, length.out = 100)
common_recall <- seq(0, 1, length.out = 100)

# Create "figures" directory if it doesn't exist
if (!dir.exists("figures")) {
  dir.create("figures")
}

# File patterns for different plots (updated to ensure both 1-group and 2-group data are captured)
tdr_files <- list.files(pattern = "sim_tdr_results_.*\\.rda")
roc_files <- list.files(pattern = "sim_tpr_results_.*\\.rda")
pr_files <- list.files(pattern = "sim_precision_results_.*\\.rda")

# Function to determine the DM ratio and group scenario based on the filename
get_dm_ratio_and_group <- function(file) {
  # Extract the DM ratio from the filename (e.g., "0.05", "0.1", "0.15")
  dm_ratio <- as.numeric(sub(".*_(0\\.\\d+)_.*", "\\1", file))
  
  # Determine if it's a 1-group or 2-group scenario based on the filename
  scenario <- if (grepl("2group", file, ignore.case = TRUE)) "2-group" else "1-group"
  
  return(list(dm_ratio = dm_ratio, scenario = scenario))
}

# Function to create TDR plots
create_tdr_plot <- function(file) {
  load(file) # Load the RDA file
  
  # Get the DM ratio and group scenario from the filename
  info <- get_dm_ratio_and_group(file)
  dm_ratio <- info$dm_ratio
  scenario <- info$scenario
  
  # Define cutoffs based on the DM ratio
  cutoffs <- seq(100, n_gene * dm_ratio, by = 100)
  n_cutoffs <- length(cutoffs)
  
  # Convert TDR results to data frames for plotting
  tdr_df <- data.frame(
    Iteration = rep(1:n_iterations, times = 3 * n_cutoffs),
    Model = rep(rep(c('Bayesian', 'GAM', 'Polynomial'), each = n_iterations), each = n_cutoffs),
    Cutoff = rep(rep(cutoffs, each = n_iterations), times = 3),
    TDR = c(tdr_results$Bayesian, tdr_results$GAM, tdr_results$Linear)
  )
  
  # Generate PDF filename based on the scenario (1-group or 2-group)
  pdf_filename <- file.path("figures", sub("v[69]_", "", sub(".rda", paste0("_", scenario, "_tdr.pdf"), file)))
  
  tdr_df$Model <- factor(tdr_df$Model, levels = c('Bayesian', 'Polynomial', 'GAM')) 
  # Plot TDR results
  pdf(file = pdf_filename, width = 9)
  print(ggplot(tdr_df, aes(x = Cutoff, y = TDR, color = Model, group = Model)) +
          stat_summary(fun = mean, geom = "line", size = 2) +
          stat_summary(fun = mean, geom = "point", size = 1.2) +
          labs(title = paste('True Discovery Rate (TDR) Across Simulations (', scenario, ')', sep = ''),
               x = 'Cutoff', y = 'TDR', color = 'Model') +
          theme_classic() +
          theme(axis.text.x = element_text(size = 30, vjust = 0.85, color = "black"),
                axis.text.y = element_text(hjust = 1, size = 30, color = "black"),
                axis.title.x = element_text(size = 25, color = "black"),
                axis.title.y = element_text(size = 25, color = "black"),
                legend.position = "bottom",
                axis.line = element_line(colour = "black", size = 0.8), 
                axis.ticks = element_line(colour = "black", size = 0.8),
                axis.ticks.length = unit(0.3, "cm"))+
          scale_color_manual(values = c("Bayesian" = "#A9082C", "GAM" = "#4B61A8", "Polynomial" = "#FDB76D"),
                             breaks = c("Bayesian", "GAM", "Polynomial")) +
          ylim(0, 1))
  dev.off()
}

# Function to create ROC plots
create_roc_plot <- function(file) {
  load(file) # Load the RDA file
  
  # Get the DM ratio and group scenario from the filename
  info <- get_dm_ratio_and_group(file)
  scenario <- info$scenario
  
  # Calculate median TPR at each common FPR point across iterations
  median_tpr <- list(
    Bayesian = apply(interpolated_tpr$Bayesian, 2, mean),
    GAM = apply(interpolated_tpr$GAM, 2, mean),
    Polynomial = apply(interpolated_tpr$Linear, 2, mean)
  )
  
  # Combine the median ROC data into a single data frame for plotting
  median_roc_combined <- rbind(
    data.frame(FPR = common_fpr, TPR = median_tpr$Bayesian, Method = rep("Bayesian", length(common_fpr))),
    data.frame(FPR = common_fpr, TPR = median_tpr$GAM, Method = rep("GAM", length(common_fpr))),
    data.frame(FPR = common_fpr, TPR = median_tpr$Polynomial, Method = rep("Polynomial", length(common_fpr)))
  )
  
  # Generate PDF filename based on the scenario (1-group or 2-group)
  pdf_filename <- file.path("figures", sub("v[69]_", "", sub(".rda", paste0("_", scenario, "_roc.pdf"), file)))
  
  median_roc_combined$Method <- factor(median_roc_combined$Method, levels = c('Bayesian', 'Polynomial', 'GAM')) 
  # Plot ROC results
  pdf(file = pdf_filename, width = 9)
  print(ggplot(median_roc_combined, aes(x = FPR, y = TPR, color = Method)) +
          geom_line(size = 2) +
          geom_abline(linetype = "dashed", color = "gray") +
          labs(title = paste('ROC Curve (', scenario, ')', sep = ''),
               x = "1 - Specificity", y = "Sensitivity") +
          theme_classic() +
          ylim(0,1) +
          theme(axis.text.x = element_text(size = 30, vjust = 0.85, color = "black"),
                axis.text.y = element_text(hjust = 1, size = 30, color = "black"),
                axis.title.x = element_text(size = 25, color = "black"),
                axis.title.y = element_text(size = 25, color = "black"),
                legend.position = "bottom",
                axis.line = element_line(colour = "black", size = 0.8), 
                axis.ticks = element_line(colour = "black", size = 0.8),
                axis.ticks.length = unit(0.3, "cm"))+
          scale_color_manual(values = c("Bayesian" = "#A9082C", "GAM" = "#4B61A8", "Polynomial" = "#FDB76D"),
                             breaks = c("Bayesian", "GAM", "Polynomial")))
  dev.off()
}

# Function to create Precision-Recall plots
create_pr_plot <- function(file) {
  load(file) # Load the RDA file
  
  # Get the DM ratio and group scenario from the filename
  info <- get_dm_ratio_and_group(file)
  scenario <- info$scenario
  
  # Calculate median precision at each common recall point across iterations
  median_precision <- list(
    Bayesian = apply(interpolated_precision$Bayesian, 2, mean),
    GAM = apply(interpolated_precision$GAM, 2, mean),
    Polynomial = apply(interpolated_precision$Linear, 2, mean)
  )
  
  # Combine the median PR data into a single data frame for plotting
  median_pr_combined <- rbind(
    data.frame(recall = common_recall, precision = median_precision$Bayesian, Method = rep("Bayesian", length(common_recall))),
    data.frame(recall = common_recall, precision = median_precision$GAM, Method = rep("GAM", length(common_recall))),
    data.frame(recall = common_recall, precision = median_precision$Polynomial, Method = rep("Polynomial", length(common_recall)))
  )
  
  # Generate PDF filename based on the scenario (1-group or 2-group)
  pdf_filename <- file.path("figures", sub("v[69]_", "", sub(".rda", paste0("_", scenario, "_pr.pdf"), file)))
  
  
  median_pr_combined$Method <- factor(median_pr_combined$Method, levels = c('Bayesian', 'Polynomial', 'GAM')) 
  # Plot Precision-Recall results
  pdf(file = pdf_filename, width = 9)
  print(ggplot(median_pr_combined, aes(x = recall, y = precision, color = Method)) +
          geom_line(size = 2) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
          labs(title = paste('Precision-Recall Curve (', scenario, ')', sep = ''),
               x = "Recall", y = "Precision") +
          theme_classic() +
          ylim(0,1) +
          theme(axis.text.x = element_text(size = 30, vjust = 0.85, color = "black"),
                axis.text.y = element_text(hjust = 1, size = 30, color = "black"),
                axis.title.x = element_text(size = 25, color = "black"),
                axis.title.y = element_text(size = 25, color = "black"),
                legend.position = "bottom",
                axis.line = element_line(colour = "black", size = 0.8), 
                axis.ticks = element_line(colour = "black", size = 0.8),
                axis.ticks.length = unit(0.3, "cm"))+
          scale_color_manual(values = c("Bayesian" = "#A9082C", "GAM" = "#4B61A8", "Polynomial" = "#FDB76D"),
                             breaks = c("Bayesian", "GAM", "Polynomial")))
  dev.off()
}

# Loop through TDR files and create plots
for (file in tdr_files) {
  create_tdr_plot(file)
}

# Loop through ROC files and create plots
for (file in roc_files) {
  create_roc_plot(file)
}

# Loop through Precision-Recall files and create plots
for (file in pr_files) {
  create_pr_plot(file)
}
#####################

# Load necessary libraries
library(ComplexHeatmap)
library(RColorBrewer)
library(viridis)
# Define the top regions and cell numbers
topregions <- c(100, 400, 700, 1000)
cell_counts <- c(1000, 2000, 4000, 8000)

# Create an empty matrix to store the TDR averages
tdr_avg_mat <- matrix(nrow = length(topregions) * length(cell_counts), ncol = 3)
colnames(tdr_avg_mat) <- c("Bayesian", "GAM", "Linear")

# Create row names based on the combinations of topregions and cell counts
rownames(tdr_avg_mat) <- c(paste0(topregions[1], ":", cell_counts),
                           paste0(topregions[2], ":", cell_counts),
                           paste0(topregions[3], ":", cell_counts),
                           paste0(topregions[4], ":", cell_counts))

# Index to store the average TDR values
index <- 1

# Loop over the combinations of topregions and cell counts
for (topregion in topregions) {
  for (cell_count in cell_counts) {
    
    # Load the TDR results for the specific topregion and cell count
    load(paste0("v9_sim_tdr_results_0.1_", cell_count, ".rda"))  # Adjust file path
    
    # Calculate the average TDR for each method (Bayesian, GAM, Linear)
    for (method_idx in seq_along(tdr_results)) {
      tdr_avg_mat[index, method_idx] <- colMeans(tdr_results[[method_idx]], na.rm = TRUE)[topregion/100]
    }
    
    # Move to the next row for the next topregion/cell_count combination
    index <- index + 1
  }
}

# Convert the resulting list to a matrix
tdr_avg_mat <- as.matrix(tdr_avg_mat)

rownames(tdr_avg_mat)<- 1:nrow(tdr_avg_mat)
#mns<- colMeans(tdr_avg_mat, na.rm=TRUE)
#tdr_avg_mat<- tdr_avg_mat[,order(mns, decreasing = T)]
annotdf <- data.frame(row.names = rownames(tdr_avg_mat), TopRegion = c(rep(c("100", "400", "700", "1000"), each = length(cell_counts))))
row_ha = rowAnnotation(TopGene = factor(c(rep(c("100", "400", "700", "1000"), each = length(cell_counts))),levels = c("100", "400", "700", "1000")), col = list(TopGene = c(
  "100" = "wheat",            # Warm, light beige
  "400" = "tan2",             # Medium tan
  "700" = "chocolate2",       # Rich chocolate brown
  "1000" = "sienna4"          # Deep, dark brown
)),
annotation_name_gp= gpar(fontsize = 18))
#tdr_avg_mat[tdr_avg_mat < 0.5]<- 0.5
names(annotdf) <- "TopGene"
png(paste0("1group_TDR_heatmap_0.1.png"), width = 1100*4.166667, height = 960*4.166667, res=300)
Heatmap(tdr_avg_mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        rect_gp = gpar(col = "grey", lwd = 2),
        row_labels = rep(cell_counts, length(topregions)),
        split=factor(annotdf$TopGene, levels = c("1000","700", "400", "100")),
        left_annotation = row_ha,
        gap = unit(3, "mm"),
        column_names_rot = 360,
        column_names_centered = TRUE,
        col = viridis(100),
        heatmap_legend_param = list(
          legend_gp = gpar(fontsize = 18),
          title = "TDR",
          at = seq(0.4, 1, by = 0.2),
          labels = seq(0.4, 1, by = 0.2)
        )
        ,row_names_gp = gpar(fontsize = 18),
        column_names_gp = gpar(fontsize = 18)
)

dev.off()

#################
# Load necessary libraries
library(ComplexHeatmap)
library(RColorBrewer)

# Define the top regions and cell numbers
topregions <- c(100, 400, 700, 1000)
cell_counts <- c(1000, 2000, 4000, 8000)

# Create an empty matrix to store the TDR averages
tdr_avg_mat <- matrix(nrow = length(topregions) * length(cell_counts), ncol = 3)
colnames(tdr_avg_mat) <- c("Bayesian", "GAM", "Linear")

# Create row names based on the combinations of topregions and cell counts
rownames(tdr_avg_mat) <- c(paste0(topregions[1], ":", cell_counts),
                           paste0(topregions[2], ":", cell_counts),
                           paste0(topregions[3], ":", cell_counts),
                           paste0(topregions[4], ":", cell_counts))

# Index to store the average TDR values
index <- 1

# Loop over the combinations of topregions and cell counts
for (topregion in topregions) {
  for (cell_count in cell_counts) {
    
    # Load the TDR results for the specific topregion and cell count
    load(paste0("v6_2group_sim_tdr_results_0.1_", cell_count, ".rda"))  # Adjust file path
    
    # Calculate the average TDR for each method (Bayesian, GAM, Linear)
    for (method_idx in seq_along(tdr_results)) {
      tdr_avg_mat[index, method_idx] <- colMeans(tdr_results[[method_idx]], na.rm = TRUE)[topregion/100]
    }
    
    # Move to the next row for the next topregion/cell_count combination
    index <- index + 1
  }
}

# Convert the resulting list to a matrix
tdr_avg_mat <- as.matrix(tdr_avg_mat)

rownames(tdr_avg_mat)<- 1:nrow(tdr_avg_mat)
#mns<- colMeans(tdr_avg_mat, na.rm=TRUE)
#tdr_avg_mat<- tdr_avg_mat[,order(mns, decreasing = T)]
annotdf <- data.frame(row.names = rownames(tdr_avg_mat), TopRegion = c(rep(c("100", "400", "700", "1000"), each = length(cell_counts))))
row_ha = rowAnnotation(TopRegion = factor(c(rep(c("100", "400", "700", "1000"), each = length(cell_counts))),levels = c("100", "400", "700", "1000")), col = list(TopRegion = c(
  "100" = "wheat",            # Warm, light beige
  "400" = "tan2",             # Medium tan
  "700" = "chocolate2",       # Rich chocolate brown
  "1000" = "sienna4"          # Deep, dark brown
)),
annotation_name_gp= gpar(fontsize = 18))
#tdr_avg_mat[tdr_avg_mat < 0.5]<- 0.5
png(paste0("2group_TDR_heatmap_0.1.png"), width = 1100*4.166667, height = 960*4.166667, res=300)
Heatmap(tdr_avg_mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        rect_gp = gpar(col = "grey", lwd = 2),
        row_labels = rep(cell_counts, length(topregions)),
        split=factor(annotdf$TopRegion, levels = c("1000","700", "400", "100")),
        left_annotation = row_ha,
        gap = unit(3, "mm"),
        column_names_rot = 360,
        column_names_centered = TRUE,
        col = viridis(100),
        heatmap_legend_param = list(
          legend_gp = gpar(fontsize = 18),
          title = "TDR"
        )
        ,row_names_gp = gpar(fontsize = 18),
        column_names_gp = gpar(fontsize = 18)
)

dev.off()