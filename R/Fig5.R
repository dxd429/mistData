rm(list = ls())
load("dat_cg.rda")
# Extract gene names from the first column of mat_cg
genes <- mat_cg$gene
rownames(mat_cg) <- genes
mat_cg <- mat_cg[, -1]  # Remove the "gene" column

# Ensure rownames in meta correspond to column names in mat_cg
colnames(mat_cg) <- sub("^X", "", colnames(mat_cg))  # Remove "X" prefix from colnames in mat_cg
rownames(mat_cg) <- gsub("_CG", "", rownames(mat_cg))
# Handle NA or NaN in the expression data by setting them to zero
mat_cg[is.na(mat_cg)] <- 0

get_earliest_principal_node <- function(cds){
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex))))]
  
  root_pr_nodes
}
# Function to process a subset with feature selection and save results
process_subset <- function(subset_meta, subset_expr, label) {
  # Round and scale expression data
  subset_expr <- round(100 * subset_expr)
  
  # Create a Cell Data Set (CDS)
  cds <- new_cell_data_set(
    as.matrix(subset_expr),
    cell_metadata = subset_meta,
    gene_metadata = data.frame(gene_short_name = rownames(subset_expr), row.names = rownames(subset_expr))
  )
  
  # Feature selection: Retain genes expressed in >10% of cells
  data <- counts(cds)  # Extract the counts matrix from the CDS
  remain_idx <- which(rowSums(data > 0) > (ncol(data) * 0.1))  # Genes expressed in >10% of cells
  cds <- cds[remain_idx, ]
  
  # Preprocess the data
  cds <- preprocess_cds(cds, num_dim = 100)
  
  # Reduce dimensions using UMAP
  cds <- reduce_dimension(cds, reduction_method = "UMAP")
  
  # Cluster cells using UMAP coordinates
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  
  # Learn the trajectory graph
  cds <- learn_graph(cds, use_partition = FALSE)
  
  # Choose a root cell for pseudotime ordering
  #root_cell <- colnames(subset_expr)[1]  # Adjust as needed
  cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds[, colData(cds)$L2 == "RG"]))
  
  # Extract pseudotime
  # pseudotime_values <- pseudotime(cds)
  # 
  # # Combine pseudotime with cell names
  # result <- data.frame(
  #   Cell = names(pseudotime_values),
  #   Pseudotime = pseudotime_values
  # )
  # cds_filtered <- cds[, is.finite(result$Pseudotime)]
  # cds_filtered <- preprocess_cds(cds_filtered, num_dim = 50)
  # cds_filtered <- reduce_dimension(cds_filtered, reduction_method = "UMAP")
  # cds_filtered <- cluster_cells(cds_filtered, reduction_method = "UMAP")
  # cds_filtered <- learn_graph(cds_filtered)
  # cds_filtered <- order_cells(cds_filtered, root_pr_nodes=get_earliest_principal_node(cds_filtered))
  
  # Save the result
  #write.csv(result, paste0(label, "_pseudotime.csv"), row.names = FALSE)
  #pseudotime(cds_filtered)[is.finite(pseudotime(cds_filtered))] <- pseudotime(cds_filtered)[is.finite(pseudotime(cds_filtered))]/max(pseudotime(cds_filtered)[is.finite(pseudotime(cds_filtered))], na.omit=TRUE)
  # Free up memory
  #rm(cds, cds_filtered, subset_expr, data, remain_idx)
  #gc()
  cds
}
library(dplyr)
library(monocle3)
require(Seurat)
require(data.table)
library(tibble)
library(ggplot2)
meta <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
rownames(meta) <- gsub("\\.", "-", rownames(meta))
table(meta$Region)
table(meta$Age_groups)
meta <- meta[,-1] %>%
  rownames_to_column(var = "Sample")
rownames(meta) <- meta[,"Sample"]
cds_filtered <- process_subset(meta[, c("L2", "Age_groups"), drop = FALSE], mat_cg, paste0("All"))
save(cds_filtered, file = 'cds_filtered.rda')
colData(cds_filtered)$Age_groups <- meta$Age_groups
plot_cells(cds_filtered,
           color_cells_by = "Age_groups",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
ggsave(file.path(paste0("All_UMAP_Age_group.png")))

plot_cells(cds_filtered,
           color_cells_by = "L2",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
ggsave(file.path(paste0("All_UMAP_L2.png")))
############
# For Region
region_counts <- table(meta$Region)
print(region_counts)

# For Age_groups
age_counts <- table(meta$Age_groups)
print(age_counts)

# For a cross-tabulation of Region and Age_groups
cross_tab <- table(meta$Region, meta$Age_groups)
print(cross_tab)
##############
library(ggplot2)

# # Bar plot for Region distribution
# pdf("demo_region.pdf")
# ggplot(meta, aes(x = Region)) +
#   geom_bar(fill = "skyblue") +
#   theme_classic() +
#   labs(title = "Distribution of Regions", x = "Region", y = "Count")
# dev.off()
# # Bar plot for Age_groups distribution
# pdf("demo_age.pdf")
# ggplot(meta, aes(x = Age_groups)) +
#   geom_bar(fill = "salmon") +
#   theme_classic() +
#   labs(title = "Distribution of Age Groups", x = "Age Group", y = "Count")
# dev.off()
meta$Age_groups <- factor(meta$Age_groups, levels = c("2T", "3T", "infant", "adult"))
pdf("demo_all.pdf", width = 10)
ggplot(meta, aes(x = Age_groups, fill = Age_groups)) +
  geom_bar() +
  facet_wrap(~ Region) +
  # Option 1: Use a Brewer palette (change "Set2" to "Dark2", "Paired", etc. as desired)
  scale_fill_brewer(palette = "Dark2") +
  # Option 2 (alternative): Use the viridis package for a colorblind-friendly palette
  # scale_fill_viridis_d(option = "plasma") +
  theme_classic(base_size = 14) +  # This sets a base font size for many elements
  theme(
    axis.text = element_text(size = 14),        # Axis tick labels
    axis.title = element_text(size = 16),       # Axis titles
    plot.title = element_text(size = 18, face = "bold"),  # Plot title
    legend.title = element_text(size = 14),     # Legend title
    legend.text = element_text(size = 12),      # Legend items
    strip.text = element_text(size = 20)        # Facet labels
  ) +
  labs(
    title = "Age Groups Distribution by Region", 
    x = "Age Group", 
    y = "Count"
  )

dev.off()
###############


#####
library(ggplot2)
library(dplyr)

# Prepare the data for L1
l1_df <- meta %>%
  group_by(L1) %>%
  summarise(Count = n())

# Create the pie chart for L1 with count labels on the slices
# pdf("Pie_L1.pdf")
# ggplot(l1_df, aes(x = "", y = Count, fill = L1)) +
#   geom_bar(stat = "identity", width = 1) +
#   coord_polar(theta = "y") +
#   geom_text(aes(label = Count), 
#             position = position_stack(vjust = 0.5), 
#             color = "white", size = 4) +
#   theme_void() +
#   ggtitle("Cell Identity Distribution: L1") +
#   theme(legend.title = element_blank())
# dev.off()


# Prepare the data for L2
# Prepare the data for L2
# Prepare the data for L2
library(viridis)
l2_df <- meta %>%
  group_by(L2) %>%
  summarise(Count = n())

# Create the pie chart for L2 with count labels on the slices
pdf("Pie_L2.pdf", width = 10)
ggplot(l2_df, aes(x = "", y = Count, fill = L2)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  geom_text(aes(label = Count), 
            position = position_stack(vjust = 0.5), 
            color = "black", size = 4) +
  theme_void() +
  ggtitle("Cell Identity Distribution: L2") +
  theme(legend.title = element_blank()) +
  scale_fill_brewer(palette = "Set3")

dev.off()
############

rm(list = ls())
gc()
library(SingleCellExperiment)

library(ggplot2)
library(pheatmap)
poly_function <- function(x, coeffs) {
  return(1/(1+exp(-coeffs[5]*x^4 - coeffs[4]*x^3 - coeffs[3]*x^2 - coeffs[2]*x - coeffs[1])))
}

load("int_CTX_HPC.rda")
# load("int_infant_adult.rda")
# load("int_infant_CTX.rda")
# load("int_infant_HPC.rda")
# load("int_CTX_HPC.rda")

load("processed_sce_CTX_1.rda")
dat_sce_CTX <- Dat_sce_new
rm(Dat_sce_new)

load("processed_sce_HPC_1.rda")
dat_sce_HPC <- Dat_sce_new
rm(Dat_sce_new)


mat_CTX <- assay(dat_sce_CTX, "methylation")
mat_HPC <- assay(dat_sce_HPC, "methylation")
rownames(mat_CTX) <- gsub("_CG", "", rownames(mat_CTX))
rownames(mat_HPC) <- gsub("_CG", "", rownames(mat_HPC))

ptime_CTX <- colData(dat_sce_CTX)[["Pseudotime"]]
ptime_HPC <- colData(dat_sce_HPC)[["Pseudotime"]]


####################order by fitted mean
library(pheatmap)
library(viridis)
num_gene <- 50
# Subset the rows for clustering
subset_CTX <- mat_CTX[names(int_CTX_HPC)[1:num_gene], ]

# Perform hierarchical clustering on the rows
row_clusters <- hclust(dist(subset_CTX))  # Default: Euclidean distance and complete linkage
row_order <- row_clusters$order  # Extract the order of clustered rows

# Get the reordered row names
reordered_row_names <- rownames(subset_CTX)[row_order]

# Heatmap for mat_CTX (with clustered rows)
tiff(paste0("real_heat_CTX_", num_gene, ".tiff"), width = 1920, height = 1440, res = 300)
pheatmap(
  mat_CTX[reordered_row_names, order(ptime_CTX)],  # Apply clustered row order
  color = plasma(256),
  scale = "none",
  cluster_rows = FALSE,  # Disable row clustering (already ordered)
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = FALSE,
  angle_col = 0,
  breaks = seq(0, 1, length.out = 256),
  main = "Cerebral Cortex"
)
dev.off()

# Heatmap for mat_HPC (using the same row order as mat_CTX)
tiff(paste0("real_heat_HPC_", num_gene, ".tiff"), width = 1920, height = 1440, res = 300)
pheatmap(
  mat_HPC[reordered_row_names, order(ptime_HPC)],  # Apply the same row order
  color = plasma(256),
  scale = "none",
  cluster_rows = FALSE,  # Ensure rows are not re-clustered
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = FALSE,
  angle_col = 0,
  breaks = seq(0, 1, length.out = 256),
  main = "Hippocampus"
)
dev.off()
##############

rm(list = ls())
gc()
load("processed_sce_CTX_1.rda")
dat_sce_HPC <- Dat_sce_new
rm(Dat_sce_new)

load("processed_sce_HPC_1.rda")
dat_sce_CTX <- Dat_sce_new
rm(Dat_sce_new)
poly_function <- function(x, coeffs) {
  return(1/(1+exp(-coeffs[5]*x^4 - coeffs[4]*x^3 - coeffs[3]*x^2 - coeffs[2]*x - coeffs[1])))
}

load("int_CTX_HPC.rda")
library(SingleCellExperiment)

pars_list <- list()
pars_list[["TOLLIP_CTX"]] <- rowData(dat_sce_CTX)$mist_pars["TOLLIP", ]
pars_list[["TOLLIP_HPC"]] <- rowData(dat_sce_HPC)$mist_pars["TOLLIP", ]
pars_list[["SOX8_CTX"]] <- rowData(dat_sce_CTX)$mist_pars["SOX8", ]
pars_list[["SOX8_HPC"]] <- rowData(dat_sce_HPC)$mist_pars["SOX8", ]
pars_list[["FGFR3_CTX"]] <- rowData(dat_sce_CTX)$mist_pars["FGFR3", ]
pars_list[["FGFR3_HPC"]] <- rowData(dat_sce_HPC)$mist_pars["FGFR3", ]
pars_list[["CXXC5_CTX"]] <- rowData(dat_sce_CTX)$mist_pars["CXXC5", ]
pars_list[["CXXC5_HPC"]] <- rowData(dat_sce_HPC)$mist_pars["CXXC5", ]
pars_list[["CAMK2A_CTX"]] <- rowData(dat_sce_CTX)$mist_pars["CAMK2A", ]
pars_list[["CAMK2A_HPC"]] <- rowData(dat_sce_HPC)$mist_pars["CAMK2A", ]
pars_list[["KLF9_CTX"]] <- rowData(dat_sce_CTX)$mist_pars["KLF9", ]
pars_list[["KLF9_HPC"]] <- rowData(dat_sce_HPC)$mist_pars["KLF9", ]
# pars_list[["SEMA7A_Infant"]] <- rowData(dat_sce_infant)$mist_pars["SEMA7A_CG", ]
# pars_list[["SEMA7A_Adult"]] <- rowData(dat_sce_adult)$mist_pars["SEMA7A_CG", ]
# pars_list[["SLIT1_Infant"]] <- rowData(dat_sce_infant)$mist_pars["SLIT1_CG", ]
# pars_list[["SLIT1_Adult"]] <- rowData(dat_sce_adult)$mist_pars["SLIT1_CG", ]
# pars_list[["FOXP1_Infant"]] <- rowData(dat_sce_infant)$mist_pars["FOXP1_CG", ]
# pars_list[["FOXP1_Adult"]] <- rowData(dat_sce_adult)$mist_pars["FOXP1_CG", ]

dat_list <- list()
dat_list[["TOLLIP_CTX"]] <- assay(dat_sce_CTX, "methylation")["TOLLIP", ]
dat_list[["TOLLIP_HPC"]] <- assay(dat_sce_HPC, "methylation")["TOLLIP", ]
dat_list[["SOX8_CTX"]] <- assay(dat_sce_CTX, "methylation")["SOX8", ]
dat_list[["SOX8_HPC"]] <- assay(dat_sce_HPC, "methylation")["SOX8", ]
dat_list[["FGFR3_CTX"]] <- assay(dat_sce_CTX, "methylation")["FGFR3", ]
dat_list[["FGFR3_HPC"]] <- assay(dat_sce_HPC, "methylation")["FGFR3", ]
dat_list[["CXXC5_CTX"]] <- assay(dat_sce_CTX, "methylation")["CXXC5", ]
dat_list[["CXXC5_HPC"]] <- assay(dat_sce_HPC, "methylation")["CXXC5", ]
dat_list[["CAMK2A_CTX"]] <- assay(dat_sce_CTX, "methylation")["CAMK2A", ]
dat_list[["CAMK2A_HPC"]] <- assay(dat_sce_HPC, "methylation")["CAMK2A", ]
dat_list[["KLF9_CTX"]] <- assay(dat_sce_CTX, "methylation")["KLF9", ]
dat_list[["KLF9_HPC"]] <- assay(dat_sce_HPC, "methylation")["KLF9", ]

# dat_list[["SEMA7A_Infant"]] <- assay(dat_sce_infant, "methylation")["SEMA7A_CG", ]
# dat_list[["SEMA7A_Adult"]] <- assay(dat_sce_adult, "methylation")["SEMA7A_CG", ]
# dat_list[["SLIT1_Infant"]] <- assay(dat_sce_infant, "methylation")["SLIT1_CG", ]
# dat_list[["SLIT1_Adult"]] <- assay(dat_sce_adult, "methylation")["SLIT1_CG", ]
# dat_list[["FOXP1_Infant"]] <- assay(dat_sce_infant, "methylation")["FOXP1_CG", ]
# dat_list[["FOXP1_Adult"]] <- assay(dat_sce_adult, "methylation")["FOXP1_CG", ]

ptime_list <- list()
ptime_list[["TOLLIP_CTX"]] <- colData(dat_sce_CTX)[["Pseudotime"]]
ptime_list[["TOLLIP_HPC"]] <- colData(dat_sce_HPC)[["Pseudotime"]]
ptime_list[["SOX8_CTX"]] <- colData(dat_sce_CTX)[["Pseudotime"]]
ptime_list[["SOX8_HPC"]] <- colData(dat_sce_HPC)[["Pseudotime"]]
ptime_list[["FGFR3_CTX"]] <- colData(dat_sce_CTX)[["Pseudotime"]]
ptime_list[["FGFR3_HPC"]] <- colData(dat_sce_HPC)[["Pseudotime"]]
ptime_list[["CXXC5_CTX"]] <- colData(dat_sce_CTX)[["Pseudotime"]]
ptime_list[["CXXC5_HPC"]] <- colData(dat_sce_HPC)[["Pseudotime"]]
ptime_list[["CAMK2A_CTX"]] <- colData(dat_sce_CTX)[["Pseudotime"]]
ptime_list[["CAMK2A_HPC"]] <- colData(dat_sce_HPC)[["Pseudotime"]]
ptime_list[["KLF9_CTX"]] <- colData(dat_sce_CTX)[["Pseudotime"]]
ptime_list[["KLF9_HPC"]] <- colData(dat_sce_HPC)[["Pseudotime"]]

# ptime_list[["SEMA7A_Infant"]] <- colData(dat_sce_infant)[["Pseudotime"]]
# ptime_list[["SEMA7A_Adult"]] <- colData(dat_sce_adult)[["Pseudotime"]]
# ptime_list[["SLIT1_Infant"]] <- colData(dat_sce_infant)[["Pseudotime"]]
# ptime_list[["SLIT1_Adult"]] <- colData(dat_sce_adult)[["Pseudotime"]]
# ptime_list[["FOXP1_Infant"]] <- colData(dat_sce_infant)[["Pseudotime"]]
# ptime_list[["FOXP1_Adult"]] <- colData(dat_sce_adult)[["Pseudotime"]]

# genenames <- c("NFIX_CTX", "NFIX_HPC", "FGFR3_CTX", 
#                "FGFR3_HPC", "SEMA7A_Infant", "SEMA7A_Adult", 
#                "SLIT1_Infant", "SLIT1_Adult",
#                "FOXP1_Infant", "FOXP1_Adult")
genenames <- c("TOLLIP_CTX", "TOLLIP_HPC", "SOX8_CTX", "SOX8_HPC", "FGFR3_CTX", "FGFR3_HPC",
               "CXXC5_CTX", "CXXC5_HPC", "CAMK2A_CTX", "CAMK2A_HPC",
               "KLF9_CTX", "KLF9_HPC")
for (i in 1:length(genenames)) {
  genename <- genenames[i]
  
  dat <- dat_list[[genename]]
  ptime <- ptime_list[[genename]]
  
  # Remove NaN values
  ind_remove <- is.nan(dat)
  dat <- dat[!ind_remove]
  ptime <- ptime[!ind_remove]
  
  # Keep only finite, non-NA pseudotime values and normalize them
  ptime <- ptime[is.finite(ptime) & !is.na(ptime)]
  ptime <- ptime / max(ptime)
  
  # Correction: Replace 0 and 1 values with small or high random values, respectively
  dat[dat == 0] <- runif(length(dat[dat == 0]), 0.00001, 0.01)
  dat[dat == 1] <- runif(length(dat[dat == 1]), 0.95, 0.99999)
  
  df <- data.frame(Pseudotime = ptime, Methylation = dat)
  df <- df[df$Methylation >= 0.05, ]
  
  beta_mu_mean1 <- pars_list[[genename]][1:5]
  
  # Set the line color based on the gene name suffix
  if (grepl("_CTX", genename)) {
    line_color <- "#00FFFF"  # Cyan for _CTX
  } else if (grepl("_HPC", genename)) {
    line_color <- "#00FF00"  # Magenta for _HPC
  } else {
    line_color <- "#000000"  # Fallback color if needed
  }
  df <- df[df$Methylation >= 0.05, ]
  library(ggplot2)
  png(file = paste0("trend_", genename, ".png"), width = 2770, height = 1500, res = 300)
  
  print(
    ggplot(df, aes(x = Pseudotime, y = Methylation, color = Methylation)) +
      geom_point(size = 1) +
      geom_line(aes(y = poly_function(Pseudotime, beta_mu_mean1)),
                color = line_color, linewidth = 1.8) +
      geom_vline(xintercept = 0.25, linetype = "dashed", color = "black") +
      geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") +
      geom_vline(xintercept = 0.75, linetype = "dashed", color = "black") +
      labs(title = paste0(genename), x = "Pseudotime", y = "Methylation Level") +
      scale_color_gradient(low = "blue", high = "red", limits = c(0, 1)) +
      theme_classic() +
      theme(
        plot.title   = element_text(size = 25, face = "bold"),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x  = element_text(size = 20),
        axis.text.y  = element_text(size = 20)
      )
  )
  
  dev.off()
}


