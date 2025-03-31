rm(list = ls())
gc()
library(SingleCellExperiment)
library(pheatmap)
poly_function <- function(x, coeffs) {
  return(1/(1+exp(-coeffs[5]*x^4 - coeffs[4]*x^3 - coeffs[3]*x^2 - coeffs[2]*x - coeffs[1])))
}

load('m_reduced_clean.rda')
load('ptime_all.rda')
int_res <- read.csv('Integral_res.csv')

library(pheatmap)
library(viridis)
num_gene <- 100
# Subset the rows for clustering
subset_dat <- m_reduced_clean[int_res$gene[1:num_gene], ]

# Perform hierarchical clustering on the rows
row_clusters <- hclust(dist(subset_dat))  # Default: Euclidean distance and complete linkage
row_order <- row_clusters$order  # Extract the order of clustered rows

# Get the reordered row names
reordered_row_names <- rownames(subset_dat)[row_order]

# Heatmap for mat_2T (with clustered rows)
tiff(paste0("BR_real_heat_g1_", num_gene, ".tiff"), width = 1920, height = 1440, res = 300)
pheatmap(
  m_reduced_clean[reordered_row_names, order(ptime_all)],  # Apply clustered row order
  color = plasma(256),
  scale = "none",
  cluster_rows = FALSE,  # Disable row clustering (already ordered)
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = FALSE,
  angle_col = 0,
  breaks = seq(0, 1, length.out = 256),
  main = " "
)
dev.off()
############################UMAP
library(monocle3)
require(Seurat)
library(ggplot2)
load('m_reduced_clean.rda')
m_reduced_clean[is.na(m_reduced_clean)] <- 0
genebody_SCE <- readRDS("./Data/genebody_SCE.rds")
meta_dat <- as.matrix(as.data.frame(colData(genebody_SCE)@listData))
rownames(meta_dat) <- meta_dat[,'id_met']

get_earliest_principal_node <- function(cds, time_bin="E4.5"){
  cell_ids <- which(colData(cds)[, "stage"] == time_bin)
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
# Function to process a subset with feature selection and save results
process_subset <- function(subset_meta, subset_expr) {
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
  cds <- reduce_dimension(cds)
  
  # Cluster cells using UMAP coordinates
  cds <- cluster_cells(cds, resolution=1e-5)
  
  # Learn the trajectory graph
  cds <- learn_graph(cds, use_partition = FALSE)
  
  # Choose a root cell for pseudotime ordering
  #root_cell <- colnames(subset_expr)[1]  # Adjust as needed
  cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
  
  plot_cells(cds,
             cell_size=1,
             color_cells_by = "pseudotime",
             label_cell_groups=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE,
             graph_label_size=1.5)
  ggsave(file.path(paste0("UMAP_pseudotime.png")), width = 10, height = 7)
  
  plot_cells(cds,
             cell_size=1,
             color_cells_by = "stage",
             label_cell_groups=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE,
             graph_label_size=1.5)
  ggsave(file.path(paste0("UMAP_stage.png")), width = 10, height = 7)
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
  rm(cds, subset_expr, data, remain_idx)
  gc()
}

process_subset(meta_dat[, c("stage", "lineage10x_2"), drop = FALSE], m_reduced_clean)

################
genebody_SCE <- readRDS("C:/Users/dd284/OneDrive/harry/pseudo/genebody_SCE.rds")
meta_dat <- as.data.frame(colData(genebody_SCE)@listData)
library(ggplot2)
library(dplyr)

library(viridis)
CT_df <- meta_dat %>%
  group_by(lineage10x_2) %>%
  summarise(Count = n()) %>%
  na.omit()

# Create the pie chart for L2 with count labels on the slices
pdf("Pie_CT.pdf", width = 10)
ggplot(CT_df, aes(x = "", y = Count, fill = lineage10x_2)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  geom_text(aes(label = Count), 
            position = position_stack(vjust = 0.5), 
            color = "black", size = 4) +
  theme_void() +
  ggtitle("Cell Identity Distribution") +
  theme(legend.title = element_blank()) +
  scale_fill_brewer(palette = "Set3")

dev.off()
################
load('m_reduced_clean.rda')
gene_mat <- m_reduced_clean
genebody_SCE <- readRDS("C:/Users/dd284/OneDrive/harry/pseudo/genebody_SCE.rds")
pseudotime <- genebody_SCE$monocle3_level_pseudotime

n_genes <- nrow(gene_mat)  # number of genes
pval_gam <- numeric(n_genes)
pval_lm  <- numeric(n_genes)

library(mgcv)

for (i in seq_len(n_genes)) {
  print(i)
  # 1) Extract response for gene i
  response_i <- as.numeric(gene_mat[i, ])  # row i
  
  # 2) Create data frame
  df_i <- data.frame(
    time = pseudotime,
    response = response_i
  )
  
  # ---- GAM part ----
  # 3) Fit GAM and null GAM
  tryCatch({
    model_gam <- gam(
      response ~ s(time, k = 5, bs = "cr"),
      knots = list(time = c(0, 0.25, 0.5, 0.75, 1)),
      data = df_i,
      family = gaussian(link = "identity")
    )
    model_gam_null <- gam(
      response ~ 1,
      data = df_i,
      family = gaussian(link = "identity")
    )
    
    # 4) LRT for GAM
    lrt_gam <- anova(model_gam_null, model_gam, test = "Chisq")
    pval_gam[i] <- lrt_gam$`Pr(>Chi)`[2]  # second row = actual GAM vs. null
    
    # ---- LM part ----
    # We'll fit a polynomial of degree 4 as in your example,
    # you can adjust the degree if you like
    df_i$poly1 <- poly(df_i$time, degree = 4)[,1]
    df_i$poly2 <- poly(df_i$time, degree = 4)[,2]
    df_i$poly3 <- poly(df_i$time, degree = 4)[,3]
    df_i$poly4 <- poly(df_i$time, degree = 4)[,4]
    
    model_lm <- lm(response ~ poly1 + poly2 + poly3 + poly4, data = df_i)
    model_lm_null <- lm(response ~ 1, data = df_i)
    
    # 5) LRT (ANOVA) for LM
    lrt_lm <- anova(model_lm_null, model_lm)
    pval_lm[i] <- lrt_lm$`Pr(>F)`[2]
  }, error = function(e) {
    message(sprintf("Skipping gene %d due to error: %s", i, e$message))
    # pval_gam[i] and pval_lm[i] remain NA
  })
}

# Transform p-values to -log10
neglog_p_gam <- -log10(pval_gam)
neglog_p_lm  <- -log10(pval_lm)

# Order from largest to smallest (-log10 p-value)
ord_gam <- order(neglog_p_gam, decreasing = TRUE)
ord_lm  <- order(neglog_p_lm,  decreasing = TRUE)

gene_names <- rownames(gene_mat)  # if row names are set
gam_sorted_gene_names <- gene_names[ord_gam]
lm_sorted_gene_names  <- gene_names[ord_lm]

# Inspect top 10
head(gam_sorted_gene_names, 10)
head(lm_sorted_gene_names, 10)

df_gam <- data.frame(
  gene       = gene_names,
  test_stats = neglog_p_gam,
  stringsAsFactors = FALSE
)

df_lm <- data.frame(
  gene       = gene_names,
  test_stats = neglog_p_lm,
  stringsAsFactors = FALSE
)

# 3) Sort (rank) each data frame from largest to smallest
df_gam_ranked <- df_gam[order(df_gam$test_stats, decreasing = TRUE), ]
df_lm_ranked  <- df_lm[order(df_lm$test_stats, decreasing = TRUE), ]

Integral_res <- read.csv("Integral_res.csv")
df_bayes_ranked <- Integral_res[,2:3]
names(df_bayes_ranked) <- c("gene", "test_stats")
###########
# Assuming each data frame has columns:
#   "gene"       (character)  
#   "test_stats" (numeric, already sorted from largest to smallest -log10 p-values)

top_bayes <- df_bayes_ranked$gene[1:10000]
top_gam   <- df_gam_ranked$gene[1:10000]
top_lm    <- df_lm_ranked$gene[1:10000]

# Named list of sets (top 10,000 genes for each method)
my_sets <- list(
  Bayesian   = top_bayes,
  GAM        = top_gam,
  Polynomial = top_lm
)
library(VennDiagram)
venn.plot <- venn.diagram(
  x = my_sets,
  filename = NULL, # don't write to a file; just capture for plotting
  main = "Overlap among Top 10,000 Genes",
  category.names = c("Bayesian", "GAM", "Polynomial"),
  
  # Specify your desired fill colors here
  fill  = c("#F39C12", "#E74C3C", "#2ECC71"), # example: gold, red, green
  alpha = c(0.5, 0.5, 0.5),   # adjust transparency if needed
  
  # Outline and text appearance
  lty         = "blank", 
  cex         = 1.2,
  main.cex    = 1.4,
  cat.cex     = 1.2,
  cat.fontface= "bold",
  margin      = 0.1
)

# Render the Venn diagram
library(grid)
png(file = "venn_diag.png")
grid.newpage()
grid.draw(venn.plot)
dev.off()

pdf(file = "venn_diag.pdf")
grid.newpage()
grid.draw(venn.plot)
dev.off()
##### Signa Distribution
load("inter_Beta_sigma_all_0.025.rda")

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggridges)

# 1. Extract the last 4 positions into a data frame
df_sigma <- do.call(rbind, lapply(beta_mu_mean, function(x) {
  c(sigma2_1 = x[6],
    sigma2_2 = x[7],
    sigma2_3 = x[8],
    sigma2_4 = x[9])
})) %>%
  as.data.frame()

# 2. Reshape into long format
df_long <- df_sigma %>%
  pivot_longer(
    cols = everything(),
    names_to = "Parameter",
    values_to = "Value"
  ) %>%
  mutate(
    Parameter = factor(
      Parameter,
      levels = rev(c("sigma2_1", "sigma2_2", "sigma2_3", "sigma2_4")),
      labels = rev(c(
        expression(sigma^2 * "[1]"),
        expression(sigma^2 * "[2]"),
        expression(sigma^2 * "[3]"),
        expression(sigma^2 * "[4]")
      ))
    )
  )

# 3. Filter out infinite and/or missing values
df_long_clean <- df_long %>%
  filter(!is.infinite(Value), !is.na(Value))

# 4. Plot (no x-range restriction, removing Inf values)
pdf("ridge_sigma.pdf", width = 9, height = 6)
ggplot(subset(df_long_clean, Value < 30), aes(x = Value, y = Parameter, fill = Parameter)) +
  geom_density_ridges(from = 0) +
  #geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  theme_classic() +
  labs(
    title = "Distribution of Sigma^2 Parameters",
    x = "Value",
    y = "Parameter"
  ) +
  theme(
    plot.title   = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x  = element_text(size = 15),
    axis.text.y  = element_text(size = 15),
    axis.line    = element_line(colour = "black", size = 0.8),
    axis.ticks   = element_line(colour = "black", size = 0.8),
    axis.ticks.length = unit(0.2, "cm")
  ) +
  scale_fill_manual(values = c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4")) +
  guides(fill = "none")
dev.off()

