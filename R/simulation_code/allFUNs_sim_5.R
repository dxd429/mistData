.libPaths("/mnt/pan/SOM_PQHS_HXF155/daoyu/migrate/library")
setwd("/mnt/pan/SOM_PQHS_HXF155/daoyu/Pseudotime/scripts/pseu_sim_1")

source('functions.R')
# Define the function to be parallelized
process_gene <- function(ind, n_gene, z_t_0, beta_mu_sigmatop1800, Num_t, stage_num, Num_cell, dm_ratio) {
  source('functions.R')
  library(mgcv)
  if(ind <= round(dm_ratio * n_gene)){
    b <- sample(1000:3000, 1)
    sigma2_set <- beta_mu_sigmatop1800[[b]][6:9]
    mu_t <- z_t_0 %*% beta_mu_sigmatop1800[[b]][1:5] # Compute mu_t
    n_t1 <- sample_cells(Num_cell = Num_cell, Num_t = Num_t)
    u <- vector("list", length(n_t1)) 
    breaks <- seq(1, Num_t, by = Num_t/stage_num)
    for (k in seq_along(n_t1)) {
      j <- findInterval(k, breaks)
      u[[k]] <- rnorm(n_t1[k], mean = mu_t[k], sd = (sqrt(sigma2_set[j]))) 
    }
  } else {
    b <- sample(3000:6000, 1)
    sigma2_set <- beta_mu_sigmatop1800[[b]][6:9]
    mu_t <- z_t_0 %*% beta_mu_sigmatop1800[[b]][1:5] # Compute mu_t
    n_t1 <- sample_cells(Num_cell = Num_cell, Num_t = Num_t)
    u <- vector("list", length(n_t1)) 
    breaks <- seq(1, Num_t, by = Num_t/stage_num)
    for (k in seq_along(n_t1)) {
      j <- findInterval(k, breaks)
      u[[k]] <- rnorm(n_t1[k], mean = mu_t[k], sd = (sqrt(sigma2_set[j]))) 
    }
  }
  int_res <- Integral_compare(u)
  return(int_res)
}
############Simulate u
load('inter_Beta_sigma_all_0.025.rda')
load('m_reduced_clean.rda')# Function to check for Inf or NA
#####
Num_t <- 40
Num_cell <- 1000
stage_num <- 4
poly_d <- 4
poly_scale <- 0.025
stage_num<- 4 # Set number of stages
z_t<- poly(poly_scale * (1:Num_t), poly_d, raw = T, simple = T) # Generate polynomial base
z_t_0<- cbind(rep(1, nrow(z_t)), z_t)
#####
names(beta_mu_mean) <- rownames(m_reduced_clean)
rm(m_reduced_clean)

integral_ori<- read.csv("Integral_res.csv")
gene_byIntegral<- integral_ori$gene
gene_top1800<- gene_byIntegral[1:10000] # asumme 10% DM
# Function to check for Inf or NA
contains_inf_or_na <- function(x) {
  any(is.infinite(x) | is.na(x))
}

# Apply this function to each element of the list using lapply
bad_elements <- lapply(beta_mu_mean, contains_inf_or_na)
beta_mu_mean <- beta_mu_mean[!unlist(bad_elements)]
beta_mu_sigmatop1800 <- beta_mu_mean[gene_top1800]
rm(beta_mu_mean)
#####simulate logit methylation level
set.seed(429)
n_iterations <- 20
n_gene <- 10000
dm_ratio <- 0.05
cutoffs <- seq(100, n_gene * dm_ratio, by = 100)

# Initialize list to store results
tdr_results <- list(Bayesian = matrix(0, n_iterations, length(cutoffs)),
                    GAM = matrix(0, n_iterations, length(cutoffs)),
                    Linear = matrix(0, n_iterations, length(cutoffs)))
#metrics_results <- list(Bayesian = matrix(0, n_iterations, 8),
#                        GAM = matrix(0, n_iterations, 8),
#                        Linear = matrix(0, n_iterations, 8))
# Initialize a list to store the ROC data for each method and iteration
roc_results <- list(
  Bayesian = list(),
  GAM = list(),
  Linear = list()
)

pr_results <- list(
  Bayesian = list(),
  GAM = list(),
  Linear = list()
)

# Thresholds to evaluate
thresholds <- c(seq(100, 1000, by = 100), seq(1500, 10000, by = 500))

common_fpr <- seq(0, 1, length.out = 100)
common_recall <- seq(0, 1, length.out = 100)

interpolated_tpr <- list(
  Bayesian = matrix(NA, nrow = n_iterations, ncol = length(common_fpr)),
  GAM = matrix(NA, nrow = n_iterations, ncol = length(common_fpr)),
  Linear = matrix(NA, nrow = n_iterations, ncol = length(common_fpr))
)

interpolated_precision <- list(
  Bayesian = matrix(NA, nrow = n_iterations, ncol = length(common_recall)),
  GAM = matrix(NA, nrow = n_iterations, ncol = length(common_recall)),
  Linear = matrix(NA, nrow = n_iterations, ncol = length(common_recall))
)

# True signal indices (rows 1 to 100)
true_signals <- 1:round(dm_ratio * n_gene)
int_list <- vector(mode = 'list', length = n_iterations)
for (i in 1:n_iterations) {
  print(i)
  library(BiocParallel)
  
  # Run the process in parallel
  results <- bplapply(1:n_gene, process_gene, n_gene = n_gene, z_t_0 = z_t_0, beta_mu_sigmatop1800 = beta_mu_sigmatop1800, Num_t = Num_t, stage_num = stage_num, Num_cell = Num_cell, dm_ratio = dm_ratio, BPPARAM = MulticoreParam())
  
  # Combine results into a data matrix
  data_matrix <- do.call(rbind, results)
  
  colnames(data_matrix) <- c('Bayesian', 'GAM', 'Linear')
  int_list[[i]] <- data_matrix
   # Calculate TDR for each column with cutoffs
  tdr_results$Bayesian[i, ] <- calculate_tdr(data_matrix[, 'Bayesian'], true_signals, cutoffs)
  tdr_results$GAM[i, ] <- calculate_tdr(data_matrix[, 'GAM'], true_signals, cutoffs)
  tdr_results$Linear[i, ] <- calculate_tdr(data_matrix[, 'Linear'], true_signals, cutoffs)
  
  roc_results$Bayesian[[i]] <- calculate_roc_data(data_matrix[, 'Bayesian'], true_signals, n_gene, thresholds)
  roc_results$GAM[[i]] <- calculate_roc_data(data_matrix[, 'GAM'], true_signals, n_gene, thresholds)
  roc_results$Linear[[i]] <- calculate_roc_data(data_matrix[, 'Linear'], true_signals, n_gene, thresholds)
  
  interpolated_tpr$Bayesian[i, ] <- interpolate_tpr(roc_results$Bayesian[[i]], common_fpr)
  interpolated_tpr$GAM[i, ] <- interpolate_tpr(roc_results$GAM[[i]], common_fpr)
  interpolated_tpr$Linear[i, ] <- interpolate_tpr(roc_results$Linear[[i]], common_fpr)
  
  pr_results$Bayesian[[i]] <- calculate_pr_data(data_matrix[, 'Bayesian'], true_signals, n_gene, thresholds)
  pr_results$GAM[[i]] <- calculate_pr_data(data_matrix[, 'GAM'], true_signals, n_gene, thresholds)
  pr_results$Linear[[i]] <- calculate_pr_data(data_matrix[, 'Linear'], true_signals, n_gene, thresholds)
  
  interpolated_precision$Bayesian[i, ] <- interpolate_precision(pr_results$Bayesian[[i]], common_recall)
  interpolated_precision$GAM[i, ] <- interpolate_precision(pr_results$GAM[[i]], common_recall)
  interpolated_precision$Linear[i, ] <- interpolate_precision(pr_results$Linear[[i]], common_recall)
  # Calculate general metrics for each column
  #metrics_results$Bayesian[i, ] <- calculate_metrics(data_matrix[, 'Bayesian'], true_signals, n_gene)
  #metrics_results$GAM[i, ] <- calculate_metrics(data_matrix[, 'GAM'], true_signals, n_gene)
  #metrics_results$Linear[i, ] <- calculate_metrics(data_matrix[, 'Linear'], true_signals, n_gene)

  save(tdr_results, file = paste0("v9_sim_tdr_results_", dm_ratio, '_', Num_cell, '.rda'))
  save(roc_results, file = paste0("v9_sim_roc_results_", dm_ratio, '_', Num_cell, '.rda'))
  save(interpolated_tpr, file = paste0("v9_sim_tpr_results_", dm_ratio, '_', Num_cell, '.rda'))
  save(pr_results, file = paste0("v9_sim_pr_results_", dm_ratio, '_', Num_cell, '.rda'))
  save(interpolated_precision, file = paste0("v9_sim_precision_results_", dm_ratio, '_', Num_cell, '.rda'))
  save(int_list, file = paste0("v9_sim_int_list_", dm_ratio, '_', Num_cell, '.rda'))
  #save(metrics_results, file = paste0("sim_metrics_results_", dm_ratio, '.rda'))
}

combined_int <- do.call(rbind, int_list)
save(combined_int, file = paste0("v9_sim_int_lcomp_", dm_ratio, '_', Num_cell, '.rda'))


###################plot
# 
# # Convert results to data frames for plotting
# metrics_df <- data.frame(
#   Iteration = rep(1:n_iterations, times = 3 * length(cutoffs)),
#   Column = rep(rep(c('Bayesian', 'GAM', 'Linear'), each = n_iterations), times = length(cutoffs)),
#   Cutoff = rep(cutoffs, each = n_iterations * 3),
#   TP = c(metrics_results$Bayesian[, 1, ], metrics_results$GAM[, 1, ], metrics_results$Linear[, 1, ]),
#   FP = c(metrics_results$Bayesian[, 2, ], metrics_results$GAM[, 2, ], metrics_results$Linear[, 2, ]),
#   TN = c(metrics_results$Bayesian[, 3, ], metrics_results$GAM[, 3, ], metrics_results$Linear[, 3, ]),
#   FN = c(metrics_results$Bayesian[, 4, ], metrics_results$GAM[, 4, ], metrics_results$Linear[, 4, ]),
#   Precision = c(metrics_results$Bayesian[, 5, ], metrics_results$GAM[, 5, ], metrics_results$Linear[, 5, ]),
#   Recall = c(metrics_results$Bayesian[, 6, ], metrics_results$GAM[, 6, ], metrics_results$Linear[, 6, ]),
#   F1_Score = c(metrics_results$Bayesian[, 7, ], metrics_results$GAM[, 7, ], metrics_results$Linear[, 7, ]),
#   Specificity = c(metrics_results$Bayesian[, 8, ], metrics_results$GAM[, 8, ], metrics_results$Linear[, 8, ]),
#   TDR = c(metrics_results$Bayesian[, 9, ], metrics_results$GAM[, 9, ], metrics_results$Linear[, 9, ])
# )
# 
# # Plot TDR results
# ggplot(metrics_df, aes(x = factor(Cutoff), y = TDR, color = Column)) +
#   geom_violin() +
#   labs(title = 'True Discovery Rate (TDR) Across Simulations',
#        x = 'Cutoff',
#        y = 'TDR',
#        color = 'Column') +
#   theme_minimal()
# 
# # Plot Precision results
# ggplot(metrics_df, aes(x = factor(Cutoff), y = Precision, color = Column)) +
#   geom_violin() +
#   labs(title = 'Precision Across Simulations',
#        x = 'Cutoff',
#        y = 'Precision',
#        color = 'Column') +
#   theme_minimal()
# 
# # Plot Recall results
# ggplot(metrics_df, aes(x = factor(Cutoff), y = Recall, color = Column)) +
#   geom_violin() +
#   labs(title = 'Recall Across Simulations',
#        x = 'Cutoff',
#        y = 'Recall',
#        color = 'Column') +
#   theme_minimal()
# 
# # Plot F1 Score results
# ggplot(metrics_df, aes(x = factor(Cutoff), y = F1_Score, color = Column)) +
#   geom_violin() +
#   labs(title = 'F1 Score Across Simulations',
#        x = 'Cutoff',
#        y = 'F1 Score',
#        color = 'Column') +
#   theme_minimal()
# 
# # Plot Specificity results
# ggplot(metrics_df, aes(x = factor(Cutoff), y = Specificity, color = Column)) +
#   geom_violin() +
#   labs(title = 'Specificity Across Simulations',
#        x = 'Cutoff',
#        y = 'Specificity',
#        color = 'Column') +
#   theme_minimal()


