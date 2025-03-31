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
    betas_set <- beta_mu_sigmatop1800[[b]][1:5]
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
    betas_set <- beta_mu_sigmatop1800[[b]][1:5]
    mu_t <- z_t_0 %*% beta_mu_sigmatop1800[[b]][1:5] # Compute mu_t
    n_t1 <- sample_cells(Num_cell = Num_cell, Num_t = Num_t)
    u <- vector("list", length(n_t1)) 
    breaks <- seq(1, Num_t, by = Num_t/stage_num)
    for (k in seq_along(n_t1)) {
      j <- findInterval(k, breaks)
      u[[k]] <- rnorm(n_t1[k], mean = mu_t[k], sd = (sqrt(sigma2_set[j]))) 
    }
  }
  bias_res <- Bias_compare(betas_set, sigma2_set, u)
  addi_coef <- addi_coef(u)
  bias_ress <- c(bias_res, addi_coef)
  names(bias_ress) <- c("Beta0_set", "Beta1_set", "Beta2_set", "Beta3_set", "Beta4_set", 
                        "Sigma2_1_set", "Sigma2_2_set", "Sigma2_3_set", "Sigma2_4_set", 
                        "Beta0_bayes", "Beta1_bayes", "Beta2_bayes", "Beta3_bayes", "Beta4_bayes",
                        "Sigma2_1_bayes", "Sigma2_2_bayes", "Sigma2_3_bayes", "Sigma2_4_bayes",
                        "Sigma2_gam", "Sigma2_lm", "Beta0_lm", "Beta1_lm", "Beta2_lm", "Beta3_lm", "Beta4_lm")
  return(bias_ress)
}
############Simulate u
load('inter_Beta_sigma_all_0.025.rda')
load('m_reduced_clean.rda')# Function to check for Inf or NA
#####
Num_t <- 40
Num_cell <- 4000
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
n_gene <- 10000
dm_ratio <- 0.1


library(BiocParallel)
 
all_results_list <- vector("list", 20)  
for (i in 1:20) {
# Run the process in parallel
bias_results <- bplapply(1:n_gene, process_gene, n_gene = n_gene, z_t_0 = z_t_0, beta_mu_sigmatop1800 = beta_mu_sigmatop1800, Num_t = Num_t, stage_num = stage_num, Num_cell = Num_cell, dm_ratio = dm_ratio, BPPARAM = MulticoreParam())
  
# Combine results into a data matrix
bias_res <- do.call(rbind, bias_results)

all_results_list[[i]] <- bias_res

}
save(all_results_list, file = paste0("sim20_bias_res_", dm_ratio, '_', Num_cell, '.rda'))


