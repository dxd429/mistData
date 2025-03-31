poly_function <- function(x, coeffs) {
  return(1/(1+exp(-coeffs[5]*x^4 - coeffs[4]*x^3 - coeffs[3]*x^2 - coeffs[2]*x - coeffs[1])))
}
library(SingleCellExperiment)
load('inter_Beta_sigma_all_0.025.rda')
load('m_reduced_clean.rda')
names(beta_mu_mean)<- rownames(m_reduced_clean)

m_reduced<- assay(genebody_SCE, "methylation")/assay(genebody_SCE, "coverage")

genenames <- c("ENSMUSG00000002033", "ENSMUSG00000087166", "ENSMUSG00000079042", "ENSMUSG00000087166", "ENSMUSG00000053914", "ENSMUSG00000029635", "ENSMUSG00000021953")

for ( i in 1:length(genenames)){
  genename<- genenames[i]
  
  dat <- m_reduced[genename,]
  ptime<- colData(genebody_SCE)$monocle3_level_pseudotime
  #remove 0 total count
  ind_remove<- is.nan(dat)
  dat<- dat[!ind_remove]
  ptime<- ptime[!ind_remove]
  ###correction
  dat[dat == 0]<- runif(length(dat[dat == 0]),0.00001, 0.01)
  dat[dat == 1]<- runif(length(dat[dat == 1]),0.95, 0.99999)
  df <- data.frame(Pseudotime = ptime, Methylation = dat)
  #df <- df[df$Methylation >= 0.05, ]
  beta_mu_mean1<- beta_mu_mean[[genename]][1:5]
  library(ggplot2)
  pdf(file = paste0("trend1_", genename, ".pdf"), width = 10, height = 5)
print(ggplot(df, aes(x = Pseudotime, y = Methylation, color = Methylation)) +
    geom_point(size=3) +
    geom_line(aes(y = poly_function(Pseudotime, beta_mu_mean1)), color = "darkgreen", linewidth = 1.2) +
    geom_vline(xintercept = 0.25, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0.75, linetype = "dashed", color = "black") +
    labs(title = paste0(genename), x = "Pseudotime", y = "Methylation Level") +
    scale_color_gradient(low="blue", high="red", limits = c(0,1)) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 25, face = "bold"),
      axis.title.x = element_text(size = 25),
      axis.title.y = element_text(size = 25),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20)
    ))
dev.off()
}

