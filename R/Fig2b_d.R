load("./Data/sim20_bias_res_0.1_4000.rda")
bias_res <- do.call(rbind, all_results_list)
bias_res <- as.data.frame(bias_res)
colnames(bias_res) <- c("Beta0_set", "Beta1_set", "Beta2_set", "Beta3_set", "Beta4_set", 
                        "Sigma2_1_set", "Sigma2_2_set", "Sigma2_3_set", "Sigma2_4_set", 
                        "Beta0_bayes", "Beta1_bayes", "Beta2_bayes", "Beta3_bayes", "Beta4_bayes",
                        "Sigma2_1_bayes", "Sigma2_2_bayes", "Sigma2_3_bayes", "Sigma2_4_bayes",
                        "Sigma2_gam", "Sigma2_lm", "Beta0_lm", "Beta1_lm", "Beta2_lm", "Beta3_lm", "Beta4_lm")
# Libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Step 1: Repeat Sigma2_gam and Sigma2_lm for all Sigma2 parameters
bias_res <- bias_res %>%
  mutate(Sigma2_1_gam = Sigma2_gam,
         Sigma2_2_gam = Sigma2_gam,
         Sigma2_3_gam = Sigma2_gam,
         Sigma2_4_gam = Sigma2_gam,
         Sigma2_1_polynomial = Sigma2_lm,
         Sigma2_2_polynomial = Sigma2_lm,
         Sigma2_3_polynomial = Sigma2_lm,
         Sigma2_4_polynomial = Sigma2_lm)

# Step 2: Calculate Relative Bias for Sigma2 (in percentages)
bias_res <- bias_res %>%
  mutate(
    RB_Sigma2_1_bayes = ((Sigma2_1_bayes - Sigma2_1_set) / Sigma2_1_set) * 100,
    RB_Sigma2_1_gam = ((Sigma2_1_gam - Sigma2_1_set) / Sigma2_1_set) * 100,
    RB_Sigma2_1_polynomial = ((Sigma2_1_polynomial - Sigma2_1_set) / Sigma2_1_set) * 100,
    
    RB_Sigma2_2_bayes = ((Sigma2_2_bayes - Sigma2_2_set) / Sigma2_2_set) * 100,
    RB_Sigma2_2_gam = ((Sigma2_2_gam - Sigma2_2_set) / Sigma2_2_set) * 100,
    RB_Sigma2_2_polynomial = ((Sigma2_2_polynomial - Sigma2_2_set) / Sigma2_2_set) * 100,
    
    RB_Sigma2_3_bayes = ((Sigma2_3_bayes - Sigma2_3_set) / Sigma2_3_set) * 100,
    RB_Sigma2_3_gam = ((Sigma2_3_gam - Sigma2_3_set) / Sigma2_3_set) * 100,
    RB_Sigma2_3_polynomial = ((Sigma2_3_polynomial - Sigma2_3_set) / Sigma2_3_set) * 100,
    
    RB_Sigma2_4_bayes = ((Sigma2_4_bayes - Sigma2_4_set) / Sigma2_4_set) * 100,
    RB_Sigma2_4_gam = ((Sigma2_4_gam - Sigma2_4_set) / Sigma2_4_set) * 100,
    RB_Sigma2_4_polynomial = ((Sigma2_4_polynomial - Sigma2_4_set) / Sigma2_4_set) * 100
  )

# Step 3: Calculate Relative Bias for Beta (in percentages)
bias_res <- bias_res %>%
  mutate(
    RB_Beta0_bayes = ((Beta0_bayes - Beta0_set) / Beta0_set) * 100,
    RB_Beta0_polynomial = ((Beta0_lm - Beta0_set) / Beta0_set) * 100,
    
    RB_Beta1_bayes = ((Beta1_bayes - Beta1_set) / Beta1_set) * 100,
    RB_Beta1_polynomial = ((Beta1_lm - Beta1_set) / Beta1_set) * 100,
    
    RB_Beta2_bayes = ((Beta2_bayes - Beta2_set) / Beta2_set) * 100,
    RB_Beta2_polynomial = ((Beta2_lm - Beta2_set) / Beta2_set) * 100,
    
    RB_Beta3_bayes = ((Beta3_bayes - Beta3_set) / Beta3_set) * 100,
    RB_Beta3_polynomial = ((Beta3_lm - Beta3_set) / Beta3_set) * 100,
    
    RB_Beta4_bayes = ((Beta4_bayes - Beta4_set) / Beta4_set) * 100,
    RB_Beta4_polynomial = ((Beta4_lm - Beta4_set) / Beta4_set) * 100
  )

# Add Squared Error for Sigma2
bias_res <- bias_res %>%
  mutate(
    SE_Sigma2_1_bayes = (Sigma2_1_bayes - Sigma2_1_set)^2,
    SE_Sigma2_1_gam = (Sigma2_1_gam - Sigma2_1_set)^2,
    SE_Sigma2_1_polynomial = (Sigma2_1_polynomial - Sigma2_1_set)^2,
    
    SE_Sigma2_2_bayes = (Sigma2_2_bayes - Sigma2_2_set)^2,
    SE_Sigma2_2_gam = (Sigma2_2_gam - Sigma2_2_set)^2,
    SE_Sigma2_2_polynomial = (Sigma2_2_polynomial - Sigma2_2_set)^2,
    
    SE_Sigma2_3_bayes = (Sigma2_3_bayes - Sigma2_3_set)^2,
    SE_Sigma2_3_gam = (Sigma2_3_gam - Sigma2_3_set)^2,
    SE_Sigma2_3_polynomial = (Sigma2_3_polynomial - Sigma2_3_set)^2,
    
    SE_Sigma2_4_bayes = (Sigma2_4_bayes - Sigma2_4_set)^2,
    SE_Sigma2_4_gam = (Sigma2_4_gam - Sigma2_4_set)^2,
    SE_Sigma2_4_polynomial = (Sigma2_4_polynomial - Sigma2_4_set)^2
  )

# Add Squared Error for Beta
bias_res <- bias_res %>%
  mutate(
    SE_Beta0_bayes = (Beta0_bayes - Beta0_set)^2,
    SE_Beta0_polynomial = (Beta0_lm - Beta0_set)^2,
    
    SE_Beta1_bayes = (Beta1_bayes - Beta1_set)^2,
    SE_Beta1_polynomial = (Beta1_lm - Beta1_set)^2,
    
    SE_Beta2_bayes = (Beta2_bayes - Beta2_set)^2,
    SE_Beta2_polynomial = (Beta2_lm - Beta2_set)^2,
    
    SE_Beta3_bayes = (Beta3_bayes - Beta3_set)^2,
    SE_Beta3_polynomial = (Beta3_lm - Beta3_set)^2,
    
    SE_Beta4_bayes = (Beta4_bayes - Beta4_set)^2,
    SE_Beta4_polynomial = (Beta4_lm - Beta4_set)^2
  )




# Reshape Sigma2 Relative Bias columns for plotting
rb_sigma2_long <- data.frame(
  Parameter = rep(c("Sigma2_1", "Sigma2_2", "Sigma2_3", "Sigma2_4"), each = 3 * nrow(bias_res)),
  Method = rep(rep(c("bayes", "gam", "polynomial"), each = nrow(bias_res)), times = 4),
  Relative_Bias = c(
    bias_res$RB_Sigma2_1_bayes, bias_res$RB_Sigma2_1_gam, bias_res$RB_Sigma2_1_polynomial,
    bias_res$RB_Sigma2_2_bayes, bias_res$RB_Sigma2_2_gam, bias_res$RB_Sigma2_2_polynomial,
    bias_res$RB_Sigma2_3_bayes, bias_res$RB_Sigma2_3_gam, bias_res$RB_Sigma2_3_polynomial,
    bias_res$RB_Sigma2_4_bayes, bias_res$RB_Sigma2_4_gam, bias_res$RB_Sigma2_4_polynomial
  )
)

custom_colors <- c("Bayesian" = "#A9082C", "GAM" = "#4B61A8", "Polynomial" = "#FDB76D")

# Update the Method names in rb_sigma2_long for better readability
rb_sigma2_long$Method <- recode(rb_sigma2_long$Method,
                                "bayes" = "Bayesian",
                                "gam" = "GAM",
                                "polynomial" = "Polynomial")

# Plot
pdf(file = "Sigma2_RB.pdf", width = 9)
ggplot(subset(rb_sigma2_long, Relative_Bias < 400), aes(x = Parameter, y = Relative_Bias, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors
  labs(
    title = "Relative Bias (%) of Sigma2 Parameters by Method",
    x = "Sigma2 Parameter",
    y = "Relative Bias (%)",
    fill = "Method"
  ) +
  theme(axis.text.x = element_text(size = 30, vjust = 0.85, color = "black"),
        axis.text.y = element_text(hjust = 1, size = 30, color = "black"),
        axis.title.x = element_text(size = 25, color = "black"),
        axis.title.y = element_text(size = 25, color = "black"),
        legend.position = "bottom",
        axis.line = element_line(colour = "black", size = 0.8), 
        axis.ticks = element_line(colour = "black", size = 0.8),
        axis.ticks.length = unit(0.3, "cm")) +
  theme_classic()

dev.off()
# Reshape Sigma2 Squared Error columns for plotting
se_sigma2_long <- data.frame(
  Parameter = rep(c("Sigma2_1", "Sigma2_2", "Sigma2_3", "Sigma2_4"), each = 3 * nrow(bias_res)),
  Method = rep(rep(c("bayes", "gam", "polynomial"), each = nrow(bias_res)), times = 4),
  Squared_Error = c(
    bias_res$SE_Sigma2_1_bayes, bias_res$SE_Sigma2_1_gam, bias_res$SE_Sigma2_1_polynomial,
    bias_res$SE_Sigma2_2_bayes, bias_res$SE_Sigma2_2_gam, bias_res$SE_Sigma2_2_polynomial,
    bias_res$SE_Sigma2_3_bayes, bias_res$SE_Sigma2_3_gam, bias_res$SE_Sigma2_3_polynomial,
    bias_res$SE_Sigma2_4_bayes, bias_res$SE_Sigma2_4_gam, bias_res$SE_Sigma2_4_polynomial
  )
)

# Update Method names for better readability
se_sigma2_long$Method <- recode(se_sigma2_long$Method,
                                "bayes" = "Bayesian",
                                "gam" = "GAM",
                                "polynomial" = "Polynomial")

# Plot
pdf(file = "Sigma2_SE.pdf", width = 5, height = 3)
ggplot(subset(se_sigma2_long, Squared_Error < 2), aes(x = Parameter, y = Squared_Error, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +  # Remove outliers
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors
  labs(
    title = "Squared Error of Sigma2 Parameters by Method",
    x = "Sigma2 Parameter",
    y = "Squared Error",
    fill = "Method"
  ) +
  theme(axis.text.x = element_text(size = 50, vjust = 0.85, color = "black"),
        axis.text.y = element_text(hjust = 1, size = 50, color = "black"),
        axis.title.x = element_text(size = 45, color = "black"),
        axis.title.y = element_text(size = 45, color = "black"),
        legend.position = "bottom",
        axis.line = element_line(colour = "black", size = 0.8), 
        axis.ticks = element_line(colour = "black", size = 0.8),
        axis.ticks.length = unit(0.5, "cm")) +
  theme_classic()
dev.off()
# Reshape Beta Squared Error columns for plotting
se_beta_long <- data.frame(
  Parameter = rep(c("Beta0", "Beta1", "Beta2", "Beta3", "Beta4"), each = 2 * nrow(bias_res)),
  Method = rep(rep(c("bayes", "polynomial"), each = nrow(bias_res)), times = 5),
  Squared_Error = c(
    bias_res$SE_Beta0_bayes, bias_res$SE_Beta0_polynomial,
    bias_res$SE_Beta1_bayes, bias_res$SE_Beta1_polynomial,
    bias_res$SE_Beta2_bayes, bias_res$SE_Beta2_polynomial,
    bias_res$SE_Beta3_bayes, bias_res$SE_Beta3_polynomial,
    bias_res$SE_Beta4_bayes, bias_res$SE_Beta4_polynomial
  )
)

# Update Method names for better readability
se_beta_long$Method <- recode(se_beta_long$Method,
                              "bayes" = "Bayesian",
                              "polynomial" = "Polynomial")

# Plot
pdf(file = "Beta_SE.pdf", width = 5, height = 3)
ggplot(subset(se_beta_long, Squared_Error < 30), aes(x = Parameter, y = Squared_Error, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +  # Remove outliers
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors
  labs(
    title = "Squared Error of Beta Parameters by Method",
    x = "Beta Parameter",
    y = "Squared Error",
    fill = "Method"
  ) +
  theme(axis.text.x = element_text(size = 50, vjust = 0.85, color = "black"),
        axis.text.y = element_text(hjust = 1, size = 50, color = "black"),
        axis.title.x = element_text(size = 45, color = "black"),
        axis.title.y = element_text(size = 45, color = "black"),
        legend.position = "bottom",
        axis.line = element_line(colour = "black", size = 0.8), 
        axis.ticks = element_line(colour = "black", size = 0.8),
        axis.ticks.length = unit(0.5, "cm")) +
  theme_classic()
dev.off()

# Calculate relative bias for Bayesian method
relative_bias_bayes <- data.frame(
  beta_0 = ((bias_res$Beta0_bayes - bias_res$Beta0_set) / bias_res$Beta0_set) * 100,
  beta_1 = ((bias_res$Beta1_bayes - bias_res$Beta1_set) / bias_res$Beta1_set) * 100,
  beta_2 = ((bias_res$Beta2_bayes - bias_res$Beta2_set) / bias_res$Beta2_set) * 100,
  beta_3 = ((bias_res$Beta3_bayes - bias_res$Beta3_set) / bias_res$Beta3_set) * 100,
  beta_4 = ((bias_res$Beta4_bayes - bias_res$Beta4_set) / bias_res$Beta4_set) * 100,
  sigma2_1 = ((bias_res$Sigma2_1_bayes - bias_res$Sigma2_1_set) / bias_res$Sigma2_1_set) * 100,
  sigma2_2 = ((bias_res$Sigma2_2_bayes - bias_res$Sigma2_2_set) / bias_res$Sigma2_2_set) * 100,
  sigma2_3 = ((bias_res$Sigma2_3_bayes - bias_res$Sigma2_3_set) / bias_res$Sigma2_3_set) * 100,
  sigma2_4 = ((bias_res$Sigma2_4_bayes - bias_res$Sigma2_4_set) / bias_res$Sigma2_4_set) * 100
)

# Check the dimensions of the resulting data frame
print(dim(relative_bias_bayes))  # Should return 10000 x 9

df <- as.data.frame(bias_res[,c("RB_Beta0_bayes", "RB_Beta1_bayes", "RB_Beta2_bayes", "RB_Beta3_bayes", "RB_Beta4_bayes",
                                "RB_Sigma2_1_bayes", "RB_Sigma2_2_bayes", "RB_Sigma2_3_bayes", "RB_Sigma2_4_bayes")])
colnames(df) <- c("beta_0", "beta_1", "beta_2", "beta_3", "beta_4", 
                  "sigma2_1", "sigma2_2", "sigma2_3", "sigma2_4")

# Convert data to long format and multiply by 100 for percentage
df_long <- df %>%
  mutate(across(everything(), ~ . * 1)) %>%  # Scale bias values to percentages
  pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Bias")

# Create the boxplot
ggplot(subset(df_long, Bias <= 100), aes(x = Parameter, y = Bias)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Bias of Estimated Parameters",
       x = "Parameter",
       y = "Bias (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_classic()
##################################

library(ggridges)


df_long <- df_long %>%
  mutate(
    Parameter = factor(Parameter, 
                       levels = rev(c("beta_0", "beta_1", "beta_2", "beta_3", "beta_4", 
                                      "sigma2_1", "sigma2_2", "sigma2_3", "sigma2_4")),
                       labels = rev(c(expression(beta[0]), expression(beta[1]), expression(beta[2]), 
                                      expression(beta[3]), expression(beta[4]), 
                                      expression(sigma^2 * "[1]"), expression(sigma^2 * "[2]"), 
                                      expression(sigma^2 * "[3]"), expression(sigma^2 * "[4]")))
    )
  )

# Create the ridgeline plot with visual clipping at x = 0 and custom colors
pdf("ridge_bias.pdf", width = 9, height = 6)
ggplot(subset(df_long, abs(Bias) <= 50), aes(x = Bias, y = Parameter, fill = Parameter)) +
  geom_density_ridges(stat = "binline", binwidth=1,
                      draw_baseline = F) + # Adjust bandwidth
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") + # Add vertical line at x = 0
  scale_x_continuous(breaks = c(-50, -30, -20, -10, 0, 10, 20, 30, 50)) + # Set custom x-axis ticks
  scale_y_discrete(
    labels = rev(c(expression(beta[0]), expression(beta[1]), expression(beta[2]), 
                   expression(beta[3]), expression(beta[4]), 
                   expression(sigma[1]^2), expression(sigma[2]^2), 
                   expression(sigma[3]^2), expression(sigma[4]^2)))
  ) + # Set math labels
  coord_cartesian(xlim = c(-50, NA)) + # Clip the x-axis visually at 0
  theme_classic() +
  labs(title = "Relative Bias of Parameters",
       x = "Absolute Bias (%)",
       y = "Parameter") +
  theme(
    plot.title = element_text(size = 0, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    axis.line = element_line(colour = "black", size = 0.8), 
    axis.ticks = element_line(colour = "black", size = 0.8),
    axis.ticks.length = unit(0.2, "cm")
  ) +
  scale_fill_manual(
    values = c(rep("seagreen", 4), rep("peru", 5)) 
  ) +
  guides(fill = "none") # Remove the legend
dev.off()
