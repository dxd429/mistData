sample_cells <- function(Num_cell, Num_t) {
  w <- 1 + 2 * rexp(Num_t, rate = 0.5)
  p <- w / sum(w)
  as.vector(rmultinom(1, Num_cell, prob = p))
}

calculate_integral <- function(params_A, params_B = c(0, 0, 0, 0, 0)) {
  t_points <- seq(0, 1, length.out = 100)  # Generate 100 points between 0 and 1
  differences <- numeric(length(t_points))  # Initialize a vector to store differences
  coef_A <- params_A[-1]  # Exclude intercept from curve A
  coef_B <- params_B[-1]  # Exclude intercept from curve B
  
  for (i in 1:length(t_points)) {
    t <- t_points[i]
    z_t <- c(t, t^2, t^3, t^4)  # Polynomial terms, exclude '1' for the intercept term
    
    # Compute inverse logit transformed values without intercepts
    curve_A <- 1 / (1 + exp(-(sum(z_t * coef_A) + params_A[1])))
    curve_B <- 1 / (1 + exp(-(sum(z_t * coef_B) + params_B[1])))
    
    # Store the difference between the curves
    differences[i] <- curve_A - curve_B
  }
  
  # Calculate the optimal shift C
  C_optimal <- mean(differences)
  
  # Adjust the differences by subtracting C
  adjusted_differences <- abs(differences - C_optimal)
  # Compute the integral using the trapezoidal rule
  integral <- sum((adjusted_differences[-1] + adjusted_differences[-length(adjusted_differences)]) / 2 * diff(t_points))
  
  return(integral)  # Return the computed integral
}
####
calculate_integral_gam <- function(gam_model1, gam_model2 = NULL) {
  # Generate a sequence of values from 0 to 1
  t_points <- seq(0, 1, length.out = 100)
  
  # Predict using the first GAM model on these values
  predictions1 <- predict(gam_model1, newdata = data.frame(time = t_points), type = "response")
  predictions1 <- 1 / (1 + exp(-predictions1))
  
  # Check if a second GAM model is provided or if we use a constant line
  if (is.null(gam_model2)) {
    # Use a constant line for the predictions of the second model
    predictions2 <- rep(0.5, length(t_points))
  } else {
    # Predict using the second GAM model on these values
    predictions2 <- predict(gam_model2, newdata = data.frame(time = t_points), type = "response")
    predictions2 <- 1 / (1 + exp(-predictions2))
  }
  
  # Calculate differences between the two curves
  differences <- predictions1 - predictions2
  
  # Calculate the optimal shift C
  C_optimal <- mean(differences)
  
  # Adjust the differences by subtracting C
  adjusted_differences <- abs(differences - C_optimal)
  
  # Compute the integral using the trapezoidal rule
  integral <- sum((adjusted_differences[-1] + adjusted_differences[-length(adjusted_differences)]) / 2 * diff(t_points))
  
  return(integral)  # Return the computed integral and optimal C
}


########
T_gamma<- function(shape, scale, lower, upper) {
  p_lower <- pgamma(lower, shape = shape, scale = scale)
  p_upper <- pgamma(upper, shape = shape, scale = scale)
  smallvalue = 1e-08
  if (((p_lower > 1 - smallvalue) & (p_upper > 1 - smallvalue)) | 
      ((p_lower < smallvalue) & (p_upper < smallvalue))) {
    return(lower* ((p_lower > 1 - smallvalue) & (p_upper  > 1 - smallvalue)) + 
             upper * ((p_lower < smallvalue) & (p_upper  < smallvalue)))
  }else{
    random_p <- runif(1, p_lower, p_upper)
    
    return(qgamma(random_p, shape, scale = scale))
  }
  
  
}
####
# Function to run Bayesian estimation
run_bayesian_estimation_sigma <- function(u, poly_d = 4, poly_scale = 0.025, stage_num = 4, Num_t = 40, Num_sim = 5000) {
  library(MCMCpack)
  # Bayesian beta estimation (need specifics of the Bayesian model setup)
  ###only keep unique methylation values
  #dat_unique_list<- lapply(u, unique)
  #u<- dat_unique_list
  n_t<- lengths(u)
  
  #### Set initial values
  z_t<- poly(poly_scale * (1:Num_t), poly_d, raw = T, simple = T) # Generate polynomial base
  beta_0 = rnorm(1,mean = 0,sd=3)
  lambda2_0 = 1
  tau2_0 = 1
  
  sigma2<- rep(2, stage_num) 
  #### Matrices to store values in iterations
  
  u_eachT<- unlist(lapply(u, sum))
  ###compute shape parameters
  stages<- (seq(Num_t) - 1) %/% (Num_t/stage_num) + 1
  shape_ig<- tapply(n_t, stages, sum)/2
  lambda2<- rep(lambda2_0, poly_d)
  tau2<- tau2_0
  beta_mu_m<- matrix(nrow = Num_sim, ncol = poly_d)
  #beta_mu_m<- matrix(nrow = Num_sim, ncol = poly_d)
  #lambda2_m<- matrix(nrow = Num_sim, ncol = poly_d)
  #tau2_v<- NULL
  sigma2_m<- matrix(nrow = Num_sim, ncol = stage_num)
  beta0_v<- NULL
  for(i in 1:Num_sim){
    #print(i)
    u_eachT_beta0<- unlist(lapply(u, function(x) (sum(x - beta_0))))
    s1<- lapply(u, function(x) (x - beta_0))
    v_beta_mu<- chol2inv(chol(t(z_t) %*% diag(n_t) %*% diag(1/rep(sigma2, each = Num_t/stage_num)) %*% z_t +
                                diag(1/(lambda2 * tau2)))) ## V_beta_mu
    
    m_beta_mu<- v_beta_mu %*% colSums(diag(u_eachT_beta0) %*% (diag(1/rep(sigma2, each = Num_t/stage_num)) %*% z_t)) ## m_beta_mu
    
    # v_beta_phi<- chol2inv(chol(t(z_t) %*% diag(1/(rep(sigma2, 50))) %*% z_t + 
    #                              diag(1/(v2 * w2)))) ## V_beta_phi
    # 
    # m_beta_phi<- v_beta_phi %*% colSums((1/sigma2) * (diag(log(phi_t)) %*% z_t)) ## m_beta_phi
    
    ###Beta_u, beta_phi
    # beta_mu<- m_beta_mu+t(Rfast::cholesky(v_beta_mu))%*%rnorm(poly_d)
    # beta_phi<- m_beta_phi+t(Rfast::cholesky(v_beta_phi))%*%rnorm(poly_d)
    # 
    beta_mu<- as.vector(mvtnorm::rmvnorm(1, m_beta_mu, (t(v_beta_mu) + v_beta_mu)/2))
    #beta_phi<- as.vector(mvtnorm::rmvnorm(1, m_beta_phi, (t(v_beta_phi) + v_beta_phi)/2))
    
    ###compute shape parameters
    
    #shapes[shapes < 0]<- 0
    ###compute scale parameters
    #s1<- u
    s2<- z_t %*% beta_mu
    s3<- mapply(function(s1, s2) s1 - s2, s1, s2, SIMPLIFY = FALSE)
    s4<- sapply(s3, function(x) sum((x^2)/2))
    scale_ig<- tapply(s4, stages, sum)
    
    sigma2<- MCMCpack::rinvgamma(4, shape = shape_ig, scale = scale_ig)
    
    ###horseshoe for lambda
    eta<- 1/lambda2
    upsi<- runif(poly_d, 0, 1/(1+eta))
    ub<- (1-upsi)/upsi
    Fub<- 1 - exp(-(beta_mu^2/(2 * tau2)) * ub) 
    up<- runif(poly_d, 0,Fub) 
    eta<- -log(1-up)/(beta_mu^2/(2 * tau2))
    lambda2<- 1/eta
    
    ###horseshoe for tau
    et<-  1/tau2 
    utau<- runif(1, 0, 1/(1+et))
    ubt<- (1-utau)/utau
    scale_tau<- sum((beta_mu^2/(2 * lambda2)))
    shape_tau<- (poly_d+1)/2
    upper_tau<- ubt
    et<- T_gamma(shape = shape_tau, scale = 1/scale_tau, lower = 0, upper = upper_tau)
    tau2<- 1/et
    ###beta0
    v_beta0 <- 1/((1/9) + as.numeric(n_t %*% (1/rep(sigma2, each = Num_t/stage_num))))
    m_beta0 <- v_beta0 * as.numeric(t(u_eachT - s2) %*% (1/rep(sigma2, each = Num_t/stage_num))) ## m_beta_mu
    beta_0 <- rnorm(1,mean = m_beta0, sd = sqrt(v_beta0))
    ###Record Betas
    ###Record Betas
    beta_mu_m[i,]<- m_beta_mu
    sigma2_m[i,]<- sigma2
    beta0_v<- c(beta0_v, m_beta0)
    
  }
  #########Bayesian Betas
  beta_mu_mean<- apply(beta_mu_m[round(0.8*Num_sim):Num_sim,], 2, mean)
  sigma2_mean<- apply(sigma2_m[round(0.8*Num_sim):Num_sim,], 2, mean)
  beta0_mean<- mean(beta0_v[round(0.8*Num_sim):Num_sim])
  rm(beta_mu_m, sigma2_m, s1, s2, s3, z_t, beta0_v)
  return(c(beta0_mean, beta_mu_mean, sigma2_mean))
}

run_bayesian_estimation <- function(u, poly_d = 4, poly_scale = 0.025, stage_num = 4, Num_t = 40, Num_sim = 5000) {
  library(MCMCpack)
  # Bayesian beta estimation (need specifics of the Bayesian model setup)
  ###only keep unique methylation values
  #dat_unique_list<- lapply(u, unique)
  #u<- dat_unique_list
  n_t<- lengths(u)
  
  #### Set initial values
  z_t<- poly(poly_scale * (1:Num_t), poly_d, raw = T, simple = T) # Generate polynomial base
  beta_0 = rnorm(1,mean = 0,sd=3)
  lambda2_0 = 1
  tau2_0 = 1
  
  sigma2<- rep(2, stage_num) 
  #### Matrices to store values in iterations
  
  u_eachT<- unlist(lapply(u, sum))
  ###compute shape parameters
  stages<- (seq(Num_t) - 1) %/% (Num_t/stage_num) + 1
  shape_ig<- tapply(n_t, stages, sum)/2
  lambda2<- rep(lambda2_0, poly_d)
  tau2<- tau2_0
  beta_mu_m<- matrix(nrow = Num_sim, ncol = poly_d)
  #beta_mu_m<- matrix(nrow = Num_sim, ncol = poly_d)
  #lambda2_m<- matrix(nrow = Num_sim, ncol = poly_d)
  #tau2_v<- NULL
  sigma2_m<- matrix(nrow = Num_sim, ncol = stage_num)
  beta0_v<- NULL
  for(i in 1:Num_sim){
    #print(i)
    u_eachT_beta0<- unlist(lapply(u, function(x) (sum(x - beta_0))))
    s1<- lapply(u, function(x) (x - beta_0))
    v_beta_mu<- chol2inv(chol(t(z_t) %*% diag(n_t) %*% diag(1/rep(sigma2, each = Num_t/stage_num)) %*% z_t +
                                diag(1/(lambda2 * tau2)))) ## V_beta_mu
    
    m_beta_mu<- v_beta_mu %*% colSums(diag(u_eachT_beta0) %*% (diag(1/rep(sigma2, each = Num_t/stage_num)) %*% z_t)) ## m_beta_mu
    
    # v_beta_phi<- chol2inv(chol(t(z_t) %*% diag(1/(rep(sigma2, 50))) %*% z_t + 
    #                              diag(1/(v2 * w2)))) ## V_beta_phi
    # 
    # m_beta_phi<- v_beta_phi %*% colSums((1/sigma2) * (diag(log(phi_t)) %*% z_t)) ## m_beta_phi
    
    ###Beta_u, beta_phi
    # beta_mu<- m_beta_mu+t(Rfast::cholesky(v_beta_mu))%*%rnorm(poly_d)
    # beta_phi<- m_beta_phi+t(Rfast::cholesky(v_beta_phi))%*%rnorm(poly_d)
    # 
    beta_mu<- as.vector(mvtnorm::rmvnorm(1, m_beta_mu, (t(v_beta_mu) + v_beta_mu)/2))
    #beta_phi<- as.vector(mvtnorm::rmvnorm(1, m_beta_phi, (t(v_beta_phi) + v_beta_phi)/2))
    
    ###compute shape parameters
    
    #shapes[shapes < 0]<- 0
    ###compute scale parameters
    #s1<- u
    s2<- z_t %*% beta_mu
    s3<- mapply(function(s1, s2) s1 - s2, s1, s2, SIMPLIFY = FALSE)
    s4<- sapply(s3, function(x) sum((x^2)/2))
    scale_ig<- tapply(s4, stages, sum)
    
    sigma2<- MCMCpack::rinvgamma(4, shape = shape_ig, scale = scale_ig)
    
    ###horseshoe for lambda
    eta<- 1/lambda2
    upsi<- runif(poly_d, 0, 1/(1+eta))
    ub<- (1-upsi)/upsi
    Fub<- 1 - exp(-(beta_mu^2/(2 * tau2)) * ub) 
    up<- runif(poly_d, 0,Fub) 
    eta<- -log(1-up)/(beta_mu^2/(2 * tau2))
    lambda2<- 1/eta
    
    ###horseshoe for tau
    et<-  1/tau2 
    utau<- runif(1, 0, 1/(1+et))
    ubt<- (1-utau)/utau
    scale_tau<- sum((beta_mu^2/(2 * lambda2)))
    shape_tau<- (poly_d+1)/2
    upper_tau<- ubt
    et<- T_gamma(shape = shape_tau, scale = 1/scale_tau, lower = 0, upper = upper_tau)
    tau2<- 1/et
    ###beta0
    v_beta0 <- 1/((1/9) + as.numeric(n_t %*% (1/rep(sigma2, each = Num_t/stage_num))))
    m_beta0 <- v_beta0 * as.numeric(t(u_eachT - s2) %*% (1/rep(sigma2, each = Num_t/stage_num))) ## m_beta_mu
    beta_0 <- rnorm(1,mean = m_beta0, sd = sqrt(v_beta0))
    ###Record Betas
    ###Record Betas
    beta_mu_m[i,]<- m_beta_mu
    sigma2_m[i,]<- sigma2
    beta0_v<- c(beta0_v, m_beta0)
    
  }
  #########Bayesian Betas
  beta_mu_mean<- apply(beta_mu_m[round(0.8*Num_sim):Num_sim,], 2, mean)
  sigma2_mean<- apply(sigma2_m[round(0.8*Num_sim):Num_sim,], 2, mean)
  beta0_mean<- mean(beta0_v[round(0.8*Num_sim):Num_sim])
  rm(beta_mu_m, sigma2_m, s1, s2, s3, z_t, beta0_v)
  return(c(beta0_mean, beta_mu_mean))
}

######################
# Function to perform the entire analysis workflow
Bias_compare_top1000 <- function(betas_set, sigma2_set, u, poly_d = 4, poly_scale = 0.025, Num_t = 40) {
  #library(mgcv)
  Beta_bayes <- run_bayesian_estimation_sigma(u)
  
  bias_res <- 100 * (Beta_bayes - c(betas_set, sigma2_set))/(c(betas_set, sigma2_set))
  
  return(bias_res)
}
Bias_compare <- function(betas_set, sigma2_set, u, poly_d = 4, poly_scale = 0.025, Num_t = 40) {
  #library(mgcv)
  Beta_bayes <- run_bayesian_estimation_sigma(u)
  
  bias_res <- (Beta_bayes - c(betas_set, sigma2_set))
  
  return(c(betas_set, sigma2_set, Beta_bayes))
}

addi_coef <- function(u, poly_d = 4, poly_scale = 0.025, Num_t = 40) {

  data <- data.frame(response = unlist(u), time = 0.025 * (rep(1:length(u), sapply(u, length))))
  z_t<- poly(poly_scale * (1:Num_t), poly_d, raw = T, simple = T)
  polynomial_terms <- z_t[rep(1:nrow(z_t), sapply(u, length)), ]
  colnames(polynomial_terms) <- c("poly1", "poly2", "poly3", "poly4")
  data <- cbind(data, polynomial_terms)
  
  model_gam <- gam(response ~ s(time, k = 5, bs = 'cr'), 
                   knots = list(time = c(0, 0.25, 0.5, 0.75, 1)), 
                   data = data, family = gaussian(link = "identity"))
  gam_var <- model_gam$sig2
  
  model_gam_null <- gam(response ~ 1, data = data, family = gaussian(link = "identity"))
  
  
  model_lm <- lm(response ~ poly1 + poly2 + poly3 + poly4, data = data)
  beta_lm <- coef(model_lm)
  lm_summary <- summary(model_lm)
  lm_var <- lm_summary$sigma^2
  model_lm_null <- lm(response ~ 1, data = data)

  
  # Return the calculated integrals
  return(c(gam_var = gam_var, lm_var = lm_var, beta_lm))
}


Integral_compare <- function(u, poly_d = 4, poly_scale = 0.025, Num_t = 40) {
  #library(mgcv)
  Beta_bayes <- run_bayesian_estimation(u)
  
  # Fit the GAM model
  data <- data.frame(response = unlist(u), time = 0.025 * (rep(1:length(u), sapply(u, length))))
  z_t<- poly(poly_scale * (1:Num_t), poly_d, raw = T, simple = T)
  polynomial_terms <- z_t[rep(1:nrow(z_t), sapply(u, length)), ]
  colnames(polynomial_terms) <- c("poly1", "poly2", "poly3", "poly4")
  data <- cbind(data, polynomial_terms)

model_gam <- gam(response ~ s(time, k = 5, bs = 'cr'), 
                 knots = list(time = c(0, 0.25, 0.5, 0.75, 1)), 
                 data = data, family = gaussian(link = "identity"))


model_gam_null <- gam(response ~ 1, data = data, family = gaussian(link = "identity"))


lrt_gam <- anova(model_gam_null, model_gam, test = "Chisq")
p_value_gam <- lrt_gam$`Pr(>Chi)`[2] 

model_lm <- lm(response ~ poly1 + poly2 + poly3 + poly4, data = data)

model_lm_null <- lm(response ~ 1, data = data)

lrt_lm <- anova(model_lm_null, model_lm)
p_value_lm <- lrt_lm$`Pr(>F)`[2] 

  
  # Calculate integrals using the adjusted `calculate_integral` function
  integral_bayesian <- calculate_integral(Beta_bayes)
  #integral_gam <- calculate_integral_gam(model_gam)
  #integral_lm <- calculate_integral(coef(model_lm))
  
  # Return the calculated integrals
  c(integral_bayesian = integral_bayesian, p_value_gam = -log(p_value_gam), p_value_lm = -log(p_value_lm))
}
######################
# Function to perform the entire analysis workflow
Integral_compare_2group <- function(u1, u2, poly_d = 4, poly_scale = 0.025, Num_t = 40) {
  #library(mgcv)
  Beta_bayes1 <- run_bayesian_estimation(u1)
  Beta_bayes2 <- run_bayesian_estimation(u2)
  # Fit the GAM model
  data1 <- data.frame(response = unlist(u1), time = 0.025 * (rep(1:length(u1), sapply(u1, length))))
  z_t<- poly(poly_scale * (1:Num_t), poly_d, raw = T, simple = T)
  polynomial_terms <- z_t[rep(1:nrow(z_t), sapply(u1, length)), ]
  colnames(polynomial_terms) <- c("poly1", "poly2", "poly3", "poly4")
  data1 <- cbind(data1, polynomial_terms)
  data1$group <- factor("group1")
  
  data2 <- data.frame(response = unlist(u2), time = 0.025 * (rep(1:length(u2), sapply(u2, length))))
  z_t<- poly(poly_scale * (1:Num_t), poly_d, raw = T, simple = T)
  polynomial_terms <- z_t[rep(1:nrow(z_t), sapply(u2, length)), ]
  colnames(polynomial_terms) <- c("poly1", "poly2", "poly3", "poly4")
  data2 <- cbind(data2, polynomial_terms)
  data2$group <- factor("group2")
  combined_data <- rbind(data1, data2)
  model_gam_combined <- gam(response ~ s(time, by = group, k = 5, bs = 'cr') + group, 
                          knots = list(time = c(0, 0.25, 0.5, 0.75, 1)), 
                          data = combined_data, 
                          family = gaussian(link = "identity"))

  model_gam_null <- gam(response ~ s(time, k = 5, bs = 'cr') + group, 
                          knots = list(time = c(0, 0.25, 0.5, 0.75, 1)), 
                          data = combined_data, 
                          family = gaussian(link = "identity"))

  lrt_gam <- anova(model_gam_null, model_gam_combined, test = "Chisq")
  p_value_gam <- lrt_gam$`Pr(>Chi)`[2]

model_lm_combined <- lm(response ~ poly1 * group + poly2 * group + poly3 * group + poly4 * group, data = combined_data)

model_lm_null <- lm(response ~ poly1 + poly2 + poly3 + poly4 + group, data = combined_data)
lrt_lm <- anova(model_lm_null, model_lm_combined)
  p_value_lm <- lrt_lm$`Pr(>F)`[2]
  # Calculate integrals using the adjusted `calculate_integral` function
  integral_bayesian <- calculate_integral(Beta_bayes1, Beta_bayes2)
  #integral_gam <- calculate_integral_gam(model_gam1, model_gam2)
  #integral_lm <- calculate_integral(coef(model_lm1), coef(model_lm2))
  
  # Return the calculated integrals
  c(integral_bayesian = integral_bayesian, p_value_gam = -log(p_value_gam), p_value_lm = -log(p_value_lm))
}
##############Function to compute the matrices
# Load necessary libraries
library(ggplot2)
library(dplyr)

# Function to calculate evaluation metrics
calculate_metrics <- function(column_data, true_signals, total_rows) {

  top_indices <- order(column_data, decreasing = TRUE)[true_signals]
  TP <- length(intersect(top_indices, true_signals))
  FP <- length(setdiff(top_indices, true_signals))
  FN <- length(setdiff(true_signals, top_indices))
  TN <- total_rows - (TP + FP + FN)

  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  f1_score <- 2 * (precision * recall) / (precision + recall)
  specificity <- TN / (TN + FP)

  return(c(TP = TP, FP = FP, TN = TN, FN = FN, 
           precision = precision, recall = recall, 
           f1_score = f1_score, specificity = specificity))

}
# Function to calculate True Discovery Rate (TDR) at various cutoffs
calculate_tdr <- function(column_data, true_signals, cutoffs) {
  tdr_list <- numeric(length(cutoffs))
  
  for (i in seq_along(cutoffs)) {
    cutoff <- cutoffs[i]
    top_indices <- order(column_data, decreasing = TRUE)[1:cutoff]
    tp <- length(intersect(top_indices, true_signals))
    tdr <- tp / cutoff
    tdr_list[i] <- tdr
  }
  
  return(tdr_list)
}

# Function to calculate TPR and FPR for ROC
calculate_roc_data <- function(column_data, true_signals, total_rows, thresholds) {
  
  roc_data <- data.frame(threshold = integer(), TPR = numeric(), FPR = numeric(), stringsAsFactors = FALSE)
  
  for (threshold in thresholds) {
    
    # Get the top indices based on the current threshold
    top_indices <- order(column_data, decreasing = TRUE)[1:threshold]
    
    # Calculate True Positives (TP), False Positives (FP), False Negatives (FN), True Negatives (TN)
    TP <- length(intersect(top_indices, true_signals))  # True Positives
    FP <- length(setdiff(top_indices, true_signals))    # False Positives
    FN <- length(setdiff(true_signals, top_indices))    # False Negatives
    TN <- total_rows - (TP + FP + FN)                   # True Negatives
    
    # Calculate True Positive Rate (TPR) and False Positive Rate (FPR)
    TPR <- ifelse(TP + FN == 0, 0, TP / (TP + FN))  # TPR is the same as recall
    FPR <- ifelse(FP + TN == 0, 0, FP / (FP + TN))  # FPR
    
    # Store TPR and FPR at this threshold
    roc_data <- rbind(roc_data, data.frame(threshold = threshold, TPR = TPR, FPR = FPR))
  }
  
  return(roc_data)
}


interpolate_tpr <- function(roc_data, common_fpr) {
  approx(roc_data$FPR, roc_data$TPR, xout = common_fpr, rule = 2)$y  # Interpolate TPR at common FPR
}

calculate_pr_data <- function(column_data, true_signals, total_rows, thresholds) {
  
  pr_data <- data.frame(threshold = integer(), precision = numeric(), recall = numeric(), stringsAsFactors = FALSE)
  
  for (threshold in thresholds) {
    
    # Get the top indices based on the current threshold
    top_indices <- order(column_data, decreasing = TRUE)[1:threshold]
    
    # Calculate True Positives (TP), False Positives (FP), False Negatives (FN)
    TP <- length(intersect(top_indices, true_signals))  # True Positives
    FP <- length(setdiff(top_indices, true_signals))    # False Positives
    FN <- length(setdiff(true_signals, top_indices))    # False Negatives
    TN <- total_rows - (TP + FP + FN)                   # True Negatives, but not used for PR curve
    
    # Calculate Precision and Recall
    precision <- ifelse(TP + FP == 0, 0, TP / (TP + FP))  # Precision
    recall <- ifelse(TP + FN == 0, 0, TP / (TP + FN))     # Recall
    
    # Store Precision and Recall at this threshold
    pr_data <- rbind(pr_data, data.frame(threshold = threshold, precision = precision, recall = recall))
  }
  
  return(pr_data)
}

interpolate_precision <- function(pr_data, common_recall) {
  approx(pr_data$recall, pr_data$precision, xout = common_recall, rule = 2)$y  # Interpolate precision at common recall
}


