library(rms)
library(ggplot2)
library(dplyr)
library(mvtnorm)
library(tidyverse)
library(scales)
library(shiny)
if (!requireNamespace("shinyBS", quietly = TRUE)) {
  install.packages("shinyBS")
}
library(shinyBS)
library(readr)
library(shinyWidgets)
library(rsconnect)
library(DT)
library(pROC)
library(bslib)

find_intermediate_cor <- function(x, cor_obj, frac_positive, N){
  R <- matrix(data = c(1, x, x, 1), ncol = 2)
  set.seed(4)
  d <- as.data.frame(rmvnorm(N, mean = c(rep(0, nrow(R))), sigma = R))
  # Dichotomize
  d$V2 <- ifelse(d$V2 > sort(d$V2)[nrow(d)*(1-frac_positive)], 1, 0)
  # Minimize to the input correlation
  (as.numeric(cor(d)[1,2])-cor_obj)^2
}
Generate_testscores <- function(N, cor_tests, frac_positive, type, seed=4){
  # Our final matrix:
  R <- matrix(data = c(1, cor_tests, cor_tests, 1), ncol = 2)
  # Generate data
  set.seed(seed)
  d <- as.data.frame(rmvnorm(N, mean = c(rep(0, nrow(R))), sigma = R))
  colnames(d) <- c("Reference", "Newtest")
  # Clone data set and add chemotherapy:
  df <- rbind(d,d)
  df$Treatment <- c(rep(c(0,1), each = nrow(df)/2))
  
  if (type == "binary"){
    # Store cutoff value for estimation analytic solution
    cutoff <- sort(df$Newtest)[nrow(df)*(1-frac_positive)]
    df$Newtest <- ifelse(df$Newtest > cutoff, 1, 0)
    # Output
    return(list(cutoff, df))
  } else {
    # Output
    return(list(0, df))
  }
}
find_h0_surv <- function(x, obs, lp){
  Y <- obs
  Yhat <- 1-mean(exp(-x*exp(lp)))
  # Minimize
  -sum(Y*log(Yhat) + (1-Y)*log((1-Yhat)))
}
find_b0_logistic <- function(x, obs, lp){
  Y <- obs
  Yhat <- mean(exp(x + lp)/(1+exp(x + lp)))
  # Minimize
  -sum(Y*log(Yhat) + (1-Y)*log((1-Yhat)))
}
finding_b1_b2 <- function(event_rate, cor_est, AUC_pref, linkfunction, 
                          frac_t1, typet1, AUC_preft1, sens_preft1, 
                          spec_preft1, b1est, b2est){
  
  # Set counter to track iterations
  counter <- new.env()
  counter$iter_count <- 0
  # Finetune balancing precision and computation time
  counter$max_iter <- 15
  
  b1res = b1est
  b2res = b2est
  # Finetune k balancing precision and computation time
  k=25
  b1=b2=numeric(k)
  for (i in 1:k) {
    # Finetune N here balancing precision and computation time
    d <- Generate_testscores(N = 1500, cor_tests = cor_est, 
                             frac_positive = frac_t1, 
                             type = typet1, seed = i)[[2]]
    
    # Reset iteration counter:
    counter$iter_count <- 0
    # The optimization function:
    ress <- optim(par = c(b1res,b2res), fn = optimization_b0_b1_b2,
                  df = d, counter = counter, event_rate = event_rate, 
                  linkfunction = linkfunction, AUC_ref_pref = AUC_pref, 
                  AUC_preft1 = AUC_preft1, sens_preft1 = sens_preft1, 
                  spec_preft1 = spec_preft1, typet1 = typet1, i=i)
    # Store results:
    b1res = b1[i] = as.numeric(ress$par[1])
    b2res = b2[i] = as.numeric(ress$par[2])
    
    print(i)
    print(c(b1[i],b2[i]))
  }
  print(colMeans(cbind(b1,b2)))
  return(cbind(b1,b2))
}
optimization_b0_b1_b2 <- function(x, df, event_rate, AUC_ref_pref, AUC_preft1, sens_preft1, spec_preft1, typet1, counter, i, linkfunction){
  
  # Set iteration counter, to prevent using too many iterations:
  counter$iter_count <- counter$iter_count + 1
  if(counter$iter_count > counter$max_iter) return(1e9)
  if(any(x<=0)){x=x+5;print("error")}
  
  # Prepare data set:
  d = df
  # Find these coefficients:
  lp <- log(x[1])*d$Reference + log(x[2])*d$Newtest
  
  if (linkfunction == "logistic"){
    # Update b0:
    findb0 = optimize(find_b0_logistic, interval = c(-19,20), 
                      obs = event_rate, lp = lp)
    b0 = as.numeric(findb0$minimum)
    # Simulate Y
    probs_Y <- exp(b0+lp)/(1+exp(b0+lp))
    set.seed(i)
    d$Y <- ifelse(runif(n=nrow(d)) < probs_Y, 1, 0)
  } else if (linkfunction == "survival"){
    # Update h0:
    findh0 = optimize(find_h0_surv, interval = c(0.0001,20), 
                      obs = event_rate, lp = lp)
    h0 = as.numeric(findh0$minimum)
    # Generate outcome to find coefficients:
    hazs <- h0*exp(lp)
    set.seed(i)
    times <- rexp(n=nrow(d), rate = hazs)
    d$Y <- (times < 1)*1
  }
  ########################################
  # Find coeffcients with:
  # To find B1: Always use:
  AUC_est <- as.numeric(roc(d$Y, d$Reference, quiet = TRUE)$auc)
  est0 <- -sum(AUC_ref_pref * log(AUC_est) + (1 - AUC_ref_pref) * log(1 - AUC_est))
  
  # Depending on the type of test:
  if (typet1 == "binary") {
    # test1 sens/spec
    sens_est <- sum(d$Newtest == 1 & d$Y == 1) / sum(d$Y == 1)
    spec_est <- sum(d$Newtest == 0 & d$Y == 0) / sum(d$Y == 0)
    est1 <- -sum(sens_preft1 * log(sens_est) + (1 - sens_preft1) * log(1 - sens_est))
    est2 <- -sum(spec_preft1 * log(spec_est) + (1 - spec_preft1) * log(1 - spec_est))
    # Minimize
    total <- sum(c(est0, est1, est2))
  } else {
    # test1 AUC
    AUC_estt1 <- as.numeric(roc(d$Y, d$Newtest, quiet = TRUE)$auc)
    est1 <- -sum(AUC_preft1 * log(AUC_estt1) + (1 - AUC_preft1) * log(1 - AUC_estt1))
    # Minimize
    total <- sum(c(est0, est1))
  }
  # Minimize
  total
}
probX2givenX1 <- function(X_values, thres_v2, rho){
  mu_x = 0;mu_y = 0;var_x = 1;var_y = 1
  musY_X = mu_y + rho * (sqrt(var_y)/sqrt(var_x)) * (X_values - mu_x)
  varY_X = var_y*(1-rho^2)
  # Compute probabilities for all x in the vector
  probY_t = 1-pnorm(thres_v2, mean = musY_X, sd = sqrt(varY_X))
  return(probY_t)
}
ref_model_probs_log_bin <- function(d, cutoff, rho, b0, b1, b2, b3){
  # P(Y|X_2=0,X1)
  lp <- b0 + log(b1)*d$Reference + log(b2)*0 + log(b3)*d$Treatment
  probX2_0 <- exp(lp)/(1+exp(lp))
  # P(Y|X_2=1,X1)
  lp <- b0 + log(b1)*d$Reference + log(b2)*1 + log(b3)*d$Treatment
  probX2_1 <- exp(lp)/(1+exp(lp))
  # P(X_2=1|X_1)
  x2_1_givenx1 <- probX2givenX1(d$Reference, thres_v2 = cutoff, rho = rho)
  # P(X_2=0|X_1)
  x2_0_givenx1 <- 1-x2_1_givenx1
  
  # Analytical solution reference model:
  probs_ref <- probX2_0*x2_0_givenx1 + probX2_1*x2_1_givenx1
  
  # Output
  return(probs_ref)
}
ref_model_probs_surv_bin <- function(cutoff, rho, d, h0, b1, b2, b3){
  # P(Y|X_2=0,X1)
  probX2_0 <- 1-exp(-h0*exp(log(b1)*d$Reference + log(b2)*0 + log(b3)*d$Treatment))
  # P(Y|X_2=1,X1)
  probX2_1 <- 1-exp(-h0*exp(log(b1)*d$Reference + log(b2)*1 + log(b3)*d$Treatment))
  # P(X_2=1|X_1)
  x2_1_givenx1 <- probX2givenX1(d$Reference, thres_v2 = cutoff, rho = rho)
  # P(X_2=0|X_1)
  x2_0_givenx1 <- 1-x2_1_givenx1
  
  # Analytical solution reference model:
  probs_ref <- probX2_0*x2_0_givenx1 + probX2_1*x2_1_givenx1
  
  # Output
  return(probs_ref)
}
ref_model_probs_log_cont <- function(d, rho, b0, b1, b2, b3){
  # Define bins and midpoints
  values <- seq(-5, 5, by = 0.5)
  test_values <- (values[-1] + values[-length(values)]) / 2
  
  # Conditional mean & sd for all patients
  muY_X <- rho*d$Reference
  sdY_X <- sqrt(1 - rho^2)
  
  # Compute probability of each bin for each patient
  probY_t_mat <- t(sapply(muY_X, function(mu) {
    pnorm(values[-1], mean = mu, sd = sdY_X) - 
      pnorm(values[-length(values)], mean = mu, sd = sdY_X)}))
  
  # Get the risk for every mid point of the bins
  lp_mat <- outer(b0 + log(b1) * d$Reference + log(b3) * d$Treatment, 
                  log(b2) * test_values, "+")
  ests_mat <- exp(lp_mat)/(1+exp(lp_mat))
  # Weighted sum over bins
  probs_ref <- rowSums(ests_mat * probY_t_mat)
  # Output
  return(probs_ref)
}
ref_model_probs_surv_cont <- function(d, rho, h0, b1, b2, b3){
  values <- seq(-5, 5, by = 0.5)
  test_values <- (values[-1] + values[-length(values)]) / 2
  
  # Conditional mean & sd for all patients
  muY_X <- rho*d$Reference
  sdY_X <- sqrt(1 - rho^2)
  
  # Compute probability of each bin for each patient
  probY_t_mat <- t(sapply(muY_X, function(mu) {
    pnorm(values[-1], mean = mu, sd = sdY_X) - 
      pnorm(values[-length(values)], mean = mu, sd = sdY_X)}))
  
  # Get the risk for every mid point of the bins
  lp_mat <- outer(log(b1) * d$Reference + log(b3) * d$Treatment, 
                  log(b2) * test_values, "+")
  ests_mat <- 1-exp(-h0*exp(lp_mat))
  # Weighted sum over bins
  probs_ref <- rowSums(ests_mat * probY_t_mat)
  # Output
  return(probs_ref)
}
generate_risks <- function(dataset, coefs, linkfunction, typet1, cutoff, rho){
  d <- dataset
  b0 <- coefs[1]
  b1 <- coefs[2]
  b2 <- coefs[3]
  b3 <- coefs[4]
  
  # Get the correct outcome and model fits depending on link function:
  if (linkfunction == "logistic"){
    # Set the linear predictors to generate probabilities:
    lp <- b0 + log(b1)*d$Reference + log(b2)*d$Newtest + log(b3)*d$Treatment
    d$risk_bothtest <- exp(lp)/(1+exp(lp))
    set.seed(4)
    d$Y <- ifelse(runif(n=nrow(d)) < d$risk_bothtest, 1, 0)
    
    if (typet1 == "binary"){
      # Generate risks based on reference model only:
      d$risk_reftest <- ref_model_probs_log_bin(cutoff = cutoff,
                                                rho = rho, 
                                                d = d, 
                                                b0 = b0, 
                                                b1 = b1, 
                                                b2 = b2, 
                                                b3 = b3)
    } else {
      d$risk_reftest <- ref_model_probs_log_cont(rho = rho, 
                                                 d = d, 
                                                 b0 = b0, 
                                                 b1 = b1, 
                                                 b2 = b2, 
                                                 b3 = b3)
    }
    
    
  } else if (linkfunction == "survival"){
    h0=b0
    # Make LP:
    lp <- log(b1)*d$Reference + log(b2)*d$Newtest + log(b3)*d$Treatment
    d$risk_bothtest <- 1-exp(-h0*exp(lp))
    # Set the linear predictors to generate probabilities:
    hazs <- h0*exp(lp)
    set.seed(4)
    times <- rexp(n=nrow(d), rate = hazs)
    d$Y <- (times < 1)*1
    
    if (typet1 == "binary"){
      # Generate risks based on reference model only:
      d$risk_reftest <- ref_model_probs_surv_bin(cutoff = cutoff,
                                                 rho = rho, 
                                                 d = d, 
                                                 h0 = h0, 
                                                 b1 = b1, 
                                                 b2 = b2, 
                                                 b3 = b3)
    } else {
      d$risk_reftest <- ref_model_probs_surv_cont(rho = rho, 
                                                  d = d, 
                                                  h0 = h0, 
                                                  b1 = b1, 
                                                  b2 = b2, 
                                                  b3 = b3)
    }
  }
  
  # Make data set with only the unique individuals:
  d2 <- d[d$Treatment == 0,]
  d2$risk_bothtestchemo <- d$risk_bothtest[d$Treatment==1]
  d2$risk_reftestchemo <- d$risk_reftest[d$Treatment==1]
  
  # Generate expected treatment benefit, used to guide treatment decisions
  d2$expben_bothtest <- d2$risk_bothtest - d2$risk_bothtestchemo
  d2$expben_reftest <- d2$risk_reftest - d2$risk_reftestchemo
  
  # Format Outcomes: Y under treatment and no treatment
  colnames(d2)[which(colnames(d2)=="Y")] <- "Y_tr0"
  d2$Y_tr1 <- d$Y[d$Treatment == 1]
  
  return(as.data.frame(d2))
}
Opt_function_all <- function(Nunique, linkfunction, event_rate, treatmenteffect,
                             AUCrefmodel, type_biomarker, AUC_biomarker=0,
                             sens_biomarker=0, spec_biomarker=0, cortestbiomarker){
  
  # Store the correct values:
  N <- 25000 # Nunique is used later
  linkfunction <- linkfunction # surv or logistic
  ev_rate <- event_rate
  b3 <- treatmenteffect   
  AUC_ref <- AUCrefmodel
  # Types
  typemarker <- type_biomarker
  # Test 1 inputs
  AUC_t1  <- AUC_biomarker
  sens_t1 <- sens_biomarker
  spec_t1 <- spec_biomarker
  # Correlations
  cor1 <- cortestbiomarker
  
  ######################################################################
  # Step 1 estimate the correlations:
  if (typemarker == "binary") {
    # Correlation b0 & t1
    frac_t1 = ev_rate * sens_t1 + (1-ev_rate) * (1-spec_t1)
    est_cor = optimize(find_intermediate_cor, 
                       interval = c(-0.5,0.98), 
                       cor_obj = cor1, N = N,
                       frac_positive = frac_t1)
    cor_est = as.numeric(est_cor$minimum)
  } else {
    frac_t1 = 0 # just to make functions run smoothly
    cor_est = cor1
  }
  
  # Step 2a: Get initial values for b1 b2 and b3
  est_b1 = exp(sqrt(2) * qnorm(AUC_ref))*(1-cor_est)/2+0.5
  if (typemarker == "continuous"){
    est_b2 = exp(sqrt(2) * qnorm(AUC_t1))*(1-cor_est)/2+0.5
  } else {
    eb2 = (sens_t1 * spec_t1) / ((1 - sens_t1) * (1 - spec_t1))
    est_b2 = exp(log(eb2)*cor_est)/2+0.5
  }
  
  # Step 2b: Estimate coefficients:
  system.time(
    ests_b1b2 <- finding_b1_b2(event_rate = ev_rate, 
                               AUC_pref = AUC_ref, 
                               linkfunction = linkfunction,
                               typet1 = typemarker, 
                               AUC_preft1 = AUC_t1, 
                               sens_preft1 = sens_t1, 
                               spec_preft1 = spec_t1, 
                               cor_est = cor_est, 
                               frac_t1 = frac_t1,
                               b1est = est_b1, 
                               b2est = est_b2))
  b1 = colMeans(ests_b1b2)[1]
  b2 = colMeans(ests_b1b2)[2]
  
  ##########################################################################
  ##########################################################################
  ##########################################################################
  # After estimation:
  # Generate data with Nunique -> finetune balancing precision and computation time
  N = Nunique
  ress <- Generate_testscores(N = N, 
                              cor_tests = cor_est, 
                              frac_positive = frac_t1,
                              type = typemarker)
  d <- ress[[2]]
  lp <- log(b1)*d$Reference + log(b2)*d$Newtest
  
  # Finally, find the correct b0 or h0:
  if (linkfunction == "logistic"){
    b0_h0 = as.numeric(optimize(find_b0_logistic, interval = c(-20,20), 
                                obs = event_rate, lp = lp)$minimum)
  } else if (linkfunction == "survival"){
    b0_h0 = as.numeric(optimize(find_h0_surv, interval = c(0.0001,20), 
                                obs = event_rate, lp = lp)$minimum)
  }
  # Fit the correct model: (chooses the b0_h0 value accordingly)
  dataset <- generate_risks(dataset = d, 
                            coefs = c(b0_h0,b1,b2,b3), 
                            linkfunction = linkfunction,
                            type = typemarker,
                            cutoff = ress[[1]],
                            rho = cor_est)
  
  ########################################################################
  ########################################################################
  # Perform checks:
  # Event rate should equal:
  EV_check = mean(dataset$Y_tr0)
  # AUC ref and increment in AUC by extra testing:
  AUC_ref = roc(dataset$Y_tr0, dataset$risk_reftest, quiet = TRUE)$auc
  AUCadd = roc(dataset$Y_tr0, dataset$risk_bothtest, quiet = TRUE)$auc
  # Input correlation should equal:
  corchecks <- cor(dataset[,1:2])
  cor_check1 = corchecks[1,2]
  
  # Check treatment effect
  if (linkfunction == "logistic"){
    or_treat = unique(round((dataset$risk_bothtestchemo / 
                               (1 - dataset$risk_bothtestchemo)) /
                              (dataset$risk_bothtest / (1 - dataset$risk_bothtest)),3))[1]
  } else if (linkfunction == "survival"){
    or_treat = unique(round(log(1 - dataset$risk_bothtestchemo) / 
                              log(1 - dataset$risk_bothtest),3))[1]
  }
  
  # Store results:
  res_df <- data.frame(N              = Nunique,
                       event_rate     = EV_check,
                       or_treat       = or_treat,
                       AUC_ref_model  = AUC_ref,
                       Incr_AUC       = AUCadd,
                       cortests       = cor_check1)
  
  # Depending on which test, add to results:
  if (typemarker == "continuous"){
    AUC_test1 = roc(dataset$Y_tr0, dataset$Newtest, quiet = TRUE)$auc
    t1_perf <- data.frame(AUC_test1 = AUC_test1)
  } else {
    sens_checkt1 <- sum(dataset$Newtest == 1 & dataset$Y_tr0 == 1) / 
      sum(dataset$Y_tr0 == 1)
    spec_checkt1 <- sum(dataset$Newtest == 0 & dataset$Y_tr0 == 0) / 
      sum(dataset$Y_tr0 == 0)
    t1_perf <- data.frame(Sens_test1 = sens_checkt1,
                          Spec_test1 = spec_checkt1)
  }
  
  # Add all the results:
  all_res <- cbind(res_df, t1_perf)
  
  # Return inputs and outputs:
  return(list(dataset, all_res))
}
model_pop_net_benefit <- function(d, thresholds) {
  # Put the risk models into a matrix
  ben_mat <- cbind(d$expben_reftest, d$expben_bothtest)
  true_ben <- d$expben_bothtest
  
  # Treat all strategy
  all <- sapply(thresholds, function(th) {
    (mean(true_ben - th)) * 1000})
  
  # Model strategies
  df <- t(vapply(thresholds, function(th) {
    if (any(ben_mat > th)) {
      colMeans((true_ben - th) * (ben_mat > th)) * 1000
    } else {c(0,0)}}, numeric(2)))
  
  # Treat none, treat all, risk models
  resdf <- as.data.frame(cbind(thresholds,
                               numeric(length(thresholds)),
                               all,
                               df,
                               df[,2] - df[,1]))
  colnames(resdf) <- c("thresholds","NB0","NB1","NB2","NB3","NBadd")
  return(resdf)
}
Decision_Analysis_plots <- function(da_data_output, plot_width = 800) {
  combined <- da_data_output
  # Prepare data frames as you have
  x = combined$thresholds
  sdf1 <- data.frame(x, y = combined$NB0, Type = "Treat none", Method = "1")
  sdf2 <- data.frame(x, y = combined$NB1, Type = "Treat all", Method = "1")
  sdf3 <- data.frame(x, y = combined$NB2, Type = "M1", Method = "2")
  sdf4 <- data.frame(x, y = combined$NB3, Type = "M2", Method = "2")
  
  df.merged <- rbind(sdf1, sdf2, sdf3, sdf4)
  df.merged$Group <- factor(paste(df.merged$Type, df.merged$Method))
  levels(df.merged$Group)[levels(df.merged$Group) == "Treat none 1"] <- "Treat none"
  levels(df.merged$Group)[levels(df.merged$Group) == "Treat all 1"] <- "Treat all"
  levels(df.merged$Group)[levels(df.merged$Group) == "M1 2"] <- "Reference model"
  levels(df.merged$Group)[levels(df.merged$Group) == "M2 2"] <- "Reference model + biomarker"
  
  # Set adaptive legend sizes based on plot width
  legend_text_size <- if (plot_width < 150) 4 else if (plot_width < 300) 8 else if (plot_width < 700) 14 else 20
  legend_key_size <- if (plot_width < 150) 0.25 else if (plot_width < 300) 0.5 else if (plot_width < 700) 1 else 1.5
  
  g <- ggplot(df.merged, aes(x, y, group = Group)) +
    geom_line(aes(colour = Group, linetype = Group)) +
    xlab('Treatment Decision Threshold') +
    ylab('Net Benefit') +
    scale_x_continuous(expand = c(0, 0), labels = scales::label_percent(scale = 100)) +
    scale_y_continuous(expand = c(0, 0), limits = c(-(0.3 * max(df.merged$y)), 1.01 * max(df.merged$y))) +
    theme_bw(base_size = 20) +
    theme(
      legend.position = c(0.99, 0.99),
      legend.justification = c("right", "top"),
      legend.text = element_text(size = legend_text_size),
      legend.title = element_text(size = legend_text_size + 2),
      legend.key.size = grid::unit(legend_key_size, "lines")
    ) +
    scale_linetype_manual(
      "Decision-making Strategy",
      values = c(1, 1, 2, 3),
      guide = guide_legend(override.aes = list(colour = c("#B4464B", "#4682B4", "black", "black")))
    ) +
    scale_colour_manual(
      "Decision-making Strategy",
      values = c("#B4464B", "#4682B4", "black", "black")
    )
  
  # Second plot with same approach
  df.merged2 <- data.frame(x, y = combined$NBadd,
                           Type = as.factor("Reference model + biomarker"))
  
  g2 <- ggplot(df.merged2, aes(x, y)) +
    geom_line(aes(colour = Type, linetype = Type)) + 
    xlab('Treatment Decision Threshold') + 
    ylab('Delta Net Benefit') +
    scale_x_continuous(expand = c(0, 0), labels = scales::label_percent(scale = 100)) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(min(-0.5,-(max(df.merged2$y) * 0.3)), 
                                  max(5,max(df.merged2$y, na.rm = TRUE) * 1.2))) +
    geom_hline(aes(yintercept = 0, colour = "Reference model", linetype = "Reference model")) + 
    scale_colour_manual(
      "Decision-making Strategy",
      breaks = c("Reference model + biomarker", "Reference model"), 
      values = c("#4682B4", "black")) +
    scale_linetype_manual(
      "Decision-making Strategy",
      breaks = c("Reference model + biomarker", "Reference model"), 
      values = c("solid", "dashed")) +
    theme_bw(base_size = 20) +
    theme(
      legend.position = c(0.99, 0.99),
      legend.justification = c("right", "top"),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.text = element_text(size = legend_text_size),
      legend.title = element_text(size = legend_text_size + 2),
      legend.key.size = grid::unit(legend_key_size, "lines")
    )
  
  return(list(a = g, b = g2))
}
make_table <- function(threshold, riskdf){
  threshold = threshold
  DFres = riskdf
  Npatients = 1000
  
  Comparison_results <- data.frame(Strategy=character(),
                                   Event=double(),
                                   Prevented=double(),
                                   Fraction_chemo=double(), 
                                   Fraction_tested=double(),
                                   Net_Benefit=double(),
                                   stringsAsFactors=FALSE)
  ###############################################################################
  # Reference strategies:
  # Stats:
  Strategy <- "Treat none"
  refEvent <- Event <- round(Npatients*mean(DFres$risk_bothtest[DFres$Treatment==0]))
  Prevented <- 0
  Frac_chemo <- 0
  Frac_test <- 0
  NB <- 0
  Comparison_results[1,] <- c(Strategy, Event, Prevented, Frac_chemo, Frac_test, NB)
  ###################################################################
  # Stats:
  Strategy <- "Treat all"
  Prevented <- round(Npatients*mean(DFres$expben_bothtest))
  Event <- refEvent-Prevented
  Frac_chemo <- Npatients
  Frac_test <- 0
  NB <- ifelse(threshold==0, round(Npatients*(mean(DFres$expben_bothtest) - threshold))
               ,round(Npatients*(mean(DFres$expben_bothtest) - threshold),1))
  Comparison_results[2,] <- c(Strategy, Event, Prevented, Frac_chemo, Frac_test, NB)
  ###################################################################
  # Decision analytic models
  # Baseline model
  chemo_ind = ifelse(DFres$expben_reftest > threshold, 1, 0)
  # Stats:
  Strategy <- "Reference model"
  Prevented <- ifelse(all(chemo_ind==F),0,round(Npatients*mean(DFres$expben_reftest[chemo_ind==1])*mean(chemo_ind==1)))
  Event <- refEvent-Prevented
  Frac_test <- 0
  Frac_chemo <- round(Npatients*mean(chemo_ind == 1))
  NB <- ifelse(all(chemo_ind == 0), 0, 
               round(Npatients*(mean(DFres$expben_reftest[chemo_ind==1])*mean(chemo_ind==1) - mean(chemo_ind==1)*threshold),1))
  Comparison_results[3,] <- c(Strategy, Event, Prevented, Frac_chemo, Frac_test, NB)
  ######################################################################################
  # MammaPrint model
  chemo_ind = ifelse(DFres$expben_bothtest > threshold, 1, 0)
  # Stats prep:
  lowbound <- min(DFres$Reference[DFres$expben_bothtest > threshold])
  uppbound <- min(DFres$Reference[DFres$expben_bothtest > threshold & DFres$Newtest == 0])
  # Stats:
  Strategy <- "Reference model + Biomarker"
  Prevented <- ifelse(all(chemo_ind==F),0,round(Npatients*mean(DFres$expben_bothtest[chemo_ind==1])*mean(chemo_ind==1)))
  Event <- refEvent-Prevented
  Frac_test <- round(Npatients*length(which(DFres$Reference>lowbound &
                                              DFres$Reference<uppbound))/nrow(DFres))
  Frac_chemo <- round(Npatients*mean(chemo_ind == 1))
  NB <- ifelse(all(chemo_ind == 0), 0, 
               round(Npatients*(mean(DFres$expben_bothtest[chemo_ind==1])*mean(chemo_ind==1) - mean(chemo_ind==1)*threshold),1))
  Comparison_results[4,] <- c(Strategy, Event, Prevented, Frac_chemo, Frac_test, NB)
  #######################################################################################
  Comparison_results[,7] <- c(0,0,0,round(1000*(sum(DFres$expben_reftest < threshold & 
                                                      DFres$expben_bothtest > threshold)/nrow(DFres) +
                                                  sum(DFres$expben_reftest > threshold & 
                                                        DFres$expben_bothtest < threshold)/nrow(DFres))))
  
  colnames(Comparison_results) <- c(
    "Strategy",
    "Expected events",
    "Prevented events",
    "Treatments given",
    "Tests done",
    "Net Benefit",
    "Reclassification"
  )
  return(Comparison_results)
}

ui <- fluidPage(
  
  # Custom CSS to remove default Bootstrap width limit
  tags$head(
    tags$style(HTML("
      .container-fluid {
        max-width: 100% !important;
        padding-left: 30px;
        padding-right: 30px;
      }
    "))
  ),
  
  # Title and subtitle
  h1("The Clinical Utility of additional Biomarker Testing",
     style = "text-align: center; color: #2c3e50; font-family: Arial, sans-serif; font-weight: bold;"),
  
  h3("A tool that translates diagnostic accuracy into impact on treatment decisions and patient outcomes",
     style = "text-align: center; color: #7f8c8d; font-family: Arial, sans-serif; font-style: italic;"),
  
  # Introductory HTML text (full width)
  HTML("
    <p>New biomarkers are increasingly identified that can improve risk predictions, 
    but evidence of their impact on patient outcomes is far more limited. Ideally, by improving 
    risk prediction, we can better identify which patients benefit from treatment, supporting 
    better treatment allocation and improving patient outcomes.Clinical usefulness of a new 
    biomarker requires that incremental value is shown over and above an established reference 
    model, ideally in a validation study. Strong empirical evidence can eventually come from 
    impact studies that assess decisions and outcomes with and without using information from 
    the biomarker. Such studies can be both costly and time consuming. Before committing to 
    such efforts, it is relevant to evaluate the potential that the biomarker will positively 
    impact patient management and outcomes.
    
    <p>This tool supports <strong>early assessment</strong> of a biomarker to appraise whether 
    a biomarker adds value in its intended clinical context. It provides a model-based evaluation 
    of how a biomarker can affect treatment decisions and patient outcomes based on inputs on its 
    initial diagnostic accuracy and context factors. The tool facilitates exploration how validated 
    diagnostic accuracy may translate into impact on treatment decisions and patient outcomes.</p>
    
    <h4>The tool requires six inputs:</h4>
    
    <ul>
      <li>The <strong>event rate</strong> in the target population without treatment</li>
      <li><strong>Reference model (AUC)</strong> in the target population without treatment</li>
      <li><strong>Biomarker AUC or sensitivity and specificity</strong> in the target population 
      without treatment</li>
      <li><strong>Correlation</strong> between the reference model and biomarker</li>
      <li><strong>Treatment effectiveness</strong></li>
    </ul>

    <p>The tool simulates directly compares treating based on a reference model alone vs. integrating a 
    biomarker. It uses the <strong>Net Benefit</strong> framework to evaluate the decisions strategies, 
    a metric that weighs the number of events prevented against the number of treatments given, across a 
    range of treatment thresholds. <strong>Delta Net Benefit</strong> refers to the improvement in net
    benefit by integrating the biomarker relative to without the biomarker.</p>

    <h4>The outputs:</h4>
    
    <ul>
      <li>Checks whether the inputs are matched in the simulation</li>
      <li>A <strong>decision curve analysis</strong> and the <strong>Delta Net Benefit</strong> of the 
      integrating the biomarker versus without the biomarker across decision thresholds</li>
      <li>A table summarizes a treat none or treat all strategy and the strategies with and without 
      the biomarker, for one specific treatment threshold</li>
    </ul>

    <h4>How to use the tool</h4>

    <p>The input values can be adjusted by clicking the <strong>Run simulation</strong> button. The tool 
    estimates all the necessary parameters and generates a data set (in ~10 seconds).<strong>Note</strong> 
    that with a mathematically impossible combination of inputs, <strong>the tool will still try to estimate 
    the parameters and give outputs</strong>. Consequently, the checks might show substantial deviation from
    the input values and interpretation of the results should be done with extra caution. For further 
    details and appropriate citations, please refer to the accompanying article and supplement.</p>
  "),
  
  h4("More information:"),
  
  actionButton("show_text1", "Net Benefit formula"),
  uiOutput("text_content1"),
  br(),
  
  actionButton("show_text2", "Simulation details"),
  uiOutput("text_content2"),
  br(),
  
  actionButton("show_text3", "Interpretation of the plots"),
  uiOutput("text_content3"),
  br(),
  
  # Inputs and Outputs
  sidebarLayout(
    
    sidebarPanel(
      h4(strong("Inputs:")),
      br(),
      
      selectInput("modelType", "Model", choices = c("Logistic" = "logistic", 
                                                    "Survival (Cox)" = "survival"), selected = "logistic"),
      
      numericInput("Eventrate_adv", "Event rate", 
                   value = 0.1, min = 0.01, max = 0.5, step = 0.01),
      shinyBS::bsTooltip("Eventrate_adv", title = "The event rate in the population without treatment", 
                         placement = "top"),
      
      numericInput("Refmodel_adv", "Reference model AUC", 
                   value = 0.7, min = 0.5, max = 0.9, step = 0.01),
      shinyBS::bsTooltip("Refmodel_adv", title = "Discrimination (AUC) of the reference model", 
                         placement = "top"),
      
      selectInput("test1_type", "Biomarker type:",
                  choices = c("Continuous" = "continuous", "Binary" = "binary")),
      
      # Continuous Test 1 (AUC)
      conditionalPanel(
        condition = "input.test1_type == 'continuous'",
        numericInput("test1_auc", "Biomarker AUC", value = 0.70, min = 0.5, max = 0.95, step = 0.01)
      ),
      
      # Binary Test 1 (Sensitivity & Specificity)
      conditionalPanel(
        condition = "input.test1_type == 'binary'",
        numericInput("Sensitivity_adv", "Biomarker Sensitivity", 
                     value = 0.65, min = 0.01, max = 0.99, step = 0.01),
        shinyBS::bsTooltip("Sensitivity_adv", title = "True positive rate for the biomarker", 
                           placement = "top"),
        numericInput("Specificity_adv", "Biomarker Specificity", 
                     value = 0.65, min = 0.01, max = 0.99, step = 0.01),
        shinyBS::bsTooltip("Specificity_adv", title = "True negative rate for the biomarker", 
                           placement = "top")
      ),
      
      numericInput("Treatment_adv", "Treatment effectiveness", 
                   value = 0.7, min = 0.5, max = 1, step = 0.01),
      shinyBS::bsTooltip("Treatment_adv", title = "Odds or hazard ratio for treatment effectiveness", 
                         placement = "top"),
      
      numericInput("Correlation_adv", "Correlation between reference model and biomarker", 
                   value = 0.5, min = 0, max = 0.7, step = 0.01),
      shinyBS::bsTooltip("Correlation_adv", title = "Correlation between reference model and biomarker", placement = "top"),
      
      numericInput("threshold_adv", "Treatment threshold (Affects only the table)", 
                   value = 0.05, min = 0, max = 0.5, step = 0.01),
      shinyBS::bsTooltip("threshold_adv", title = "Threshold for treatment benefit", placement = "top"),
      
      br(),
      actionButton("simulate", "Run simulation")
    ),
    
    mainPanel(
      uiOutput("inputcheck"),
      plotOutput("plot1", width = "100%"),
      plotOutput("plot2", width = "100%"),
      dataTableOutput("thresholdTable")
    )
  )
)

server <- function(input, output, session){
  thematic::thematic_shiny()
  
  # For the first one
  contentVisible1 <- reactiveVal(FALSE)
  observeEvent(input$show_text1, {
    contentVisible1(!contentVisible1()) # Toggle the visibility
    
    output$text_content1 <- renderUI({
      if (contentVisible1()) {
        tagList(
          HTML("<h4>What is Net Benefit?</h4>
          
          <p>
          Net Benefit is a decision-analytic measure that incorporates the clinical consequences of false positives and false negatives<sup>1</sup>.           The idea to incorporate consequences of decisions dates back from 1884 <sup>2</sup>, however, in recent year this has been formalized and            developed for diagnostic settings<sup>1</sup>.
          In this tool, we use an adapted version for treatment decisions<sup>3</sup>, where:
          </p>
          
          <ul>
            <li><strong>Benefit</strong> = number of events prevented by treatment</li>
            <li><strong>Harm</strong> = number of treatments given</li>
          </ul>
          
          <p>
          A key component is the weighting factor <strong>w</strong>, which reflects how many treatments are acceptable to prevent one event.
          This corresponds to the treatment threshold<sup>3</sup>. For example, if treatment is only justified when the absolute benefit is at 
          least 5%, then <strong>w = 0.05</strong>, meaning a willingness to treat 20 patients to prevent one event.
          </p>
          
          <p>
          Net Benefit is then calculated as: Net Benefit = Prevented events − w × Treatments given
          </p>
          
          <p>
          It is interpreted as the net number of events prevented per 1,000 patients, after adjusting for treatment burden.
          </p>
          
          <p>
          In our model, outcomes are perfectly calibrated. Therefore, the number of prevented events equals the expected absolute benefit (B) among
          treated patients. This allows Net Benefit to be expressed as:
          </p>
          
          <pre>
          Net Benefit = (E[B | B &gt; w] − w) × P(B &gt; w)
          </pre>
          
          <p>
          The number of prevented events can be expressed as the expectation of the absolute treatment benefit among treated patients (E[B ∣ B>w]
          multiplied by the probability of being treated (P(B>w)), which also reflects the number of treatments given.
          </p>
          
          <p>
          For further explanation, see the published article and supplement.
          </p>
          
          <p>
          [1] Vickers AJ, Elkin EB. Decision Curve Analysis: A Novel Method for Evaluating Prediction Models. Med Decis Making. 2006 Nov;26(6):565–574.
          </p>
          <p>
          [2] Peirce CS. The numerical measure of the success of predictions. Science. 1884 Nov 14;ns-4(93):453–454.
          </p>
          <p>
          [3] Vickers AJ, Kattan MW, Sargent DJ. Method for evaluating prediction models that apply the results of randomized trials to individual patients. Trials. 2007 Dec;8(1):14.
          </p>
          ")
        )
      } else {
        NULL  # Return nothing if the content is hidden
      }
    })
  })
  
  contentVisible2 <- reactiveVal(FALSE)
  observeEvent(input$show_text2, {
    contentVisible2(!contentVisible2()) # Toggle the visibility
    
    output$text_content2 <- renderUI({
      if (contentVisible2()) {
        tagList(
          HTML("
                <h4>What happens in the background?</h4>
                
                <p>
                This tool runs a large simulated study in the background to estimate how a new biomarker might affect treatment decisions and 
                patient outcomes. Here's an overview of how it works:
                </p>
                
                <h5>1. Simulating Patients</h5>
                <ul>
                  <li>We generate 500,000 unique patients to ensure stable results while keeping the tool fast.</li>
                  <li>Each patient gets a <strong>reference test score</strong> (e.g., clinical risk) and a <strong>biomarker result</strong> 
                  (high or low risk), based on the correlation you specify.</li>
                  <li>Each patient is <strong>cloned</strong>: one version receives treatment, the other does not.</li>
                </ul>
                
                <h5>2. Modeling Risk and Outcomes</h5>
                <ul>
                  <li>We assume a <strong>binary outcome</strong> (e.g., survival vs. death).</li>
                  <li>One can choose between two common medical models: <strong>logistic regression</strong> or a 
                  <strong>survival model</strong>, to model the link between predictors and the outcome.</li>
                  <li>The risk for each patient depends on the reference test, the biomarker result, and whether they were treated.</li>
                  <li>The model esimates its parameters to match the inputs: AUC for the reference test, sensitivity and specificity for the
                  biomarker, the event rate, and the treatment effect.</li>
                </ul>
                
                <h5>3. Making Treatment Decisions</h5>
                <ul>
                  <li>For each unique patient, we calculate the <strong>individual treatment benefit</strong>, 
                  i.e. the absolute risk reduction by treatment.</li>
                  <li>If the benefit is higher than the chosen treatment threshold, the model assumes the patient would be treated.</li>
                  <li>This is done twice: once using the reference test alone, and once with the biomarker added.</li>
                </ul>
                
                <h5>4. Deciding Who Needs Biomarker Testing</h5>
                <ul>
                  <li>The biomarker should only be used if it might change the treatment decision.</li>
                  <li>This happens when the patient’s estimated benefit (based on the reference test alone) is near the treatment threshold.</li>
                  <li>In other words, the biomarker is only used if the biomarker can change the decision, otherwise testing is futile.</li>
                </ul>
                
                <h5>⚠️ A Note on Warnings</h5>
                <p>
                In most real-world cases, the tool runs smoothly. But in some combinations of inputs the model may show a warning. 
                These combinations may be at or over the boundary of what is theoretical possible. For more technical details, 
                see the article and supplement.
                </p>
                ")
        )
      } else {
        NULL
      }
    })
  })
  
  contentVisible3 <- reactiveVal(FALSE)
  observeEvent(input$show_text3, {
    contentVisible3(!contentVisible3()) # Toggle the visibility
    
    output$text_content3 <- renderUI({
      if (contentVisible3()) {
        tagList(
          HTML("
              <h3>How to interpret the output</h3>
              
              <p>
              The <strong>upper plot</strong> shows a decision curve analysis comparing four strategies:
              <ul>
                <li><em>Treat None</em>: no patients receive treatment</li>
                <li><em>Treat All</em>: all patients are treated</li>
                <li><em>Reference Test</em>: treatment is guided by the reference test alone</li>
                <li><em>Reference + New Test</em>: treatment is guided by the reference test plus the new biomarker</li>
              </ul>
              </p>
              
              <p>
              The Y-axis shows <strong>Net Benefit</strong>, representing the number of events prevented per 1,000 patients after adjusting for 
              the harms of treatment. The X-axis is the <strong>treatment threshold</strong>: the minimum benefit needed to justify treatment in
              terms of absolute risk reduction.
              </p>
              
              <p>
              The <em>Treat None</em> line is always zero as no treatments means no harms or benefits. The <em>Treat All</em> line starts high 
              when the threshold is low (because harms are accepted), but drops below zero when the threshold increases, as harms begin to be 
              weighed more heavily, eventually becoming a harmful strategy.
              </p>
              
              <p>
              The lines for the two test-based strategies show how treatment guided by either the reference test or both tests performs across
              different thresholds. These also decline with increasing threshold, but typically outperform the two reference strategies across 
              (most; unless miscalibrated) of the range
              </p>
              
              <p>
              The <strong>lower plot</strong> focuses on the added value of the new test. It shows the <strong>Delta Net Benefit</strong>: the
              difference in Net Benefit between using both tests versus using the reference test alone. This directly reflects the clinical
              utility of the new biomarker Higher values mean more net events prevented due to better treatment targeting.
              </p>
              
              <h4>Table with strategy outcomes at a specific threshold</h4>
              
              <p>
              Below the plots, a table summarizes results for the four strategies at a specific treatment threshold that can be chosen by the user.
              Changing the input for the treatment threshold performs a quick update of only the table. The table helps understand how many patients
              benefit, how often the test is needed, and whether the new biomarker improves outcomes enough to justify its use. Values are shown 
              per 1,000 patients.
              </p>
              
              <ul>
                <li><strong>Expected events</strong>: the number of events expected under each strategy</li>
                <li><strong>Prevented events</strong>: the number of events avoided by treatment</li>
                <li><strong>Treatments given</strong>: how many patients received treatment under each strategy</li>
                <li><strong>Tests done</strong>: how often the new test is used, only in cases where it could affect the decision</li>
                <li><strong>Net Benefit</strong>: the weighted balance between events prevented and the burden of treatment, weighted by the
                treatment threshold</li>
              </ul>
              ")
          
        )
      } else {
        NULL  # Return nothing if the content is hidden
      }
    })
  })
  
  sim_data <- reactiveVal(NULL)
  make_plots <- eventReactive(input$simulate, {
    
    sensitivity <- input$Sensitivity_adv 
    specificity <- input$Specificity_adv
    dOR <- (sensitivity / (1 - sensitivity)) / ((1 - specificity) / specificity)
    if ((dOR >= 81) | (dOR < 1)) {
      return(NULL)  # Stop further execution
    }
    
    if (input$test1_type == "binary"){
      d <- Opt_function_all(Nunique = 50000,
                            linkfunction = input$modelType,
                            event_rate = input$Eventrate_adv,
                            treatmenteffect = input$Treatment_adv,
                            AUCrefmodel = input$Refmodel_adv,
                            type_biomarker = input$test1_type,
                            sens_biomarker = input$Sensitivity_adv,
                            spec_biomarker = input$Specificity_adv, 
                            cortestbiomarker = input$Correlation_adv)
    } else {
      d <- Opt_function_all(Nunique = 50000,
                            linkfunction = input$modelType,
                            event_rate = input$Eventrate_adv,
                            treatmenteffect = input$Treatment_adv,
                            AUCrefmodel = input$Refmodel_adv,
                            type_biomarker = input$test1_type,
                            AUC_biomarker = input$test1_auc,
                            cortestbiomarker = input$Correlation_adv)
    }
    
    # # Example inputs
    # Nunique         <- 50000            # number of simulated patients
    # linkfunction    <- "survival"       # type of outcome model
    # event_rate      <- 0.12             # baseline event rate
    # treatmenteffect <- 0.6              # relative risk reduction by treatment
    # AUCrefmodel     <- 0.70             # discrimination of the reference model
    # type_biomarker  <- "binary"         # type of biomarker
    # AUC_biomarker   <- 0.65             # AUC of the biomarker
    # sens_biomarker  <- 0.63             # sensitivity
    # spec_biomarker  <- 0.68             # specificity
    # cortestbiomarker<- 0.45             # correlation with reference model
    # 
    # d <- Opt_function_all(
    #   Nunique = Nunique,
    #   linkfunction = linkfunction,
    #   event_rate = event_rate,
    #   treatmenteffect = treatmenteffect,
    #   AUCrefmodel = AUCrefmodel,
    #   type_biomarker = type_biomarker,
    #   AUC_biomarker = AUC_biomarker,
    #   sens_biomarker = sens_biomarker,
    #   spec_biomarker = spec_biomarker,
    #   cortestbiomarker = cortestbiomarker
    # )
    
    df <- d[[1]]
    res <- d[[2]]
    threshold = 0.05#input$threshold_adv
    range_thresholds = seq(0, max(df$expben_bothtest), length.out = 50)
    # Make the plots:
    tb_nb <- make_table(threshold = threshold, riskdf = df)
    nb_df <- model_pop_net_benefit(d = df, thresholds = range_thresholds)
    plots <- Decision_Analysis_plots(nb_df)
    
    all <- list(a = plots$a, b = plots$b, c = tb_nb, d = res)
    sim_data(df)
    all
  }, ignoreNULL = FALSE)
  
  table_data <- reactive({
    req(sim_data())
    threshold <- input$threshold_adv
    make_table(threshold = threshold, riskdf = sim_data())
  })
  
  # Outputs:
  output$inputcheck <- renderUI({
    all <- make_plots()
    tbl <- all$d
    # Test 1 checks:
    if (input$test1_type == "continuous") {
      perf1 <- paste0("Biomarker AUC: ", round(tbl$AUC_test1, 2))
    } else {
      perf1 <- paste0("Biomarker Sensitivity: ", round(tbl$Sens_test1, 2),
                      ", Specificity: ", round(tbl$Spec_test1, 2))
    }
    # Message:
    HTML(paste0("<b>Check these values with your inputs:</b><br>",
                "Event rate: ", round(tbl$event_rate, 2), "<br>",
                "Treatment OR: ", round(tbl$or_treat, 2), "<br>",
                "AUC reference model: ", round(tbl$AUC_ref_model, 2), "<br>",
                "Correlation reference model & biomarker: ", round(tbl$cortests, 2), "<br>",
                perf1, "<br>"))
  })
  
  output$plot1 <- renderPlot({
    all <- make_plots()
    if (is.null(all)) {
      sensitivity <- input$Sensitivity_adv 
      specificity <- input$Specificity_adv
      dOR <- (sensitivity / (1 - sensitivity)) / ((1 - specificity) / specificity)
      plot(NA, xlim = c(0, 1), ylim = c(0, 1), xaxt = 'n', yaxt = 'n', bty = 'n', 
           xlab = "", ylab = "", main = paste0("⚠️ No calculations performed: Sensitivity and Specificity are too high or low. \nDiagnostic odds ratio (dOR) can not exceed 81 or lower than 1. Now the dOR (dOR = (sens / (1-sens)) / ((1-spec) / spec)) equals: ", dOR))
    } else {
      all$a
    }
  })
  
  output$plot2 <- renderPlot({
    all <- make_plots()
    if (is.null(all)) {
      sensitivity <- input$Sensitivity_adv 
      specificity <- input$Specificity_adv
      dOR <- (sensitivity / (1 - sensitivity)) / ((1 - specificity) / specificity)
      plot(NA, xlim = c(0, 1), ylim = c(0, 1), xaxt = 'n', yaxt = 'n', bty = 'n', 
           xlab = "", ylab = "", main = paste0("⚠️ No calculations performed: Sensitivity and Specificity are too high or low. \nDiagnostic odds ratio (dOR) can not exceed 81 or lower than 1. Now the dOR (dOR = (sens / (1-sens)) / ((1-spec) / spec)) equals: ", dOR))
    } else {
      all$b
    }
  })
  
  output$thresholdTable <- renderDataTable({
    datatable(
      table_data(),
      options = list(
        dom = 't',
        paging = FALSE
      )
    )
  })
}

shinyApp(ui = ui, server = server)



