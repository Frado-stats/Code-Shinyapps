library(rms)
library(ggplot2)
library(dplyr)
library(mvtnorm)
library(tidyverse)
library(scales)
library(shiny)
library(shinyBS)
library(readr)
library(shinyWidgets)
library(rsconnect)
library(DT)
library(pROC)


find_intermediate_cor <- function(x, cor_obj, frac_positive, N){
  R <- matrix(data = c(1, x, x, 1), ncol = 2)
  set.seed(4)
  d <- as.data.frame(rmvnorm(N, mean = c(rep(0, nrow(R))), sigma = R))
  # Dichotomize
  d$V2 <- ifelse(d$V2 > sort(d$V2)[nrow(d)*(1-frac_positive)], 1, 0)
  # Minimize to the input correlation
  (as.numeric(cor(d)[1,2])-cor_obj)^2
}
find_intermediate_cor_bb <- function(x, cor_obj, frac_positive1, frac_positive2, N){
  R <- matrix(data = c(1, x, x, 1), ncol = 2)
  set.seed(4)
  d <- as.data.frame(rmvnorm(N, mean = c(rep(0, nrow(R))), sigma = R))
  # Dichotomize
  d$V1 <- ifelse(d$V1 > sort(d$V1)[nrow(d)*(1-frac_positive1)], 1, 0)
  d$V2 <- ifelse(d$V2 > sort(d$V2)[nrow(d)*(1-frac_positive2)], 1, 0)
  # Minimize to the input correlation
  (as.numeric(cor(d)[1,2])-cor_obj)^2
}
Generate_scores_opt <- function(N, cor_bt1, cor_bt2, cort1t2, frac_pos1, frac_pos2, typet1, typet2, i){
  # Our final matrix:
  R <- matrix(c(1,cor_bt1,cor_bt2,
                cor_bt1,1,cort1t2,
                cor_bt2,cort1t2,1), ncol = 3, byrow = T)
  # Generate data
  set.seed(i)
  df <- as.data.frame(rmvnorm(N, mean = c(rep(0, nrow(R))), sigma = R))
  colnames(df) <- c("Reference", "test1", "test2")
  df <- rbind(df,df)
  df$Chemo <- c(rep(c(0,1), each = nrow(df)/2))
  # Adjust accordingly
  if (typet1 == "binary"){
    df$test1 <- ifelse(df$test1 > sort(df$test1)[nrow(df)*(1-frac_pos1)], 1, 0)}
  if (typet2 == "binary"){
    df$test2 <- ifelse(df$test2 > sort(df$test2)[nrow(df)*(1-frac_pos2)], 1, 0)}
  # Output:
  return(df)
}
Generate_testscores <- function(N, cor_bt1, cor_bt2, cort1t2, frac_pos1, frac_pos2, typet1, typet2){
  # Our final matrix:
  R <- matrix(c(1,cor_bt1,cor_bt2,
                cor_bt1,1,cort1t2,
                cor_bt2,cort1t2,1), ncol = 3, byrow = T)
  # Generate data
  set.seed(4)
  df <- as.data.frame(rmvnorm(N, mean = c(rep(0, nrow(R))), sigma = R))
  colnames(df) <- c("Reference", "test1", "test2")
  df <- rbind(df,df)
  df$Chemo <- c(rep(c(0,1), each = nrow(df)/2))
  # Adjust accordingly
  if (typet1 == "binary"){
    df$test1 <- ifelse(df$test1 > sort(df$test1)[nrow(df)*(1-frac_pos1)], 1, 0)}
  if (typet2 == "binary"){
    df$test2 <- ifelse(df$test2 > sort(df$test2)[nrow(df)*(1-frac_pos2)], 1, 0)}
  # Output:
  return(df)
}
###############################################################################
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
######################################################################
generate_risks <- function(dataset, coefs, linkfunction){
  d <- dataset
  b0 <- coefs[1]
  b1 <- coefs[2]
  b2 <- coefs[3]
  b3 <- coefs[4]
  b4 <- coefs[5]
  
  # Get the correct outcome and model fits depending on link function:
  if (linkfunction == "logistic"){
    # Set the linear predictors to generate probabilities:
    lp <- b0 + log(b1)*d$Reference + log(b2)*d$test1 + log(b3)*d$test2 + log(b4)*d$Chemo
    preds <- exp(lp)/(1+exp(lp))
    set.seed(4)
    d$Y <- ifelse(runif(n=nrow(d)) < preds, 1, 0)
    
    # Store fits:
    fit_ref <- glm(Y ~ Reference + Chemo, data = d, family = binomial)
    fit_both <- glm(Y ~ Reference + test1 + test2 + Chemo, data = d, family = binomial)
    fit_test1 <- glm(Y ~ Reference + test1 + Chemo, data = d, family = binomial)
    fit_test2 <- glm(Y ~ Reference + test2 + Chemo, data = d, family = binomial)
    # Get risks:
    d$risks_ref <- predict(fit_ref, newdata = d, type = "response")
    d$risks_both <- predict(fit_both, newdata = d, type = "response")
    d$risks_test1 <- predict(fit_test1, newdata = d, type = "response")
    d$risks_test2 <- predict(fit_test2, newdata = d, type = "response")
    
  } else if (linkfunction == "survival"){
    h0=b0
    lp <- log(b1)*d$Reference + log(b2)*d$test1 + log(b3)*d$test2 + log(b4)*d$Chemo
    hazs <- h0*exp(lp)
    set.seed(4)
    times <- rexp(n=nrow(d), rate = hazs)
    d$Y <- (times < 1)*1
    d$S <- Surv(time = times, event = d$Y)
    
    # Store fits:
    fit_ref <- cph(S ~ Reference + Chemo, x=T,y=T, 
                   linear.predictors=T, surv=T, data = d)
    fit_both <- cph(S ~ Reference + test1 + test2 + Chemo, x=T,y=T, 
                    linear.predictors=T, surv=T, data = d)
    fit_test1 <- cph(S ~ Reference + test1 + Chemo, x=T,y=T, 
                     linear.predictors=T, surv=T, data = d)
    fit_test2 <- cph(S ~ Reference + test2 + Chemo, x=T,y=T, 
                     linear.predictors=T, surv=T, data = d)
    # Get risks:
    d$risks_ref <- as.numeric(1-survest(fit_ref, d, se.fit=F, times=1)$surv)
    d$risks_both <- as.numeric(1-survest(fit_both, d, se.fit=F, times=1)$surv)
    d$risks_test1 <- as.numeric(1-survest(fit_test1, d, se.fit=F, times=1)$surv)
    d$risks_test2 <- as.numeric(1-survest(fit_test2, d, se.fit=F, times=1)$surv)
  }
  
  # Get risks with chemo:
  d2 <- d[d$Chemo == 0,]
  d2$risks_ref_chemo   <- d$risks_ref[d$Chemo == 1]
  d2$risks_both_chemo  <- d$risks_both[d$Chemo == 1]
  d2$risks_test1_chemo <- d$risks_test1[d$Chemo == 1]
  d2$risks_test2_chemo <- d$risks_test2[d$Chemo == 1]
  # Expected treatment benefit from each model
  d2$expben_ref   <- d2$risks_ref   - d2$risks_ref_chemo
  d2$expben_both  <- d2$risks_both  - d2$risks_both_chemo
  d2$expben_test1 <- d2$risks_test1 - d2$risks_test1_chemo
  d2$expben_test2 <- d2$risks_test2 - d2$risks_test2_chemo
  # Format Outcomes: Y under treatment and no treatment
  colnames(d2)[which(colnames(d2)=="Y")] <- "Y_tr0"
  d2$Y_tr1 <- d$Y[d$Chemo == 1]
  
  return(as.data.frame(d2))
}

optimization_b0_b1_b2_b3 <- function(x, df, event_rate, AUC_pref0, AUC_preft1, AUC_preft2, sens_preft1, spec_preft1, sens_preft2, spec_preft2, typet1, typet2, counter, i, linkfunction){
  
  # Set iteration counter, to prevent using too many iterations:
  counter$iter_count <- counter$iter_count + 1
  if(counter$iter_count > counter$max_iter) return(1e9)
  if(any(x<=0)){x=x+5;print("error")}
  
  # Prepare data set:
  d = df
  # Find these coefficients:
  lp <- log(x[1])*d$Reference + log(x[2])*d$test1 + log(x[3])*d$test2
  
  if (linkfunction == "logistic"){
    # Update b0:
    findb0 = optimize(find_b0_logistic, interval = c(-20,20), obs = event_rate, lp = lp)
    b0 = as.numeric(findb0$minimum)
    # Simulate Y
    probs_Y <- exp(b0+lp)/(1+exp(b0+lp))
    set.seed(i)
    d$Y <- ifelse(runif(n=nrow(d)) < probs_Y, 1, 0)
  } else if (linkfunction == "survival"){
    # Update h0:
    findh0 = optimize(find_h0_surv, interval = c(0.0001,20), obs = event_rate, lp = lp)
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
  AUC_est0 <- as.numeric(roc(d$Y, d$Reference, quiet = TRUE)$auc)
  est0 <- -sum(AUC_pref0 * log(AUC_est0) + (1 - AUC_pref0) * log(1 - AUC_est0))
  
  # Depending on the type of test:
  if (typet1 == "binary" & typet2 == "binary") {
    # test1 sens/spec
    sens_est <- sum(d$test1 == 1 & d$Y == 1) / sum(d$Y == 1)
    spec_est <- sum(d$test1 == 0 & d$Y == 0) / sum(d$Y == 0)
    est1 <- -sum(sens_preft1 * log(sens_est) + (1 - sens_preft1) * log(1 - sens_est))
    est2 <- -sum(spec_preft1 * log(spec_est) + (1 - spec_preft1) * log(1 - spec_est))
    # test2 sens/spec
    sens_est2 <- sum(d$test2 == 1 & d$Y == 1) / sum(d$Y == 1)
    spec_est2 <- sum(d$test2 == 0 & d$Y == 0) / sum(d$Y == 0)
    est3 <- -sum(sens_preft2 * log(sens_est2) + (1 - sens_preft2) * log(1 - sens_est2))
    est4 <- -sum(spec_preft2 * log(spec_est2) + (1 - spec_preft2) * log(1 - spec_est2))
    # Minimize
    total <- sum(c(est0, est1, est2, est3, est4))
  } else if (typet1 == "binary" & typet2 == "continuous") {
    # test1 sens/spec
    sens_est <- sum(d$test1 == 1 & d$Y == 1) / sum(d$Y == 1)
    spec_est <- sum(d$test1 == 0 & d$Y == 0) / sum(d$Y == 0)
    est1 <- -sum(sens_pref * log(sens_est) + (1 - sens_pref) * log(1 - sens_est))
    est2 <- -sum(spec_pref * log(spec_est) + (1 - spec_pref) * log(1 - spec_est))
    # test2 AUC
    AUC_est2 <- as.numeric(roc(d$Y, d$test2, quiet = TRUE)$auc)
    est3 <- -sum(AUC_preft2 * log(AUC_est2) + (1 - AUC_preft2) * log(1 - AUC_est2))
    # Minimize
    total <- sum(c(est0, est1, est2, est3))
  } else if (typet1 == "continuous" & typet2 == "binary") {
    # test1 AUC
    AUC_est1 <- as.numeric(roc(d$Y, d$test1, quiet = TRUE)$auc)
    est1 <- -sum(AUC_preft1 * log(AUC_est1) + (1 - AUC_preft1) * log(1 - AUC_est1))
    # test2 sens/spec
    sens_est2 <- sum(d$test2 == 1 & d$Y == 1) / sum(d$Y == 1)
    spec_est2 <- sum(d$test2 == 0 & d$Y == 0) / sum(d$Y == 0)
    est2 <- -sum(sens_pref2 * log(sens_est2) + (1 - sens_pref2) * log(1 - sens_est2))
    est3 <- -sum(spec_pref2 * log(spec_est2) + (1 - spec_pref2) * log(1 - spec_est2))
    # Minimize
    total <- sum(c(est0, est1, est2, est3))
  } else if (typet1 == "continuous" & typet2 == "continuous") {
    # test1 AUC
    AUC_est1 <- as.numeric(roc(d$Y, d$test1, quiet = TRUE)$auc)
    est1 <- -sum(AUC_preft1 * log(AUC_est1) + (1 - AUC_preft1) * log(1 - AUC_est1))
    # test2 AUC
    AUC_est2 <- as.numeric(roc(d$Y, d$test2, quiet = TRUE)$auc)
    est2 <- -sum(AUC_preft2 * log(AUC_est2) + (1 - AUC_preft2) * log(1 - AUC_est2))
    # Minimize
    total <- sum(c(est0, est1, est2))
  }
  # Minimize
  total
}
finding_b1_b2_b3 <- function(event_rate, AUC_pref, typet1, typet2, 
                             AUC_preft1, AUC_preft2, c_b0t1, c_b0t2, c_t1t2,
                             frac_t1, frac_t2, linkfunction,
                             sens_preft1, spec_preft1, sens_preft2, spec_preft2,
                             b1est, b2est, b3est){
  
  # Set counter to track iterations
  counter <- new.env()
  counter$iter_count <- 0
  # Finetune balancing precision and computation time
  counter$max_iter <- 15
  
  b1res = b1est
  b2res = b2est
  b3res = b3est
  # Finetune k balancing precision and computation time
  k=25
  b1=b2=b3=numeric(k)
  for (i in 1:k) {
    # Finetune N here balancing precision and computation time
    d <- Generate_scores_opt(N = 1500, cor_bt1 = c_b0t1,
                             cor_bt2 = c_b0t2, cort1t2 = c_t1t2,
                             frac_pos1 = frac_t1, frac_pos2 = frac_t2,
                             typet1 = typet1, typet2 = typet2, i = i)
    
    # Reset iteration counter:
    counter$iter_count <- 0
    # The optimization function:
    ress <- optim(par = c(b1res,b2res,b3res), fn = optimization_b0_b1_b2_b3,
                  df = d, counter = counter, event_rate = event_rate, linkfunction = linkfunction,
                  AUC_pref0 = AUC_pref, AUC_preft1 = AUC_preft1, AUC_preft2 = AUC_preft2,
                  sens_pref = sens_preft1, spec_pref = spec_preft1, 
                  sens_preft2 = sens_preft2, spec_preft2 = spec_preft2,
                  typet1 = typet1, typet2 = typet2, i=i)
    # Store results:
    b1res = as.numeric(ress$par[1])
    b2res = as.numeric(ress$par[2])
    b3res = as.numeric(ress$par[3])
    
    # Repeat once for more optimal starting points of temporary b1b2b3 results:
    counter$iter_count <- 0
    ress <- optim(par = c(b1res,b2res,b3res), fn = optimization_b0_b1_b2_b3,
                  df = d, counter = counter, event_rate = event_rate, linkfunction = linkfunction,
                  AUC_pref0 = AUC_pref, AUC_preft1 = AUC_preft1, AUC_preft2 = AUC_preft2,
                  sens_preft1 = sens_preft1, spec_preft1 = spec_preft1, 
                  sens_preft2 = sens_preft2, spec_preft2 = spec_preft2,
                  typet1 = typet1, typet2 = typet2, i=i)
    
    b1[i] = as.numeric(ress$par[1])
    b2[i] = as.numeric(ress$par[2])
    b3[i] = as.numeric(ress$par[3])
    
    
    
    # Prepare data set:
    # Find these coefficients:
    lp <- log(b1[i])*d$Reference + log(b2[i])*d$test1 + log(b3[i])*d$test2
    
    if (linkfunction == "logistic"){
      # Update b0:
      findb0 = optimize(find_b0_logistic, interval = c(-20,20), obs = event_rate, lp = lp)
      b0 = as.numeric(findb0$minimum)
      # Simulate Y
      probs_Y <- exp(b0+lp)/(1+exp(b0+lp))
      set.seed(i)
      d$Y <- ifelse(runif(n=nrow(d)) < probs_Y, 1, 0)
    } else if (linkfunction == "survival"){
      # Update h0:
      findh0 = optimize(find_h0_surv, interval = c(0.0001,20), obs = event_rate, lp = lp)
      h0 = as.numeric(findh0$minimum)
      # Generate outcome to find coefficients:
      hazs <- h0*exp(lp)
      set.seed(i)
      times <- rexp(n=nrow(d), rate = hazs)
      d$Y <- (times < 1)*1
    }
    
    print(i)
    print(c(b1[i],b2[i],b3[i]))
    ########################################
    AUC_est0 <- as.numeric(roc(d$Y, d$Reference, quiet = TRUE)$auc)
    # Depending on the type of test:
    if (typet1 == "binary" & typet2 == "binary") {
      sens_est <- sum(d$test1 == 1 & d$Y == 1) / sum(d$Y == 1)
      spec_est <- sum(d$test1 == 0 & d$Y == 0) / sum(d$Y == 0)
      sens_est2 <- sum(d$test2 == 1 & d$Y == 1) / sum(d$Y == 1)
      spec_est2 <- sum(d$test2 == 0 & d$Y == 0) / sum(d$Y == 0)
      print(c(AUC_est0, sens_est, spec_est, sens_est2, spec_est2))
    } else if (typet1 == "continuous" & typet2 == "continuous") {
      # test1 AUC
      AUC_est1 <- as.numeric(roc(d$Y, d$test1, quiet = TRUE)$auc)
      AUC_est2 <- as.numeric(roc(d$Y, d$test2, quiet = TRUE)$auc)
      print(c(AUC_est0, AUC_est1, AUC_est2))
    }
  }
  print(colMeans(cbind(b1,b2,b3)))
  return(cbind(b1,b2,b3))
}
Opt_function_all <- function(Nunique, linkfunction, event_rate, treatmenteffect, AUCrefmodel, 
                             typet1, typet2, AUC_preft1=0, AUC_preft2=0, 
                             sens_preft1=0, spec_preft1=0, sens_preft2=0, spec_preft2=0, 
                             cor_b0t1, cor_b0t2, cor_t1t2){
  
  # Store the correct values:
  N <- 25000 # Nunique is used later, only to estimate correlations:
  linkfunction <- linkfunction # logistic or survival
  b4 <- treatmenteffect   
  AUC_ref_model <- AUCrefmodel
  # Test 1 inputs
  AUC_t1  <- AUC_preft1
  sens_t1 <- sens_preft1
  spec_t1 <- spec_preft1
  # Test 2 inputs
  AUC_t2  <- AUC_preft2
  sens_t2 <- sens_preft2
  spec_t2 <- spec_preft2
  # Correlations
  cor_b0t1 <- cor_b0t1
  cor_b0t2 <- cor_b0t2
  cor_t1t2 <- cor_t1t2
  # Types
  typet1 <- typet1
  typet2 <- typet2
  
  ######################################################################
  # Step 1 estimate the correlations:
  if (typet1 == "continuous" & typet2 == "binary"){
    # Correlation b0 & t1
    c_b0t1 = cor_b0t1
    # Correlation b0 & t2
    frac_t2 = event_rate * sens_t2 + (1-event_rate) * (1-spec_t2)
    est_cor = optimize(find_intermediate_cor, 
                       interval = c(-0.5,0.98), 
                       cor_obj = cor_b0t2, N = N,
                       frac_positive = frac_t2)
    c_b0t2 = as.numeric(est_cor$minimum)
    # Correlation t1 & t2
    est_cor = optimize(find_intermediate_cor, 
                       interval = c(-0.5,0.98), 
                       cor_obj = cor_t1t2, N = N,
                       frac_positive = frac_t2)
    c_t1t2 = as.numeric(est_cor$minimum)
  } else if (typet1 == "binary" & typet2 == "continuous") {
    # Correlation b0 & t1
    frac_t1 = event_rate * sens_t1 + (1-event_rate) * (1-spec_t1)
    est_cor = optimize(find_intermediate_cor, 
                       interval = c(-0.5,0.98), 
                       cor_obj = cor_b0t1, N = N,
                       frac_positive = frac_t1)
    c_b0t1 = as.numeric(est_cor$minimum)
    # Correlation b0 & t2
    c_b0t2 = cor_b0t2
    # Correlation t1 & t2
    est_cor = optimize(find_intermediate_cor, 
                       interval = c(-0.5,0.98), 
                       cor_obj = cor_t1t2, N = N,
                       frac_positive = frac_t1)
    c_t1t2 = as.numeric(est_cor$minimum)
  } else if (typet1 == "binary" & typet2 == "binary"){
    # Correlation b0 & t1
    frac_t1 = event_rate * sens_t1 + (1-event_rate) * (1-spec_t1)
    est_cor = optimize(find_intermediate_cor, 
                       interval = c(-0.5,0.98), 
                       cor_obj = cor_b0t1, N = N,
                       frac_positive = frac_t1)
    c_b0t1 = as.numeric(est_cor$minimum)
    # Correlation b0 & t2
    frac_t2 = event_rate * sens_t2 + (1-event_rate) * (1-spec_t2)
    est_cor = optimize(find_intermediate_cor, 
                       interval = c(-0.5,0.98), 
                       cor_obj = cor_b0t2, N = N,
                       frac_positive = frac_t2)
    c_b0t2 = as.numeric(est_cor$minimum)
    # Correlation t1 & t2
    est_cor = optimize(find_intermediate_cor_bb, 
                       interval = c(-0.5,0.98), 
                       cor_obj = cor_t1t2, N = N,
                       frac_positive1 = frac_t1,
                       frac_positive2 = frac_t2)
    c_t1t2 = as.numeric(est_cor$minimum)
  } else {
    c_b0t1 = cor_b0t1
    c_b0t2 = cor_b0t2
    c_t1t2 = cor_t1t2
  }
  
  # Step 2a: Get initial values for b1 b2 and b3
  est_b1 = exp(sqrt(2) * qnorm(AUC_ref_model))*(1-mean(c(c_b0t1, c_b0t2)))/2+0.5
  if (typet1 == "continuous"){
    est_b2 = exp(sqrt(2) * qnorm(AUC_t1))*(1-mean(c(c_b0t1, c_t1t2)))/2+0.5
  } else {
    eb2 = (sens_t1 * spec_t1) / ((1 - sens_t1) * (1 - spec_t1))
    est_b2 = exp(log(eb2)*((1-mean(c(c_b0t1, c_t1t2)))/2+0.5))
  }
  if (typet2 == "continuous"){
    est_b3 = exp(sqrt(2) * qnorm(AUC_t2))*(1-mean(c(c_b0t2, c_t1t2)))/2+0.5
  } else {
    eb3 = (sens_t2 * spec_t2) / ((1 - sens_t2) * (1 - spec_t2))
    est_b3 = exp(log(eb3)*((1-mean(c(c_b0t2, c_t1t2)))/2+0.5))
  }
  
  # Step 2b: Estimate coefficients:
  system.time(
    ests_b1b2b3 <- finding_b1_b2_b3(event_rate = event_rate, 
                                    AUC_pref = AUC_ref_model, linkfunction = linkfunction,
                                    typet1 = typet1, typet2 = typet2, 
                                    AUC_preft1 = AUC_t1, AUC_preft2 = AUC_t2, 
                                    sens_preft1 = sens_t1, spec_preft1 = spec_t1, 
                                    sens_preft2 = sens_t2, spec_preft2 = spec_t2,
                                    c_b0t1 = c_b0t1, c_b0t2 = c_b0t2, c_t1t2 = c_t1t2,
                                    frac_t1 = frac_t1, frac_t2 = frac_t2,
                                    b1est = est_b1, b2est = est_b2, b3est = est_b3))
  b1 = colMeans(ests_b1b2b3)[1]
  b2 = colMeans(ests_b1b2b3)[2]
  b3 = colMeans(ests_b1b2b3)[3]
  
  ##########################################################################
  ##########################################################################
  ##########################################################################
  # After estimation:
  # Generate data with Nunique -> finetune balancing precision and computation time
  N = Nunique
  d <- Generate_testscores(N = N, cor_bt1 = c_b0t1,
                           cor_bt2 = c_b0t2, cort1t2 = c_t1t2,
                           frac_pos1 = frac_t1, frac_pos2 = frac_t2,
                           typet1 = typet1, typet2 = typet2)
  lp <- log(b1)*d$Reference + log(b2)*d$test1 + log(b3)*d$test2
  
  # Finally, find the correct b0 or h0:
  if (linkfunction == "logistic"){
    b0_h0 = as.numeric(optimize(find_b0_logistic, interval = c(-20,20), 
                                obs = event_rate, lp = lp)$minimum)
  } else if (linkfunction == "survival"){
    b0_h0 = as.numeric(optimize(find_h0_surv, interval = c(0.0001,20), 
                                obs = event_rate, lp = lp)$minimum)
  }
  # Fit the correct model: (chooses the b0_h0 value accordingly)
  dataset <- generate_risks(dataset = d, coefs = c(b0_h0,b1,b2,b3,b4), linkfunction = linkfunction)
  
  ########################################################################
  ########################################################################
  # Perform checks:
  # Event rate should equal:
  EV_check = mean(dataset$Y_tr0)
  # AUC ref and increment in AUC by extra testing:
  AUC_ref = roc(dataset$Y_tr0, dataset$risks_ref, quiet = TRUE)$auc
  AUCadd1 <- roc(dataset$Y_tr0, dataset$risks_test1, quiet = TRUE)$auc
  AUCadd2 <- roc(dataset$Y_tr0, dataset$risks_test2, quiet = TRUE)$auc
  # Input correlation should equal:
  corchecks <- cor(dataset[,1:3])
  cor_check1 = corchecks[1,2]
  cor_check2 = corchecks[1,3]
  cor_check3 = corchecks[2,3]
  # Check treatment effect
  if (linkfunction == "logistic"){
    or_treat = unique(round((dataset$risks_both_chemo / (1 - dataset$risks_both_chemo)) /
                              (dataset$risks_both / (1 - dataset$risks_both)),3))[1]
  } else if (linkfunction == "survival"){
    or_treat = unique(round(log(1 - dataset$risks_both_chemo) / log(1 - dataset$risks_both),3))[1]
  }
  
  # Store results:
  res_df <- data.frame(N              = Nunique,
                       event_rate       = EV_check,
                       or_treat         = or_treat,
                       AUC_ref_model    = AUC_ref,
                       Incr_AUC_t1      = AUCadd1,
                       Incr_AUC_t2      = AUCadd2,
                       cor_b0t1         = cor_check1,
                       cor_b0t2         = cor_check2,
                       cor_t1t2         = cor_check3)
  
  # Depending on which test, add to results:
  if (typet1 == "continuous"){
    AUC_test1 = roc(dataset$Y_tr0, dataset$test1, quiet = TRUE)$auc
    t1_perf <- data.frame(AUC_test1 = AUC_test1)
  } else {
    sens_checkt1 <- sum(dataset$test1 == 1 & dataset$Y_tr0 == 1) / sum(dataset$Y_tr0 == 1)
    spec_checkt1 <- sum(dataset$test1 == 0 & dataset$Y_tr0 == 0) / sum(dataset$Y_tr0 == 0)
    t1_perf <- data.frame(Sens_test1 = sens_checkt1,
                          Spec_test1 = spec_checkt1)
  }
  if (typet2 == "continuous"){
    AUC_test2 = roc(dataset$Y_tr0, dataset$test2, quiet = TRUE)$auc
    t2_perf <- data.frame(AUC_test2 = AUC_test2)
  } else {
    sens_checkt2 <- sum(dataset$test2 == 1 & dataset$Y_tr0 == 1) / sum(dataset$Y_tr0 == 1)
    spec_checkt2 <- sum(dataset$test2 == 0 & dataset$Y_tr0 == 0) / sum(dataset$Y_tr0 == 0)
    t2_perf <- data.frame(Sens_test2 = sens_checkt2,
                          Spec_test2 = spec_checkt2)
  }
  # Add all the results:
  all_res <- cbind(res_df, t1_perf, t2_perf)
  
  # Return inputs and outputs:
  return(list(dataset, all_res))
}
diagnostic_net_benefit <- function(d, risk_thresholds) {
  k=20
  NB_models <- array(0, dim = c(length(risk_thresholds), 3, k))
  NB_all <- array(0, dim = c(length(risk_thresholds), 1, k))
  risk_mat <- cbind(d$risks_ref, d$risks_test1, d$risks_test2)
  n <- nrow(d)
  
  for (s in 1:k) {
    set.seed(100 + s) # loop over different seeds to generate Y
    Y <- ifelse(runif(n=n) < d$risks_both, 1, 0)
    
    NBall <- sapply(risk_thresholds, function(th) {
      p <- th / (1 - th)
      ev <- sum(Y)
      (ev - (n-ev) * p)/n * 1000
    })
    NB_mat <- t(sapply(risk_thresholds, function(th) {
      p <- th / (1 - th)
      class_mat <- risk_mat > th
      # TP & FP by multiplying by event indicator:
      e_t   <- colSums(class_mat * Y) #TPs
      e_not <- colSums(class_mat * abs(Y-1)) #FPs
      # Estimate NB: Vickers 2006
      (e_t - e_not * p)/n * 1000
    }))
    NB_all[ , , s] <- NBall
    NB_models[ , , s] <- NB_mat
  }
  # Average across seeds
  NBa <- apply(NB_all, c(1,2), mean)
  NBm <- apply(NB_models, c(1,2), mean)
  
  return(cbind(risk_thresholds,numeric(length(risk_thresholds)),NBa,NBm))
}
prognostic_net_benefit <- function(d, thresholds) {
  k=20
  NB_models <- array(0, dim = c(length(thresholds), 3, k))
  NB_all <- array(0, dim = c(length(thresholds), 1, k))
  ben_mat <- cbind(d$expben_ref, d$expben_test1, d$expben_test2)
  n <- nrow(d)
  
  for (s in 1:k) {
    set.seed(100 + s)
    
    Ynot <- ifelse(runif(n=nrow(d)) < d$risks_both, 1, 0)
    Yt   <- ifelse(runif(n=nrow(d)) < d$risks_both_chemo, 1, 0)
    T_none <- mean(Ynot)
    
    NBall <- sapply(thresholds, function(th) {
      ((sum(Ynot) - sum(Yt))/n) - th}) * 1000
    NB_mat <- t(sapply(thresholds, function(th) {
      # Following the procedure Vickers 2007
      # We can skip the step to check if model and tr indicator coincide:
      # subdf <- d[which(d$expben_..>th & d$Chemo==1 |
      #                  d$expben_..<=th & d$Chemo==0),]
      # We replicated all patients thus need to select the correct outcome: 
      # ifelse(d$expben_.. > th, Yt, Ynot)
      
      class_mat <- ben_mat > th
      n_treat <- colMeans(class_mat)
      
      ev_mat <- ifelse(class_mat, Yt, Ynot)
      ev <- colSums(ev_mat)/n
      
      # Estimate NB: Vickers 2007
      ((T_none - ev) - th * n_treat) * 1000
    }))
    NB_all[ , , s] <- NBall
    NB_models[ , , s] <- NB_mat
  }
  # Average across seeds
  NBa <- apply(NB_all, c(1,2), mean)
  NBm <- apply(NB_models, c(1,2), mean)
  
  return(cbind(thresholds,numeric(length(thresholds)),NBa,NBm))
}
model_pop_net_benefit <- function(d, thresholds) {
  # Put the 3 risk models into a matrix (n × 3)
  ben_mat <- cbind(d$expben_ref, d$expben_test1, d$expben_test2)
  true_ben <- d$expben_both
  
  all <- sapply(thresholds, function(th) {
    (mean(true_ben - th)) * 1000})
  df <- t(sapply(thresholds, function(th) {
    if (any(ben_mat > th)){
      colMeans((true_ben - th)*(ben_mat > th)) * 1000
    } else {c(0,0,0)}
  }))
  return(cbind(thresholds,numeric(length(thresholds)),all,df))
}
individual_net_benefit <- function(d, c_scores, threshold) {
  custom_round <- function(x, grid) {grid[which.min(abs(grid - x))]}
  d$expben_ref_round <- sapply(d$expben_ref, custom_round, grid = c_scores)
  t = threshold
  df <- t(sapply(c_scores, function(score) {
    subd <- d[which(d$expben_ref_round==score),]
    ben_mat <- cbind(subd$expben_ref, subd$expben_test1, subd$expben_test2)
    true_ben <- subd$expben_both
    if (any(ben_mat > t)){
      colMeans((true_ben - t)*(ben_mat > t)) * 1000
    } else {c(0,0,0)}
  }))
  return(cbind(c_scores, df))
}
Decision_Analysis_plots_diag <- function(d){
  combined <- as.data.frame(d)
  colnames(combined) <- c("thresholds","NB0","NB1","NB2","NB3","NB4")
  x <- combined$thresholds
  
  # Build long data frame
  df.merged <- rbind(
    data.frame(x, y = combined$NB0, Type = "Treat none", Method = "1"),
    data.frame(x, y = combined$NB1, Type = "Treat all", Method = "1"),
    data.frame(x, y = combined$NB2, Type = "M1", Method = "2"),
    data.frame(x, y = combined$NB3, Type = "M2", Method = "2"),
    data.frame(x, y = combined$NB4, Type = "M3", Method = "2")
  )
  
  # Smoothing:
  df.merged <- as.data.frame(df.merged %>%
                               group_by(Type, Method) %>%
                               mutate(y = loess(y ~ x, span = 0.25)$fitted) %>% ungroup())
  
  # Set legend names
  df.merged$Group <- factor(paste(df.merged$Type, df.merged$Method))
  levels(df.merged$Group)[levels(df.merged$Group) == "Treat none 1"] <- "Treat none"
  levels(df.merged$Group)[levels(df.merged$Group) == "Treat all 1"] <- "Treat all"
  levels(df.merged$Group)[levels(df.merged$Group) == "M1 2"] <- "Reference model"
  levels(df.merged$Group)[levels(df.merged$Group) == "M2 2"] <- "Reference model + Test 1"
  levels(df.merged$Group)[levels(df.merged$Group) == "M3 2"] <- "Reference model + Test 2"
  
  # Plot
  g <- ggplot(df.merged, aes(x, y, group = Group)) +
    geom_line(aes(colour = Group, linetype = Group)) +
    xlab('Risk Threshold') +
    ylab('Diagnostic Net Benefit') +
    scale_x_continuous(expand = c(0,0), labels = scales::label_percent(scale = 100)) +
    scale_y_continuous(expand = c(0,0),
                       limits = c(-0.3 * max(df.merged$y), 1.02 * max(df.merged$y))) +
    theme_bw(base_size = 20) +
    theme(legend.position = c(0.99, 0.99),
          legend.justification = c("right","top"))
  
  g <- g +
    scale_linetype_manual("Decision-making Strategy",
                          values = c("solid","solid","solid","dotted","dotdash"),
                          guide = guide_legend(
                            override.aes = list(colour = c("#3BAA8F", "#B4464B", "#4682B4",
                                                           "black", "black"))
                          )) +
    scale_colour_manual("Decision-making Strategy",
                        values = c("#3BAA8F", "#B4464B", "#4682B4",
                                   "black","black"))
  ############################################################
  ############################################################
  ############################################################
  
  # Generate dataframe:
  x <- combined$thresholds
  df.merged2 <- rbind(data.frame(x, y = combined$NB3-combined$NB2, Type = "Test 1"),
                      data.frame(x, y = combined$NB4-combined$NB2, Type = "Test 2"))
  # Smoothing:
  df.merged2 <- as.data.frame(df.merged2 %>% group_by(Type) %>%
                                mutate(y = loess(y ~ x, span = 0.25)$fitted) %>% ungroup())
  
  g2 <- ggplot(df.merged2, aes(x, y, colour = Type, linetype = Type)) +
    geom_line() +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    xlab("Risk Threshold") +
    ylab("Net Benefit of additional testing") +
    scale_x_continuous(expand = c(0, 0),
                       labels = scales::label_percent(scale = 100)) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(min(-0.5, min(df.merged2$y)*1.2),
                                  max(5, max(df.merged2$y)*1.2))) +
    scale_colour_manual("Additional test", values = c("#B4464B", "#4682B4")) +
    scale_linetype_manual("Additional test", values = c("solid", "solid")) +
    theme_bw(base_size = 20) +
    theme(legend.position = c(0.99, 0.99),
          legend.justification = c("right","top"),
          legend.background = element_rect(fill = "transparent", color = NA))
  
  return(list(a = g, b = g2))
}
Decision_Analysis_plots <- function(d, title){
  combined <- as.data.frame(d)
  colnames(combined) <- c("thresholds","NB0","NB1","NB2","NB3","NB4")
  x <- combined$thresholds
  
  # Building long data frame
  df.merged <- rbind(
    data.frame(x, y = combined$NB0, Type = "Treat none", Method = "1"),
    data.frame(x, y = combined$NB1, Type = "Treat all", Method = "1"),
    data.frame(x, y = combined$NB2, Type = "M1", Method = "2"),
    data.frame(x, y = combined$NB3, Type = "M2", Method = "2"),
    data.frame(x, y = combined$NB4, Type = "M3", Method = "2")
  )
  
  # Smoothing:
  df.merged <- as.data.frame(df.merged %>%
                               group_by(Type, Method) %>%
                               mutate(y = loess(y ~ x, span = 0.25)$fitted) %>% ungroup())
  
  # Set names for legend
  df.merged$Group <- factor(paste(df.merged$Type, df.merged$Method))
  levels(df.merged$Group)[levels(df.merged$Group) == "Treat none 1"] <- "Treat none"
  levels(df.merged$Group)[levels(df.merged$Group) == "Treat all 1"] <- "Treat all"
  levels(df.merged$Group)[levels(df.merged$Group) == "M1 2"] <- "Reference model"
  levels(df.merged$Group)[levels(df.merged$Group) == "M2 2"] <- "Reference model + Test 1"
  levels(df.merged$Group)[levels(df.merged$Group) == "M3 2"] <- "Reference model + Test 2"
  
  # Plot
  g <- ggplot(df.merged, aes(x, y, group = Group)) +
    geom_line(aes(colour = Group, linetype = Group)) +
    xlab('Treatment Threshold') +
    ylab(title) +
    scale_x_continuous(expand = c(0,0), labels = scales::label_percent(scale = 100)) +
    scale_y_continuous(expand = c(0,0),
                       limits = c(-0.3 * max(df.merged$y), 1.02 * max(df.merged$y))) +
    theme_bw(base_size = 20) +
    theme(legend.position = c(0.99, 0.99),
          legend.justification = c("right","top"))
  
  g <- g +
    scale_linetype_manual("Decision-making Strategy",
                          values = c("solid","solid","solid","dotted","dotdash"),
                          guide = guide_legend(
                            override.aes = list(colour = c("#3BAA8F", "#B4464B", "#4682B4",
                                                           "black", "black"))
                          )) +
    scale_colour_manual("Decision-making Strategy",
                        values = c("#3BAA8F", "#B4464B", "#4682B4",
                                   "black","black"))
  ############################################################
  # Make data frame:
  x <- combined$thresholds
  df.merged2 <- rbind(data.frame(x, y = combined$NB3-combined$NB2, Type = "Test 1"),
                      data.frame(x, y = combined$NB4-combined$NB2, Type = "Test 2"))
  
  df.merged2 <- as.data.frame(df.merged2 %>%
                                group_by(Type) %>% mutate(y = loess(y ~ x, span = 0.25)$fitted) %>% ungroup())
  
  g2 <- ggplot(df.merged2, aes(x, y, colour = Type, linetype = Type)) +
    geom_line() +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    xlab("Treatment Threshold") +
    ylab("Net Benefit of additional testing") +
    scale_x_continuous(expand = c(0, 0),
                       labels = scales::label_percent(scale = 100)) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(min(-0.5, min(df.merged2$y)*1.2),
                                  max(5, max(df.merged2$y)*1.2))) +
    scale_colour_manual("Additional test", values = c("#B4464B", "#4682B4")) +
    scale_linetype_manual("Additional test", values = c("solid", "solid")) +
    theme_bw(base_size = 20) +
    theme(legend.position = c(0.99, 0.99),
          legend.justification = c("right","top"),
          legend.background = element_rect(fill = "transparent", color = NA))
  
  return(list(a = g, b = g2))
}
individual_NBplot <- function(d, treatthreshold){
  t <- treatthreshold
  results_df <- as.data.frame(d)
  colnames(results_df) <- c("Clin_scores","Reference model", 
                            "Reference model + Test 1", "Reference model + Test 2")
  
  results_df$`Test 1` <- results_df$`Reference model + Test 1` - results_df$`Reference model`
  results_df$`Test 2` <- results_df$`Reference model + Test 2` - results_df$`Reference model`
  
  # Reshape
  long_df <- pivot_longer(results_df, 
                          cols = -Clin_scores, 
                          names_to = "Strategy", 
                          values_to = "NetBenefit")
  
  # Split into two datasets
  long_df_add <- long_df %>% filter(!grepl("Reference", Strategy))
  long_df_main <- long_df %>% filter(grepl("Reference", Strategy))
  
  # Common theme with dynamic legend
  custom_theme <- theme_minimal(base_size = 14) +
    theme(
      legend.position = c(0.02, 0.98),
      legend.justification = c(0.02, 0.98),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
  
  # Color palette for strategies
  strategy_colors_a <- c(
    "Reference model"   = "#3BAA8F",
    "Reference model + Test 1" = "#B4464B",
    "Reference model + Test 2" = "#4682B4")
  strategy_colors_b <- c(
    "Test 1" = "#B4464B",
    "Test 2" = "#4682B4")
  
  legend_position <- if (t > 0.025) c(0.02, 0.98) else c(0.98, 0.98)
  custom_theme <- theme_minimal(base_size = 14) +
    theme(
      legend.position = legend_position,
      legend.justification = legend_position,
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
  
  # Plot 1: Main strategies
  p1 <- ggplot(long_df_main, aes(x = Clin_scores, y = NetBenefit, color = Strategy)) +
    geom_line(size = 1.2) +
    scale_color_manual(values = strategy_colors_a) +
    labs(title = "The individual Net Benefit of testing",
         subtitle = paste0("Assuming treatment threshold of ", t*100," %"),
         x = "Estimated treatment benefit based on the reference model only",
         y = "Individual Net Benefit per strategy",
         color = "Strategy") +
    custom_theme
  
  # Plot 2: Additional benefit
  p2 <- ggplot(long_df_add, aes(x = Clin_scores, y = NetBenefit, color = Strategy)) +
    geom_line(size = 1.2) +
    scale_color_manual(values = strategy_colors_b) +
    labs(title = "The individual Net Benefit of the additional test",
         x = "Estimated treatment benefit based on the reference model only",
         y = "Individual Net Benefit of additional testing",
         color = "Additional test") +
    custom_theme
  
  return(list(a = p1, b = p2))
}
make_table <- function(d, threshold){
  threshold = threshold
  DFres = d
  Npatients = 1000
  
  Comparison_results <- data.frame(Strategy=character(),
                                   Event=double(),
                                   Prevented=double(),
                                   Fraction_chemo=double(), 
                                   Fraction_tested=double(),
                                   Net_Benefit1=double(),
                                   Net_Benefit2=double(),
                                   Net_Benefit3=double(),
                                   Reclassified=double(),
                                   stringsAsFactors=FALSE)
  ###############################################################################
  # Reference strategies:
  # Stats:
  Strategy <- "Treat none"
  refEvent <- Event <- round(Npatients*mean(DFres$Y_tr0))
  Prevented <- 0
  Frac_chemo <- 0
  Frac_test <- 0
  NB1 <- 0
  NB2 <- 0
  NB3 <- 0
  Reclassified <- 0
  Comparison_results[1,] <- c(Strategy, Event, Prevented, Frac_chemo, Frac_test, NB1, NB2, NB3, Reclassified)
  ###################################################################
  # Stats:
  Strategy <- "Treat all"
  Prevented <- round(Npatients*(mean(DFres$Y_tr0)-mean(DFres$Y_tr1)))
  Event <- refEvent-Prevented
  Frac_chemo <- Npatients
  Frac_test <- 0
  w <- DFres$risks_both[which(abs(DFres$expben_both-threshold)==min(abs(DFres$expben_both-threshold)))]
  NB1 <- round(Npatients * ((sum(DFres$Y_tr0)-(nrow(DFres)-sum(DFres$Y_tr0))*(w/(1-w)))/nrow(DFres)),1)
  NB2 <- round(Npatients * (((sum(DFres$Y_tr0)-sum(DFres$Y_tr1))/nrow(DFres)) - threshold),1)
  NB3 <- round(Npatients * (mean(DFres$expben_both) - threshold),1)
  Reclassified <- 0
  Comparison_results[2,] <- c(Strategy, Event, Prevented, Frac_chemo, Frac_test, NB1, NB2, NB3, Reclassified)
  ###################################################################
  # Decision analytic models
  # Baseline model
  chemo_ind = chemo_ind_ref = ifelse(DFres$expben_ref > threshold, 1, 0)
  # Stats:
  mean(DFres$expben_ref > threshold)
  Strategy <- "Reference model"
  Event <- ifelse(all(chemo_ind==F),0,round(Npatients*mean(ifelse(chemo_ind==1, DFres$Y_tr1, DFres$Y_tr0))))
  Prevented <- refEvent-Event
  Frac_test <- 0
  Frac_chemo <- round(Npatients*mean(chemo_ind == 1))
  NB1 <- round(diagnostic_net_benefit(d = DFres, risk_thresholds = w)[4],1)
  NB2 <- round(prognostic_net_benefit(d = DFres, thresholds = threshold)[4],1)
  NB3 <- round(model_pop_net_benefit(d = DFres, thresholds = threshold)[4],1)
  Reclassified <- 0
  Comparison_results[3,] <- c(Strategy, Event, Prevented, Frac_chemo, Frac_test, NB1, NB2, NB3, Reclassified)
  ######################################################################################
  # test 1
  chemo_ind = ifelse(DFres$expben_test1 > threshold, 1, 0)
  x1=max(DFres$expben_ref[DFres$expben_ref > threshold & DFres$expben_test1 < threshold])
  x2=min(DFres$expben_ref[DFres$expben_ref < threshold & DFres$expben_test1 > threshold])
  # Stats:
  Strategy <- "Reference model + test 1"
  Event <- ifelse(all(chemo_ind==F),0,round(Npatients*mean(ifelse(chemo_ind==1, DFres$Y_tr1, DFres$Y_tr0))))
  Prevented <- refEvent-Event
  Frac_test <- round(mean(DFres$expben_ref > x2 & DFres$expben_ref < x1)*100)
  Frac_chemo <- round(Npatients*mean(chemo_ind == 1))
  NB1 <- round(diagnostic_net_benefit(d = DFres, risk_thresholds = w)[5],1)
  NB2 <- round(prognostic_net_benefit(d = DFres, thresholds = threshold)[5],1)
  NB3 <- round(model_pop_net_benefit(d = DFres, thresholds = threshold)[5],1)
  Reclassified <- round(100*(nrow(DFres) - sum(chemo_ind == chemo_ind_ref)) / nrow(DFres))
  Comparison_results[4,] <- c(Strategy, Event, Prevented, Frac_chemo, Frac_test, NB1, NB2, NB3, Reclassified)
  ######################################################################################
  # test 2
  chemo_ind = ifelse(DFres$expben_test2 > threshold, 1, 0)
  x1=max(DFres$expben_ref[DFres$expben_ref > threshold & DFres$expben_test2 < threshold])
  x2=min(DFres$expben_ref[DFres$expben_ref < threshold & DFres$expben_test2 > threshold])
  # Stats:
  Strategy <- "Reference model + test 2"
  Event <- ifelse(all(chemo_ind==F),0,round(Npatients*mean(ifelse(chemo_ind==1, DFres$Y_tr1, DFres$Y_tr0))))
  Prevented <- refEvent-Event
  Frac_test <- round(mean(DFres$expben_ref > x2 & DFres$expben_ref < x1)*100)
  Frac_chemo <- round(Npatients*mean(chemo_ind == 1))
  NB1 <- round(diagnostic_net_benefit(d = DFres, risk_thresholds = w)[6],1)
  NB2 <- round(prognostic_net_benefit(d = DFres, thresholds = threshold)[6],1)
  NB3 <- round(model_pop_net_benefit(d = DFres, thresholds = threshold)[6],1)
  Reclassified <- round(100*(nrow(DFres) - sum(chemo_ind == chemo_ind_ref)) / nrow(DFres))
  Comparison_results[5,] <- c(Strategy, Event, Prevented, Frac_chemo, Frac_test, NB1, NB2, NB3, Reclassified)
  
  colnames(Comparison_results) <- c("Strategy",
                                    "Expected events",
                                    "Prevented events",
                                    "Treatments given",
                                    "Tests done",
                                    "Diagnostic Net Benefit",
                                    "Prognostic Net Benefit",
                                    "Model based Net Benefit",
                                    "Reclassified by additional testing")
  return(Comparison_results)
}

ui <- fluidPage(
  
  tags$head(
    tags$style(HTML("
    .main-container, 
    .container, 
    .container-fluid {
      max-width: 100% !important;
      width: 100% !important;
      padding-left: 0px !important;
      padding-right: 0px !important;
    }
    body {
      margin: 0 !important;
      padding: 0 !important;
    }
  "))
  ),
  
  # Title:
  h1("The Clinical Usefulness of Additional Testing", 
     style = "text-align: center; color: #2c3e50; font-family: Arial, sans-serif; font-weight: bold;"),
  # Subtitle:
  h3("From diagnostic performance to impact on individual and population level treatment allocation and patient outcomes", 
     style = "text-align: center; color: #7f8c8d; font-family: Arial, sans-serif; font-style: italic;"),
  # Introductory text
  
  HTML("<p>Predictive biomarkers are increasingly developed to refine risk estimates and guide treatment decisions. 
        However, not all biomarkers that improve prediction necessarily improve patient outcomes. Even if a test might improve patient outcomes, it does not necessarily do so for every patient. To justify clinical implementation, it is important to evaluate whether a biomarker is likely to improve clinical decision-making, referred to as its 
        <strong>clinical utility</strong>. This can be evaluated for populations or individual level.</p>

        <p>This tool supports <strong>early health technology assessment (HTA)</strong> of new diagnostic tests in its intended clinical context, or a <strong>didactic</strong> evaluation of the relation between <strong>diagnostic performance</strong>, the increment in <strong>AUC</strong>, the amount of <strong>reclassification</strong> by the extra test, and how it influences different versions of <strong>Net Benefit</strong>. To use the model, several inputs are required.</p>
      
        <h4>The inputs:</h4>
        
        <ul>
          <li>The <strong>link function</strong>: Either <strong>logistic</strong> or <strong>survival</strong> at given timepoint</li>
          <li>The <strong>event rate</strong> in the target population without treatment</li>
          <li><strong>Reference model performance (AUC)</strong>  in the target population without treatment</li>
          <li><strong>Treatment effectiveness</strong> in terms of an (equal for all) odds or hazard ratio</li>
          <li><strong>Correlation</strong> between the reference model and the additional tests, as well as between the additional tests</li>
          <li>The characteristics of the additional tests, either in terms of <strong>AUC</strong> when continuous or <strong>Sensitivity</strong> and <strong>Specificity</strong> when binary</li>
          <li>The <strong>Treatment threshold</strong> which is defined in terms of absolute treatment benefit (i.e. absolute risk reduction by treatment). This can be used on population and individual level</li>
        </ul>

        <h4>The outputs:</h4>
        <ul>
          <li>A check if the inputs are accurately replicated (some small variations can be expected due to randomness)</li>
          <li>Different versions of Net Benefit that evaluates <strong>population level</strong> performance: <strong>Diagnostic, Prognostic, Model-based Net Benefit</strong></li>
          <li>A table with the performance in terms of the <strong>Net Benefit</strong>, <strong>reclassification</strong> and <strong>incremental AUC</strong>, at a <strong>specific</strong> treatment threshold.</li>
          <li>The <strong>individual Net Benefit</strong> of performing additional tests. This shows the expected value of performing the additional tests for individual patients.</li>
        </ul>

        <h4>How to use the tool</h4>
        
        <p>Input values can be given and by clicking run simulation, the tool starts estimating the necessary parameters and generating a synthetic data set (~15seconds). The user can also only change the treatment threshold, which results in an instant update for only the table and individual Net Benefit plot, without waiting.</p>

        <p><strong>Note:</strong> When considering more extreme theoretical cases, the model does try to match inputs and does not give a warning or error when impossible combinations are inputed. For example, a combination that can be mathematically impossible is adding two highly predictive tests with a very low correlation. This results in poor convergence in the optimization steps to match the inputs. Therefore, the check is added if the values match closely with the inputs. There can be some slight deviations due to randomness by relatively small sample sizes to speed up the algorithm. Most practical scenarios can encounter some deviation, but will not encounter problems in convergence. For further explanation, and approiate citations, please refer to the accompanying article and the detailed supplement.</p>
      "),
  
  h4("More information:"),
  
  actionButton("show_text1", "AUC, reclassification, and Net Benefit variants"),
  uiOutput("text_content1"),
  br(),
  
  actionButton("show_text2", "Simulation details"),
  uiOutput("text_content2"),
  br(),
  
  actionButton("show_text3", "Interpretation of the plots"),
  uiOutput("text_content3"),
  br(),
  
  sidebarLayout(
    sidebarPanel(
      h4(strong("Inputs for the model:")),
      #####################################################################
      h3("Context parameters"),
      
      selectInput("modelType", "Model",
                  choices = c("Logistic" = "logistic",
                              "Survival (Cox)" = "survival"),
                  selected = "logistic"),
      
      numericInput("Eventrate0", "Event rate",
                   value = 0.10, min = 0.001, max = 0.5, step = 0.001),
      
      numericInput("Refmodel0", "Reference model AUC",
                   value = 0.70, min = 0.5, max = 0.9, step = 0.01),
      
      numericInput("Treatment0", "Treatment effectiveness",
                   value = 0.70, min = 0.5, max = 1, step = 0.01),
      
      numericInput("Correlation1", "Correlation (reference model vs test 1)",
                   value = 0.50, min = 0, max = 0.7, step = 0.01),
      
      numericInput("Correlation2", "Correlation (reference model vs test 2)",
                   value = 0.50, min = 0, max = 0.7, step = 0.01),
      
      numericInput("Correlation3", "Correlation (test 1 vs test 2)",
                   value = 0.50, min = 0, max = 0.7, step = 0.01),
      ## ---------------------------------------------------------------
      ##  Test 1
      ## ---------------------------------------------------------------
      h3("Additional Test 1"),
      selectInput("test1_type", "Test 1 type:",
                  choices = c("Continuous" = "continuous", "Binary" = "binary")),
      
      # Continuous Test 1 (AUC)
      conditionalPanel(
        condition = "input.test1_type == 'continuous'",
        numericInput("test1_auc", "AUC (Test 1)", value = 0.70, min = 0.5, max = 0.95, step = 0.01)
      ),
      
      # Binary Test 1 (Sensitivity & Specificity)
      conditionalPanel(
        condition = "input.test1_type == 'binary'",
        numericInput("test1_sens", "Sensitivity (Test 1)", value = 0.65, min = 0.01, max = 0.99, step = 0.01),
        numericInput("test1_spec", "Specificity (Test 1)", value = 0.65, min = 0.01, max = 0.99, step = 0.01)
      ),
      
      ## ---------------------------------------------------------------
      ##  Test 2
      ## ---------------------------------------------------------------
      h3("Additional Test 2"),
      selectInput("test2_type", "Test 2 type:",
                  choices = c("Continuous" = "continuous", "Binary" = "binary")),
      
      # Continuous Test 2 (AUC)
      conditionalPanel(
        condition = "input.test2_type == 'continuous'",
        numericInput("test2_auc", "AUC (Test 2)", value = 0.70, min = 0.5, max = 0.95, step = 0.01)
      ),
      
      # Binary Test 2 (Sensitivity & Specificity)
      conditionalPanel(
        condition = "input.test2_type == 'binary'",
        numericInput("test2_sens", "Sensitivity (Test 2)", value = 0.65, min = 0.01, max = 0.99, step = 0.01),
        numericInput("test2_spec", "Specificity (Test 2)", value = 0.65, min = 0.01, max = 0.99, step = 0.01)
      ),
      
      ## ---------------------------------------------------------------
      h3("Only updates the table and individual Net Benefit"),
      numericInput("threshold0", "Treatment threshold",
                   value = 0.05, min = 0.0, max = 0.15, step = 0.01),
      #####################################################################
      br(),
      actionButton("simulate", "Run simulation")
    ),
    
    mainPanel(
      uiOutput("inputcheck"),
      plotOutput("plot1"), 
      plotOutput("plot2"),
      plotOutput("plot3"),
      plotOutput("plot4"), 
      plotOutput("plot5"),
      plotOutput("plot6"),
      dataTableOutput("thresholdTable"),
      textOutput("tableNote"),
      plotOutput("plot7"),
      plotOutput("plot8")
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
          HTML("
          
          <h4>The evaluation of new tests</h4>
          
          <p>The tool tries to provide an overview of different methods to evaluate the performance of a new test in their intended clinical context. The tool provided an overview of their marginal diagnostic performance (as inputs), and how that translates in improvements in AUC, reclassification, and variants of Net Benefit. It is not meant as a complete overview of the performance of a test. Different overviews can be found elsewhere [cite ewout & cook].</p>
          
          <h4>What is the area under the receiver operator curve (AUC)?</h4>
          
          <p>The AUC is a famous measure that informs about the discriminative performance of a model. It asks how well we can discriminate between events and non-events [cite ewout]. For survival analysis it is often refered to as C-statistic or time-dependent AUC. The different between those is that the C-statistic compares predictions with survival times, whereas the time-dependent AUC at the status indicator at time point t [cite?].</p>
          
          <h4>What is reclassification?</h4>
          
          <p>The idea of reclassification is that added a new predictor changes the risk estimates in such a way that they reclassify a decision, which was made based on a decision threshold. In other words, for how many patients did the decision change by integrating the additional test (i.e. the risk estimate is updated across the decision threshold) [cite hunink]. A(n) (in)famous is the net reclassification index (NRI) [cite pencina]. This measure is widely used in clinical papers, but caution is advised given its criticisms [cite kerr]. Reclassification here refers to the amount of patients whose treatment decision did change by integrating the additional test (per 1000 patients).</p>
          
          <h4>What is Net Benefit?</h4>
          
          <p>Net Benefit is a decision-analytic measure that incorporates the clinical consequences of false positives and false negatives<sup>1</sup>. The idea to incorporate consequences of decisions dates back from 1884<sup>2</sup>, however, in recent years this idea has been formalized and developed for diagnostic settings<sup>1</sup>. A less famous variant to the diagnostic Net Benefit is the prognostic Net Benefit that integrates treatment effectiveness<sup>3</sup>. We included two extra variants. The first is a model-based Net Benefit, which is the same as the prognostic Net Benefit under very large numbers and strong calibration (i.e. predicted risks align with actual risks). The second is the individual Net Benefit, which shows the Net Benefit of additional testing on a more individual level.</p>
          
          <p>The diagnostic Net Benefit is estimated by:</p>
          
          <pre>
            NB = (TP − FP * (p/(1-p))) / n
          </pre>
          
          <p>The prognostic Net Benefit is estimated by:</p>
          
          <pre>
            NB = (Σᵢ=1ⁿ⁰ x₀ᵢ / n₀) − (Σᵢ=1ⁿ¹ x₁ᵢ / n₁) × (n₁ × w / (n₀ + n₁))
          </pre>
          
          <p>The model-based Net Benefit is estimated by:</p>
          
          <pre>
            NB = (E[B | B &gt; w] − w) × P(B &gt; w)
          </pre>
          
          <p>The individual Net Benefit is estimated by:</p>
          
          <pre>
            iNB = E[(B − w) * 1(B &gt; w) | X]
          </pre>
          
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
        tagList(HTML("
                <h4>What happens in the background?</h4>
                
                <p> After pressing 'run simulation', the tool starts to generate a synthetic data set. The sample size of 50000 is chosen by balancing accuracy and computational speed. This is meant as a general overview, more details can be found in the code.</p>
                
                <h5>1. Simulating predictors</h5>
                <ul>
                  <li>The three tests scores are generated with a multivariate standard normal distribution with given correlation structure. Depending on the choice of binary or continuous variants of the additional tests, their underlying continuous correlation is estimated first to simulate the data and later the variables are dichotomized.</li>
                  <li>Dichotomization of the test scores is done by defining the proportion of high first. This is defined by: P(test = high risk) = sensitivity * event rate + specificity * (1-event rate). The cutoff score is found by ordering the continuous scores and find the cutoff that makes this proportion.</li>
                  <li>Each patient is cloned: one with, and one without treatment.</li>
                </ul>
                
                <h5>2. Finding coefficients and simulating the outcome</h5>
                <ul>
                  <li>Two link functions can be chosen: logistic or survival.</li>
                  <li>Coefficients are iterively estimated by making sure they align with the input values.</li>
                  <li>Risks are estimated with and without treatment, the difference between those is the absolute risk reduction by treatment, i.e. treatment benefit.</li>
                </ul>
                
                <h5>3. Evaluating performance</h5>
                <ul>
                  <li>Checks whether the estimated values align with input values is incorporated. It is crucial to make sure the data set aligns with what is expected.</li>
                  <li>The AUC is defined for all models, the increment by adding the additional test is evaluated.</li>
                  <li>Treatment decisions are guided by all models, and are evaluated. The reclassification is evaluated at a specific treatment threshold and counts for how many patients (per 1000) the decision changed relative to using only the reference model.</li>
                  <li>For the diagnostic Net Benefit, the treatment threshold is transformed to the risk scale. For the other variants the treatment threshold is used as given by the input. The results are shown across a very wide range of treatment thresholds. Not all thresholds might be relevant. Some additional information on how the Net Benefit variants are estimated is given under 'AUC, reclassification, and Net Benefit variants'.</li>
                </ul>"))
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
              <h3>How to interpret the outputs</h3>
              
              <p>Central in the outputs are the variants of Net Benefit. All variants have two figures. The <strong>upper figures</strong> illustrates a traditional decision curve analysis. This figure compares five strategies:
              <ul>
                <li><em>Treat None</em>: no patients receive treatment</li>
                <li><em>Treat All</em>: all patients are treated</li>
                <li><em>Reference model</em>: treatments are guided by using only the reference model</li>
                <li><em>Reference model + Test 1</em>: treatments are guided by the reference model plus test 1</li>
                <li><em>Reference model + Test 2</em>: treatments are guided by the reference model plus test 2</li>
              </ul></p>
              
              <p>The Y-axis shows the <strong>Net Benefit</strong> estimated in their respective ways (see 'AUC, reclassification, and Net Benefit variants'). For the diagnostic Net Benefit this means the net true positive classifications; number of true positives penalized for the false positives. For the prognostic and model-based Net Benefit this means the net events prevented: the number of events prevented penalized for the number of treatments given. For the individual Net Benefit this means the net events prevented: the number of events prevented minus the minimal events prevented to accept treatment for the patient.</p>
              
              <p>The X-axis is the <strong>treatment threshold</strong>. For the diagnostic variant, this means the minimum risk to classify a patient as high risk. For the other variants, the treatment threshold refers to the the minimum benefit needed to justify treatment in terms of absolute risk reduction.</p>
              
              <p>The <em>Treat None</em> line is always zero as no treatments means no harms or benefits. The <em>Treat All</em> line starts high when the threshold is low (because harms are accepted), but drops below zero when the threshold increases, as harms begin to be weighed more heavily, eventually becoming a harmful strategy.</p>
              
              <p>The lines for the test-based strategies show how treatment guided by either the reference model or with additional testing across different thresholds. These also decline with increasing threshold, but typically outperform the two reference strategies across (unless miscalibrated, i.e. predicted risks do not align with true risks).</p>
              
              <p>The <strong>lower plot</strong> offsets the increase in Net Benefit (<strong>Delta Net Benefit</strong>) of additional testing against the reference model performance. The higher the delta Net Benefit, the higher the clinical usefulness of the new test.</p>
              
              <h4>Table with strategy outcomes at a specific threshold</h4>
              
              <p>Below the plots, a table summarizes results for the five strategies at a specific treatment threshold. Values are shown per 1,000 patients. Below the table the increment in AUC is given for both additional tests.</p>
              <ul>
                <li><strong>Expected events</strong>: the number of events expected under each strategy</li>
                <li><strong>Prevented events</strong>: the number of events avoided by treatment</li>
                <li><strong>Treatments given</strong>: how many patients received treatment under each strategy</li>
                <li><strong>Tests done</strong>: how many times the new test needs to be performed to get the increase in Net Benefit. This means the number of patients for which the additional test can reclassify the decision</li>
                <li><strong>Net Benefit variants</strong></li>
                <li><strong>Reclassified by additional testing</strong>: the number of reclassifications of the decision by updating the reference model estimate with the new test.</li>
              </ul>
              "))
      } else {
        NULL  # Return nothing if the content is hidden
      }
    })
  })
  # Store reactive data set
  sim_data <- reactiveVal(NULL)
  # Make plots if something changes in and hit simulation:
  make_plots <- eventReactive(input$simulate, {
    
    # Run optimization based on test type:
    if(input$test1_type == "continuous" & input$test2_type == "continuous"){
      res <- Opt_function_all(Nunique         = 25000,
                              linkfunction    = input$modelType,
                              event_rate      = input$Eventrate0,
                              treatmenteffect = input$Treatment0,
                              AUCrefmodel     = input$Refmodel0,
                              typet1          = input$test1_type,
                              typet2          = input$test2_type,
                              AUC_preft1      = input$test1_auc,
                              AUC_preft2      = input$test2_auc,
                              cor_b0t1        = input$Correlation1,
                              cor_b0t2        = input$Correlation2,
                              cor_t1t2        = input$Correlation3)
    } else if (input$test1_type == "continuous" & input$test2_type == "binary"){
      res <- Opt_function_all(Nunique         = 25000,
                              linkfunction    = input$modelType,
                              event_rate      = input$Eventrate0,
                              treatmenteffect = input$Treatment0,
                              AUCrefmodel     = input$Refmodel0,
                              typet1          = input$test1_type,
                              typet2          = input$test2_type,
                              AUC_preft1      = input$test1_auc,
                              sens_preft2     = input$test2_sens, 
                              spec_preft2     = input$test2_spec,
                              cor_b0t1        = input$Correlation1,
                              cor_b0t2        = input$Correlation2,
                              cor_t1t2        = input$Correlation3)
    } else if (input$test1_type == "binary" & input$test2_type == "continuous"){
      res <- Opt_function_all(Nunique         = 25000,
                              linkfunction    = input$modelType,
                              event_rate      = input$Eventrate0,
                              treatmenteffect = input$Treatment0,
                              AUCrefmodel     = input$Refmodel0,
                              typet1          = input$test1_type,
                              typet2          = input$test2_type,
                              sens_preft1     = input$test1_sens, 
                              spec_preft1     = input$test1_spec,
                              AUC_preft2      = input$test2_auc,
                              cor_b0t1        = input$Correlation1,
                              cor_b0t2        = input$Correlation2,
                              cor_t1t2        = input$Correlation3)
    } else if (input$test1_type == "binary" & input$test2_type == "binary"){
      res <- Opt_function_all(Nunique         = 25000,
                              linkfunction    = input$modelType,
                              event_rate      = input$Eventrate0,
                              treatmenteffect = input$Treatment0,
                              AUCrefmodel     = input$Refmodel0,
                              typet1          = input$test1_type,
                              typet2          = input$test2_type,
                              sens_preft1     = input$test1_sens, 
                              spec_preft1     = input$test1_spec,
                              sens_preft2     = input$test2_sens, 
                              spec_preft2     = input$test2_spec,
                              cor_b0t1        = input$Correlation1,
                              cor_b0t2        = input$Correlation2,
                              cor_t1t2        = input$Correlation3)
    } # end of analysis if something above changes + simulation hit
    
    # Get datasets:
    d <- res[[1]]
    r <- res[[2]]
    # If threshold changes, the plots do not change
    th = input$threshold0
    
    # Net Benefit variants:
    dd1 <- diagnostic_net_benefit(d=d, risk_thresholds = seq(0,min(0.95,max(d$risks_both)), length.out = 25))
    dd2 <- prognostic_net_benefit(d=d, thresholds = seq(0, max(d$expben_both), length.out = 25))
    dd3 <- model_pop_net_benefit(d=d, thresholds = seq(0, max(d$expben_both), length.out = 25))
    
    # Some extra steps for the individual Net Benefit
    c_scores = seq(0.005, round(max(d$expben_ref)*200)/200, by = 0.005)
    custom_round <- function(x, grid) {grid[which.min(abs(grid - x))]}
    d$expben_ref_round <- sapply(d$expben_ref, custom_round, grid = c_scores)
    max_score = as.numeric(max(names(table(d$expben_ref_round)/nrow(d) > 0.01)))
    c_scores = seq(0.005, max_score, by = 0.005)
    
    dd4 <- individual_net_benefit(d=d, c_scores = seq(0.005, max_score, by = 0.005), 
                                  threshold = th)
    
    # Make Net Benefit plots and table
    aa <- Decision_Analysis_plots_diag(d = dd1)
    bb <- Decision_Analysis_plots(d = dd2, title = "Prognostic Net Benefit")
    cc <- Decision_Analysis_plots(d = dd3, title = "Model-based Net Benefit")
    dd <- make_table(d = d, threshold = th)
    ee <- individual_NBplot(d = dd4, treatthreshold = th)
    
    all <- list(a = aa$a, a2 = aa$b, 
                b = bb$a, b2 = bb$b, 
                c = cc$a, c2 = cc$b, 
                d = dd, e = r,
                f = ee$a, f2 = ee$b)
    sim_data(d)
    all
  }, ignoreNULL = FALSE)
  
  # Three functions reactive of the threshold input:
  table_data <- reactive({
    req(sim_data())
    threshold <- input$threshold0
    make_table(d = sim_data(), threshold = threshold)
  })
  
  i_plot <- reactive({
    req(sim_data())
    threshold <- input$threshold0
    d2 <- sim_data()
    
    c_scores = seq(0.005, round(max(d2$expben_ref)*200)/200, by = 0.005)
    custom_round <- function(x, grid) {grid[which.min(abs(grid - x))]}
    d2$expben_ref_round <- sapply(d2$expben_ref, custom_round, grid = c_scores)
    max_score = as.numeric(max(names(table(d2$expben_ref_round)/nrow(d2) > 0.01)))
    c_scores = seq(0.005, max_score, by = 0.005)
    
    dd4 <- individual_net_benefit(d=d2,
                                  c_scores = c_scores,
                                  threshold = threshold)
    pp <- individual_NBplot(d = dd4, treatthreshold = threshold)
    pp$a
  })
  
  i_plot2 <- reactive({
    req(sim_data())
    threshold <- input$threshold0
    d2 <- sim_data()
    
    c_scores = seq(0.005, round(max(d2$expben_ref)*200)/200, by = 0.005)
    custom_round <- function(x, grid) {grid[which.min(abs(grid - x))]}
    d2$expben_ref_round <- sapply(d2$expben_ref, custom_round, grid = c_scores)
    max_score = as.numeric(max(names(table(d2$expben_ref_round)/nrow(d2) > 0.01)))
    c_scores = seq(0.005, max_score, by = 0.005)
    
    dd4 <- individual_net_benefit(d=d2,
                                  c_scores = c_scores,
                                  threshold = threshold)
    pp <- individual_NBplot(d = dd4, treatthreshold = threshold)
    pp$b
  })  
  ###################################################
  # Outputs:
  output$inputcheck <- renderUI({
    all <- make_plots()
    tbl <- all$e
    # Test 1 and Test 2 checks:
    if (input$test1_type == "continuous") {
      perf1 <- paste0("Test 1 AUC: ", round(tbl$AUC_test1, 2))
    } else {
      perf1 <- paste0("Test 1 Sensitivity: ", round(tbl$Sens_test1, 2),
                      ", Specificity: ", round(tbl$Spec_test1, 2))
    }
    if (input$test2_type == "continuous") {
      perf2 <- paste0("Test 2 AUC: ", round(tbl$AUC_test2, 2))
    } else {
      perf2 <- paste0("Test 2 Sensitivity: ", round(tbl$Sens_test2, 2),
                      ", Specificity: ", round(tbl$Spec_test2, 2))
    }
    # Message:
    HTML(paste0("<b>Check these values with your inputs:</b><br>",
                "Event rate: ", round(tbl$event_rate, 3), "<br>",
                "Treatment OR: ", round(tbl$or_treat,2), "<br>",
                "AUC reference model: ", round(tbl$AUC_ref_model, 2), "<br>",
                "Correlation reference model & Test 1: ", round(tbl$cor_b0t1, 2), "<br>",
                "Correlation reference model & Test 2: ", round(tbl$cor_b0t2, 2), "<br>",
                "Correlation Test 1 & Test 2: ", round(tbl$cor_t1t2, 2), "<br>",
                perf1, "<br>",
                perf2, "<br>"))
  })
  
  output$plot1 <- renderPlot({
    all <- make_plots()
    all$a
  })
  
  output$plot2 <- renderPlot({
    all <- make_plots()
    all$a2
  })
  
  output$plot3 <- renderPlot({
    all <- make_plots()
    all$b
  })
  
  output$plot4 <- renderPlot({
    all <- make_plots()
    all$b2
  })
  
  output$plot5 <- renderPlot({
    all <- make_plots()
    all$c
  })
  
  output$plot6 <- renderPlot({
    all <- make_plots()
    all$c2
  })
  
  output$thresholdTable <- renderDataTable({
    datatable(table_data(),options = list(dom = 't',paging = FALSE))
  })
  
  output$tableNote <- renderText({
    all <- make_plots()
    tbl <- all$e
    auc1 <- round(tbl$Incr_AUC_t1, 3)
    auc2 <- round(tbl$Incr_AUC_t2, 3)
    paste0("Incremental AUC by adding: Test 1 = ", auc1, ", Test 2 = ", auc2)
  })
  
  output$plot7 <- renderPlot({
    i_plot()
  })
  
  output$plot8 <- renderPlot({
    i_plot2()
  })
  ########################################################
}

shinyApp(ui = ui, server = server)