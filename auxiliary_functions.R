# Function that fits a Elastic-Net or Lasso
enet_fitting <- function(dt, # The data with a Skipped variable and all the other independent ones
                         grid, fittingControl){ # The grid of parameters for caret to search the optimal set
  
  # Adjusting the model using the dt and grid arguments
  fitObj <- caret::train(Skipped ~ ., data = dt, 
                         # Using the glmnet package
                         method = "glmnet", 
                         # Binomial family to perform logistic regression
                         family = "binomial", 
                         trControl = fittingControl, 
                         # Using ROC to optimization, given by twoClassSummary function in fitControl
                         metric = "ROC", 
                         tuneGrid = grid)
  return(fitObj)
}

# Function that extracts the  coefficients already with significance marks
extr <- function(modobj){
  sum <- summary(modobj$finalModel)
  coef <- sum$coefficients
  coef[, 1] <- round(coef[, 1], 4)
  coefs <- ifelse(coef[, 4] <=0.01, 
                  paste0(coef[, 1], "**"), 
                  ifelse(coef[, 4] <= 0.05, 
                         paste0(coef[, 1], "*"), 
                         ifelse(coef[, 4] <= 0.1, 
                                paste0(coef[, 1], "."), paste0(coef[, 1], ""))))
  return(coefs)
}

# Function that adds predicted probabilities of Skipped to a data frame
pred <- function(modelObj, newData){
  newData %>% mutate(prob = predict(modelObj, newData, type = "prob")[, "y"])
}

# Defining function that fits the stepwise chosen model and returns the coefficients
return_coefficients <- function(data, # Dataset
                                idx, # Indexes of resample
                                fitControlObj){ # control object for caret modelling
  
  bootdf <- data[idx, ]
  
  # Training the model
  model <- train(Skipped ~ AvDBP + LDL + Trig + Chol + Sodium, data = bootdf, 
                 method = "glm",
                 family = "binomial",
                 metric = "ROC",
                 trControl = fitControlObj)
  
  cfs <- coefficients(model$finalModel)
  
  coefs <- tibble(Intercep = cfs[1], AvDBP = cfs[2],
                  LDL = cfs[3], Trig = cfs[4], 
                  Chol = cfs[5], Sodium = cfs[6])
  
  return(coefs)
}

# Function that performs the bootstrap resampling
boot_resamp <- function(r, # Number of resamples
                        resamp_func, # Function that performs the resample measurement
                        df, parallelize = FALSE){ # Should furrr prallel processin be used?
  
  # Creating indexes
  index2boot <- map(.x = 1:r, 
                    .f = ~sample(x = 1:nrow(df), 
                                 size = nrow(df), 
                                 replace = T))
  
  # Aplying the function to the indexes and returning a df
  if(parallelize){
    # Needs plan(multisession) to be run beforehand
    plan(multisession)
    return(future_map_dfr(.x = index2boot, 
                          .f = ~resamp_func(., df)))
  }else{
    return(map_dfr(.x = index2boot, 
                   .f = ~resamp_func(., df)))
  }
}

# Function that obtains the coefficients confidence intervals
boot_ci <- function(df, conf = 0.95){
  return(df %>% 
           gather(key = "key", val = "val") %>% 
           group_by(key) %>% 
           summarise(Mean = mean(val),
                     IL = quantile(val, (1-conf)/2),
                     SL = quantile(val, 1-(1-conf)/2)))
}

# Function that fits the selected model and returns the measures
boot632Step <- function(ind, df){
  
  # Fitting the model to resample
  mod <- train(Skipped ~ LDL + Trig + Chol + Sodium, 
               data = df[ind, ], 
               method = "glm",family = "binomial")
  
  # Obtaining the coefficients
  cof <- coefficients(mod$finalModel)
  
  # Obtaining predictions on resample (train)
  tru_trn <- df[ind, ]$Skipped
  prb_trn <- predict(mod$finalModel, df[ind, ], type = "response")
  
  # Obtaining best cutoff on 
  modRoc <- roc(tru_trn, prb_trn)
  modCut <- coords(modRoc, "best", 
                   input="threshold", best.method="youden")
  
  prd_trn <- ifelse(prb_trn > modCut[1], "y", "n")
  
  # Obtaining predictions on un-sampled (testing)
  tru_tst <- df[-ind, ]$Skipped
  prb_tst <- predict(mod$finalModel, df[-ind, ], type = "response")
  prd_tst <- ifelse(prb_tst > modCut[1], "y", "n")
  
  # Obtaining AUC, ACC, SENS and SPEC
  ## Training (resample)
  AUC_trn <- modRoc$auc
  ACC_trn <- mean(prd_trn == tru_trn)
  SNS_trn <- mean(prd_trn[tru_trn == "y"] == "y")
  SPC_trn <- mean(prd_trn[tru_trn == "n"] == "n")
  ## Test (unsampled)
  AUC_tst <- roc(tru_tst, prb_tst)$auc
  ACC_tst <- mean(prd_tst == tru_tst)
  SNS_tst <- mean(prd_tst[tru_tst == "y"] == "y")
  SPC_tst <- mean(prd_tst[tru_tst == "n"] == "n")
  
  # Returning results
  return(as.data.frame(t(c(cof,
                           "cutoff" = modCut[[1]],
                           "AUC" = 0.632*AUC_tst + 0.368*AUC_trn, 
                           "ACC" = 0.632*ACC_tst + 0.368*ACC_trn, 
                           "SNS" = 0.632*SNS_tst + 0.368*SNS_trn,
                           "SPC" = 0.632*SPC_tst + 0.368*SPC_trn))))
}

# Function that fits the selected model and returns the measures
boot632Lasso <- function(ind, df){
  
  # Fitting the model
  mod <- caret::train(Skipped ~ HDL + Chol + Sodium, data = df[ind, ], 
                      # Using the glmnet package
                      method = "glmnet", 
                      # Binomial family to perform logistic regression
                      family = "binomial", 
                      trControl = fitReduced, 
                      # Using ROC to optimization, given by twoClassSummary function in fitControl
                      metric = "ROC", 
                      tuneGrid = expand.grid(
                        .alpha = 1, # alpha = 1 defines the lasso model
                        .lambda = lassoImp$bestTune$lambda))
  
  # Obtaining the coefficients
  cof <- coefficients(mod$finalModel, 
                      mod$finalModel$lambdaOpt)[,1]
  
  # Obtaining predictions on resample (train)
  tru_trn <- df[ind, ]$Skipped
  prb_trn <- pred(mod, df[ind, ])$prob
  
  # Obtaining best cutoff
  modRoc <- roc(tru_trn, prb_trn)
  modCut <- coords(modRoc, "best", 
                   input="threshold", best.method="youden")
  
  # Predicting class
  prd_trn <- ifelse(prb_trn > modCut[1], "y", "n")
  
  # Obtaining predictions on non-sampled (testing)
  tru_tst <- df[-ind, ]$Skipped
  prb_tst <- pred(mod, df[-ind, ])$prob
  prd_tst <- ifelse(prb_tst > modCut[1], "y", "n")
  
  # Obtaining AUC, ACC, SENS and SPEC
  ## Training (resample)
  AUC_trn <- modRoc$auc
  ACC_trn <- mean(prd_trn == tru_trn)
  SNS_trn <- mean(prd_trn[tru_trn == "y"] == "y")
  SPC_trn <- mean(prd_trn[tru_trn == "n"] == "n")
  ## Test (unsampled)
  AUC_tst <- roc(tru_tst, prb_tst)$auc
  ACC_tst <- mean(prd_tst == tru_tst)
  SNS_tst <- mean(prd_tst[tru_tst == "y"] == "y")
  SPC_tst <- mean(prd_tst[tru_tst == "n"] == "n")
  
  # Returning results
  return(as.data.frame(t(c(cof,
                           "cutoff" = modCut[[1]],
                           "AUC" = 0.632*AUC_tst + 0.368*AUC_trn, 
                           "ACC" = 0.632*ACC_tst + 0.368*ACC_trn, 
                           "SNS" = 0.632*SNS_tst + 0.368*SNS_trn,
                           "SPC" = 0.632*SPC_tst + 0.368*SPC_trn))))
}

# Function that fits the models and returns the measures
boot632Elnet <- function(ind, df){
  
  # Fitting the model
  mod <- caret::train(Skipped ~ HDL + Chol + Sodium + Potassium, data = df[ind, ], 
                      # Using the glmnet package
                      method = "glmnet", 
                      # Binomial family to perform logistic regression
                      family = "binomial", 
                      trControl = fitReduced, 
                      # Using ROC to optimization, given by twoClassSummary function in fitControl
                      metric = "ROC", 
                      tuneGrid = expand.grid(
                        .alpha = elnetImp$bestTune$alpha, # alpha = 1 defines the lasso model
                        .lambda = elnetImp$bestTune$lambda))
  
  # Obtaining the coefficients
  cof <- coefficients(mod$finalModel, 
                      mod$finalModel$lambdaOpt)[,1]
  
  # Obtaining predictions on resample (train)
  tru_trn <- df[ind, ]$Skipped
  prb_trn <- pred(mod, df[ind, ])$prob
  
  # Obtaining best cutoff
  modRoc <- roc(tru_trn, prb_trn)
  modCut <- coords(modRoc, "best", 
                   input="threshold", best.method="youden")
  
  # Predicting class
  prd_trn <- ifelse(prb_trn > modCut[1], "y", "n")
  
  # Obtaining predictions on non-sampled (testing)
  tru_tst <- df[-ind, ]$Skipped
  prb_tst <- pred(mod, df[-ind, ])$prob
  prd_tst <- ifelse(prb_tst > modCut[1], "y", "n")
  
  # Obtaining AUC, ACC, SENS and SPEC
  ## Training (resample)
  AUC_trn <- modRoc$auc
  ACC_trn <- mean(prd_trn == tru_trn)
  SNS_trn <- mean(prd_trn[tru_trn == "y"] == "y")
  SPC_trn <- mean(prd_trn[tru_trn == "n"] == "n")
  ## Test (unsampled)
  AUC_tst <- roc(tru_tst, prb_tst)$auc
  ACC_tst <- mean(prd_tst == tru_tst)
  SNS_tst <- mean(prd_tst[tru_tst == "y"] == "y")
  SPC_tst <- mean(prd_tst[tru_tst == "n"] == "n")
  
  # Returning results
  return(as.data.frame(t(c(cof,
                           "cutoff" = modCut[[1]],
                           "AUC" = 0.632*AUC_tst + 0.368*AUC_trn, 
                           "ACC" = 0.632*ACC_tst + 0.368*ACC_trn, 
                           "SNS" = 0.632*SNS_tst + 0.368*SNS_trn,
                           "SPC" = 0.632*SPC_tst + 0.368*SPC_trn))))
}
