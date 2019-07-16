# Importing needed packages
library(tidyverse) #Data wrangling
library(caret) # Pre processing
library(VIM) # Missing data analysis
library(missForest)
library(pROC) # ggplot addins for ROC curves
library(furrr) # For parallel processing

# Sourcing the auxiliary functions
source("Auxiliary_functions.R")

# Reading data
Data <-  as_tibble(read.csv2("WorkersData.csv")) #Original data
Data$Sex <- factor(Data$Sex)
Data$Skipped <- factor(Data$Skipped)
Data$Skipped <- factor(Data$Skipped, levels = levels(Data$Skipped)[c(2, 1)]) #Reordering levels

# Center and scale the data
preproc1 <- preProcess(Data, method=c("center", "scale"))
stdData <- predict(preproc1, Data)

# Removing variables with high missing values percentage (>15%)
stdData <- dplyr::select(stdData, -c(PTH, VitD, Phosphorus, NMissing))

# Performing missForest imputation
impData <- missForest(xmis = as.data.frame(stdData), verbose = T)
impData <- as.tibble(impData$ximp)

# Removing weight and height, which were only useful for imputation
# Both variables are fairly explained by BMI
impData <- dplyr::select(impData, -c(Weight, Height))
ccsData <- dplyr::select(stdData, -c(Weight, Height)) %>% # Complete-Cases data
  drop_na()

# Recoding factors for tunning requirements
impData$Skipped <- factor(ifelse(impData$Skipped == 0, "n", "y"), levels = c("n", "y"))
ccsData$Skipped <- factor(ifelse(ccsData$Skipped == 0, "n", "y"), levels = c("n", "y"))
impData$Sex <- factor(ifelse(impData$Sex == 0, "f", "m"), levels = c("f", "m"))
ccsData$Sex <- factor(ifelse(ccsData$Sex == 0, "f", "m"), levels = c("f", "m"))

# Seed for random number generation - For reproducibility purposes
set.seed(123456)

# Creating folds for cross validation
folds <- createMultiFolds(impData$Skipped, k = 10, times = 10)
# It was created with Imputed data cuz the modelos will be selected at it

# Creating fit control for caret to use in all the models
fitControl <- trainControl(
  # Stratified Indexes that will be used for CV
  index = folds, 
  ## 10-fold CV
  method = "repeatedcv", 
  number = 10,
  ## Each fold repeated 10 times
  repeats = 10,
  ## fun that returns sensitivity, specificity and AUROC
  summaryFunction = twoClassSummary,
  classProb = T, returnData = T, returnResamp = "all") 

# Fit control for cases where the model has already been selected
fitReduced <- trainControl(
  method = "none",
  ## fun that returns sensitivity, specificity and AUROC
  summaryFunction = twoClassSummary,
  classProb = T, returnData = T, returnResamp = "all")



# STEPWISE AND FULL LOGISTIC MODELS
# Fitting the full model in both datasets
fullImp <- train(Skipped ~ ., data = impData, 
                 method = "glm",
                 family = "binomial",
                 metric = "ROC",
                 trControl = fitControl)
fullCcs <- train(Skipped ~ ., data = ccsData, 
                 method = "glm",
                 family = "binomial",
                 metric = "ROC",
                 trControl = fitControl)
# Fitting the Stepwise models
stepImp <- train(Skipped ~ ., data = impData, 
                 method = "glmStepAIC",
                 family = "binomial",
                 metric = "ROC",
                 trControl = fitControl)
stepCcs <- train(Skipped ~ LDL + Trig + Chol + Sodium, data = ccsData, 
                 method = "glm",
                 family = "binomial",
                 metric = "ROC",
                 trControl = fitControl)

# Obtaining the coefficients
full_coefs <- tibble(Names = names(extr(fullImp)),
                     FullIMP = extr(fullImp),
                     FullIMPClean = coef(fullImp$finalModel),
                     FullCCS = extr(fullCcs),
                     FullCCSClean = coef(fullCcs$finalModel))

step1_coefs <- tibble(Names = names(extr(stepImp)), 
                      StepIMP = extr(stepImp))
step1_coefs_clean <- tibble(Names = names(coef(stepImp$finalModel)), 
                            StepIMPClean = coef(stepImp$finalModel))
step2_coefs <- tibble(Names = names(extr(stepCcs)), 
                      StepCCS = extr(stepCcs))
step2_coefs_clean <- tibble(Names = names(coef(stepCcs$finalModel)), 
                            StepCCSClean = coef(stepCcs$finalModel))

logist_coefs <- full_coefs %>% 
  full_join(step1_coefs, by = "Names") %>% 
  full_join(step1_coefs_clean, by = "Names") %>% 
  full_join(step2_coefs, by = "Names") %>% 
  full_join(step2_coefs_clean, by = "Names") %>% 
  replace_na(list("StepIMP" = " - ", "StepCCS" = " - ", 
                  "StepIMPClean" = " - ", "StepCCSClean" = " - "))

# Obtaining the predictions for the training data
imp_pred <- impData[, 1]
imp_pred$StepProb <- predict(stepImp$finalModel, impData, type = "response")
ccs_pred <- ccsData[, 1]
ccs_pred$StepProb <- predict(stepCcs$finalModel, ccsData, type = "response")

# Computing the ROC curve TRAINING
stepRocImp <- roc(imp_pred$Skipped, imp_pred$StepProb)
stepRocCcs <- roc(ccs_pred$Skipped, ccs_pred$StepProb)

# Performing bootstrap (r = 10000)
# Both boot_resamp and boot632Step are 
bootStepImp <- boot_resamp(r = 10000, 
                           resamp_func = boot632Step, 
                           df = impData, 
                           parallelize = T)

bootStepCcs <- boot_resamp(r = 10000, 
                           resamp_func = boot632Step, 
                           df = ccsData, 
                           parallelize = T)



# LASSO MODEL
# Fitting Lasso using the custom enet_fitting function
lassoImp <- enet_fitting(dt = impData, fittingControl = fitControl,
                         grid = expand.grid(
                           .alpha = 1, # alpha = 1 defines the lasso model
                           .lambda = seq(from = 0.015, to = 0.04, 
                                         length.out = 100)))
lassoCcs <- enet_fitting(dt = ccsData, fittingControl = fitReduced,
                         grid = expand.grid(
                           .alpha = 1, # alpha = 1 defines the lasso model
                           .lambda = lassoImp$bestTune$lambda)) 

#Lasso Fit coefficients
lasso_coefs <- tibble(Names = names(coefficients(lassoImp$finalModel, 
                                                 lassoImp$finalModel$lambdaOpt)[, 1]), 
                      LassoIMP = coefficients(lassoImp$finalModel, 
                                              lassoImp$finalModel$lambdaOpt)[, 1],
                      LassoCCS = coefficients(lassoCcs$finalModel, 
                                              lassoCcs$finalModel$lambdaOpt)[, 1])

# Obtaining the predictions for the training data
imp_pred$LassoProb <- pred(lassoImp, impData)$prob
ccs_pred$LassoProb <- pred(lassoCcs, ccsData)$prob

# Computing the ROC curve TRAINING
lassoRocImp <- roc(imp_pred$Skipped, imp_pred$LassoProb)
lassoRocCcs <- roc(ccs_pred$Skipped, ccs_pred$LassoProb)

# Performing bootstrap (r = 10000) for Imputed data
bootLassoImp <- boot_resamp(r = 10000, 
                                resamp_func = boot632Lasso, 
                                df = impData, 
                                parallelize = T)

bootLassoCcs <- boot_resamp(r = 10000, 
                            resamp_func = boot632Lasso, 
                            df = ccsData, 
                            parallelize = T)


# ELASTIC-NET MODEL
# Fitting the model with the elnet custom function
elnetImp <- enet_fitting(dt = impData, fittingControl = fitControl,
                         grid = expand.grid(
                           .alpha = seq(from = 0, to = 0.20, length.out = 100), # A grid to estimate Elastic-Net
                           .lambda = seq(from = 0.80, to = 1, length.out = 100)))
elnetCcs <- enet_fitting(dt = ccsData, fittingControl = fitReduced, 
                         grid = expand.grid(
                           .alpha = elnetImp$bestTune$alpha, # A grid to estimate Elastic-Net
                           .lambda = elnetImp$bestTune$lambda))

#Elastic-Net Fit coefficients
elnet_coefs <- cbind(ElNetIMP = coefficients(elnetImp$finalModel, 
                                             elnetImp$finalModel$lambdaOpt)[, 1],
                     ElNetCCS = coefficients(elnetCcs$finalModel, 
                                             elnetCcs$finalModel$lambdaOpt)[, 1])

# Obtaining the predictions for the training data
imp_pred$ElnetProb <- pred(elnetImp, impData)$prob
ccs_pred$ElnetProb <- pred(elnetCcs, ccsData)$prob

# Computing the ROC curve TRAINING
elnetRocImp <- roc(imp_pred$Skipped, imp_pred$ElnetProb)
elnetRocCcs <- roc(ccs_pred$Skipped, ccs_pred$ElnetProb)

# Performing bootstrap (r = 10000) for Imputed data
bootElnetImp_all <- boot_resamp(r = 10000, 
                                resamp_func = boot632Elnet, 
                                df = impData, 
                                parallelize = T)
bootElnetCcs_all <- boot_resamp(r = 10000, 
                                resamp_func = boot632Elnet, 
                                df = ccsData, 
                                parallelize = T)