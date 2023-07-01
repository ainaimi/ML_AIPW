pacman::p_load(
  rio,          
  here,         
  skimr,        
  tidyverse,     
  gtsummary,
  parallel,
  rstatix,      
  janitor,      
  scales,        
  dplyr,
  tidyverse,
  haven,
  lmtest,
  sandwich,
  MASS,
  grf,
  SuperLearner,
  pROC,
  gridExtra,
  vip,
  broom,
  xtable,
  GGally,
  ggpubr
)

## Zhaohua, this is for you!!

mean_learner <- "SL.mean"
glm_learner <- "SL.glm"

# ranger learner
ranger_learner <- create.Learner("SL.ranger",
                                 params = list(num.trees = 500, 
                                               min.node.size = 50),
                                 tune = list(mtry = c(3,4,5)))

# glmnet learner
glmnet_learner <- create.Learner("SL.glmnet",
                                 tune = list(alpha = seq(0,1,.25)))

# xgboost learner
xgboost_learner <- create.Learner("SL.xgboost",
                                  params = list(nrounds = 500),
                                  tune = list(max_depth = c(4,6),
                                              eta = c(.1,.2)))

# earth learner
earth_learner <- create.Learner("SL.earth",
                                tune = list(degree = c(3,4,5)))

.SL.require <- function(package, message = paste('loading required package (', package, ') failed', sep = '')) {
  if(!requireNamespace(package, quietly = FALSE)) {
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}
screen.glmnet1 <- function (Y, X, family, alpha = 1, minscreen = 20, nfolds = 10, 
                            nlambda = 100, ...) 
{
  .SL.require("glmnet")
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, lambda = NULL, type.measure = "deviance", 
                             nfolds = nfolds, family = family$family, alpha = alpha, 
                             nlambda = nlambda)
  whichVariable <- (as.numeric(coef(fitCV$glmnet.fit, s = fitCV$lambda.min))[-1] != 
                      0)
  if (sum(whichVariable) < minscreen) {
    warning("fewer than minscreen variables passed the glmnet screen, increased lambda to allow minscreen variables")
    sumCoef <- apply(as.matrix(fitCV$glmnet.fit$beta), 2, 
                     function(x) sum((x != 0)))
    newCut <- which.max(sumCoef >= minscreen)
    whichVariable <- (as.matrix(fitCV$glmnet.fit$beta)[, 
                                                       newCut] != 0)
  }
  return(whichVariable)
}

sl_lib <- c(ranger_learner$names,
            glmnet_learner$names, 
            xgboost_learner$names, 
            earth_learner$names,
            list(mean_learner,
                 glm_learner, 
                 c(ranger_learner$names[1], "screen.glmnet1"),
                 c(ranger_learner$names[2], "screen.glmnet1"),
                 c(ranger_learner$names[3], "screen.glmnet1"),
                 c(glmnet_learner$names[1], "screen.glmnet1"),
                 c(glmnet_learner$names[2], "screen.glmnet1"),
                 c(glmnet_learner$names[3], "screen.glmnet1"),
                 c(glmnet_learner$names[4], "screen.glmnet1"),
                 c(glmnet_learner$names[5], "screen.glmnet1"),
                 c(xgboost_learner$names[1], "screen.glmnet1"),
                 c(xgboost_learner$names[2], "screen.glmnet1"),
                 c(xgboost_learner$names[3], "screen.glmnet1"),
                 c(xgboost_learner$names[4], "screen.glmnet1"),
                 c(earth_learner$names[1], "screen.glmnet1"),
                 c(earth_learner$names[2], "screen.glmnet1"),
                 c(earth_learner$names[3], "screen.glmnet1")))

# Specify the number of folds for V-fold cross-validation
# Use same folds as used for causal_forest function
# Doing cross-validation this way automatically deploys cross-fitting
fold_dat <- tibble(id = 1:n,folds)
fold_index <- split(fold_dat$id,fold_dat$folds)

# we'll use parallel processing to speed things up
options(mc.cores = detectCores() - 2)

getOption("mc.cores")

augment_data <- covariates_matrix_w[,env_vars]*covariates_matrix_w$DOR
names(augment_data) <- paste0(names(augment_data),"_DOR")

covariates_matrix_w_augment <- cbind(covariates_matrix_w, augment_data)

fit_mu <- CV.SuperLearner(Y = outcome,
                          X = covariates_matrix_w_augment, 
                          method = "method.NNLS", 
                          family = gaussian,
                          SL.library = sl_lib,
                          cvControl = list(V = num.folds, validRows = fold_index),
                          control = list(saveCVFitLibrary = T),
                          parallel = "multicore",
                          verbose = T)

fit_pi <- CV.SuperLearner(Y = exposure,
                          X = covariates_matrix,
                          method = "method.NNLS", 
                          family = binomial,
                          SL.library = sl_lib,
                          cvControl = list(V = num.folds, validRows = fold_index),#, stratifyCV = TRUE),
                          control = list(saveCVFitLibrary = T),
                          parallel = "multicore",
                          verbose = T)

## cross-fit predictions

pscore <- as.matrix(fit_pi$SL.predict)

mu_hat <- as.matrix(fit_mu$SL.predict)

mu_hat1 <- NULL
for(i in 1:num.folds){
  mu_hat1 <- rbind(mu_hat1, 
                   predict(fit_mu$AllSL[[i]],
                           newdata = base::transform(
                             covariates_matrix_w_augment[fold_index[[i]],], DOR = 1), 
                           onlySL=T)$pred)
}

mu_hat0 <- NULL
for(i in 1:num.folds){
  mu_hat0 <- rbind(mu_hat0, 
                   predict(fit_mu$AllSL[[i]],
                           newdata = base::transform(
                             covariates_matrix_w_augment[fold_index[[i]],], DOR = 0), 
                           onlySL=T)$pred)
}

## aipw
aipw_func <- function(exposure, outcome, pscore, mu_hat, mu_hat0, mu_hat1){
  aipw_score <- ((2*exposure - 1)*(outcome - mu_hat))/((2*exposure - 1)*pscore + (1 - exposure)) + (mu_hat1 - mu_hat0)
  return(aipw_score)
}

aipw_score <- aipw_func(exposure, 
                        outcome, 
                        pscore, 
                        mu_hat, 
                        mu_hat0, 
                        mu_hat1)
colnames(aipw_score) <- NULL

aipw_psi <- mean(aipw_score)

aipw_se <- sd(aipw_score)/sqrt(n)

aipw_ate <- c(aipw_psi, aipw_se)

aipw_ate <- cbind(t(aipw_ate), 
                  aipw_ate[1] - 1.96*aipw_ate[2],
                  aipw_ate[1] + 1.96*aipw_ate[2])