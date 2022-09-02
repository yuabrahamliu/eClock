
#Single elastic net, single ensemble, single balancing####
#cross elastic net, cross ensemble, cross balancing, ensemble prediction

unregister_dopar <- function(){
  
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
  
}

singlenet <- function(trainvar, trainres, 
                      threads = 1, alphas = c(0.5), seednum = 1, 
                      foldnum = 10, errortype = 'min'){
  
  #Elastic net on bootstrapped/SMOTE train dataset
  #library(glmnet)
  #library(doParallel)
  #library(foreach)
  
  a <- alphas
  #Tune the value of alpha through a line search with the parallelism
  
  for(k in 1:length(a)){
    l <- a[k]
    
    if(threads > 1){
      
      parallelval <- TRUE
      cores <- parallel::detectCores()
      cl <- parallel::makeCluster(min(threads, cores))
      
      doParallel::registerDoParallel(cl)
      
    }else{
      parallelval <- FALSE
    }
    
    set.seed(seednum)
    cv <- glmnet::cv.glmnet(x = trainvar, y = trainres, family = "gaussian", 
                            nfold = foldnum, type.measure = "deviance", 
                            paralle = parallelval, alpha = l)
    
    if(threads > 1){
      parallel::stopCluster(cl)
      unregister_dopar()
    }
    
    
    if(k == 1){
      cvms.1ses <- c(cv$cvm[cv$lambda == cv$lambda.1se])
      cvms.mins <- c(cv$cvm[cv$lambda == cv$lambda.min])
      lambda.1ses <- c(cv$lambda.1se)
      lambda.mins <- c(cv$lambda.min)
      alphas <- c(l)
    }else{
      cvms.1ses <- c(cvms.1ses, cv$cvm[cv$lambda == cv$lambda.1se])
      cvms.mins <- c(cvms.mins, cv$cvm[cv$lambda == cv$lambda.min])
      lambda.1ses <- c(lambda.1ses, cv$lambda.1se)
      lambda.mins <- c(lambda.mins, cv$lambda.min)
      alphas <- c(alphas, l)
    }
    
  }
  
  search <- data.frame(cvms.1ses = cvms.1ses, cvms.mins = cvms.mins, 
                       lambda.1ses = lambda.1ses, lambda.mins = lambda.mins, 
                       alphas = alphas)
  
  if(errortype == 'min'){
    parameters <- search[search$cvms.mins == min(search$cvms.mins),]
    cat(paste0('Choose the elastic net mixing parameter (alpha) as ', 
               parameters$alphas, '\n'))
    cat(paste0('Choose the regularization constant (lambda) as ', 
               signif(parameters$lambda.mins, 3), '\n'))
    
    Alpha <- parameters$alphas
    Lambda <- parameters$lambda.mins
  }else{
    parameters <- search[search$cvms.1ses == min(search$cvms.1ses),]
    cat(paste0('Choose the elastic net mixing parameter (alpha) as ', 
               parameters$alphas, '\n'))
    cat(paste0('Choose the regularization constant (lambda) as ', 
               signif(parameters$lambda.1ses, 3), '\n'))
    Alpha <- parameters$alphas
    Lambda <- parameters$lambda.1ses
  }
  
  #Elastic net##########
  
  elasticnetmodel <- glmnet::glmnet(x = trainvar, y = trainres, 
                                    family = "gaussian", 
                                    lambda = Lambda, 
                                    alpha = Alpha)
  modelcoefs <- coef(elasticnetmodel)
  modelcoefs <- as.matrix(modelcoefs)
  
  reslist <- list(elasticnetmodel = elasticnetmodel, modelcoefs = modelcoefs)
  
  return(reslist)
  
}

coefselection <- function(modelcoeflist){
  
  modelcoefs <- do.call(cbind, modelcoeflist)
  coefmean <- rowMeans(modelcoefs)
  coefplus <- rowSums(modelcoefs > 0)
  coefminus <- rowSums(modelcoefs < 0)
  coefsign <- coefplus - coefminus
  coefscore <- coefsign*coefmean
  coefscoresd <- sd(coefscore)
  coefscoremean <- mean(coefscore)
  coefscoreupper <- coefscoremean + 2*coefscoresd
  coefscorelower <- coefscoremean - 2*coefscoresd
  
  finalfeatures <- coefscore[coefscore < coefscorelower | coefscore > coefscoreupper]
  
  finalfeatures <- as.matrix(finalfeatures)
  finalfeatures <- finalfeatures[-1,]
  
  finalfeatures <- finalfeatures[order(-finalfeatures)]
  
  return(finalfeatures)
  
}

singlelinear <- function(finalfeatures, 
                         resamplevar, 
                         resampleresponse){
  
  #Model training on resampled data
  subresamplevar <- resamplevar[,finalfeatures]
  
  subresamplevar <- cbind(subresamplevar, resampleresponse$Response)
  colnames(subresamplevar)[ncol(subresamplevar)] <- 'Response'
  subresamplevar <- as.data.frame(subresamplevar)
  
  fit <- lm(Response ~ ., data = subresamplevar)
  
  return(fit)
  
}

weightensemble <- function(innertrainvar, 
                           innertrainres, 
                           innertestvar, 
                           innertestres, 
                           modellist){
  
  featurenames <- names(modellist[[1]]$coefficients)[-1]
  featurenames <- gsub(pattern = "`", replacement = '', x = featurenames)
  subtrain <- innertrainvar[,featurenames]
  subtest <- innertestvar[,featurenames]
  
  subtrain <- as.data.frame(subtrain, stringsAsFactors = FALSE)
  subtest <- as.data.frame(subtest, stringsAsFactors = FALSE)
  
  trainprelist <- list()
  testprelist <- list()
  for(m in 1:length(modellist)){
    subfit <- modellist[[m]]
    
    trainpre <- predict(subfit, subtrain)
    testpre <- predict(subfit, subtest)
    
    trainprelist[[m]] <- trainpre
    testprelist[[m]] <- testpre
    
  }
  
  trainpreres <- do.call(cbind, trainprelist)
  testpreres <- do.call(cbind, testprelist)
  
  colnames(trainpreres) <- paste0('model_', 1:ncol(trainpreres))
  colnames(testpreres) <- paste0('model_', 1:ncol(testpreres))
  
  trainpccs <- cor(trainpreres, innertrainres)
  testpccs <- cor(testpreres, innertestres)
  modelpccs <- cbind(trainpccs, testpccs)
  colnames(modelpccs) <- c('training', 'testing')
  
  #Only keep base learners with an R^2 > 0.5 in training data
  Rsquare <- 0.5
  
  if(sum(modelpccs[,1] > sqrt(Rsquare)) == 0){
    return(NULL)
  }
  
  modellist <- modellist[modelpccs[,1] > sqrt(Rsquare)]
  trainpreres <- trainpreres[,modelpccs[,1] > sqrt(Rsquare)]
  testpreres <- testpreres[,modelpccs[,1] > sqrt(Rsquare)]
  
  modelpccs <- modelpccs[modelpccs[,1] > sqrt(Rsquare),]
  
  if(!is.matrix(trainpreres)){
    
    trainpreres <- as.matrix(trainpreres)
    testpreres <- as.matrix(testpreres)
    colnames(trainpreres) <- colnames(testpreres) <- 'model_1'
    
    modelpccs <- t(as.matrix(modelpccs))
    rownames(modelpccs) <- 'model_1'
  }
  
  #Calcuate the weight for each base learner only from training data
  weights <- 0.5 * log(modelpccs[,1]^2/(1 - modelpccs[,1]^2))
  weights <- weights/sum(weights)
  weights <- as.matrix(weights)
  weights <- t(weights)
  
  
  
  #Ensemble the prediction from all base learners
  trainpreres <- trainpreres %*% t(weights)
  testpreres <- testpreres %*% t(weights)
  rownames(weights) <- 'weights'
  
  if(is.null(colnames(weights))){
    colnames(weights) <- 'model_1'
  }
  
  colnames(trainpreres) <- colnames(testpreres) <- 'Prediction'
  
  trainpcc <- cor(trainpreres[,1], innertrainres)
  testpcc <- cor(testpreres[,1], innertestres)
  
  traincomp <- cbind(trainpreres, innertrainres)
  testcomp <- cbind(testpreres, innertestres)
  colnames(traincomp)[2] <- colnames(testcomp)[2] <- 'True'
  
  reslist <- list(baselearners = modellist, 
                  baseweights = weights, 
                  traincomp = traincomp, 
                  testcomp = testcomp)
  
  return(reslist)
  
}

prescreenfun <- function(oripd, orivar){
  
  orivar <- orivar[oripd$sampleid, , drop = FALSE]
  
  cor.pval <- function(x, y){
    
    pvalue <- cor.test(x = x, y = y)$p.value
    
    return(pvalue)
  }
  
  pccpvals <- apply(X = orivar, MARGIN = 2, 
                    FUN = cor.pval, y = oripd$Response)
  
  pccpvals <- pccpvals[!is.na(pccpvals)]
  
  if(length(pccpvals) == 0){
    return(orivar)
  }
  
  screenfeatures <- pccpvals[pccpvals < 0.05]
  
  if(length(screenfeatures) == 0){
    return(orivar)
  }
  
  screenvar <- orivar[, names(screenfeatures), drop = FALSE]
  
  return(screenvar)
}

#'Train an elastic net model to fit the response and do feature selection
#'
#'Divide the original dataset into training and testing for one time and then 
#'  train an elastic net model and do feature selection.
#'
#'@param orivar A matrix recording the values of candidate features for all 
#'  samples (training samples + testing samples). Each row represents one 
#'  sample and each column represents one feature. The column names are the 
#'  feature names while the row names should be sample IDs. Each column 
#'  records a series of numeric values.
#'@param oripd A data.frame indicating the meta information of all samples 
#'  (training samples + testing samples). Each row represents one sample and 
#'  each column represents one feature. Two columns must be included in this 
#'  data.frame, one is a column named "sampleid", which records the unique IDs 
#'  of samples, and the other is a column named "Response", which records the 
#'  values of the response variable that will be used to construct the 
#'  downstream elastic net model, and this column can be a factor, or record 
#'  numeric values. If it is a factor, it will be converted to numeric values.
#'@param fold A float number between 0 and 1, indicating the proportion of 
#'  total samples that will be included in the training dataset, so that the 
#'  proportion of total samples that will be included in the testing dataset 
#'  is \code{1 - fold}. Default is 0.9.
#'@param seednum Random seed need to be set for dividing the dataset into 
#'  training and testing randomly. Default if 1.
#'@param cores Number of threads need to be used to parallelize the computing. 
#'  Default is 1. 
#'@param alphas A vector with the candidate elastic net mixing parameters, 
#'  which should be between 0 and 1. If the mixing parameter (alpha) is 1, the 
#'  regularization is actually a LASSO regularization. If alpha is 0, the 
#'  regularization is actually a ridge regularization. Default is 0.5.
#'@param errortype The method to automatically choose the regularization 
#'  constant lambda and the elastic net mixing parameter alpha. If it is 
#'  "min", the lambda and alpha combination giving minimum mean 10-fold cross-
#'  validation error on the training set will be used, while if it is "1ses", 
#'  the lambda and alpha combination giving the cross-validation error within 
#'  one standard error of the minimum will be used. Default is "min".
#'@param localrsquare Whether the local R squares of several small sample 
#'  groups within specific response value intervals need to be calculated. 
#'  Default is FALSE.
#'@param plotting Whether several plots need to be generated to show the 
#'  performance of the model, including the scatter plots, heatmap plots, and 
#'  residual plots on both training and testing datasets. Default is FALSE.
#'@param balanceidx Whether or not to calculate the balance index of the 
#'  original training data response distribution.
#'@param binwith If the parameter code{balanceidx} is set as TRUE, this 
#'  parameter will be needed to set the start value of the bin width to divide 
#'  the response variable and to measure this distribution. Its default value 
#'  is 1.
#'@param prescreen If this parameter is set as TRUE, before the elastic net 
#'  model construction, a preliminary screen on the features included in 
#'  \code{orivar} will be performed first according to the p-value of their 
#'  Pearson correlation coefficients to the response variable, and the ones 
#'  with such a p-value less than 0.05 will be selected to train the elastic 
#'  net model. Default is FALSE.
#'@return A list with several slots. The slot named "elasticmodel" contains 
#'  an S3 object for the final model with the optimal lambda and alpha values. 
#'  The slot "modelcoefs" records the regression coefficients for all the 
#'  candidate features as well as the intercept of the model. The slots 
#'  "traincomp" and "testcomp" are data.frames with the predicted and true 
#'  response values for training and testing datasets. The column names are 
#'  "Prediction" and "True", the row names are the sample IDs. If the 
#'  parameter \code{plotting} is set as TRUE, several plots reflecting the 
#'  model performance will also be returned.
#'@export
singleselection <- function(orivar, 
                            oripd, 
                            fold = 0.9, 
                            seednum = 1, 
                            cores = 1, 
                            alphas = c(0.5), 
                            errortype = 'min', 
                            localrsquare = FALSE, 
                            plotting = FALSE, 
                            balanceidx = FALSE, 
                            binwith = 1, 
                            prescreen = FALSE){
  
  if(prescreen == TRUE){
    
    orivar <- prescreenfun(oripd = oripd, orivar = orivar)
  }
  
  #Part1, divide the original dataset into train and test sets
  sampledivide <- dividedat(i = seednum, fold = fold, totalpd = oripd, totalvar = orivar)
  if(is.null(sampledivide)){
    return(NULL)
  }
  
  
  if(balanceidx == TRUE){
    
    balanceindex <- resamplebin(resrange = sampledivide$trainset$innertrainpd, 
                                binwith = binwith, 
                                sampleseed = j, 
                                samplevar = sampledivide$trainset$innertrainvar, 
                                upsampling = 'SMOTE', 
                                plotting = FALSE)
    
  }
  
  #Part2, elastic net on train dataset
  foldnum <- max(min(floor(nrow(sampledivide$trainset$innertrainvar)/10), 10), 3)
  
  elasticres <- singlenet(trainvar = sampledivide$trainset$innertrainvar, 
                          trainres = sampledivide$trainset$innertrainpd$Response, 
                          threads = cores, alphas = alphas, seednum = seednum, 
                          foldnum = foldnum, errortype = errortype)
  
  #Part3, summarize prediction on train and test datasets
  trainpre <-predict(elasticres$elasticnetmodel, 
                     s = elasticres$elasticnetmodel$lambda, 
                     newx = sampledivide$trainset$innertrainvar)
  
  testpre <- predict(elasticres$elasticnetmodel, 
                     s = elasticres$elasticnetmodel$lambda, 
                     newx = sampledivide$testset$innertestvar)
  
  traincomp <- cbind(trainpre, sampledivide$trainset$innertrainresponse)
  colnames(traincomp) <- c('Prediction', 'True')
  
  testcomp <- cbind(testpre, sampledivide$testset$innertestpd$Response)
  colnames(testcomp) <- c('Prediction', 'True')
  
  res <- list(elasticmodel = elasticres$elasticnetmodel, 
              modelcoeffs = elasticres$modelcoefs, 
              traincomp = traincomp, 
              testcomp = testcomp)
  
  if(balanceidx == TRUE){
    
    res$balanceidx <- balanceindex$balanceidx
    
  }
  
  if(plotting == TRUE){
    featurenum <- sum(res$modelcoeffs[-1,] != 0)
    
    scatterplot(comptab = traincomp, 
                featurenum = featurenum, 
                title = 'Training data', colorful = TRUE)
    
    scatterplot(comptab = testcomp, 
                featurenum = featurenum, 
                title = 'Testing data', colorful = TRUE)
    
    features <- res$modelcoeffs[-1,1]
    features <- features[features != 0]
    features <- features[order(features)]
    
    heatmapplot(orivar = sampledivide$trainset$innertrainvar, 
                comptab = traincomp, 
                featurenames = names(features), 
                title = 'Training data')
    
    heatmapplot(orivar = sampledivide$testset$innertestvar, 
                comptab = testcomp, 
                featurenames = names(features), 
                title = 'Testing data')
    
    residualplot(comptab = traincomp, 
                 featurenum = featurenum, 
                 title = 'Training data', colorful = TRUE)
    
    residualplot(comptab = testcomp, 
                 featurenum = featurenum, 
                 title = 'Testing data', colorful = TRUE)
    
  }
  
  if(localrsquare == TRUE){
    
    trainlocalr <- tryCatch({
      localrplot(comptab = traincomp, 
                 featurenum = featurenum, 
                 title = 'Training data', 
                 stepsize = NULL, 
                 averagepre = FALSE, 
                 plotting = plotting)
    }, error = function(err){
      NULL
    })
    
    testlocalr <- tryCatch({
      localrplot(comptab = testcomp, 
                 featurenum = featurenum, 
                 title = 'Testing data', 
                 stepsize = NULL, 
                 averagepre = FALSE, 
                 plotting = plotting)
    }, error = function(err){
      NULL
    })
    
    res$trainlocalrsquare <- trainlocalr
    res$testlocalrsquare <- testlocalr
    
  }
  
  return(res)
  
}

#'Calculate model predicted response values for the input data
#'
#'Predict response values for the input data using a trained elastic net model 
#'  or ensemble model.
#'
#'@param balanceres The model trained. The result generated by the function 
#'  \code{singlenet} or \code{singlebalance} can be directly used for this 
#'  parameter, while for the function \code{crosstraining}, the element in its 
#'  sub-list named "modellist" could also be used for this parameter. 
#'@param vardat The input data whose response values need to be predicted by 
#'  the model. Each row is a sample and each column is a feature (e.g. a 
#'  methylation probe). The values in \code{vardat} should be numeric values. 
#'@return A one-column matrix containing the model predicted response values. 
#'  Column name of this matrix is "Prediction" and each row of it represents a 
#'  sample.
#'@export
ensemblepredict <- function(balanceres, 
                            vardat){
  
  vardattmp <- vardat
  #vardat <- vardat[,rownames(balanceres$modelcoeffs)[-1]]
  
  
  if('elasticmodel' %in% names(balanceres)){
    vardat <- vardattmp[,rownames(balanceres$modelcoeffs)[-1]]
    
    vardat <- as.matrix(vardat)
    
    baselearner <- balanceres$elasticmodel
    baselearnerpre <- glmnet::predict.glmnet(object = baselearner, 
                                             newx = vardat)
    baselearnerpre <- as.matrix(baselearnerpre)
    
    ensembleres <- baselearnerpre
    
  }else{
    
    num <- length(balanceres$baselearners)
    for(o in 1:num){
      
      baselearner <- balanceres$baselearners[[o]]
      
      featurenames <- gsub(pattern = '`', replacement = '',
                           x = names(baselearner$coefficients)[-1])

      vardat <- vardattmp[,featurenames]
      vardat <- as.data.frame(vardat, stringsAsFactors = FALSE)
      
      baselearnerpre <- predict(object = baselearner, newdata = vardat)
      baselearnerpre <- as.matrix(baselearnerpre)
      if(o == 1){
        baselearnerpres <- baselearnerpre
      }else{
        baselearnerpres <- cbind(baselearnerpres, baselearnerpre)
      }
      
    }
    colnames(baselearnerpres) <- colnames(balanceres$baseweights)
    
    ensembleres <- baselearnerpres %*% t(balanceres$baseweights)
    
  }
  
  colnames(ensembleres)[1] <- 'Prediction'
  return(ensembleres)
  
}

#'Train an ensemble elastic net model to fit the response and do feature 
#'  selection
#'
#'Divide the original dataset into training and testing for one time and then 
#'  train an ensemble elastic net model and do feature selection.
#'
#'@param orivar A matrix recording the values of candidate features for all 
#'  samples (training samples + testing samples). Each row represents one 
#'  sample and each column represents one feature. The column names are the 
#'  feature names while the row names should be sample IDs. Each column 
#'  records a series of numeric values.
#'@param oripd A data.frame indicating the meta information of all samples 
#'  (training samples + testing samples). Each row represents one sample and 
#'  each column represents one feature. Two columns must be included in this 
#'  data.frame, one is a column named "sampleid", which records the unique IDs 
#'  of samples, and the other is a column named "Response", which records the 
#'  values of the response variable that will be used to construct the 
#'  downstream ensemble model, and this column can be a factor, or record 
#'  numeric values. If it is a factor, it will be converted to numeric values.
#'@param fold A float number between 0 and 1, indicating the proportion of 
#'  total samples that will be included in the training dataset, so that the 
#'  proportion of total samples that will be included in the testing dataset 
#'  is \code{1 - fold}. Default is 0.9.
#'@param seednum Random seed need to be set for dividing the dataset into 
#'  training and testing randomly. Default is 1.
#'@param cores Number of threads need to be used to parallelize the computing. 
#'  Default is 1.
#'@param balancing When taking samples to construct base learners, whether a 
#'  balancing sampling strategy, or a normal bootstrapping sampling strategy 
#'  should be used. If it is TRUE, a balancing sampling strategy will be used, 
#'  which means the original response distribution of the training samples 
#'  will be adjusted to a balanced one via over-sampling on the samples with 
#'  low distribution density and under-sampling on the ones with high 
#'  distribution density. If it is FALSE, a normal bootstrapping sampling will 
#'  be used and the original sample distribution will not be adjusted. Default 
#'  is TRUE.
#'@param samplenum The number of sampling rounds need to be taken to construct 
#'  the ensemble model. It can also be understood as the number of base 
#'  learners included in the model. Default is 10.
#'@param binwith If the parameter \code{balancing} or \code{balanceidx} is set 
#'  as TRUE, this parameter will be needed to set the start value of the bin 
#'  width to divide the response variable and to measure this distribution. 
#'  Its default value is 1.
#'@param upsampling If the parameter \code{balancing} is set as TRUE, this 
#'  parameter will be needed to set the method to do over-sampling for the 
#'  samples with a low response distribution density. The over-sampling will 
#'  make the density of these samples increase to the level of other samples. 
#'  Its default value is "SMOTE", which means that the SMOTE (Synthetic 
#'  Minority Over-sampling Technique) method will be used. The value of this 
#'  parameter can also be "boot", which means bootstrapping will be used, but 
#'  it is not recommended because it is more likely to bring over-fitting 
#'  problems.
#'@param alphas A vector with the candidate elastic net mixing parameters, 
#'  which should be between 0 and 1. If the mixing parameter (alpha) is 1, the 
#'  regularization is actually a LASSO regularization. If alpha is 0, the 
#'  regularization is actually a ridge regularization. Default is 0.5.
#'@param errortype The method to automatically choose the regularization 
#'  constant lambda and the elastic net mixing parameter alpha for each 
#'  elastic net base learner. If it is "min", the lambda and alpha combination 
#'  giving minimum mean 10-fold cross-validation error on the sampled training 
#'  set of the base learner will be used, while if it is "1ses", the lambda 
#'  and alpha combination giving the cross-validation error within one 
#'  standard error of the minimum will be used. Default is "min".
#'@param localrsquare Whether the local R squares of several small sample 
#'  groups within specific response value intervals need to be calculated. 
#'  Default is FALSE.
#'@param plotting Whether several plots need to be generated to show the 
#'  performance of the model, including the scatter plots, heatmap plots, and 
#'  residual plots on both training and testing datasets. Default is FALSE.
#'@param balanceidx Whether or not to calculate the balance index of the 
#'  original training data response distribution.
#'@param prescreen If this parameter is set as TRUE, before the model 
#'  construction, a preliminary screen on the features included in 
#'  \code{orivar} will be performed first according to the p-value of their 
#'  Pearson correlation coefficients to the response variable, and the ones 
#'  with such a p-value less than 0.05 will be selected to train the model. 
#'  Default is FALSE.
#'@return A list with several slots. The slot named "baselearners" contains S3 
#'  objects for the base learners. The slot "baseweights" records the weights 
#'  of the base learners when using them to construct the ensemble model. The 
#'  slots "traincomp" and "testcomp" are data.frames with the predicted and 
#'  true response values for training and testing datasets. The column names 
#'  are "Prediction" and "True", the row names are the sample IDs. The slot 
#'  "modelscores" records the features finally selected and their scores 
#'  reflecting their importance. The slot "modelcoefs" records the regression 
#'  coefficients for all the selected features as well as the intercepts in 
#'  the base learners. The slot "balanceidx" is the balance index of the 
#'  response distribution of the original training dataset. If the parameter 
#'  \code{plotting} is set as TRUE, several plots reflecting the ensemble 
#'  model performance will also be returned.
#'@export
singlebalance <- function(orivar, 
                          oripd, 
                          fold = 0.9, 
                          seednum = 1, 
                          cores = 1, 
                          balancing = TRUE, 
                          samplenum = 10, 
                          binwith = 1, 
                          upsampling = 'SMOTE', 
                          alphas = c(0.5), 
                          errortype = 'min', 
                          #errortype can be '1ses' or 'min'
                          localrsquare = FALSE, 
                          plotting = FALSE, 
                          balanceidx = FALSE, 
                          prescreen = FALSE){
  
  if(prescreen == TRUE){
    
    orivar <- prescreenfun(oripd = oripd, orivar = orivar)
  }
  
  
  #Part1, divide the original dataset into train and test sets
  sampledivide <- dividedat(i = seednum, fold = fold, totalpd = oripd, totalvar = orivar)
  if(is.null(sampledivide)){
    return(NULL)
  }
  
  probelist <- list()
  balancedsamplelist <- list()
  balancedpdlist <- list()
  for(j in 1:samplenum){
    
    #Part2, data balancing
    if(balancing == TRUE){
      
      balancingdat <- resamplebin(resrange = sampledivide$trainset$innertrainpd, 
                                  binwith = binwith, 
                                  sampleseed = j, 
                                  samplevar = sampledivide$trainset$innertrainvar, 
                                  upsampling = upsampling, 
                                  plotting = FALSE, 
                                  balanceidx = balanceidx)
      
    }else{
      
      balancingdat <- simpleboot(resrange = sampledivide$trainset$innertrainpd, 
                                 sampleseed = j, 
                                 samplevar = sampledivide$trainset$innertrainvar, 
                                 plotting = FALSE)
      
      
      if(balanceidx == TRUE){
        
        balanceindex <- resamplebin(resrange = sampledivide$trainset$innertrainpd, 
                                    binwith = binwith, 
                                    sampleseed = j, 
                                    samplevar = sampledivide$trainset$innertrainvar, 
                                    upsampling = 'SMOTE', 
                                    plotting = FALSE, 
                                    balanceidx = balanceidx)
        
        balancingdat$balanceidx <- balanceindex$balanceidx
        
        rm(balanceindex)
        
      }
      
      
      
      
      
      
    }
    
    
    
    
    balancedsamplelist[[j]] <- balancingdat$resamplevars
    balancedpdlist[[j]] <- balancingdat$resampleresponse
    
    #Part3, Elastic net on balanced train dataset
    foldnum <- max(min(floor(nrow(sampledivide$trainset$innertrainvar)/10), 10), 3)
    
    elasticres <- singlenet(trainvar = balancingdat$resamplevars, 
                            trainres = balancingdat$resampleresponse$Response, 
                            threads = cores, alphas = alphas, seednum = j, 
                            foldnum = foldnum, errortype = errortype)
    
    probelist[[j]] <- elasticres$modelcoefs
    
    
  }
  
  
  #Part4, feature integration
  finalprobes <- coefselection(modelcoeflist = probelist)
  
  #Part5, linear base learner training
  modellist <- list()
  for(m in 1:length(balancedsamplelist)){
    
    balancedsamples <- list(resamplevars = balancedsamplelist[[m]], 
                            resampleresponse = balancedpdlist[[m]])
    linearfit <- singlelinear(finalfeatures = names(finalprobes), 
                              resamplevar = balancedsamples$resamplevars, 
                              resampleresponse = balancedsamples$resampleresponse)
    
    modellist[[m]] <- linearfit
    
  }
  
  #Part6, final ensemble model
  resampleensemble <- weightensemble(innertrainvar = sampledivide$trainset$innertrainvar, 
                                     innertrainres = sampledivide$trainset$innertrainpd$Response, 
                                     innertestvar = sampledivide$testset$innertestvar, 
                                     innertestres = sampledivide$testset$innertestpd$Response, 
                                     modellist = modellist)
  if(is.null(resampleensemble)){
    return(NULL)
  }
  
  modelscores <- as.matrix(finalprobes)
  colnames(modelscores)[1] <- 'score'
  if(length(resampleensemble$baselearners) >= 1){
    for(n in 1:length(resampleensemble$baselearners)){
      modelcoeff <- resampleensemble$baselearners[[n]]$coefficients
      names(modelcoeff) <- gsub(pattern = "`", replacement = '', 
                                x = names(modelcoeff))
      modelcoeff <- modelcoeff[c('(Intercept)', row.names(modelscores))]
      modelcoeff <- as.matrix(modelcoeff)
      if(n == 1){
        modelcoeffs <- modelcoeff
      }else{
        modelcoeffs <- cbind(modelcoeffs, modelcoeff)
      }
      
    }
    
    colnames(modelcoeffs) <- paste0('model_', seq(1, ncol(modelcoeffs), 1))
  }
  
  
  modelcoeffs <- modelcoeffs[complete.cases(modelcoeffs),]
  modelcoeffs <- as.matrix(modelcoeffs)
  colnames(modelcoeffs) <- paste0('model_', seq(1, ncol(modelcoeffs), 1))
  modelscores <- modelscores[row.names(modelcoeffs)[row.names(modelcoeffs) != 
                                                      '(Intercept)'],]
  modelscores <- as.matrix(modelscores)
  colnames(modelscores)[1] <- 'score'
  
  resampleensemble$modelscores <- modelscores
  resampleensemble$modelcoeffs <- modelcoeffs
  
  if(balanceidx == TRUE){
    
    resampleensemble$balanceidx <- balancingdat$balanceidx
    
  }
  
  
  #Part7, plotting data
  if(plotting == TRUE){
    
    bootres <- do.call(rbind, balancedpdlist)
    bootres <- bootres$Response
    realres <- sampledivide$trainset$innertrainpd$Response
    
    if(balancing == TRUE){
      if(upsampling == 'SMOTE'){
        bootlabel <- 'SMOTE data'
      }else{
        bootlabel <- 'bootstrapped data'
      }
    }else{
      bootlabel <- 'bootstrappped data'
    }
    
    bootres <- data.frame(Response = bootres, 
                          Group = bootlabel, 
                          stringsAsFactors = FALSE)
    realres <- data.frame(Response = realres, Group = 'real data', 
                          stringsAsFactors = FALSE)
    allres <- rbind(realres, bootres)
    allres$Group <- factor(allres$Group, 
                           levels = c('real data', bootlabel), ordered = TRUE)
    
    #library(ggplot2)
    
    p <- ggplot2::ggplot(allres, ggplot2::aes(x = Response))
    p <- p + ggplot2::geom_histogram(ggplot2::aes(fill = Group, y = ..density..), 
                                     position = 'identity', alpha = 0.3, color = 'black', 
                                     binwidth = 1) + 
      ggplot2::stat_density(geom = 'line', position = 'identity', size = 1.5, 
                            ggplot2::aes(color = Group)) + 
      ggplot2::ggtitle('Response Distribution') + 
      ggplot2::ylab('Density') + ggplot2::xlab('Response') + 
      ggplot2::scale_fill_discrete(guide = FALSE) + 
      ggplot2::scale_color_discrete(guide = FALSE) + 
      ggplot2::facet_grid(Group~., scales = 'free_y') + 
      ggplot2::theme_bw()
    
    print(p)
    
    
    #Scatter plot on SMOTE data
    
    simprobelist <- c()
    for(i in 1:length(probelist)){
      
      subprobelist <- probelist[[i]]
      subprobes <- subprobelist[subprobelist[,1] != 0,]
      subprobes <- names(subprobes)
      
      simprobelist <- union(simprobelist, subprobes)
      
    }
    
    simprobelist <- simprobelist[simprobelist != '(Intercept)']
    simprobelist <- unique(simprobelist)
    
    simbalancedsamplelist <- list()
    for(i in 1:length(balancedsamplelist)){
      
      subbalancedsamplelist <- balancedsamplelist[[i]]
      subbalancedsamplelist <- subbalancedsamplelist[,simprobelist]
      simbalancedsamplelist[[i]] <- subbalancedsamplelist
      
    }
    
    balancedsamplelist <- simbalancedsamplelist
    rm(simbalancedsamplelist)
    
    
    bootvar <- do.call(rbind, balancedsamplelist)
    
    bootpre <- ensemblepredict(balanceres = resampleensemble, 
                               vardat = bootvar)
    bootcomp <- cbind(bootpre, bootres$Response)
    colnames(bootcomp) <- c('Prediction', 'True')
    
    scatterplot(comptab = bootcomp, 
                featurenum = nrow(resampleensemble$modelcoeffs) - 1, 
                title = bootlabel, colorful = TRUE)
    
    #Heatmap on SMOTE data
    features <- resampleensemble$modelcoeffs
    features <- as.matrix(features)
    features <- rowSums(features)
    features <- features[-1]
    features <- features[order(features)]
    
    row.names(bootvar) <- row.names(bootcomp)
    heatmapplot(orivar = bootvar, 
                comptab = bootcomp, 
                featurenames = names(features), 
                title = bootlabel)
    
    #Residual plot on SMOTE data
    residualplot(comptab = bootcomp, 
                 featurenum = nrow(resampleensemble$modelcoeffs) - 1, 
                 title = bootlabel, colorful = TRUE)
    
    #Local R square plot on SMOTE data
    if(localrsquare == TRUE){
      
      bootlocalr <- tryCatch({
        localrplot(comptab = bootcomp, 
                   featurenum = nrow(resampleensemble$modelcoeffs) - 1, 
                   title = bootlabel, 
                   stepsize = NULL, 
                   averagepre = FALSE, 
                   plotting = plotting)
      }, error = function(err){
        NULL
      })
    }
    
    
    #Scatter plot on training and testing data
    scatterplot(comptab = resampleensemble$traincomp, 
                featurenum = nrow(resampleensemble$modelcoeffs) - 1, 
                title = 'Training data', colorful = TRUE)
    
    scatterplot(comptab = resampleensemble$testcomp, 
                featurenum = nrow(resampleensemble$modelcoeffs) - 1, 
                title = 'Testing data', colorful = TRUE)
    
    #Heatmap on training and testing data
    
    heatmapplot(orivar = sampledivide$trainset$innertrainvar, 
                comptab = resampleensemble$traincomp, 
                featurenames = names(features), 
                title = 'Training data')
    
    heatmapplot(orivar = sampledivide$testset$innertestvar, 
                comptab = resampleensemble$testcomp, 
                featurenames = names(features), 
                title = 'Testing data')
    
    #Residual plot on training and testing data
    residualplot(comptab = resampleensemble$traincomp, 
                 featurenum = nrow(resampleensemble$modelcoeffs) - 1, 
                 title = 'Training data', colorful = TRUE)
    
    residualplot(comptab = resampleensemble$testcomp, 
                 featurenum = nrow(resampleensemble$modelcoeffs) - 1, 
                 title = 'Testing data', colorful = TRUE)
    
    
  }
  
  #Local R square plot on training and testing data
  if(localrsquare == TRUE){
    
    trainlocalr <- tryCatch({
      localrplot(comptab = resampleensemble$traincomp, 
                 featurenum = nrow(resampleensemble$modelcoeffs) - 1, 
                 title = 'Training data', 
                 stepsize = NULL, 
                 averagepre = FALSE, 
                 plotting = plotting)
    }, error = function(err){
      NULL
    })
    
    testlocalr <- tryCatch({
      localrplot(comptab = resampleensemble$testcomp, 
                 featurenum = nrow(resampleensemble$modelcoeffs) - 1, 
                 title = 'Testing data', 
                 stepsize = NULL, 
                 averagepre = FALSE, 
                 plotting = plotting)
    }, error = function(err){
      NULL
    })
    
    resampleensemble$trainlocalrsquare <- trainlocalr
    resampleensemble$testlocalrsquare <- testlocalr
    
  }
  
  return(resampleensemble)
  
}

#'Cross-training for ensemble or elastic net models
#'
#'Divide the original dataset into training and testing for several rounds and 
#'  then for each round, train an ensemble or elastic net model.
#'
#'@param orivar A matrix recording the values of candidate features for all 
#'  samples (training samples + testing samples). Each row represents one 
#'  sample and each column represents one feature. The column names are the 
#'  feature names while the row names should be sample IDs. Each column 
#'  records a series of numeric values.
#'@param oripd A data.frame indicating the meta information of all samples 
#'  (training samples + testing samples). Each row represents one sample and 
#'  each column represents one feature. Two columns must be included in this 
#'  data.frame, one is a column named "sampleid", which records the unique IDs 
#'  of samples, and the other is a column named "Response", which records the 
#'  values of the response variable that will be used to construct the 
#'  downstream ensemble model, and this column can be a factor, or record 
#'  numeric values. If it is a factor, it will be converted to numeric values.
#'@param fold A float number between 0 and 1, indicating the proportion of 
#'  total samples that will be included in the training dataset for each 
#'  sample dividing round, so that the proportion of total samples that will 
#'  be included in the testing dataset is \code{1 - fold}. Default is 0.9.
#'@param roundnum Number of rounds need to divide the original dataset into 
#'  training and testing to train the models. It can also be understood as the 
#'  number of ensemble or elastic net models need to be trained. Default value 
#'  is 10.
#'@param cores Number of threads need to be used to parallelize the computing. 
#'  Default is 1.
#'@param learner What kind of models need to be trained. If it is set as 1, 
#'  the models will be ensemble models with elastic net base learners, and for 
#'  each base learner, its training sample distribution will be adjusted to a 
#'  balance status via over-sampling and under-sampling. The over-sampling 
#'  will be performed using a SMOTE method, while the under-sampling will be 
#'  performed using a bootstrapping method. If it is set as 2, also ensemble 
#'  models with distribution adjustment will be trained. Its difference from 1 
#'  is that both the over-sampling and under-sampling are performed using a 
#'  bootstrapping strategy. If it is set as 3, ensemble models without sample 
#'  distribution adjustment will be trained. If this parameter is 4, simple 
#'  elastic net models will be trained. Default is 1.
#'@param samplenum If set \code{learner} as 1, 2, or 3, i.e., the models are 
#'  set as ensemble models, for each dividing round, after getting the 
#'  training dataset, several rounds of sampling need to be further performed 
#'  on the training dataset to train several base learners and then construct 
#'  the ensemble model. It can also be understood as the number of base 
#'  learners included in one ensemble model generated from one round of 
#'  training/testing division. Default is 10.
#'@param binwith If the parameter \code{learner} is set as 1 or 2, i.e. the 
#'  training sample distribution need to be adjusted before training a base 
#'  learner, this parameter will be needed to set the start value of the bin 
#'  width to divide the response variable and to measure and adjust this 
#'  distribution. Its default value is 1.
#'@param alphas A vector with the candidate elastic net mixing parameters, 
#'  which should be between 0 and 1. If the mixing parameter (alpha) is 1, the 
#'  regularization is actually a LASSO regularization. If alpha is 0, the 
#'  regularization is actually a ridge regularization. Default is 0.5.
#'@param errortype The method to automatically choose the regularization 
#'  constant lambda and the elastic net mixing parameter alpha for each 
#'  elastic net model. If it is "min", the lambda and alpha combination 
#'  giving minimum mean 10-fold cross-validation error on the specific 
#'  training set of the model will be used, while if it is "1ses", the lambda 
#'  and alpha combination giving the cross-validation error within one 
#'  standard error of the minimum will be used. Default is "min".
#'@param localrsquare Whether the local R squares of several small sample 
#'  groups within specific response value intervals need to be calculated. 
#'  Default is FALSE.
#'@param plotting Whether several plots need to be generated to show the 
#'  performance of the models. Each model for each training/testing division 
#'  round will get a set of plots, including the scatter plots, heatmap plots, 
#'  and residual plots on both training and testing datasets. Default value 
#'  is FALSE.
#'@param prescreen If this parameter is set as TRUE, before the model 
#'  construction, a preliminary screen on the features included in 
#'  \code{orivar} will be performed first according to the p-value of their 
#'  Pearson correlation coefficients to the response variable, and the ones 
#'  with such a p-value less than 0.05 will be selected to train the model. 
#'  Default is FALSE.
#'@return A list with several slots. The slot named "modellist" is a sub-list, 
#'  which contains the model training results for all the training/testing 
#'  division rounds, and each of its elements corresponds to a single round. 
#'  The slots "traininglist" and "testinglist" record the sample IDs of the 
#'  training and testing sets for each round. The slots "trainingsummary" and 
#'  "testingsummary" record the metrics reflecting the model performance on 
#'  training and testing sets for the rounds, including the R square and the 
#'  MSE (mean squared error) of the models. If the parameter \code{plotting} 
#'  is set as TRUE, several plots reflecting the model performances will also 
#'  be returned.
#'@export
crosstraining <- function(orivar, 
                          oripd, 
                          fold = 0.9, 
                          roundnum = 10, 
                          cores = 1, 
                          learner = 1, 
                          #learner can be 1 (for balancing model with SMOTE upsampling), 
                          #2 (for balancing model with bootstrapping upsampling), 
                          #3 (for bootstrapping model), and 
                          #4 (for simple elasticnet model)
                          samplenum = 10, 
                          binwith = 1, 
                          alphas = c(0.5), 
                          errortype = 'min', 
                          #errortype can be '1ses' or 'min'
                          localrsquare = FALSE, 
                          plotting = FALSE, 
                          prescreen = FALSE){
  
  if(prescreen == TRUE){
    
    orivar <- prescreenfun(oripd = oripd, orivar = orivar)
  }
  
  
  machinelist <- list()
  trainsamplelist <- list()
  testsamplelist <- list()
  
  trainpcclist <- c()
  testpcclist <- c()
  trainresiduallist <- c()
  testresiduallist <- c()
  
  if(localrsquare == TRUE){
    
    trainlocalrlist <- c()
    testlocalrlist <- c()
    
  }
  
  
  if(learner == 1){
    balancing <- TRUE
    upsampling <- 'SMOTE'
  }else if(learner == 2){
    balancing <- TRUE
    upsampling <- 'boot'
  }else{
    balancing <- FALSE
    upsampling <- NULL
  }
  
  for(i in 1:roundnum){
    
    if(learner == 1 | learner == 2 | learner == 3){
      
      singleround <- singlebalance(orivar = orivar, 
                                   oripd = oripd, 
                                   fold = fold, 
                                   seednum = i, 
                                   cores = cores, 
                                   balancing = balancing, 
                                   samplenum = samplenum, 
                                   binwith = binwith, 
                                   upsampling = upsampling, 
                                   alphas = alphas, 
                                   errortype = errortype, 
                                   #errortype can be '1ses' or 'min'
                                   localrsquare = localrsquare, 
                                   plotting = FALSE)
      
    }else{
      
      singleround <- singleselection(orivar = orivar, 
                                     oripd = oripd, 
                                     fold = fold, 
                                     seednum = i, 
                                     cores = cores, 
                                     alphas = alphas, 
                                     errortype = errortype, 
                                     localrsquare = localrsquare, 
                                     plotting = FALSE)
      
    }
    
    if(is.null(singleround)){
      return(NULL)
    }
    
    trainpcc <- cor(singleround$traincomp[,1], 
                    singleround$traincomp[,2])
    testpcc <- cor(singleround$testcomp[,1], 
                   singleround$testcomp[,2])
    
    
    trainresidual <- mean((singleround$traincomp[,2] - 
                             singleround$traincomp[,1])^2)
    testresidual <- mean((singleround$testcomp[,2] - 
                            singleround$testcomp[,1])^2)
    
    trainpcclist <- c(trainpcclist, trainpcc)
    testpcclist <- c(testpcclist, testpcc)
    
    trainresiduallist <- c(trainresiduallist, trainresidual)
    testresiduallist <- c(testresiduallist, testresidual)
    
    if(localrsquare == TRUE){
      
      trainlocalr <- mean(singleround$trainlocalrsquare$localr)
      testlocalr <- mean(singleround$testlocalrsquare$localr)
      
      trainlocalrlist <- c(trainlocalrlist, trainlocalr)
      testlocalrlist <- c(testlocalrlist, testlocalr)
      
    }
    
    machinelist[[i]] <- singleround
    
    trainsamples <- row.names(singleround$traincomp)
    testsamples <- row.names(singleround$testcomp)
    
    trainsamplelist[[i]] <- trainsamples
    testsamplelist[[i]] <- testsamples
    
    cat(paste0('Completed round ', i, '\n'))
    
  }
  
  trainsummary <- data.frame(trainrsqare = trainpcclist^2, 
                             trainmse = trainresiduallist, 
                             stringsAsFactors = FALSE)
  row.names(trainsummary) <- paste0('model_', 1:nrow(trainsummary))
  
  testsummary <- data.frame(testrsquare = testpcclist^2, 
                            testmse = testresiduallist, 
                            stringsAsFactors = FALSE)
  row.names(testsummary) <- paste0('model_', 1:nrow(testsummary))
  
  if(localrsquare == TRUE){
    
    trainsummary$trainlocalrsquare <- trainlocalrlist
    testsummary$testlocalrsquare <- testlocalrlist
    
  }
  
  
  if(plotting == TRUE){
    
    summaryplot <- function(plotdat = trainsummary, title = 'Training data'){
      
      names(plotdat) <- NA
      
      if(ncol(plotdat) == 3){
        dat <- rbind(plotdat[1], plotdat[3], plotdat[2])
        
        dat$metric <- rep(x = c('R square', 'Local R square', 'Abs residual'), 
                          each = nrow(plotdat))
        dat$model <- rep(x = row.names(plotdat), times = 3)
        
        dat$metric <- factor(dat$metric, levels = c('R square', 
                                                    'Local R square', 
                                                    'Abs residual'), 
                             ordered = TRUE)
        
      }else{
        dat <- rbind(plotdat[1], plotdat[2])
        
        dat$metric <- rep(x = c('R square', 'Abs residual'), 
                          each = nrow(plotdat))
        dat$model <- rep(x = row.names(plotdat), times = 2)
        
        dat$metric <- factor(dat$metric, levels = c('R square', 
                                                    'Abs residual'), 
                             ordered = TRUE)
        
      }
      
      names(dat)[1] <- 'value'
      
      row.names(dat) <- 1:nrow(dat)
      
      p <- ggplot2::ggplot(dat, ggplot2::aes(x = metric, y = value, 
                                             fill = metric))
      p <- p + ggplot2::geom_boxplot(alpha = 0.3) + 
        ggplot2::geom_jitter(shape = 21) + 
        ggplot2::xlab('') + ggplot2::ylab('Value') + 
        ggplot2::theme_bw() + 
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), 
                       axis.ticks.x = ggplot2::element_blank()) + 
        ggplot2::ggtitle(title) + 
        ggplot2::scale_fill_manual(values = scales::hue_pal()(ncol(plotdat)), 
                                   guide = FALSE) + 
        ggplot2::facet_wrap(~metric, scales = 'free')
      
      print(p)
      
    }
    
    summaryplot(plotdat = trainsummary, title = 'Training data')
    summaryplot(plotdat = testsummary, title = 'Testing data')
  }
  
  res <- list(modellist = machinelist, 
              traininglist = trainsamplelist, 
              testinglist = testsamplelist, 
              trainingsummary = trainsummary, 
              testingsummary = testsummary)
  
  return(res)
  
}

#Methylation data preprocessing, probe extraction, probe annotation###
#methylation gene summary, methylation DMR###
