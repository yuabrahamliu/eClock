
#Divide training and testing#####

#'Divide a dataset into training dataset and testing dataset
#'
#'Divide a dataset into training dataset and testing dataset randomly.
#'
#'@param i Random seed need to be set for dividing the dataset randomly. 
#'  Default is 1.
#'@param fold A float number between 0 and 1, indicating the proportion of 
#'  total samples that will be included in the training dataset, so that the 
#'  proportion of total samples that will be included in the testing dataset 
#'  is \code{1 - fold}. Default is 0.9.
#'@param totalpd A data.frame indicating the meta information of all samples. 
#'  Each row represents one sample and each column represents one feature. Two 
#'  columns must be included in this data.frame, one is a column named 
#'  "sampleid", which records the unique IDs of samples, and the other is a 
#'  column named "Response", which records the values of the response variable 
#'  that will be used to construct the downstream machine learning model, and 
#'  this column can be a factor, or record numeric values. If it is a factor, 
#'  it will be converted to numeric values.
#'@param totalvar A matrix recording the values of candidate features for all 
#'  samples. Each row represents one sample and each column represents one 
#'  feature. The column names are the feature names while the row names should 
#'  be sample IDs. Each column records a series of numeric values.
#'@return A list with tow sub-lists. One is a sub-list named "trainset" 
#'  recording the information of samples that have been attributed to the 
#'  training dataset. It contains 3 slots. The one named "innertrainvar" has 
#'  the same structure as the data transferred to the parameter 
#'  \code{totalvar}, another one named "innertrainpd" including the meta 
#'  information for the training samples similar to the data transferred to 
#'  \code{totalpd}. The third one named "innertrainresponse" is a one-column 
#'  frame with the value of the response variable as its only column and the 
#'  row names are training sample IDs. The second sub-list is named "testset", 
#'  and it contains 2 slots named "innertestvar" and "innertestpd" that record 
#'  the candidate feature values and meta information for samples in the 
#'  testing dataset.
#'@export
dividedat <- function(i = 1, fold = 0.9, totalpd, totalvar){
  
  if(is.factor(totalpd$Response)){
    totalpd$Response <- as.numeric(totalpd$Response)
  }
  
  set.seed(i)
  
  dividelist  <- list()
  
  innertrainsize <- round(nrow(totalpd)*fold)
  innertestsize <- nrow(totalpd) - innertrainsize
  
  if(innertestsize <= 0){
    return(NULL)
  }
  
  innertest_idx <- sample(x = 1:nrow(totalpd), size = innertestsize, replace = FALSE)
  innertestpd <- totalpd[innertest_idx,]
  innertrainpd <- totalpd[-innertest_idx,]
  
  innertrainpd <- innertrainpd[order(innertrainpd$sampleid, innertrainpd$Response),]
  innertestpd <- innertestpd[order(innertestpd$sampleid, innertestpd$Response),]
  row.names(innertrainpd) <- 1:nrow(innertrainpd)
  row.names(innertestpd) <- 1:nrow(innertestpd)
  
  
  innertrainresponse <- innertrainpd[c('sampleid', 'Response')]
  row.names(innertrainresponse) <- innertrainresponse$sampleid
  innertrainresponse <- innertrainresponse[-1]
  innertrainresponse <- t(innertrainresponse)
  innertrainresponse <- t(innertrainresponse)
  
  
  innertrainvar <- totalvar[innertrainpd$sampleid,]
  innertestvar <- totalvar[innertestpd$sampleid,]
  
  #trainsamplelist[[i]] <- innertrainpd$sampleid
  #testsamplelist[[i]] <- innertestpd$sampleid
  
  dividelist[[1]] <- list(innertrainvar = innertrainvar, 
                          innertrainpd = innertrainpd, 
                          innertrainresponse = innertrainresponse)
  dividelist[[2]] <- list(innertestvar = innertestvar, 
                          innertestpd = innertestpd)
  
  names(dividelist) <- c('trainset', 'testset')
  
  return(dividelist)
  
}

#bootstrapping and SMOTE sampling#####

#'Generate data with balanced response distribution from an imbalanced one
#'
#'Use bootstrapping or SMOTE method to sample the original samples and 
#'  generate new data with a balanced response distribution.
#'
#'@param resrange A data.frame with meta information of the samples. Two 
#'  columns must be contained in it. One is the column "sampleid", the other 
#'  is the column "Response", which records the response variable values of 
#'  the samples. The response values should be numeric values.
#'@param binwith The start value of the bin width used to divide the response 
#'  values into different bins and calculate their distribution. Its value 
#'  will be increased a bit automatically for each recursion loop conducted by 
#'  this function. The default value of this parameter is 1.
#'@param minwith The increased value on \code{binwith} for each recursion loop 
#'  of this function. Its default value is \code{NULL}, meaning it is the same 
#'  as the original \code{binwith} value.
#'@param sampleseed The random seed used in this function to do sampling. Its 
#'  default value is 1.
#'@param samplevar A matrix recording the values of features for all samples. 
#'  Each row represents one sample and each column represents one feature. The 
#'  column names are the feature names while the row names should be sample 
#'  IDs same as the ones recorded by the "sampleid" column in \code{resrange}. 
#'  Each column records a series of numeric values.
#'@param upsampling The method to do over-sampling for the samples with a low 
#'  response distribution density. The over-sampling will make the density of 
#'  these samples increase to the level of other samples. Its default value is 
#'  "SMOTE", which means that the SMOTE (Synthetic Minority Over-sampling 
#'  Technique) method will be used. The value of this parameter can also be 
#'  "boot", which means bootstrapping will be used, but it is not recommended 
#'  because it is more likely to bring over-fitting problems.
#'@param plotting Whether or not to plot the sample response distribution 
#'  before and after the adjustment by this function. Default value is TRUE.
#'@param balanceidx Whether or not to calculate the balance index of the 
#'  sample response distribution before the adjustment, which is the variance 
#'  of the distribution density values around their mean.
#'@return A list with several elements. One is a data.frame recording the 
#'  sample meta data of the new dataset generated by this function. The 
#'  distribution of the response variable will be more balanced. The other 
#'  element is the feature value matrix for the new samples generated. Its 
#'  format is the same as the one for the parameter \code{samplevar}. If the 
#'  parameter \code{balanceidx} is set as TRUE, also a list element containing 
#'  the balance index before distribution adjustment will be returned.
#'@export
resamplebin <- function(resrange, 
                        binwith = 1, 
                        minwith = NULL, 
                        sampleseed = 1, 
                        samplevar, 
                        upsampling = 'SMOTE', 
                        #upsampling = 'SMOTE' or 'boot'
                        plotting = TRUE, 
                        balanceidx = FALSE){
  
  resrangetmp <- resrange
  binwithtmp <- binwith
  sampleseedtmp <- sampleseed
  upsamplingtmp <- upsampling
  samplevartmp <- samplevar
  
  if(is.null(minwith)){
    minwith <- binwithtmp
  }
  
  resid <- resrange$sampleid
  resrange <- resrange$Response
  resmin <- floor(min(resrange)) - minwith/2
  resmax <- ceiling(max(resrange)) + minwith/2
  resseq <- seq(resmin, resmax, binwith)
  
  if(max(resseq) < resmax){
    resseq <- c(resseq, max(resseq) + binwith)
  }
  
  tst <- cut(resrange, resseq)
  
  if(sum(table(tst) <= 1) > 0){
    
    resamplebin(resrange = resrangetmp, binwith = binwithtmp + minwith, 
                minwith = minwith, 
                sampleseed = sampleseedtmp, samplevar = samplevartmp, 
                upsampling = upsamplingtmp, plotting = plotting, 
                balanceidx = balanceidx)
    
    
  }else{
    
    groups <- data.frame(binname = as.character(tst), sampleid = resid, Response = resrange, 
                         stringsAsFactors = FALSE)
    groups <- groups[order(groups$Response, groups$binname),]
    uniqbinids <- seq(1, length(unique(groups$binname)), 1)
    
    groupdict <- data.frame(binname = unique(groups$binname), binid = uniqbinids, stringsAsFactors = FALSE)
    groups <- merge(groups, groupdict, by = c('binname'))
    groups <- groups[order(groups$binid, groups$sampleid),]
    groups <- groups[c('binid', 'sampleid', 'Response')]
    row.names(groups) <- 1:nrow(groups)
    binsamplesize <- ceiling(nrow(groups)/length(uniqbinids))
    
    
    for(j in uniqbinids){
      
      sub <- subset(groups, binid == j)
      rownames(sub) <- 1:nrow(sub)
      
      subvar <- samplevar[sub$sampleid,]
      
      if((nrow(sub) < binsamplesize) & upsampling == 'SMOTE'){
        
        dists <- dist(subvar, method = 'euclidean', diag = FALSE)
        distmin <- min(dists)
        dists <- as.matrix(dists)
        
        pairindces <- which(dists == distmin, arr.ind = TRUE)[1,]
        sample1 <- row.names(dists)[pairindces[1]]
        sample2 <- row.names(dists)[pairindces[2]]
        sample1var <- subvar[sample1,]
        sample2var <- subvar[sample2,]
        
        diff1 <- sample1var - sample2var
        diff1 <- diff1/2
        diff2 <- -diff1
        
        diffs <- rbind(diff1, diff2)
        row.names(diffs) <- 1:2
        
        
        set.seed(sampleseed)
        SMOTEressampleididx <- sample(x = 1:nrow(sub), size = binsamplesize, replace = TRUE)
        SMOTEressampleid <- sub$sampleid[SMOTEressampleididx]
        
        set.seed(sampleseed)
        SMOTEdiffid <- sample(x = c(1, 2), size = binsamplesize, replace = TRUE)
        
        set.seed(sampleseed)
        zetas <- runif(n = binsamplesize, min = 0, max = 1)
        zetas[!duplicated(SMOTEressampleididx)] <- 0
        
        SMOTEsamplevar <- samplevar[SMOTEressampleid,]
        SMOTEdiff <- diffs[SMOTEdiffid,]
        SMOTEdiff <- SMOTEdiff * zetas
        
        SMOTEsamplevar <- SMOTEsamplevar + SMOTEdiff
        
        SMOTEsub <- sub[SMOTEressampleididx,]
        
        suffix <- rep('', nrow(SMOTEsub))
        suffix[grepl(pattern = '\\.', x = row.names(SMOTEsub))] <- 
          substring(row.names(SMOTEsub), 
                    regexpr('\\.', row.names(SMOTEsub)))[grepl(pattern = '\\.', 
                                                               x = row.names(SMOTEsub))]
        
        SMOTEsub$sampleid <- paste0(SMOTEsub$sampleid, suffix)
        row.names(SMOTEsub) <- 1:nrow(SMOTEsub)
        row.names(SMOTEsamplevar) <- SMOTEsub$sampleid
        
        bootsamplevar <- SMOTEsamplevar
        rm(SMOTEsamplevar)
        bootsub <- SMOTEsub
        rm(SMOTEsub)
        
      }else{
        
        set.seed(sampleseed)
        bootressampleididx <- sample(x = 1:nrow(sub), size = binsamplesize, replace = TRUE)
        bootressampleid <- sub$sampleid[bootressampleididx]
        
        bootsamplevar <- samplevar[bootressampleid,]
        bootsub <- sub[bootressampleididx,]
        
        suffix <- rep('', nrow(bootsub))
        suffix[grepl(pattern = '\\.', x = row.names(bootsub))] <- 
          substring(row.names(bootsub), 
                    regexpr('\\.', row.names(bootsub)))[grepl(pattern = '\\.', 
                                                              x = row.names(bootsub))]
        
        bootsub$sampleid <- paste0(bootsub$sampleid, suffix)
        row.names(bootsub) <- 1:nrow(bootsub)
        row.names(bootsamplevar) <- bootsub$sampleid
        
        
      }
      
      if(j == 1){
        bootsamplevars <- bootsamplevar
        bootsubs <- bootsub
      }else{
        bootsamplevars <- rbind(bootsamplevars, bootsamplevar)
        bootsubs <- rbind(bootsubs, bootsub)
      }
      
    }
    
    names(bootsubs)[2] <- c('resampleid')
    final <- list(resampleresponse = bootsubs, 
                  resamplevars = bootsamplevars)
    
    #Balance index calculation#####
    if(balanceidx == TRUE){
      
      
      p <- ggplot2::ggplot() + 
        ggplot2::stat_density(data = groups, 
                              mapping = ggplot2::aes(x = Response, 
                                                     y = ..density..), 
                              geom = 'line', position = 'identity')
      
      q <- ggplot2::ggplot_build(p)$data[[1]]
      
      formerdis <- q$density
      
      rm(p)
      
      rm(q)
      
      formerdis <- as.vector(formerdis)
      
      formerdis <- formerdis[!is.na(formerdis)]
      
      if(length(formerdis) == 0){
        formersqr <- NULL
      }else{
        
        formersqr <- mean((formerdis - mean(formerdis))^2)
        
        formersqr <- -log10(formersqr)
        
      }
      
      
      
      final$balanceidx <- formersqr
      
    }
    
    
    if(plotting == TRUE){
      
      bootres <- bootsubs$Response
      realres <- groups$Response
      
      if(upsampling == 'SMOTE'){
        bootlabel <- 'SMOTE data'
      }else{
        bootlabel <- 'bootstrapped data'
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
      
    }
    
    return(final)
  }
  
}

#'Do bootstrapping on the samples of a dataset
#'
#'Use bootstrapping to sample the original samples and generate new 
#'  bootstrapped data.
#'
#'@param resrange A data.frame with meta information of the samples. Two 
#'  columns must be contained in it. One is the column "sampleid", the other 
#'  is the column "Response", which records the response variable values of 
#'  the samples. The response values should be numeric values.
#'@param sampleseed The random seed used in this function to do sampling. Its 
#'  default value is 1.
#'@param samplevar A matrix recording the values of features for all samples. 
#'  Each row represents one sample and each column represents one feature. The 
#'  column names are the feature names while the row names should be sample 
#'  IDs same as the ones recorded by the "sampleid" column in \code{resrange}. 
#'  Each column records a series of numeric values.
#'@param plotting Whether or not to plot the sample response distribution 
#'  before and after the bootstrapping by this function. Default value is 
#'  TRUE.
#'@return A list with 2 elements. One is a data.frame recording the sample 
#'  meta data of the bootstrapped dataset generated by this function. The 
#'  other element is the feature value matrix for the bootstrapped samples 
#'  generated. Its format is the same as the one for the parameter 
#'  \code{samplevar}.
#'@export
simpleboot <- function(resrange, 
                       sampleseed = 1, 
                       samplevar, 
                       plotting = TRUE){
  
  resrangetmp <- resrange
  sampleseedtmp <- sampleseed
  samplevartmp <- samplevar
  
  groups <- resrangetmp[c('sampleid', 'Response')]
  row.names(groups) <- 1:nrow(groups)
  samplesize <- nrow(groups)
  
  set.seed(sampleseed)
  bootressampleididx <- sample(x = 1:samplesize, size = samplesize, replace = TRUE)
  bootressampleid <- groups$sampleid[bootressampleididx]
  
  bootsamplevar <- samplevar[bootressampleid,]
  bootgroups <- groups[bootressampleididx,]
  
  suffix <- rep('', nrow(bootgroups))
  suffix[grepl(pattern = '\\.', x = row.names(bootgroups))] <- 
    substring(row.names(bootgroups), 
              regexpr('\\.', row.names(bootgroups)))[grepl(pattern = '\\.', 
                                                           x = row.names(bootgroups))]
  
  bootgroups$sampleid <- paste0(bootgroups$sampleid, suffix)
  row.names(bootgroups) <- 1:nrow(bootgroups)
  row.names(bootsamplevar) <- bootgroups$sampleid
  
  bootgroups <- bootgroups[order(bootgroups$Response),]
  bootsamplevar <- bootsamplevar[bootgroups$sampleid,]
  row.names(bootgroups) <- 1:nrow(bootgroups)
  
  
  names(bootgroups)[1] <- c('resampleid')
  final <- list(resampleresponse = bootgroups, 
                resamplevars = bootsamplevar)
  
  if(plotting == TRUE){
    
    bootres <- bootgroups$Response
    realres <- groups$Response
    
    bootlabel <- 'bootstrapped data'
    
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
    
  }
  
  return(final)
  
}
