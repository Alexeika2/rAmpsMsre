require(stringr)
require(cvAUC)
require(caret)
require(rlang)
require(dplyr)
source('/media/HEAP-EPI/akalinkin-seq/TNBC_Sigin_2020/cvauc_helper.R')

get_rocr_list <- function(train_obj,order_label){
  # order_label <- c('SD', 'PR') # control, case
  # train_obj <- glm_two
  
  # получить список список размером с выборку = length(labels)
  fr <- str_match(train_obj$pred$Resample, "^Fold(\\d+)\\.Rep(\\d+)")
  
  rocr_list <- list(
    predictions=split(train_obj$pred[[order_label[2]]], fr[,3]), 
    labels=split(train_obj$pred$obs, fr[,3]))
  rocr_list
}

get_cv_auc <- function(train_obj, order_label){
  # order_label <- c('SD', 'PR') # control, case
  # train_obj <- glm_two

    # получить список список размером с выборку = length(labels)
    fr <- str_match(train_obj$pred$Resample, "^Fold(\\d+)\\.Rep(\\d+)")
    
    rocr_list <- list(
      predictions=split(train_obj$pred[[order_label[2]]], fr[,3]), 
      labels=split(train_obj$pred$obs, fr[,3]))
    
    cv_auc <- cvAUC(predictions = rocr_list$predictions, rocr_list$labels, label.ordering = order_label)
    cv_auc
}

# PPV NPV metrics 
get_ppv_npv <- function(train_obj) {
  # order_label <- c('bad', 'good') # control, case
  # train_obj <- gene_glm
  
  # получить список список размером с выборку = length(labels)
  fr <- str_match(train_obj$pred$Resample, "^Fold(\\d+)\\.Rep(\\d+)")
  
  prediction <- split(train_obj$pred, fr[,3])
  
  mean_ppv <- mean(unlist(lapply(prediction, function(x) {posPredValue(x$pred,x$obs)})), na.rm=TRUE)
  mean_npv <- mean(unlist(lapply(prediction, function(x) {negPredValue(x$pred,x$obs)})), na.rm=TRUE)
  
  ppv_npv <- list(mean_ppv = mean_ppv, mean_npv = mean_npv)
  ppv_npv
}


# ROC pval (taking mean) 
get_roc_pval <- function(train_obj, order_label) {
  # order_label <- c('bad', 'good') # control, case
  # train_obj <- gene_glm
  
  # получить список список размером с выборку = length(labels)
  fr <- str_match(train_obj$pred$Resample, "^Fold(\\d+)\\.Rep(\\d+)")
  
  prediction <- split(train_obj$pred, fr[,3])
  
  
  p_val <- lapply(prediction, function(x) {
    if (train_obj$method == 'rpart'){
      names(train_obj$trainingData)[2:length(names(train_obj$trainingData))] <- train_obj$coefnames;
      outcome <- train_obj$trainingData$.outcome
      names_old <- names(train_obj$trainingData)
      train_obj$trainingData[,train_obj$coefnames,drop=FALSE] %>% mutate_if('is.character','as.factor') %>% mutate_if('is.factor','as.numeric') -> data_aux
      train_obj$trainingData <- cbind(outcome, data_aux) %>% as.data.frame
      names(train_obj$trainingData) <- names_old 
      val <- predict(train_obj$finalModel, train_obj$trainingData[x$rowIndex,train_obj$coefnames, drop=F], type='vector');
    } else {
      names(train_obj$trainingData)[2:length(names(train_obj$trainingData))] <- train_obj$coefnames;
      outcome <- train_obj$trainingData$.outcome
      names_old <- names(train_obj$trainingData)
      train_obj$trainingData[,train_obj$coefnames,drop=FALSE] %>% mutate_if('is.character','as.factor') %>% mutate_if('is.factor','as.numeric') -> data_aux
      train_obj$trainingData <- cbind(outcome, data_aux) %>% as.data.frame
      names(train_obj$trainingData) <- names_old 
      val <- predict(train_obj$finalModel, train_obj$trainingData[x$rowIndex,train_obj$coefnames, drop=F], type='response'); 
    }
    coded_obs <- ifelse(x$obs==order_label[2],1,0); 
    verification::roc.area(coded_obs, val)$p.value
    })
  
  p_val_mean_roc <- mean(unlist(p_val), na.rm = TRUE)
  p_val_mean_roc_adj <- mean(p.adjust(unlist(p_val)),na.rm=TRUE)
  list(pval_mean_roc = p_val_mean_roc, pval_mean_roc_adj = p_val_mean_roc_adj)
}


mean_rep_roc <- function(train_obj, order_label,method, angle=NULL){
  #method - метод для усреднения
  
  if(method == 'avg_perf_a' & is.null(angle)){stop('angle should be set up!')}
  # решили использовать усреднение ROC из пакета cvAUC
  # with(list(p=train_obj$pred),{
  #   fr <- str_match(p$Resample, "^Fold(\\d+)\\.Rep(\\d+)")
  #   roc_set <- llply(
  #     unique(fr[,3]),
  #     function(frep){
  #       # frep <- '001'
  #       rep_roc <- with(
  #         list(s=p[fr[,3]==frep,]),
  #         roc(predictor=as.integer(factor(s$pred)), response=factor(s$obs), levels = levels(s$pred)))
  #       rep_roc
  #     })
  #   list(auc=mean(sapply(roc_set, '[[', 'auc')), 
  #        spec=mean(sapply(roc_set, function(roc){ coords(roc, 'best',ret = c("threshold","specificity","sensitivity", "accuracy"),best.method='youden')['specificity']}),na.rm=T),
  #        sens=mean(sapply(roc_set, function(roc){ coords(roc, 'best',ret = c("threshold","specificity","sensitivity", "accuracy"),best.method='youden')['sensitivity']}),na.rm=T),
  #        accur=mean(sapply(roc_set, function(roc){ coords(roc, 'best',ret = c("threshold","specificity","sensitivity", "accuracy"),best.method='youden')['accuracy']}),na.rm=T))
  # })
  cv_auc <- get_cv_auc(train_obj, order_label)
  ppv_npv <- get_ppv_npv(train_obj)
  roc_pval <- tryCatch(get_roc_pval(train_obj, order_label), error=function(cond) {list(pval_mean_roc=NA,pval_mean_roc_adj=NA)})
  #Почему то при одинаковых значениях она не хочет считать теcт шапиро - будет единичка значит
  pval_shapiro = tryCatch(shapiro.test(cv_auc$fold.AUC)$p.value, error= function(cond) { 1; })
  switch(method,
         avg_perf_v ={
           avg_auc <- flip_perf(avg_perf_vertical(cv_auc$perf))
         },
         avg_perf_h = {
           avg_auc <- flip_perf(avg_perf_horizontal(cv_auc$perf))
         },
         avg_perf_a = {
           avg_auc <- flip_perf(avg_perf_angle(cv_auc$perf, angle))
         },
         stop('You must enter average method!')
         )
  #avg_auc <- avg_perf(cv_auc$perf)
  
  best <- with(
    list(best_i =which.max(avg_auc@x.values[[1]] + avg_auc@y.values[[1]])),{
      list(spec=avg_auc@x.values[[1]][best_i],
           sens=avg_auc@y.values[[1]][best_i])
    })
  
  positive = sum(train_obj$trainingData$.outcome == order_label[2])
  negative = sum(train_obj$trainingData$.outcome == order_label[1])
  accur = (positive*best$sens + negative*best$spec)/(positive + negative)
  
  list(
    auc = cv_auc$cvAUC,
    spec = best$spec,
    sens = best$sens,
    accur = accur,
    ppv = ppv_npv$mean_ppv,
    npv = ppv_npv$mean_npv,
    pval_roc = roc_pval$pval_mean_roc,
    pval_roc_adj = roc_pval$pval_mean_roc_adj,
    shapiro_pval = pval_shapiro,
    fold_auc = cv_auc$fold.AUC )
}

glm_genes <- function(data, gene_set, ans_name, seed, number, repeats) {
  #ans_name - column with answer name. 
  # gene_set <- c('ADCY8','DPYS')
  
  cat(paste0('Panel: ',paste(gene_set,collapse = ', '),'\n'))
  set.seed(seed)
  genes_form <- paste0(gene_set,collapse ='+')
  form <- as.formula(paste0(ans_name,'~',genes_form))
  data_glm <- train(form, data = data,
                    trControl=trainControl(method='repeatedcv', number=number, repeats=repeats, savePredictions = T,
                                           classProbs = T, returnResamp = 'all',
                                           summaryFunction = twoClassSummary, allowParallel = FALSE), #,sampling = 'down'),
                    method = 'glm', family = 'binomial', metric = 'ROC', na.action = na.omit)
}

glmnet_genes <- function(data, gene_set, ans_name, seed, number, repeats) {
  #ans_name - column with answer name. 
  # gene_set <- c('ADCY8','DPYS')
  
  set.seed(seed)
  genes_form <- paste0(gene_set,collapse ='+')
  form <- as.formula(paste0(ans_name,'~',genes_form))
  data_glm <- train(form, data = data,
                    trControl=trainControl(method='repeatedcv', number=number, repeats=repeats, savePredictions = T,
                                           classProbs = T, returnResamp = 'all',
                                           summaryFunction = twoClassSummary,sampling = 'down', allowParallel=FALSE),
                    method = 'glmnet', family = 'binomial', metric = 'ROC', na.action = na.omit)
}

bayesglm_genes <- function(data, gene_set, ans_name, seed, number, repeats) {
  #ans_name - column with answer name. 
  # gene_set <- c('ADCY8','DPYS')
  
  set.seed(seed)
  genes_form <- paste0(gene_set,collapse ='+')
  form <- as.formula(paste0(ans_name,'~',genes_form))
  data_glm <- train(form, data = data,
                    trControl=trainControl(method='repeatedcv', number=number, repeats=repeats, savePredictions = T,
                                           classProbs = T, returnResamp = 'all',
                                           summaryFunction = twoClassSummary,sampling = 'down',
										   allowParallel=FALSE),
                    method = 'bayesglm',family='binomial', metric = 'ROC', na.action = na.omit)
}




# build decision tree classifier from set of amplicons
tree_genes <- function(data, gene_set, ans_name, seed, number, repeats) {
  #ans_name - column with answer name. 
  # gene_set <- c('ADCY8','DPYS')
  set.seed(seed)
  genes_form <- paste0(gene_set,collapse ='+')
  form <- as.formula(paste0(ans_name,'~',genes_form))
  data_tree <- train(form, data = data,
                    trControl=trainControl(method='repeatedcv', number=number, repeats=repeats, savePredictions = T,
                                           classProbs = T, returnResamp = 'all',
                                           summaryFunction = twoClassSummary,
										   allowParallel=FALSE),
                    method = "rpart", metric='ROC', na.action = na.omit)
}

svm_genes <- function(data, gene_set, ans_name, seed, number, repeats) {
  #ans_name - column with answer name. 
  # gene_set <- c('ADCY8','DPYS')
  set.seed(seed)
  genes_form <- paste0(gene_set,collapse ='+')
  form <- as.formula(paste0(ans_name,'~',genes_form))
  data_tree <- train(form, data = data,
                     trControl=trainControl(method='repeatedcv', number=number, repeats=repeats, savePredictions = T,
                                            classProbs = T, returnResamp = 'all',
                                            summaryFunction = twoClassSummary,
											allowParallel=FALSE),
                     method = "svmLinear", metric='ROC', tuneGrid = expand.grid(C = seq(0, 2, length = 6)))
}

