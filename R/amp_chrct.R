
#' @export amp_chrct

amp_chrct <- function(t_mean_amp_imp_ftrs,method, angle = NULL, order_label, .parallel = FALSE) {

  # t_mean_amp_imp_ftrs <- t_mean_amp_ftrs_svr
  # method ='avg_perf_a'
  # angle = pi/4
  # order_label <- c('exposure_000','exposure_150')
  # .parallel = TRUE

  amp_list_rows <- as.list(unique(colnames(t_mean_amp_imp_ftrs)[1:ncol(t_mean_amp_imp_ftrs)-1]))
  to_keep_colnames <- colnames(t_mean_amp_imp_ftrs)[1:ncol(t_mean_amp_imp_ftrs)-1]
  colnames(t_mean_amp_imp_ftrs)[1:ncol(t_mean_amp_imp_ftrs)-1] <- str_match(to_keep_colnames,"([a])(\\d+)")[,1]
  amp_list_mean <- as.list(unique(colnames(t_mean_amp_imp_ftrs)[1:ncol(t_mean_amp_imp_ftrs)-1]))

  amp_set_roc <- {
    d <- plyr::ldply(
      amp_list_mean,
      function(amp){
        # amp=c('a4962')
        # method = 'avg_perf_a'
        # angle  = pi/4
        amp_glm <- glm_genes(t_mean_amp_imp_ftrs, amp,'answer', 333, 5, 100)
        r <- mean_rep_roc(amp_glm, order_label = order_label, method = method, angle = angle)
        c(auc=r$auc, sens = r$sens, spec=r$spec,accur=r$accur, shapiro_pval = r$shapiro_pval, fold_auc = r$fold_auc)
      }, .parallel = .parallel)
    rownames(d) <- sapply(amp_list_rows, paste)
    d
  }
  # amp_set_roc <- with(list(
  #   d=amp_set_roc[
  #     amp_set_roc$auc >= auc_thrs,]),
  #   d[order(d$auc,decreasing = T),])

  # amp_set_roc[complete.cases(amp_set_roc),]
  amp_set_roc
}


#' @export
amp_comb_chrct <-function(data, combs) {

  # calc_matr1 <- data_tnbc_to
  # combs <- combs_to[1:5]
  calc_matr1 <- data
  # colnames(calc_matr)[1:ncol(calc_matr)-1] <- str_match(colnames(calc_matr)[1:ncol(calc_matr)-1],"([a])(\\d+)")[,1]
  colnames(calc_matr1)[1:ncol(calc_matr1)-1] <- str_match(colnames(calc_matr1)[1:ncol(calc_matr1)-1],"([a]|[o])(\\d+)")[,1]
  comb_list <- combs
  combs_sub <- lapply(combs, function(com) { str_match(com,"([a]|[o])(\\d+)")[,1] })

  amp_comb_set_list <- {
    d <- llply(
      combs_sub,
      function(amp_set){
        # amp_set <- combs_sub[[1]]
        # amp_set <- pool_probe_list[[1]]
        tryCatch({
          amp_comb_model <- glm_genes(calc_matr1,amp_set,'answer', 333, 5, 100)
          #amp_comb_model_wo_imp <- glm_genes(calc_matr1,amp_set,'answer', 333, 5, 100)
          r <- mean_rep_roc(amp_comb_model, order_label = c('bad','good'),method='avg_perf_a', angle=pi/4)
          #r_wo_imp <- mean_rep_roc(amp_comb_model_wo_imp, order_label = c('bad','good'),method='avg_perf_a', angle=pi/4)
          fullness_panel_bad <- sum(rowSums(is.na(calc_matr1[calc_matr1$answer == 'bad', amp_set,drop=FALSE])) == 0)
          fullness_panel_good <- sum(rowSums(is.na(calc_matr1[calc_matr1$answer == 'good', amp_set,drop=FALSE])) == 0)
          list(
            amp_set = amp_set,
            metrics = c(auc=r$auc, sens = r$sens, spec=r$spec,accur=r$accur, shapiro_pval=r$shapiro_pval, panel_length=length(amp_set),
                        ppv = r$ppv, npv = r$npv, pval_roc = r$pval_roc, pval_roc_adj = r$pval_roc_adj,
                        fullness_panel_bad = fullness_panel_bad, fullness_panel_good = fullness_panel_good, median_auc = median(r$fold_auc)),
            fold_auc = r$fold_auc,
            error=NA)
        },
        error = function(err){
          list(amp_set = amp_set, error = err)
        })
      }, .parallel = TRUE)
    names(d) <- sapply(comb_list, paste , collapse = ', ')
    d
  }
amp_comb_set_list
}

#' @export
make_combs_df <- function(amps_comb_set_list, bad, good) {
  # amps_comb_set_list <- amp_panels_lumb_to

  amp_comb <- {
    # сделаем таблицу с метриками
      df = ldply(
        amps_comb_set_list,
        function(amp_comb){
        # amp_comb <- amp_comb_set_list_roc_all_new[[1]]

          amp_comb$metrics[1:13]
        })

      metrics_sorted <- df[order(df$auc,decreasing = T),]
      metrics_filter <- metrics_sorted[metrics_sorted$fullness_panel_bad >= bad & metrics_sorted$fullness_panel_good >= good, ]
      list(
        metrics=metrics_filter,
        #aucs = lapply(amps_comb_set_list, '[[', 'fold_auc')[order(df$auc, decreasing = T)])
        aucs = lapply(amps_comb_set_list, '[[', 'fold_auc')[metrics_filter$.id])
    }
    amp_comb
}


#' @export
make_combs_df_full <- function(amps_comb_set_list) {
  # amps_comb_set_list <- amp_panels_lumb_to

  amp_comb <- {
    # сделаем таблицу с метриками
    df = ldply(
      amps_comb_set_list,
      function(amp_comb){
        # amp_comb <- amp_comb_set_list_roc_all_new[[1]]

        amp_comb$metrics[1:13]
      })

    #metrics_sorted <- df[order(df$auc,decreasing = T),]
    #metrics_filter <- metrics_sorted[metrics_sorted$fullness_panel_bad >= 10 & metrics_sorted$fullness_panel_good >=7, ]
    list(
      metrics=df[order(df$auc,decreasing = T),],
      aucs = lapply(amps_comb_set_list, '[[', 'fold_auc')[order(df$auc, decreasing = T)])
      #aucs = lapply(amps_comb_set_list, '[[', 'fold_auc')[metrics_filter$.id])
  }
  amp_comb
}
