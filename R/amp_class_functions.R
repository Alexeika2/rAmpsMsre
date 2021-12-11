require(factoextra)
require(ComplexHeatmap)
require(RColorBrewer)
require(reshape2)
require(ggpubr)

#Hierachical cluster tree ----
amp_clust <- function(ampList) {
#Иерархическая кластеризация ----
  # x <- amp_roc_class1
  stopifnot(class(ampList) == 'ampList')
  clust <- hclust(dist(t(ampList@after_impt[,ampList@amp_choose]),method='manhattan'),method='ward.D2')
  plot(clust, hang=-1)
}

#PCA ----
amp_pca <- function(ampList) {
  stopifnot(class(ampList) == 'ampList')
  pca <- prcomp(t(ampList@after_impt[,ampList@amp_choose]))
  fviz_pca_ind(pca, repel = T)
}

#Merged pval and roc in class slot ----
get_roc_pval <- function(ampList) {
  #Get ROC and pval together

  stopifnot(class(ampList) == 'ampList')

  #ampList = amp_roc_class_003_07

  pval_table <- ampList@genes_pval
  roc_table <- ampList@amp_roc
  pools_table <- ampList@pools

  pools_table$gene_name <- str_split_fixed(pools_table$amp_name,'_',2)[,1]
  fixed_amp_name <- gsub(', ','_',pools_table$amp_name)
  fixed_amp_name <- gsub('-','_',fixed_amp_name)

  pools_table$fixed <- fixed_amp_name

  pools_pval <- merge(pools_table, pval_table[,c('gene','adj_p_value')], by.x='gene_name', by.y = 'gene')

  roc_table$fixed <- rownames(roc_table)

  roc_pools <- merge(roc_table,pools_pval[,c('fixed','amp_name','adj_p_value')], by.x = 'fixed', by.y='fixed')

  roc_pools <- roc_pools[!duplicated(roc_pools$fixed),]

  rownames(roc_pools) <- roc_pools$amp_name

  roc_pools$amp_name <- NULL
  roc_pools$fixed <- NULL

  new(
    'ampList',
    genes_pval=ampList@genes_pval,
    pools = ampList@pools,
    bef_impt = ampList@bef_impt,
    after_impt = ampList@after_impt,
    amp_roc = roc_pools,
    amp_choose = ampList@amp_choose,
    description = ampList@description
  )

}

#Violin plots of choosen amplicons ----
get_violins <- function(ampList, matr_type = 'bef_impt') {

  stopifnot(class(ampList) == 'ampList')

  #По какой матрице рисовать: До импутации и со всеми образцами, до импутации после фильтрации и после фильтрации?

  bef_impt <- ampList@bef_impt
  bef_impt$answer <- ampList@after_impt$answer

  bef_impt <-bef_impt[,c(ampList@amp_choose, 'answer')]

  melted_bi <- melt(bef_impt, id.vars = 'answer')

  p <- ggviolin(melted_bi, x = 'answer', y = 'value',fill='answer', add = "boxplot",
                add.params = list(fill = "white", width=0.2), facet.by = 'variable', width=0.5)
  p + stat_compare_means(aes(group = answer), label = 'p.signif', method = 'wilcox',label.x.npc = 'centre', size = 10)

  switch(matr_type,
         bef_impt = {
           matr <- ampList@bef_impt
           matr$answer <- ampList@after_impt$answer

           matr <- matr[,c(ampList@amp_choose, 'answer')]

           melted_matr <- melt(matr, id.vars = 'answer')

           p <- ggviolin(melted_matr, x = 'answer', y = 'value',fill='answer', add = "boxplot",
                         add.params = list(fill = "white", width=0.2), facet.by = 'variable', width=0.5)
           p + stat_compare_means(aes(group = answer), label = 'p.signif', method = 'wilcox',label.x.npc = 'centre', size = 10)
         },
         after_impt = {
           matr <- ampList@after_impt

           matr <- matr[,c(ampList@amp_choose, 'answer')]

           melted_matr <- melt(matr, id.vars = 'answer')

           p <- ggviolin(melted_matr, x = 'answer', y = 'value',fill='answer', add = "boxplot",
                         add.params = list(fill = "white", width=0.2), facet.by = 'variable', width=0.5)
           p + stat_compare_means(aes(group = answer), label = 'p.signif', method = 'wilcox',label.x.npc = 'centre', size = 10)
         })

}

#Heatmap on choosen amplicons ----
get_heatmap <- function(ampList, matr_type = 'bef_impt') {

  stopifnot(class(ampList) == 'ampList')

  switch(matr_type,
         bef_impt = {
           matr <- ampList@bef_impt[,ampList@amp_choose]
           annotation <- data.frame('NACT' = ampList@after_impt$answer)
           rownames(annotation) <- rownames(matr)
           Heatmap(t(matr), cluster_rows = T, cluster_columns = T,
                   col =  colorRampPalette(rev(brewer.pal(n = 7, name= "RdYlGn")))(50),
                   top_annotation = HeatmapAnnotation(NACT = annotation$NACT,
                                                      col = list(NACT = c('good'='cyan','bad'='indianred'))),
                   clustering_method_rows = 'ward.D2', clustering_distance_rows = 'manhattan', clustering_method_columns = 'ward.D2',
                   clustering_distance_columns = 'manhattan',
                   heatmap_legend_param = list(title = 'Avg. B-value'))
         },
         after_impt = {
           matr <- ampList@after_impt[,ampList@amp_choose]
           annotation <- data.frame('NACT' = ampList@after_impt$answer)
           rownames(annotation) <- rownames(matr)
           Heatmap(t(matr), cluster_rows = T, cluster_columns = T,
                   col =  colorRampPalette(rev(brewer.pal(n = 7, name= "RdYlGn")))(50),
                   top_annotation = HeatmapAnnotation(NACT = annotation$NACT,
                                                      col = list(NACT = c('good'='cyan','bad'='indianred'))),
                   clustering_method_rows = 'ward.D2', clustering_distance_rows = 'manhattan',
                   clustering_method_columns = 'ward.D2',
                   clustering_distance_columns = 'manhattan',
                   heatmap_legend_param = list(title = 'Avg. B-value'))
         })
}

#' @export
get_fold_boxplot <- function(fold_auc, to, main,multiplex=NULL) {

  # fold_auc <- amp_panels
  fold_auc$metrics$.id <- gsub('_a\\d+',"",fold_auc$metrics$.id)
  #fold_auc$metrics$.id <- gsub('_',"",fold_auc$metrics$.id)
  names(fold_auc$aucs) <- gsub('_a\\d+',"",names(fold_auc$aucs))
  #names(fold_auc$aucs) <- gsub('_',"",names(fold_auc$aucs))

  test_dataframe <- ldply(fold_auc$aucs[1:to], data.frame)
  colnames(test_dataframe)[2] <- 'x'
  test_dataframe$.id <- factor(test_dataframe$.id, levels = fold_auc$metrics$.id[1:to])
  test_dataframe$Group <- sapply(stringr::str_split(test_dataframe$.id, ','),length)
  test_dataframe$Group <- as.factor(test_dataframe$Group)
  if(is.null(multiplex)) {
    ggplot(test_dataframe,aes(x = .id, y = x,fill=Group)) + stat_boxplot(geom='errorbar') + geom_boxplot() +
      theme(legend.title=element_text(size=15),legend.text = element_text(size=15),legend.position='top',axis.line = element_line(colour = "black"),
            panel.background = element_blank(),axis.text.x=element_text(angle=77, hjust=1),
            plot.margin = margin(0.5,0.5,0.5,0.5,'cm')) + scale_fill_brewer('Number of amplicons in panel',palette="RdYlGn") + ylab('cvAUC') +xlab(NULL)+ylim(0.4,1)+ggtitle(main)
  } else {
    ggplot(test_dataframe,aes(x = .id, y = x,fill=Group)) + stat_boxplot(geom='errorbar') + geom_boxplot() +
      theme(legend.title=element_text(size=15),legend.text = element_text(size=15),legend.position='top',axis.line = element_line(colour = "black"),
            panel.background = element_blank(),axis.text.x=element_text(angle=77, hjust=1),
            plot.margin = margin(0.5,0.5,0.5,0.5,'cm')) + scale_fill_brewer('Number of amplicons in panel',palette="RdYlGn") + ylab('cvAUC') +xlab(NULL)+ylim(0.4,1)+ggtitle(main)+
      scale_x_discrete(labels=multiplex)
}

}



get_fold_boxplot_noimp <- function(fold_auc) {

  # fold_auc <- amp_panels
  fold_auc$metrics <- fold_auc$metrics[order(fold_auc$metrics$auc_wo_imp, decreasing = TRUE),]
  fold_auc$metrics$.id <- gsub('_a\\d+',"",fold_auc$metrics$.id)
  #fold_auc$metrics$.id <- gsub('RUSC1_RUSC1_AS1, ',"RUSC1/RUSC1-AS1, ",fold_auc$metrics$.id)
  #fold_auc$metrics$.id <- gsub('_',"",fold_auc$metrics$.id)
  names(fold_auc$aucs_wo_imp) <- gsub('_a\\d+',"",names(fold_auc$aucs_wo_imp))
  #names(fold_auc$aucs) <- gsub('_',"",names(fold_auc$aucs))
  aucs_to_plot <- c(fold_auc$aucs_wo_imp[1:27], fold_auc$aucs_wo_imp[filter_panel_name(names(fold_auc$aucs_wo_imp), c('RUSC1_RUSC1_AS1','MXRA5','ANKRD46'))])
  test_dataframe <- ldply(aucs_to_plot, data.frame)
  colnames(test_dataframe)[2] <- 'x'
  test_dataframe$.id <- factor(test_dataframe$.id, levels = names(aucs_to_plot))
  test_dataframe$Group <- sapply(str_split(test_dataframe$.id, ','),length)
  test_dataframe$Group <- as.factor(test_dataframe$Group)








  ggplot(test_dataframe,aes(x = .id, y = x,fill=Group)) + stat_boxplot(geom='errorbar') + geom_boxplot() +
    theme(legend.title=element_text(size=15),legend.text = element_text(size=15),legend.position='top',axis.line = element_line(colour = "black"),
          panel.background = element_blank(),axis.text.x=element_text(angle=77, hjust=1),
          plot.margin = margin(0.5,0.5,0.5,0.5,'cm')) + scale_fill_brewer('Number of amplicons in panel',palette="RdYlGn") + ylab('cvAUC') +xlab(NULL)
}


filter_panel_name <- function(panel_names, markers) {

  # panel_names <- names(amp_panels$aucs_wo_imp)

  panel_names_split <- stringr::str_split(panel_names, ', ')
  unlist(llply(panel_names_split,
        function(panel_markers, markers){
          all(
            llply(panel_markers,
                  function(panel_marker, markers){
                    panel_marker2 <- gsub('_a\\d+','', panel_marker)
                    panel_marker2 %in% markers
                  }, markers))
        }, markers))
}

