require(plyr)

# Функция для формирования комбинаций из вектора выбранных маркеров ----
# С учетом матрицы дистанций по нормализованному Манхэттанскому расстоянию

form_comb <- function(matrix, k) {
  
  # matrix - матрица с образцами и ампликонами
  # k - сколько комбинаций
  # matrix <- as.matrix(t_mean_amp_imp_comb)
  # k <- 2:8
  
  #dist_matr <- dist_pairwise(as.matrix(matrix),method = 'manhattan')
  #as_matr_dist <- as.matrix(dist_matr)
  amp_vector <- colnames(matrix[,1:ncol(matrix)])
  combn_list <- sapply(k, combn, x=amp_vector)
  comnb_list_ftrs <- {
    d <- llply(
      combn_list,
      function(comb_matr){
        # comb_matr <- combn_list[[1]]
        t_comb_matr <- t(comb_matr)
        # t_comb_matr_idx <- sapply(seq(nrow(t_comb_matr)), 
        #                           function(r){ vec_amp <- t_comb_matr[r,]; short <- as_matr_dist[vec_amp,vec_amp]; 
        #                             short <- as.dist(short); 
        #                             if(sum(short < 0.2) >= 1){r}})
        t_comb_matr_ftrs <- unlist(t_comb_matr)
      }
      
    )
    d
  }
  comnb_list_ftrs <- comnb_list_ftrs[lengths(comnb_list_ftrs)>0]
  comnb_list_ftrs <- unlist(lapply(comnb_list_ftrs, function(m) {lapply(seq(nrow(m)),function(r) m[r,])}),recursive = F)
  comnb_list_ftrs
}
