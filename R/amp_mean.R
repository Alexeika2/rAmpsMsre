require(GenomicRanges)
require(missRanger)
require(rtracklayer)
require(stringr)
require(plyr)
require(doParallel)
registerDoParallel()


#' @export
get_amps <- function(res_dist, rrbs_m_msre, cut_thrsh=100, .parallel = FALSE) {


  # res_dist <- res_dist_new
  # rrbs_m_msre <- rrbs_m_msre_new

  sprintf('Branching...\n')

  branch_by_chrom <- llply(res_dist,
                           function(dist_chr) {
                             cutree( hclust(dist_chr, method = "complete"), h = cut_thrsh)
                           }, .parallel = TRUE)

  last_num <- 0
  for (chr_i in 1:length(branch_by_chrom)){
    # chr_hc <- branch_by_chrom[[1]]
    branch_by_chrom[[chr_i]] <- branch_by_chrom[[chr_i]] + last_num
    last_num <- max(branch_by_chrom[[chr_i]])
  }
  cutted <- do.call(c, branch_by_chrom)
  names(cutted) <- gsub('.*\\.',"",names(cutted))

  sprintf('Branching is done!\n')

   # названия ампликонов ----

  sprintf('Getting gene names...\n')

  amp_names_list <- llply(
    unique(cutted),
    function(amp){
      # amp <- 1
      amp_coords_txt <- names(cutted[cutted==amp])
      m <- str_match(amp_coords_txt, "(chr[^:]+):(\\d+)-(\\d+)")
      amp_coords_gr <- GRanges(
        seqnames = m[,2], ranges = IRanges(start = as.numeric(m[,3]), end=as.numeric(m[,4])))
      gene <- paste(
        unique(unlist(rrbs_m_msre@genes[queryHits(findOverlaps(rrbs_m_msre@coords, amp_coords_gr))])),
        collapse = ',')

      sprintf("%s_a%d", gene, amp)
    }, .parallel = .parallel)

  # pools <- as.data.frame(unlist(pools))
  # colnames(pools)[1] <- "amp_name"
  # pools$n_amp <- rownames(pools)
  # pools <- merge(as.data.frame(cutted),pools, by.x='cutted',by.y='n_amp' )
  # pools$cpg<-names(cutted)

  sprintf('Getting gene names is done! \n')

  amps <- data.frame(
      as.data.frame(cutted),
      cpg = names(cutted),
      amp_name = unlist(amp_names_list)[cutted])

  #Убираем ампликоны где лишь одна CpG
  amps_flt <- subset(amps, amps$amp_name %in% names(which(table(amps$amp_name) >= 2)))
  amps_flt$amp_name <- as.character(amps_flt$amp_name)
  amps_flt$cpg <- as.character(amps_flt$cpg)

  sprintf('Amps is ready!\n')

  amps_flt
}


#' @export
get_mean_amp <- function(amps, rrbs_m_msre) {
#Подход с средними по амликону и по образцу ----
  mean_amp <- ldply(
    unique(amps$amp_name),
    function(amp_name){
    # amp_name <- unique(pools$pools)[1]
      cpg_sel <- amps$cpg[amps$amp_name==amp_name]
      apply(rrbs_m_msre@.Data[cpg_sel,], 2, mean, na.rm=T)
  })
#Заменяем - и , иначе в логит будет ошибка
  rownames(mean_amp) <- gsub(', ','_',unique(amps$amp_name))
  rownames(mean_amp) <- gsub('-','_',rownames(mean_amp))
  as.data.frame(t(mean_amp))
}


filter_miss_amp <- function(x, max_miss = 0.2) {

  # x=t_mean_amp_h110_a075
  # max_miss = 0.2
  miss_x <- sum(is.na(x))/(sum(!is.na(x))+sum(is.na(x)));

  while (
    miss_x > max_miss
  ) {
    amp_thrs <- min(colMeans(!is.na(x)))
    # s_thrs <- min(rowMeans(!is.na(x)))
    m_row <- x[!(rownames(x) %in% names(head(sort(pmin(rowSums(!is.na(x)))),1))),]
    m_col <- x[,which(colMeans(!is.na(x)) > amp_thrs)]
    miss_row <-sum(is.na(m_row))/(sum(!is.na(m_row))+sum(is.na(m_row)))
    miss_col <-sum(is.na(m_col))/(sum(!is.na(m_col))+sum(is.na(m_col)))
    if(miss_row < miss_col) {
      x <- m_row
    } else { x <- m_col}
    miss_x <- sum(is.na(x))/(sum(!is.na(x))+sum(is.na(x)));
  }
  return(x)
}

get_impt_res <- function(t_mean_amp_ftrs,k) {

#Сделал импутацию пропущенных значений ----
non_miss_mean_sample <- rowSums(!is.na(t_mean_amp_ftrs))
#non_miss_mean_amp <- colSums(!is.na(t_mean_amp_flt))

to_keep_colnames <- colnames(t_mean_amp_ftrs)
#Так как imputation затрудняется работать с некоторомы именами
#Cтобцов, то я заменил их просто на имя ампликона и затем переведу обратно
colnames(t_mean_amp_ftrs) <- str_match(to_keep_colnames,"([a])(\\d+)")[,1]
t_mean_amp_imp_ftrs <- missRanger(t_mean_amp_ftrs,
                                         pmm.k = k, seed=333,
                                         case.weights = non_miss_mean_sample, verbose = 0)
colnames(t_mean_amp_imp_ftrs) <- to_keep_colnames
t_mean_amp_imp_ftrs
}

#' @export
get_ftrs_table <- function(t_mean_amp, answer, paired=FALSE) {

  ftrs_amp_all <- {
    d <- ldply(
      with(list(g_list = colnames(t_mean_amp)[1:ncol(t_mean_amp)-1]), g_list ),
      function(amp, t_mean_amp, paired, answer){
        # amp <- "GLT1D1_a1"
        # получить два вектора чисел уровня метилирования

        bad_v <- c(t_mean_amp[t_mean_amp$answer == answer[1],colnames(t_mean_amp) == amp])
        good_v <- c(t_mean_amp[t_mean_amp$answer == answer[2],colnames(t_mean_amp) == amp])

        p_value = tryCatch(wilcox.test(bad_v, good_v, paired = paired)$p.value, error= function(cond) { 1; })

        c(amp=amp, bad = mean(bad_v, na.rm=T), good=mean(good_v, na.rm=T), p_value=p_value,
          fullness_good = mean(!is.na(good_v)), fullness_bad = mean(!is.na(bad_v)),
          fullness = mean(!is.na(c(bad_v,good_v))),
          emptyness = mean(is.na(c(bad_v,good_v))),good_sample = sum(!is.na(good_v)), bad_sample = sum(!is.na(bad_v)))
      }, t_mean_amp, paired, answer,
      .parallel=TRUE)
    d$adj_p_value <- p.adjust(as.numeric(d$p_value, method='fdr'))
    d
  }
  ftrs_amp_all
}
