#' @name amps_cutted
#' @usage amps_cutted(rrbs_m_cor)
#' @description Split matrix of rrbs by chromosome, calculate distance between CpGs with correction factor of correlation, 
#' build hierarchical clustering with method single and cut tree by height in cut_thrsh parameter.
#' @param rrbs_m_cor RRBS object (matrix with b-value in rows and samples in columns)
#' @param cut_thrsh height at which we perform cutting (default 110)
#' @return named vector with numbers (amplicon numbers) and CpGs in names

amps_cutted <- function(rrbs_m_cor, .parallel=FALSE) {
  
  
  
  # rrbs_m_cor <- rrbs_m_20
  # rrbs_m_cor <- rrbs_m_msre_ccgg
  # разделим по хромосомам
  rrbs_m_by_chrom <- {
    tmp.rrbs_m_coords <- with(list(
      val= split(rrbs_m_cor@coords, seqnames(rrbs_m_cor@coords))),
      val[lengths(val) > 0]
    )
    tmp.rrbs_m_by_chrom <- lapply(
      tmp.rrbs_m_coords, 
      function(chr){
        rrbsObj = rrbs_m_cor[queryHits(findOverlaps(rrbs_m_cor@coords, chr)),]
      })
    tmp.rrbs_m_by_chrom[lapply(tmp.rrbs_m_by_chrom, nrow) >= 2]
  }
  
  
  branch_by_chrom <- llply(
    rrbs_m_by_chrom,
    # rrbs_m_by_chrom <- rrbs_m_by_chrom[['chr20']]
    function(chr_rrbs_m){
      # chr_rrbs_m <- rrbs_m_by_chrom[['chr17_gl000204_random']]
      cat(sprintf('Starting calculate %s \n', unique(seqnames(chr_rrbs_m@coords))))
      chr_coords_meth <- llply(
        1:length(chr_rrbs_m@coords),
        function(coord_i, answer){
          
          signif = tryCatch(
            wilcox.test(
              chr_rrbs_m@.Data[coord_i, answer == 'bad'],
              chr_rrbs_m@.Data[coord_i, answer == 'good']),
            error = function(cond) {list('p.value'=1)})
          signif <- if (is.na(signif$p.value)) { list('p.value'=1) } else { signif }

          list(
            coord = chr_rrbs_m@coords[coord_i],
            # bad = chr_rrbs_m@.Data[coord_i, names(answer[answer %in% 'bad'])],
            # good = chr_rrbs_m@.Data[coord_i, names(answer[answer %in% 'good'])]
            meth = chr_rrbs_m@.Data[coord_i,],
            signif = signif
          )
        }, answer = get_answer_nact(chr_rrbs_m))
      names(chr_coords_meth) <- rownames(chr_rrbs_m)
      
      # посчитаем составное расстояние между CpG-парами сайтов,
      # чтобы потом cut ампликоны
      
      chr_dists <- dist_list3(
        chr_coords_meth,
        coords_dist = 100,
        far_value = 1000,
        fun = function(cpg_site1, cpg_site2, meth, coords){
          
          # cpg_site1 <- chr_coords_meth$`chr20:57042460-57042461/+`
          # cpg_site2 <- chr_coords_meth$`chr20:57042484-57042485/+`
          
          
          coords_dist <- IRanges::distance(cpg_site1$coord,cpg_site2$coord)
          prev_cond <- coords_dist < 100 && 
            !is.na(cpg_site1$signif$p.value) && cpg_site1$signif$p.value < 0.2 &&
            !is.na(cpg_site2$signif$p.value) && cpg_site2$signif$p.value < 0.2
          
          if (prev_cond){
            
            manh_dist <- rUtils::dist_manhattan(cpg_site1$meth, cpg_site2$meth)
            if( !is.na(manh_dist) & manh_dist < 0.15) {
              coords_cpg <- GRanges(seqnames = seqnames(cpg_site1$coord), IRanges(start = end(cpg_site1$coord), end = start(cpg_site2$coord)))
              idx <- queryHits(findOverlaps(coords, coords_cpg, type='within'))
              between_cpg <- meth[idx]
              max_pvalue <-  max(sapply(between_cpg, function(cpg) {cpg$signif$p.value}))
              if (max_pvalue < 0.2) { 0 } else { 1000 }
            } else { 1000 }
          } else { 1000 }
          
          
          #print(paste0('coord 1: ', cpg_site1$coord, ' coord 2: ',cpg_site2$coord, ' dist genom: ',coords_dist, ' dist manh: ', manh_dist, ' pval site1: ', cpg_site1$signif$p.value,
          #            ' pval site 2: ',cpg_site2$signif$p.value))
          #ifelse(coords_dist < 100 & (manh_dist < 0.15 & !is.na(manh_dist)) & (cpg_site1$signif$p.value < 0.2 & !is.na(cpg_site1$signif$p.value)) & 
          #         (cpg_site2$signif$p.value < 0.2 & !is.na(cpg_site2$signif$p.value)), 0,
          #       1000)
        }, .parallel = .parallel, chr_coords_meth, chr_rrbs_m@coords)
      
      # метод single - чтобы считать групповое расстояние по ближайшим 
      #chr_branch <- cutree(hclust(chr_dists, method = "complete"), h=cut_thrsh)
      cat(sprintf('Calculation is done for %s \n', unique(seqnames(chr_rrbs_m@coords))))
      chr_dists
    })
  branch_by_chrom
  
  # branch_by_chrom <- list(chr1 = chr_branch, chr22 =  chr_branch )
  
  # Пофиксил немного цикл, в первоначальной версии были перескоки по индексам бранчей.
}
