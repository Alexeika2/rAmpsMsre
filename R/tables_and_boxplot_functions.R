
#' @export
gather_charact <- function(ftrs_mrgd_df, amps, msre_sites) {



  ftrs_mrgd_df$num_amp <- str_match(ftrs_mrgd_df$amp,"([a])(\\d+)")[,3]
  amp_names <- ftrs_mrgd_df$amp
  ftrs_mrgd_df$amp <- NULL

  amps$chr <- str_match(amps$cpg,"(chr[^:]+):(\\d+)-(\\d+)")[,2]
  amps$start <- str_match(amps$cpg,"(chr[^:]+):(\\d+)-(\\d+)")[,3]
  amps$end <- str_match(amps$cpg,"(chr[^:]+):(\\d+)-(\\d+)")[,4]
  amps$start <- as.numeric(amps$start)
  amps$end <- as.numeric(amps$end)



  amps_first <- amps %>% group_by(amp_name) %>% arrange(start) %>% filter(row_number()==1)
  amps_last <- amps %>% group_by(amp_name) %>% arrange(start) %>% filter(row_number()==n())
  amps_merged <- data.frame('amp_name'=amps_first$amp_name,'cutted'=amps_first$cutted,'chr' = amps_first$chr,
                                'start' = amps_first$start-1, 'end' = amps_last$end+1)
  amps_merged$amp_range <- paste0(amps_merged$chr,':',amps_merged$start,'-',amps_merged$end)
  amps_merged$amp_len <- abs(amps_merged$start - amps_merged$end)+1

  # Количество рестриктаз внутри ампликона
  overlaps <- ldply(
    unique(amps_merged$amp_name),
    function(amp_name) {
      # amp_name <- unique(amps_flt_merged$amp_name)[2120:2125]
      # amps <- amps_flt_merged
      gr = GRanges(seqnames=amps_merged$chr[amps_merged$amp_name==amp_name],
                   ranges=IRanges(start=amps_merged$start[amps_merged$amp_name==amp_name],
                                  end=amps_merged$end[amps_merged$amp_name==amp_name]))
      countOverlaps(gr, msre_sites)
    }, .parallel = TRUE
  )

  amps_merged$count_restr <- overlaps$V1

  # ftrs_mrgd_df <- ftrs_mrgd_all
  # amps_merged <- amps_flt_merged

  ftrs_mrgd_aux <- merge(ftrs_mrgd_df, amps_merged[,c('amp_name','cutted','amp_range','chr','start','end','amp_len','count_restr')], by.x='num_amp', by.y='cutted', sort = FALSE)

  # Разбираемся с паддингом

  padding_left <- ftrs_mrgd_aux[,c('amp_name','chr','start')]

  colnames(padding_left)[3] <- 'end'

  padding_left$start <- padding_left$end-100

  padding_right <- ftrs_mrgd_aux[,c('amp_name','chr','end')]

  colnames(padding_right)[3] <- 'start'

  padding_right$end <- padding_right$start+100

  overlaps_left <- llply(
    unique(padding_left$amp_name),
    function(amp_name) {
      # amp_name <- unique(amps_flt_merged$amp_name)[2120:2125]
      # amps <- amps_flt_merged
      gr <- GRanges(seqnames=padding_left$chr[padding_left$amp_name==amp_name],
                   ranges=IRanges(start=padding_left$start[padding_left$amp_name==amp_name],
                                  end=padding_left$end[padding_left$amp_name==amp_name]))
      hits <- findOverlaps(msre_sites, gr, type='within')
      p <- Pairs(msre_sites,gr,hits=hits)
      inter <- pintersect(p)
      distance(inter, GRanges(seqnames=padding_left$chr[padding_left$amp_name==amp_name], IRanges(start=padding_left$end[padding_left$amp_name==amp_name],
                                                                                                  end=padding_left$end[padding_left$amp_name==amp_name]+1)))
    }, .parallel = TRUE
  )


   overlaps_right <- llply(
    unique(padding_right$amp_name),
    function(amp_name) {
      # amp_name <- unique(amps_flt_merged$amp_name)[2120:2125]
      # amps <- amps_flt_merged
      gr <- GRanges(seqnames=padding_right$chr[padding_right$amp_name==amp_name],
                    ranges=IRanges(start=padding_right$start[padding_right$amp_name==amp_name],
                                   end=padding_right$end[padding_right$amp_name==amp_name]))
      hits <- findOverlaps(msre_sites, gr, type='within')
      p <- Pairs(msre_sites,gr,hits=hits)
      inter <- pintersect(p)
      c(distance(inter, GRanges(seqnames=padding_right$chr[padding_right$amp_name==amp_name], IRanges(start=padding_right$start[padding_right$amp_name==amp_name],
                                                                                                  end=padding_right$start[padding_right$amp_name==amp_name]+1))))
    }, .parallel = TRUE
  )

   ftrs_mrgd_aux$distance_left <- data.frame(distance_left = unlist(lapply(overlaps_left, paste, collapse = ", ")))$distance_left
   ftrs_mrgd_aux$distance_right <- data.frame(distance_right = unlist(lapply(overlaps_right, paste, collapse = ", ")))$distance_right


   ftrs_mrgd_aux <- ftrs_mrgd_aux[,c(16,17,21,22,23,24,2:15)]
   ftrs_mrgd_aux


}

#' @export
plot_detailed_boxplot <- function(region, rrbs_m_msre, sites=msre_sites, pad = 100, ftrs_mrgd_tbl = ftrs_mrgd_all_test) {

  # region <- ftrs_mrgd_all_test$amp_range[1]

  region_aux <- str_match(region,"(chr[^:]+):(\\d+)-(\\d+)")
  gr <- GRanges(seqnames = region_aux[,2], IRanges(start=as.numeric(region_aux[,3]), end=as.numeric(region_aux[,4])))
  gr_pad <- GRanges(seqnames = region_aux[,2], IRanges(start=as.numeric(region_aux[,3])-pad, end=as.numeric(region_aux[,4])+pad))

  rrbs_aux <- rrbs_m_msre[queryHits(findOverlaps(rrbs_m_msre@coords, gr_pad, type='within')),]

  cpg_regions <- rownames(rrbs_m_msre[queryHits(findOverlaps(rrbs_m_msre@coords, gr, type='within')),])

  sites <- mcols(msre_sites[subjectHits(findOverlaps(rrbs_aux@coords,msre_sites, type='within'))])
  coords_sites <- rrbs_aux@coords
  coords_sites$restr <- sites
  df_sites <- as.data.frame(coords_sites)
  df_sites$rrbs_names <- paste0(df_sites$seqnames,':',df_sites$start,'-',df_sites$end,'/+')
  df_sites$new_names <- paste0(df_sites$seqnames,':',df_sites$start,':',df_sites$restr)

  aux_matr <- t(rrbs_aux@.Data)

  amplicon <- rownames(ftrs_mrgd_tbl)[ftrs_mrgd_tbl$amp_range == region]

  aux_matr <- cbind(aux_matr, get_answer_nact(rrbs_m_msre))
  colnames(aux_matr)[ncol(aux_matr)] <- 'answer'
  melted <- melt(as.data.frame(aux_matr), id.vars=c('answer'))
  melted$value <- as.numeric(melted$value)
  melted <- melted[complete.cases(melted$answer),]
  melted <- merge(melted, df_sites[,c('rrbs_names','new_names')], by.x='variable',by.y='rrbs_names')
  theme_set(theme_classic())
  gp <- ggplot(melted, aes(x=new_names,y=value, fill=answer))+geom_boxplot()+
    theme(axis.text.x = element_text(size=12,angle = 60, hjust=1), axis.title.x = element_text(size=15), axis.text.y=element_text(size=12), axis.title.y = element_text(size=15),
          legend.position = 'top')+
    geom_vline(xintercept = which(levels(melted$variable)==cpg_regions[1])-0.5,col='blue',lwd=2)+
    geom_vline(xintercept = which(levels(melted$variable)==cpg_regions[length(cpg_regions)])+0.5, col='blue',lwd=2)+
    labs(x='Координаты CpG',y='Уровень метилирования b-value', title = paste0('Локус: ',amplicon,', Область ',region,' с отступами в ', pad))
  #title = paste0('Ампликон: ',amplicon,', Region ',region,' with padding of ', pad)
  return(gp)

}


plot_boxplot_amp <- function(region, rrbs_m_msre, sites_gr=msre_sites, marker, FUN=get_answer_nact, pad = 100) {

  # region <- table_two_sites_svr_ccgg$amp_range[1]
  # rrbs_m_msre <- rrbs_svr_ccgg
  # sites_gr <- CCGG_sites

  region_aux <- str_match(region,"(chr[^:]+):(\\d+)-(\\d+)")
  gr <- GRanges(seqnames = region_aux[,2], IRanges(start=as.numeric(region_aux[,3]), end=as.numeric(region_aux[,4])))
  gr_pad <- GRanges(seqnames = region_aux[,2], IRanges(start=as.numeric(region_aux[,3])-pad, end=as.numeric(region_aux[,4])+pad))




  rrbs_aux <- rrbs_m_msre[queryHits(findOverlaps(rrbs_m_msre@coords, gr_pad, type='within')),]

  # # checking on PC
  # if (grepl('PC', region_aux[,5])) {
  #   if (all(is.na(rrbs_aux@.Data) | length(rrbs_aux@.Data)==0)) {
  #     stop('Checking of PC is ok, no data is found\n')
  #   }
  #   else {stop('You should redesign or found another PC')}
  # }



  cpg_regions <- rownames(rrbs_m_msre[queryHits(findOverlaps(rrbs_m_msre@coords, gr, type='within')),])

  #sites <- mcols(sites_gr[subjectHits(findOverlaps(rrbs_aux@coords,sites_gr, type='within'))])
  coords_sites <- rrbs_aux@coords
  #coords_sites$restr <- sites

  df_sites <- as.data.frame(coords_sites)
  df_sites$rrbs_names <- paste0(df_sites$seqnames,':',df_sites$start,'-',df_sites$end,'/+')
  all_with_pads <- as.data.frame(sites_gr[queryHits(findOverlaps(sites_gr, gr_pad, type='within')),])
  all_with_pads$with_pads <- 1:nrow(all_with_pads)
  df_sites$with_pads <-subjectHits(findOverlaps(rrbs_aux@coords,sites_gr[queryHits(findOverlaps(sites_gr, gr_pad, type='within')),], type='within'))
  df_sites <- merge(all_with_pads, df_sites, by.x='with_pads', by.y='with_pads',all.x=TRUE)
  df_sites$RRBS <- ifelse(!is.na(df_sites$rrbs_names), TRUE, FALSE)
  df_sites$rrbs_names <- ifelse(df_sites$RRBS==TRUE,  paste0(df_sites$seqnames.y,':',df_sites$start.y,'-',df_sites$end.y,'/+'),  paste0(df_sites$seqnames.x,':',df_sites$start.x,'-',df_sites$end.x,'/+'))
  df_sites$new_names <- ifelse(df_sites$RRBS==TRUE,paste0(df_sites$seqnames.y,':',df_sites$start.y,':',df_sites$restr), paste0(df_sites$seqnames.x,':',df_sites$start.x,':',df_sites$restr))

  vec_pads <- df_sites$with_pads[df_sites$RRBS==FALSE]
  aux_matr <- rrbs_aux@.Data

  if (length(vec_pads) > 0) {
    matr_with_pads <- matrix(NA, nrow=nrow(df_sites),ncol=ncol(aux_matr), dimnames = list(df_sites$new_names, colnames(aux_matr)))
    matr_with_pads[-vec_pads,] <- aux_matr
  } else {matr_with_pads <- aux_matr; rownames(matr_with_pads) <- df_sites$new_names}

  # попытаемся прицепить к этой матрице даже сайты по которым нет данных в RRBS


  aux_matr <- t(matr_with_pads)

  #amplicon <- rownames(ftrs_mrgd_tbl)[ftrs_mrgd_tbl$amp_range == region]

  aux_matr <- cbind(aux_matr, FUN(rrbs_m_msre))
  colnames(aux_matr)[ncol(aux_matr)] <- 'answer'
  melted <- melt(as.data.frame(aux_matr), id.vars=c('answer'))
  melted$value <- as.numeric(melted$value)
  #melted <- melted[complete.cases(melted$answer),]
  melted <- merge(melted, df_sites[,c('rrbs_names','new_names')], by.x='variable',by.y='new_names')
  melted$rrbs_names <- factor(melted$rrbs_names)
  theme_set(theme_classic())
  gp <- ggplot(melted, aes(x=variable,y=value, fill=answer))+geom_boxplot()+
    theme(axis.text.x = element_text(size=12,angle = 60, hjust=1), axis.title.x = element_text(size=15), axis.text.y=element_text(size=12), axis.title.y = element_text(size=15))+
    geom_vline(xintercept = which(levels(melted$rrbs_names)==cpg_regions[1])-0.5,col='blue',lwd=2)+
    geom_vline(xintercept = which(levels(melted$rrbs_names)==cpg_regions[length(cpg_regions)])+0.5, col='blue',lwd=2)+
    labs(x='Coordinates',y='b-value', title = paste0(region,' with padding of ', pad,' Amplicon: ', marker))
  return(gp)

}
