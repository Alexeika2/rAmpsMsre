require(GenomicRanges)
require(doParallel)
require(stringr)
require(plyr)
require(dplyr)
require(caret)
require(rlang)
require(cvAUC)
require(ggplot2)
require(googlesheets4)

require(rrbsData)

# Парралеризм (если нужен)
registerDoParallel()

source('mean_rep_roc.R') # Для вычисления диагностических показателей (cvAUC, PPV, NPV, ...) и построение моделей различными методами
source('amp_mean.R') # Формирование матриц ампликонов с средними значениями уровня метилирования, формирование таблиц с ампликонами
source('amp_class_function.R') # Функции для рисования боксплотов cvAUC панелей
source('amps_cutted_new.R') # Основная функция для подбора ампликонов
source('cvauc_helper.R') # вспомогательные функции для рисования и вычисления cvAUC-ROC
source('plot_cv_roc.R') # Рисование cvAUC-ROC кривых
source('form_comb.R') # Формирование комбинаций панелей
source('amp_chrct.R') # Оценка диагностических ценностей сформированных панелей
source('tables_and_boxplots_functions.R') # Рисование боксплотов по одному ампликону с разметкой границ и формирование подробной диагностической таблицы для индивидуальных ампликонов 
source('dist_list3.R') # Вспомогательная ф-ция для amps_cutted_new

# Разметка для рестриктаз HpaII и HhaI

CCGG_sites <- load_g_sites(BSgenome.Hsapiens.UCSC.hg19,'CCGG', .parallel=TRUE)
GCGC_sites <- load_g_sites(BSgenome.Hsapiens.UCSC.hg19,'GCGC', .parallel=TRUE)
CCGG_sites$restr <- 'HpaII'
GCGC_sites$restr <- 'HhaI'

CGCG_sites <- load_g_sites(BSgenome.Hsapiens.UCSC.hg19,'CGCG', .parallel=TRUE)

# Объект GRanges c координатами для двух рестриктаз
msre_sites <-  sort(c(CCGG_sites, GCGC_sites))

rm(CCGG_sites, GCGC_sites)

# Подразумевается, что уже есть объект типа rrbsMatrix, 
# формирование rrbsMatrix объекта можно посмотреть на примере в пакете rrbsData

# rrbs_m_msre - объект класса rrbsMatrix
rrbs_m_msre <- rrbs_m_tn[queryHits(findOverlaps(rrbs_m_msre@coords, msre_sites, type='within')), answer %in% c('bad','good')]

answer <- get_answer_nact(rrbs_m_msre)


# Подбор ампликонов 

res_dist_new <- amps_cutted(rrbs_m_msre, .parallel = TRUE)

# Сырая таблица с CpG парами и принадлжености к сформированным ампликонам
amps_flt <- get_amps(res_dist_new, rrbs_m_msre_new, .parallel = TRUE)

# Матрица со средними значениями уровня метилирования b-value для ампликонов
t_mean_amp <- get_mean_amp(amps_flt, rrbs_m_msre)

# Достаем ответ на лечение из rrbsMatrix и перекодируем его в bad/good
answer <- get_answer_nact(rrbs_m_msre)

# Приклеиваем ответ на лечение к матрице ампликонов
t_mean_amp$answer <- answer[names(answer) %in% rownames(t_mean_amp)]

# Сделаем таблицу с базовыми характеристиками ампликонов (средний уровень метилирования в группах, количество доступной информации для групп, p-value и adjusted p-value)
ftrs_table <- get_ftrs_table(t_mean_amp, answer)


# Отфильтруем по количеству образцов в обеих группах и зафильтруем матрицу ампликонов
ftrs_amps <- ftrs_table$amp[as.numeric(ftrs_table$good_sample) >= 10 & as.numeric(ftrs_table$bad_sample) >= 10]
t_mean_amp_ftrs <- t_mean_amp[,c(ftrs_amps, 'answer')]

# Диагностические характеристики

set.seed(333)
amp_roc_all <- amp_chrct(t_mean_amp_ftrs, method = 'avg_perf_a', angle = pi/4, .parallel = TRUE)
amp_roc_all <- amp_roc_all[,1:4]
amp_roc_all$id <- rownames(amp_roc_all)
ftrs_mrgd_all <- merge(ftrs_table, amp_roc_all, by.x='amp',by.y='id')
amp_roc_all$id <- NULL

# Сформируем подробную таблицу для всех ампликонов
ftrs_mrgd_all_two_sites <- gather_charact(ftrs_mrgd_all, amps_flt, msre_sites)

# Нарисуем районы индивидуальных ампликонов с cvAUC >= 0.75

# имена маркеров с cvAUC >= 0.75

top_markers_new <- ftrs_mrgd_all_two_sites$amp_name[ftrs_mrgd_all_two_sites$auc >= 0.75]

# Список с боксплотами индивидуальных ампликонов
list_plots_amps <- llply(top_markers_new, function(mark, ftrs_mrgd_all_two_sites, rrbs_m_msre) {
  region <- ftrs_mrgd_all_two_sites$amp_range[rownames(ftrs_mrgd_all_two_sites)==mark]
  plot_detailed_boxplot(region, rrbs_m_msre, ftrs_mrgd_tbl = ftrs_mrgd_all_two_sites, sites=msre_sites)
  }, ftrs_mrgd_all_two_sites, rrbs_m_msre, msre_sites)

# Нарисуем все ампликоны в отдельный файл
cairo_pdf(filename='plots_ind_amps.pdf', onefile = TRUE)
lapply(list_plots_amps, function(p) {print(p) } )
dev.off()

# Сформируем панели ампликонов с 2 по количество столбцов в матрице за исключением колонки с ответом на НАХТ
combs <- form_comb(as.matrix(t_mean_amp_ftrs[,-ncol(t_mean_amp_ftrs)]), k=2:ncol(t_mean_amp_imp_comb[,-ncol(t_mean_amp_ftrs)]))

set.seed(333)
# Формирование списка с панелями ампликонов и их характеристиками
amp_panels <- amp_comb_chrct(t_mean_amp_ftrs, combs)

# Формирование таблицы со всеми ампликонами и их характеристиками
panels <- make_combs_df_full(amp_panels)

# Можно так же отфильтровать таблицу на количество образцов в каждой панеле
panels_ftrs <- make_combs_df(panels, bad=10, good=10)

# Нарисуем боксплоты cvAUC панелей, to указывает какое количество панелей рисовать, в примере первые 25 штук
get_fold_boxplot(panels, to=25, main='Panels')
get_fold_boxplot(panels_ftrs, to=25, main='Panels, filtered')
