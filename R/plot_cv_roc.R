require(cvAUC)

# source('mean_rep_roc-copy2.R')

plot_cv_roc <- function(glm_model, order_label,method, angle=NULL, path=NULL,...) {
  # glm_model <- glm_model_lb1h
  # order_label <- c('bad','good')
  # order_label - нужна чтобы получать правильную кривую
  # path='2LB1h.pdf'
  # method = 'avg_perf_a'
  # angle = pi/4
  #Cначала отрицательный класс, потом положительный
  charact_cv_roc <- mean_rep_roc(glm_model, order_label, method, angle)
  model_list <- get_rocr_list(glm_model,order_label)
  cv_auc <- cvAUC(model_list$predictions, model_list$labels, label.ordering = order_label)
  ci_auc <- ci.cvAUC(model_list$predictions, model_list$labels, label.ordering = order_label)
  
  pdf(path, width = 5, height = 5)
    
  par(mar=c(5.1, 5.1, 4.1, 4.1), cex=1.0, cex.main=0.9, cex.axis=1.0, cex.lab=1.0, lwd=1.5)
  plot(c(), xlim=c(1,0), ylim=c(0,1), xlab="1 - Specificity", ylab="Sensitivity", ...)
  ROCR::plot(flip_perf(cv_auc$perf), col="grey82", lty=3, add=T, lwd=0.8)
  abline(a=1, b=-1, col="gray")
  
  switch(method,
  avg_perf_v = {
    ROCR::plot(flip_perf(avg_perf_vertical(cv_auc$perf)), col="black", add=TRUE, lwd=2.5)
  },
  avg_perf_h = {
    ROCR::plot(flip_perf(avg_perf_horizontal(cv_auc$perf)), col="black", add=TRUE, lwd=2.5)
  },
  avg_perf_a = {
    ROCR::plot(flip_perf(avg_perf_angle(cv_auc$perf, angle)), col="black", add=TRUE, lwd=2.5)
  }
  )
  
  par(cex=1.0, cex.main=1.0, cex.axis=1.0)
  text(sprintf("    AUC: %.2f", charact_cv_roc$auc), x=0.1, y=0.1, adj=1)
  text(sprintf("95%% CI: %.2f - %.2f", ci_auc$ci[1],ci_auc$ci[2]), x=0.05, y=0.05, adj=1)
  
  points(x=charact_cv_roc$spec, y=charact_cv_roc$sens, pch=16)
  text(
    paste(
      sprintf("Sensitivity: %.2f", charact_cv_roc$sens),
      sprintf("Specificity: %.2f", charact_cv_roc$spec),
      sprintf("  Accuracy: %.2f", charact_cv_roc$accur),
      sep="\n"), x=0.05, y=0.4, adj=1)
  dev.off()
}
