require(cvAUC)
require(plyr)

#' @export
flip_perf <- function(perf){
  #Function to rotate x axis

  res <- duplicate(perf, shallow = FALSE)
  res@x.values <- llply(perf@x.values, function(x) {1 - x} )
  res
}

#' @export
get_rotation_matrix_2d <- function(angle){
  matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)), nrow=2, ncol=2 )
}

#' @export
rotate_perf <- function(perf, angle){
  # perf <- test_cvAUC$perf
  # angle <- pi/4

  #perf@x.values, perf@x.values

  res <- duplicate(perf, shallow = FALSE)

  roc_lines <- mapply(function(x, y) {
    cbind(x, y)
    },
    perf@x.values,
    perf@y.values,
    SIMPLIFY = FALSE)

  rot_m <- get_rotation_matrix_2d(angle)

  roc_lines_rotated <- llply(roc_lines, '%*%', rot_m)

  res@x.values <- lapply(roc_lines_rotated, function(m) { m[,1] } )
  res@y.values <- lapply(roc_lines_rotated, function(m) { m[,2] } )
  res
}

#' @export
avg_perf_vertical <- function(perf) {
  # get single ROC averaged from ROC repeats by vertical

  perf.avg <- perf
  perf <- perf
  x.values <- seq(min(unlist(perf@x.values)), max(unlist(perf@x.values)),
                  length = max(sapply(perf@x.values, length)))
  for (i in 1:length(perf@y.values)) {
    perf.avg@y.values[[i]] <- approxfun(perf@x.values[[i]],
                                        perf@y.values[[i]], ties = mean, rule = 2)(x.values)
  }
  perf.avg@y.values <- list(rowMeans(data.frame(perf.avg@y.values)))
  perf.avg@x.values <- list(x.values)
  perf.avg@alpha.values <- list()
  perf.avg
}

#' @export
avg_perf_horizontal <- function(perf) {
  # get single ROC averaged from ROC repeats by horizontal

  perf.avg <- perf
  perf <- perf

  y.values <- seq(min(unlist(perf@y.values)), max(unlist(perf@y.values)),
                  length = max(sapply(perf@y.values, length)))
  for (i in 1:length(perf@x.values)) {
    perf.avg@x.values[[i]] <- approxfun(perf@y.values[[i]],
                                        perf@x.values[[i]], ties = mean, rule = 2)(y.values)
  }
  perf.avg@x.values <- list(rowMeans(data.frame(perf.avg@x.values)))
  perf.avg@y.values <- list(y.values)

  perf.avg@alpha.values <- list()
  perf.avg
}

#' @export
avg_perf_angle <- function(perf, angle) {
  # Get averaged ROC from ROC repeats by angle direction
  # positive angle - clockwise
  # perf <- test_set$perf
  # angle <- pi/4
  # возвращает объект ROCR::performance с одним (усредненным) Rep

  perf.rot <- rotate_perf(perf, angle) # duplicate

  x.values <- seq(min(unlist(perf.rot@x.values)), max(unlist(perf.rot@x.values)),
                  length = max(sapply(perf.rot@x.values, length)))
  for (i in 1:length(perf.rot@y.values)) {
    perf.rot@y.values[[i]] <- approxfun(perf.rot@x.values[[i]],
                                        perf.rot@y.values[[i]], ties = mean, rule = 2)(x.values)
  }

  perf.rot@y.values <- list(rowMeans(data.frame(perf.rot@y.values)))
  perf.rot@x.values <- list(x.values)
  perf.rot@alpha.values <- list()


  perf.avg <- rotate_perf(perf.rot, -angle) # duplicate

  perf.avg
}


