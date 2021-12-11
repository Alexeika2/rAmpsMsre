#' @export

dist_list3 <- function(data, fun, coords_dist, far_value, .parallel = FALSE, ...) {
  if (!is.list(data)){
    stop("Param 'data' must be a list")
  }

  n <- length(data)

  cat('Not full ij_df preparing...\n')
  ij_df <- ldply(
    1:(n-1),
    function(i){
      # i = 1
      jn <- i
      while ( jn < n && IRanges::distance(data[[i]]$coord, data[[jn+1]]$coord) < coords_dist ){
        jn <- jn+1
      }

      if ( i < jn ){
        data.frame(
          i = rep(i, (jn-i)),
          j = (i+1):jn
        )
      } else{
        data.frame(i=c(), j=c())
      }
    }
  )
  cat('Not full ij_df done!\n')


  cat('ij_df chunk preparing\n')
  ij_chunk_list <- {
    worker_count <- foreach::getDoParWorkers()
    if (.parallel == TRUE & worker_count > 1) {
      split(ij_df, 1:nrow(ij_df) %% worker_count )
    } else {
      list(ij_df)
    }
  }
  cat('ij_df chunking done!\n')

  cat('calculation of dist df...\n')
  dists_df <- ldply(
    ij_chunk_list,
    function(ij_chunk){
      ddply(
        ij_chunk,
        c('i', 'j'),
        function(row) {
          o1 <- data[[row[1, 1]]]
          o2 <- data[[row[1, 2]]]

          # значения i, j появляются в первых столбцах результата data.frame сами
          c(dist=tryCatch(
            fun(o1, o2, ...),
            error = function(cond ) { NA }
          ))
        })
    }, .parallel = .parallel)
  cat('calculation of dist df done!', nrow(dists_df),'\n')

  # В случае пустого dataframe с dist (такое может быть)

  if (nrow(dists_df) == 0) {
    dists_df <- data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("i", "j", "dist"))))
  }


  cat('calculation of ij full df...\n')
  ij_full_df <- ldply(
    1:(n-1),
    function(i){
      data.frame(
        i = rep(i, n - i),
        j = (i + 1):n,
        dist_far = rep(far_value, n-i)
      )
    })
  # ij_full_df$dist_far <- rep(far_value, nrow(ij_full_df))
  cat('calculation of ij full df done!\n')


  cat('joining of dists...\n')
  dists_join <- dplyr::left_join(ij_full_df, dists_df, by = c('i','j'))
  dists <- ifelse(is.na(dists_join$dist), dists_join$dist_far+0.001, dists_join$dist) # vector
  cat('joining of dists is done!\n')

  l <- (n/2) * (n - 1)
  if (nrow(dists_join) != l) {
    stop(sprintf("Invalid operation l=(n/2)*(n-1)=%d, nrow(dists)=%d", l, nrow(dists_join)))
  }

  structure(as.numeric(dists), Size = n, Labels = names(data),
            Diag = FALSE, Upper = FALSE, method = "user",
            class = "dist")
}
