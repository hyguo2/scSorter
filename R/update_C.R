#' Update Cluster
#'
#' Updates cluster assignments based on center estimates from \code{update_mu}
#'
#' @param dat A matrix of input data.
#' @param mu_mat Center estimates from \code{update_mu}
#' @param designmat An indicator variable matrix records specified marker genes of each cell type.
#'

update_C = function(dat, mu_mat, designmat){
  ncmk = nrow(designmat)
  ngenes <- nrow(dat)
  ncells <- ncol(dat)
  nclusters <- ncol(mu_mat)

  dat_dist_mat_mk <- dat_dist_mat_mk_cache <- matrix(0, ncells, nclusters)
  dat_dist_mat_ad <- matrix(0, ncells, nclusters)

  base_mu_vec = rep(0, ncmk)
  for(bv in 1:ncmk){
    base_mu_vec[bv] = mu_mat[bv, which(designmat[bv,] == 0)[1]]  # a tiny change here
  }

  ## I re-wrote this block to make it ~60 times faster.
  delta <- mu_mat[1:ncmk,] * designmat
  diff <- dat[1:ncmk,] - base_mu_vec
  for (j in 1 : nclusters) {
    mat1 <- (diff - delta[, j]) ^ 2
    mat2 <- diff ^ 2
    which.smaller <- (mat1 < mat2 - .Machine$double.eps ^ 0.5) * designmat[, j]
    # On the above, when you compare two float numbers, it is always a good idea to specify
    # whether you want to include or exclude the equal case
    dat_dist_mat_mk[, j] <- colSums(mat1 * which.smaller + mat2 * (!which.smaller))
    dat_dist_mat_mk_cache[, j] <- colSums(which.smaller)

    dat_dist_mat_ad[, j] <- colSums((dat[(ncmk+1):ngenes, ] - mu_mat[(ncmk+1):ngenes, j])^2)
  }

  # calculate the distance matrix
  dat_dist_mat = dat_dist_mat_mk + dat_dist_mat_ad

  # I also re-wrote the rest to make it quicker and more concise
  # typically tie is not a big problem as it happens rarely; but if you indeed
  # worry about it, my code gives an easy way to get around it.
  dat_dist_mat_rand <- dat_dist_mat + rnorm(length(dat_dist_mat)) * .Machine$double.eps ^ 0.5
  clus <- apply(dat_dist_mat_rand, 1, which.min)

  return(list(clus, dat_dist_mat_mk_cache))
}
