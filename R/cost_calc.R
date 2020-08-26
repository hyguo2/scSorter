#' Cost Function
#'
#' Calculates the cost.
#'
#' @param dat A matrix of input data.
#' @param clus A vector of predicted cell types.
#' @param mu Parameter estimates from \code{update_mu}.
#' @param designmat An indicator variable matrix records specified marker genes of each cell type.
#'

cost_func = function(dat, clus, mu, designmat){
  nmk = nrow(designmat)
  ncells = nrow(dat)
  uniclus = sort(unique(clus))
  cost = 0

  base_mu_vec = rep(0, nmk)
  for(bv in 1:nmk){
    base_mu_vec[bv] = mu[bv, which(designmat[bv,] == 0)[1]]
  }

  delta <- mu[1:nmk,]*designmat
  diff <- dat[1:nmk,] - base_mu_vec

  for(i in uniclus){
    mat1 <- (diff[, clus == i] - delta[, i]) ^ 2
    mat2 <- diff[, clus == i] ^ 2
    which.smaller <- (mat1 < mat2 - .Machine$double.eps ^ 0.5) * designmat[, i]
    cost = cost + sum(mat1 * which.smaller + mat2 * (!which.smaller))

    cost = cost + sum((dat[(nmk + 1):ncells, clus==i]-mu[(nmk + 1):ncells , i])^2)
  }

  return(cost)
}
