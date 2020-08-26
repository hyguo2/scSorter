#' Mu Update
#'
#' Solves mu and delta given sample cluster assignment.
#'
#' @param dat A matrix of input data.
#' @param designmat An indicator variable matrix records marker genes of each pre-specified cell type.
#' @param clus A vector of cluster assignment.
#'
#' @return A matrix of parameter estimates.
#'

update_mu = function(dat, designmat, clus){
  nclus = ncol(designmat)
  nsample = ncol(dat)

  mu_mat = as.matrix(designmat)
  rownames(mu_mat) = NULL
  colnames(mu_mat) = NULL

  uc = unique(clus)

  #pre determine baseline level
  n_marker_gene = nrow(designmat)
  for(i0 in 1:nrow(mu_mat)){
    pck = clus %in% which(designmat[i0,]==0)
    if(sum(pck)==0){
      mu_mat[i0,]=0
    }else{
      mu_mat[i0,] = mean(dat[i0, clus %in% which(designmat[i0,]==0), drop = F])
    }
  }


  update_delta = function(idx, mu_mat, designmat){
    marker_clus = which(designmat[idx, ]==1)
    non_empty_marker_clus = uc[uc %in% marker_clus]
    if(length(non_empty_marker_clus)==0)return(mu_mat)

    base_mu = mu_mat[idx, which(designmat[idx,]==0)[1]]
    #pick out the samples satisfy delta < 2(x-mu) and calculate delta by an iterative approach
    for(ud in non_empty_marker_clus){
      obs = dat[idx, clus==ud]
      obs = obs-base_mu
      obs = obs[obs>0]
      if(length(obs)==0){
        mu_mat[idx, ud]=0
        next
      }

      for(i in 1:20){
        delta = mean(obs)
        picker = delta < 2*obs
        if(sum(picker)==length(picker)){
          mu_mat[idx, ud]=delta
          break
        }else{
          obs = obs[picker]
        }
      }
    }

    return(mu_mat)
  }

  #this function updates baseline nu (as the symbol used in paper).
  #All samples will either be used to estimate the delta or the baseline nu.
  #This function picks out samples that are not used in delta estimation (in update_delta function) and calculates baseline nu.
  update_basemu = function(idx, mu_mat, designmat){
    marker_clus = which(designmat[idx, ]==1)
    non_empty_marker_clus = uc[uc %in% marker_clus]
    base_mu = mu_mat[idx, which(designmat[idx,]==0)[1]]
    if(length(non_empty_marker_clus)==0){
      mu_mat[idx,designmat[idx,]==0]=mean(dat[idx,])
      return(mu_mat)
    }else{
      mu_l = mean(dat[idx,])
      for(ub in non_empty_marker_clus){
        J = dat[idx, clus == ub] > base_mu + .5*mu_mat[idx, ub]
        mu_l = mu_l - sum(J)*mu_mat[idx, ub]/nsample
      }
      mu_mat[idx,designmat[idx,]==0] = mu_l
    }
    return(mu_mat)
  }

  for(i in 1:n_marker_gene){
    for(z in 1:20){
      mu_old = mu_mat[i,]
      mu_mat = update_delta(i, mu_mat, designmat)
      mu_mat = update_basemu(i, mu_mat, designmat)
      if(mean(abs(mu_old-mu_mat[i,])) < 10^-6) break
    }
  }

  #this part estimates mu for other highly variable genes which follows kmeans approach.
  n_total_gene = nrow(dat)
  mu_mat_p2 = matrix(0, n_total_gene-n_marker_gene, nclus)

  for(k in unique(clus)){
    mu_mat_p2[,k] = rowMeans(dat[(n_marker_gene+1):n_total_gene, clus == k, drop = F])
  }

  return(rbind(mu_mat, mu_mat_p2))
}
