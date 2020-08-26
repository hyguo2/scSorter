#' Update Function
#'
#' Implements the scSorter method by iteratively running \code{update_mu} and \code{update_C}.
#'
#' @param dat A matrix of input data.
#' @param design_mat An indicator variable matrix records specified marker genes of each cell type.
#' @param weightmat A matrix of weights assigned to each marker gene.
#' @param unknown_threshold1 The parameter determines undecided cells cutoff. The default value is 0.
#' @param unknown_threshold2 The parameter determines whether undecided cells are further processed. The default value is 0.05.
#' @param max_iter The maximum number of iterations for the algorithm to update parameters. The default value is 100.
#'
#' @return A list contains parameter estimates, type assignments, and the corresponding cost.
#'

update_func = function(dat, design_mat, weightmat, unknown_threshold1 = 0, unknown_threshold2 = 0.05, max_iter = 100){
  cluster = sample(1:ncol(design_mat), ncol(dat), replace = T)

  rgwmat = range(weightmat)
  if (rgwmat[1] == rgwmat[2]) {

    ldt = nrow(dat)
    lmk = nrow(design_mat)
    marker_w = rep(weightmat[1, 1], lmk)
    rest_w = rep(min(lmk/(ldt-lmk),1), ldt-lmk)
    w = c(marker_w, rest_w)
    dat = dat*sqrt(w)

    for(i in 1:max_iter){
      mu = update_mu(dat, design_mat, cluster)
      cluster_old = cluster
      cluster_ot = update_C(dat, mu, design_mat)
      cluster = cluster_ot[[1]]
      if(sum(cluster != cluster_old)==0){
        break
      }
    }
  } else {
    ldt = nrow(dat)
    lmk = nrow(design_mat)
    rest_w = rep(min(lmk/(ldt-lmk),1), ldt-lmk)

    dat_mk = dat[1:lmk, ]
    dat_rt = dat[(lmk+1):ldt, ]

    dat_rt = dat_rt * sqrt(rest_w)

    weight_decorator = function(dat_mk, weightmat, clus) {
      uniclus = unique(clus)
      dat_mk2 = dat_mk

      for (i in uniclus) {
        marker_w = weightmat[, i]
        dat_mk2[, clus == i] = dat_mk[, clus == i]*sqrt(marker_w)
      }
      return(dat_mk2)
    }

    for(i in 1:max_iter){
      dat2 = rbind(weight_decorator(dat_mk, weightmat, cluster), dat_rt)
      mu = update_mu(dat2, design_mat, cluster)
      cluster_old = cluster
      cluster_ot = update_C(dat2, mu, design_mat)
      cluster = cluster_ot[[1]]
      if(sum(cluster != cluster_old)==0){
        break
      }
    }

    dat = dat2
  }

  cache_mat = cluster_ot[[2]]
  numofmarkergenes = colSums(design_mat)
  numofhighexpmkgenes = sapply(1:nrow(cache_mat), function(x)cache_mat[x,cluster[x]]/numofmarkergenes[cluster[x]])
  pks = numofhighexpmkgenes <= unknown_threshold1

  cluster_ukn_helper = cluster[pks]
  cluster_ukn = cluster[pks]

  cluster_kn = cluster[!pks]

  dat_kn_add = dat[(lmk+1):ldt, !pks]
  dat_ukn_add = dat[(lmk+1):ldt, pks, drop = F]

  #detect unknown cells from a given cluster.
  #Now we need to determine whether to put those potential unknown cells back to the corresponding cluster.
  #The total cost of two groups of cells (cells of that cluster and unknown cells) should reach a minimum when true unknown cells are picked out while the rest are put back into the cluster.
  #Under chi-square distribution, a cutoff is set to pick out true unknown cells. Since unknown cells are expected to be far away from the samples of that cluster. The search starts from one standard deviation.
  #50 cutoffs are selected and the one leads to the minimum total cost are used as final cutoff.
  unknown_detector = function(dat_kn, dat_ukn){
    pknz = rowSums(dat_kn != 0) != 0
    df = sum(pknz)

    dat_kn = dat_kn[pknz,]
    dat_ukn = dat_ukn[pknz,]

    rm = rowMeans(dat_kn)
    sdr = apply(dat_kn, 1, sd)

    dat_kn2 = (dat_kn-rm)/sdr
    dat_ukn2 = (dat_ukn-rm)/sdr

    dsdat = colSums(dat_kn2^2)
    dsdatuk = colSums(dat_ukn2^2)

    cost_calc2 = function(x){
      mx = rowMeans(x)
      return(sum((x-mx)^2))
    }

    calc_target2 = function(p){
      pk = pchisq(dsdatuk,df) <= p
      return(cost_calc2(cbind(dat_kn, dat_ukn[,pk]))+cost_calc2(dat_ukn[,!pk, drop = F]))
    }

    pl = c(seq(pchisq(df+sqrt(2*df),df),0.9,length.out = 10),seq(0.9, 0.95, length.out = 11)[2:11], seq(0.95, 0.99, length.out = 11)[2:11], seq(0.99, 1,length.out = 21)[2:21])
    ot = sapply(pl, calc_target2)

    pl_chosen = pl[which.min(ot)[1]]

    pkfinal = pchisq(dsdatuk, df) <= pl_chosen
    return(pkfinal)
  }

  uni_clus_ukn = unique(cluster_ukn_helper)

  nc = ncol(design_mat)

  for(ucukn in uni_clus_ukn){
    #this part could be further modified to account for bad quality clusters
    if(sum(cluster_kn==ucukn) <= max(1, round(unknown_threshold2*(sum(cluster_kn==ucukn)+sum(cluster_ukn_helper==ucukn))))){
      cluster_ukn[cluster_ukn_helper==ucukn] = nc+ucukn
    }else if(sum(cluster_ukn_helper==ucukn)==1){
      cluster_ukn[cluster_ukn_helper==ucukn] = nc+ucukn
    }else{
      uknrt = unknown_detector(dat_kn_add[,cluster_kn == ucukn], dat_ukn_add[,cluster_ukn_helper == ucukn, drop = F])
      cluster_ukn[cluster_ukn_helper==ucukn][!uknrt] = nc+ucukn
    }
  }

  cluster[pks] = cluster_ukn

  uniclus = unique(cluster)

  for(allclus in 1:max(uniclus)){
    if(allclus>nc){
      design_mat[, paste('Unknown_', allclus-nc, sep = '')] = 0
      if(allclus %in% uniclus){
        mu = cbind(mu, c(rep(0, lmk), rowMeans(dat_ukn_add[, cluster_ukn == allclus, drop = F])))
      }else{
        mu = cbind(mu, rep(0, nrow(mu)))
      }
    }else{
      if(allclus %in% uniclus){
        mu[(lmk+1):ldt,allclus] = rowMeans(dat[(lmk+1):ldt, cluster == allclus, drop = F])
      }else{
        mu[(lmk+1):ldt,allclus] = rep(0, ldt-lmk)
      }
    }
  }

  cost = cost_func(dat, cluster, mu, design_mat)
  return(list(mu, cluster, cost))
}
