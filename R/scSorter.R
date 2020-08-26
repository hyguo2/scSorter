#' scSorter
#'
#' This is the main function that implements the scSorter method.
#'
#' @param expr A matrix of the input expression data. Each row represents a gene and each column represents a cell. Each row of this matrix should be named by the gene name it represents.
#' @param anno A matrix or data frame that contains marker genes specified for cell types of interest.
#' It should contain three columns named "Type", "Marker", and "Weight" that records the name and weight of marker genes specified for each cell type.
#' "Weight" column is optional. If it is not specified, the \code{default_weight} will be applied to all marker genes.
#' @param default_weight The default weight assigned to marker genes. The default value is 2.
#' @param n_start The number of possible cluster initializations. The default value is 10.
#' @param alpha The parameter determines the cutoff whether the cell type of a cell should be considered as undecided during unknown cell calling. The default value is 0.
#' @param u The parameter determines whether undecided cells are further processed. The default value is 0.05.
#' @param max_iter The maximum number of iterations for the algorithm to update parameters. The default value is 100.
#' @param setseed Random seed for cluster initialization. The default value is 0.
#'
#' @return A list contains the elements:
#'  \code{Pred_Type}: The predicted cell types.
#'  \code{Pred_param}: The parameter estimates of \code{mu} and \code{delta}.
#'
#' @export
#'
#' @examples
#' load(system.file('extdata', 'example_data.RData', package = 'scSorter'))
#' result = scSorter(expr, anno)
#' misclassification_rate = 1 - mean(result$Pred_Type == true_type)
#' table(result$Pred_Type, true_type)
#'

scSorter = function(expr, anno, default_weight = 2, n_start = 10, alpha = 0, u = 0.05, max_iter = 100, setseed = 0){
  #this is a wrapper function that implements the whole method based on the rest functions.
  #Rfast package is needed to run this method.
  anno_processed = design_matrix_builder(anno, default_weight)
  dt = data_preprocess(expr, anno_processed)

  dat = dt[[1]]
  designmat = dt[[2]]
  weightmat = dt[[3]]

  c_cost = NULL
  c_mu = list()
  c_clus = list()

  for(i in 1:n_start){
    set.seed(i+setseed)
    t1 = Sys.time()
    pred_ot = update_func(as.matrix(dat), designmat, weightmat, unknown_threshold1 = alpha, unknown_threshold2 = u, max_iter = max_iter)
    t2 = Sys.time()

    c_cost = c(c_cost, pred_ot[[3]])
    c_mu[[i]] = pred_ot[[1]]
    c_clus[[i]] = pred_ot[[2]]
  }

  pk = which.min(c_cost)

  pred_clus = c_clus[[pk]]
  pred_clus = c(colnames(designmat), rep('Unknown', ncol(designmat)))[pred_clus]
  pred_mu = c_mu[[pk]]

  return(list(Pred_Type = pred_clus, Pred_param = pred_mu))
}
