#' Design Matrix Builder
#'
#' Builds the design matrix required by \code{update_func} based on user input.
#'
#' @param anno A matrix or data frame that contains marker genes specified for cell types of interest.
#' @param weight The default weight assigned to marker genes.
#'
#' @return A list contains processed design matrix and weight matrix.
#'

design_matrix_builder = function(anno, weight){
  if(weight <= 0)stop(paste('Default weight should be positive instead of the assigned value of ', weight, '.', sep = ''))

  colnames(anno) = toupper(colnames(anno))

  if(sum(c('TYPE', 'MARKER') %in% colnames(anno)) == 2) {
    ror = c(which(colnames(anno) == 'TYPE'), which(colnames(anno) == 'MARKER'))
    lidx = c(1:ncol(anno))[!(1:ncol(anno) %in% ror)]
    ror = c(ror, lidx)
    anno = anno[, ror]
  }

  celltypes = unique(anno[,1])
  marker = unique(anno[,2])

  nc = length(celltypes)
  nm = length(marker)

  designmat = matrix(0, nm, nc)

  nc_anno = ncol(anno)

  if (nc_anno == 3) {
    anno[, 3] = as.numeric(anno[, 3])
    poscheck = (anno[, 3] <= 0)
    if(sum(poscheck) > 0) {
      stop(paste('Please assign positive weights for the following marker genes: \n',
                 paste(sapply(1:nrow(anno), function(x)paste(anno[x, 1], anno[x, 2], sep = ': '))[poscheck], collapse = ', \n'),
                 '.', sep = ''))
    }
    weight = min(anno[, 3])
  }

  weightmat = matrix(weight, nm, nc)

  for (i in 1:nrow(anno)) {
    designmat[which(marker == anno[i, 2]), which(celltypes == anno[i, 1])] = 1
    if (nc_anno == 3) weightmat[which(marker == anno[i, 2]), which(celltypes == anno[i, 1])] = anno[i, 3]
  }

  designmat = as.data.frame(designmat, stringsAsFactors = F)
  colnames(designmat) = celltypes
  rownames(designmat) = marker

  colnames(weightmat) = celltypes
  rownames(weightmat) = marker

  return(list(designmat, weightmat))
}
