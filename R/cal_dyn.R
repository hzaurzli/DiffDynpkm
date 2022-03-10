#' @title  DiffDynpkm
#
#' @description  Predic genes dynamic expression in different conditions
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#' @param design  need a dataframe,expriment matrix
#' @param exp  need a dataframe,genes expression matrix
#' @param gene  target gene
#' @param k   the dimension of the basis used to represent the smooth term
#' @param bs   a two letter character string indicating the (penalized) smoothing basis to use. (eg "tp" for thin plate regression spline, "cr" for cubic regression spline)
#'
#' @return list
#' @importFrom dplyr group_by summarise %>%
#' @importFrom stats as.formula
#' @importFrom mgcv gam
#' @export cal_dyn
#'
#' @examples
#' data('exp')
#' data('design)
#' cal_dyn(exp = exp,design = design,k=c(3,3),gene = 'gene_1',bs = 'cr')

cal_dyn = function(design,exp, gene = 'gene_1',k = c(3,3), bs = 'cr'){

  if (!is.data.frame(design)) {
    stop('design must be data.frame!')
  }

  if (!is.data.frame(exp)) {
    stop('Gene expression matrix must be data.frame!')
  }

  library(mgcv)
  library(dplyr)
  for (i in 3:ncol(design)) {
    if (!is.numeric(design[,i])) {
      design[,i] = design[,i] %>% as.factor() %>% as.numeric()
    }
  }

  expv = data.frame(t(exp[gene,]))
  all_data = cbind(expv, design)
  k = k
  variable = c()
  bs = bs
  for (i in 3:ncol(design)) {
    variable[i-2] = colnames(design)[i]
  }


  if(length(variable) == 1 & length(k) == 1 & length(variable) == length(k)){
    fit.formula = stats::as.formula(paste0("expv ~ s(",variable[1],", k = ", k[1], ",bs = '",bs,"')"))
  }else if (length(variable) == 2 & length(k) == 2 & length(variable) == length(k)) {
    condition = c(paste0("expv ~ s(",variable[1],", k = ", k[1], ",bs = '",bs,"')"),paste0("s(",variable[2],", k = ", k[2],",bs = '",bs,"')"))
    fit.formula = stats::as.formula(paste(condition[1],condition[2],sep = '+'))
  }else if (length(variable) > 3 & length(k) > 3 & length(variable) == length(k)) {
    condition = c(paste0("expv ~ s(",variable[1],", k = ", k[1], ",bs = '",bs,"')"),paste0("s(",variable[2],", k = ", k[2], ",bs = '",bs,"')"))
    for (i in 3:length(variable)) {
      condition[i] = paste0("s(",variable[i],", k = ", k[i], ",bs = '",bs,"')")
      if(i>length(variable)){
        break
      }
    }
    factor = paste(condition[1],condition[2],sep = '+')
    for (i in 3:length(condition)) {
      factor = paste(factor,condition[i],sep = '+')
      if(i>length(condition)){
        break
      }
    }
    fit.formula = stats::as.formula(factor)
  }


  dat = all_data %>% group_by(rep,tissue,time) %>% summarise(expv = mean(all_data[,1]))%>% data.frame()
  fit = gam(fit.formula ,data = dat)
  return(fit)
}
