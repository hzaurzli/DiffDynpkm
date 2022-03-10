#' @title  DiffDynpkm
#
#' @description  Calculate the differential expression of genes using fpkm/rpkm/tpm/cpm
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
#' @param treat  need a dataframe,treatment genes expression
#' @param control  need a dataframe,control genes expression
#' @param method  'glm' or 'wilcox' calculate p value
#' @param shrink   NB model shrink
#'
#' @return dataframe
#' @importFrom fitdistrplus fitdist
#' @importFrom stats rnbinom wilcox.test
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom MASS glm.nb
#' @export cal_diff
#'
#' @examples
#' data('gene_fpkm')
#' treat = gene_fpkm[,1:3]
#' control = gene_fpkm[,4:6]
#' cal_diff(treat,control)

cal_diff = function(treat,control,method = c('glm','wilcox')[1],shrink = FALSE){

  p_result = data.frame()

  if (shrink == TRUE) {
    library(fitdistrplus)
    n_treat = ncol(treat)
    n_control = ncol(control)
    counts = cbind(treat,control)
    logFC = data.frame()
    message(paste('Calculating NB model!\n'))
    pb <- txtProgressBar(style=3)
    for (i in 1:nrow(counts)) {
      setTxtProgressBar(pb, i/nrow(counts))
      treat_expression <- c(as.numeric(counts[i,colnames(treat)]))
      control_expression <- c(as.numeric(counts[i,colnames(control)]))

      if (rowSums(treat[i,]) == 0 | rowSums(control[i,]) == 0) {
        logFC[i,1] = log2(mean(as.numeric(treat[i,])) / mean(as.numeric(control[i,])))
      }
      else {
        fit_treat <- fitdist(data = treat_expression, distr = "nbinom", method = "mse")
        fit_control <- fitdist(data = control_expression, distr = "nbinom", method = "mse")
        set.seed(123)
        tst_treat = mean(rnbinom(1000, mu = summary(fit_treat)$estimate[2],size = summary(fit_treat)$estimate[1]))
        tst_control = mean(rnbinom(1000, mu = summary(fit_control)$estimate[2],size = summary(fit_control)$estimate[1]))
        logFC[i,1] = log2(tst_treat / tst_control)
      }
    }
  } else{
    n_treat = ncol(treat)
    n_control = ncol(control)
    logFC = log2(rowMeans(treat) / rowMeans(control))
  }

  message('\nCalculating p value!\n')
  pb <- txtProgressBar(style=3)
  if (method == 'glm') {
    library(MASS)
    design = c(rep('treat',ncol(treat)),rep('control',ncol(control)))
    for (i in 1:nrow(treat)) {
      setTxtProgressBar(pb, i/nrow(treat))
      expression = c(unlist(treat[i,]),unlist(control[i,]))
      result = tryCatch({
        p_val = summary(glm.nb(expression~design))$coefficients[2,4]
        p_val
      }, warning = function(w) {
        p_val = summary(glm.nb(expression~design))$coefficients[2,4]
        p_val
      }, error = function(e) {
        p_val = 1
        p_val
      })
      p_result[i,1] = result
    }
  }

  if (method == 'wilcox') {
    for (i in 1:nrow(treat)) {
      setTxtProgressBar(pb, i/nrow(treat))
      result = tryCatch({
        p_val = wilcox.test(unlist(treat[i,]),unlist(control[i,]),paired = T)
        p_val$p.value
      }, warning = function(w) {
        p_val$p.value
      }, error = function(e) {
        p_val$p.value = 1
        p_val$p.value
      })
      p_result[i,1] = result
    }
  }
  p_adj = p.adjust(p_result[,1],method = 'bonferroni')
  final = data.frame(cbind(treat,control,logFC,p_result,p_adj))
  colnames(final) = c(colnames(treat),colnames(control),'logFC','p_val','p,adj')
  return(final)
}

