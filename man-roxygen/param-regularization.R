#' @param reg_wt  the multiplicative weight(s) of the regularization terms. 
#' must be non-negative. 
#' May be a scalar or vector, but will be recycled to the length of the \code{reg_coef_idx}.
#' @param reg_zero  the \sQuote{zero} of the regularization terms. This would
#' usually be zero if you want to shrink to zero, but in some cases you may 
#' with to shrink to 1 for example.
#' May be a scalar or vector, but will be recycled to the length of the \code{reg_coef_idx}.
#' If \code{NULL} is given, defaults to zeroes for beta terms,
#' and ones for gamma terms in a Henery model fit, all zeroes
#' for Harville model fits.
#' @param reg_power  the power of the regularization terms, 2 for ridge
#' regression, 1 for lasso. 
#' @param reg_coef_idx  the index of the coefficient which the corresponding
#' regularization term is applied to. For the Harville model, the indices only refer to
#' the \eqn{\beta} coefficient vector. For the Henery model, the indices refer to the
#' \eqn{beta} coefficient vector and \eqn{\gamma} coefficient vector
#' concatenated together.
#' @param reg_standardize  if true, the \code{reg_wt} are normalized, or
#' \sQuote{standardized} with respect to the standard deviation of the
#' corresponding columns of the design matrix. That is, the weight used
#' is the given weight divided by the standard deviation of the corresponding
#' independent variable. Only terms associated with the betas are so
#' normalized.
#'
#' @note
#' The fit functions return an object of type \code{\link[maxLik]{maxLik}}
#' even when regularization penalties are applied; the statistical inference
#' functions are not valid when regularization is used. The user is warned.
#'
#' @md
#' @details 
#' # Regularization
#'
#' The regularization term is of the form
#' \deqn{\sum_i w_i |\nu_{c_i} - z_i|^{p_i},}
#' where \eqn{w_i} are the \code{reg_wt} weights,
#' \eqn{z_i} are the \code{reg_zero} zeroes,
#' \eqn{p_i} are the \code{reg_power} powers,
#' and \eqn{c_i} are the \code{reg_coef_idx} coefficient indices.
#' Note that the coefficient indices can be repeated so that the
#' regularization term can include multiple contributions from one
#' coefficient element. This allows \sQuote{elasticnet} regularization.
#' The \eqn{\nu} here refer to the regression coefficients \eqn{\beta}
#' concatenated with the \eqn{\gamma} coefficients in the Henery model.
