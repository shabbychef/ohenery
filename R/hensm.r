# /usr/bin/r
#
# Copyright 2018-2025 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav
#
# This file is part of ohenery.
#
# ohenery is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ohenery is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with ohenery.  If not, see <http://www.gnu.org/licenses/>.
#
# Created: 2018.12.20
# Copyright: Steven E. Pav, 2018-2024
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

# apparently necessary to register hensm as an S3 class
setOldClass('hensm')

.hensmlik <- function(theta, group, idx, X, wt, eta0, ...) {
  k <- ncol(X)
  beta <- theta[1:k]
  gamma <- theta[(k + 1):length(theta)]
  eta <- X %*% beta
  if (!is.null(eta0)) {
    eta <- eta + eta0
  }
  hensmlik(group, idx, eta, gamma = gamma, wt = wt) +
    .regularization_term(c(beta, gamma), ...)
}
.hensmgrad <- function(theta, group, idx, X, wt, eta0, ...) {
  k <- ncol(X)
  beta <- theta[1:k]
  gamma <- theta[(k + 1):length(theta)]
  eta <- X %*% beta
  if (!is.null(eta0)) {
    eta <- eta + eta0
  }
  hval <- hensmlik(group, idx, eta, gamma = gamma, wt = wt, deleta = X)
  grad <- c(attr(hval, 'gradient'), attr(hval, 'gradgamma'))
  grad + .regularization_grad(c(beta, gamma), ...)
}
#2FIX: why isn't there a experts version of this one?
#  @param ngamma  the number of gammas to model; we model
#        \eqn{\gamma_2} through \eqn{\gamma_n}.
.hmfit <- function(
  y,
  g,
  X,
  wt = NULL,
  eta0 = NULL,
  beta0 = NULL,
  gamma0 = NULL,
  normalize_wt = FALSE,
  ngamma = 4,
  reg_wt = NULL,
  reg_zero = NULL,
  reg_power = NULL,
  reg_coef_idx = NULL,
  reg_standardize = FALSE,
  method = c('BFGS', 'NR', 'CG', 'NM')
) {
  method <- match.arg(method)
  if (!is.null(gamma0)) {
    ngamma <- length(gamma0) + 1
  } else {
    gamma0 <- array(1, ngamma - 1)
  }
  stopifnot(ngamma >= 2)
  k <- ncol(X)
  if (is.null(beta0)) {
    beta0 <- array(0, k)
  }
  theta0 <- c(beta0, gamma0)
  reg_zero <- .regularization_default_zero(reg_zero, reg_coef_idx, num_beta = k)
  reg_wt <- .regularization_standardize(
    reg_wt,
    reg_power,
    reg_coef_idx,
    reg_standardize,
    X
  )
  .check_regularization(theta0, reg_wt, reg_zero, reg_power, reg_coef_idx)

  if (!is.null(wt) && normalize_wt) {
    wt <- wt / abs(mean(wt, na.rm = TRUE))
  } # by having the abs, negative weights still throw an error.
  # turn g into integers?
  if (is.integer(g)) {
    group <- g
  } else {
    group <- match(g, unique(g))
  }

  idx <- order(g, y, decreasing = TRUE) - 1
  rv <- maxLik(
    logLik = .hensmlik,
    grad = .hensmgrad,
    hess = NULL,
    start = theta0,
    method = method,
    group = group,
    idx = idx,
    X = X,
    wt = wt,
    eta0 = eta0,
    reg_wt = reg_wt,
    reg_zero = reg_zero,
    reg_power = reg_power,
    reg_coef_idx = reg_coef_idx
  )
  retv <- list(
    mle = rv,
    beta = rv$estimate[1:k],
    coefficients = rv,
    gammas = rv$estimate[(k + 1):length(theta0)],
    gamma2 = rv$estimate[k + 1],
    estimate = rv$estimate, # sigh
    wt = wt,
    g = g,
    y = y,
    formula = NULL,
    eta0 = eta0,
    reg_wt = reg_wt,
    reg_zero = reg_zero,
    reg_power = reg_power,
    reg_coef_idx = reg_coef_idx,
    reg_standardize = reg_standardize
  )

  gnames <- paste0('gamma', 2:ngamma)
  names(retv$mle$estimate) <- c(colnames(X), gnames)
  names(retv$beta) <- colnames(X)

  # do some summarization
  retv$deviance <- -2 * rv$maximum
  retv$deviance_df <- length(retv$estimate)
  retv$is_relative <- !is.null(eta0)
  class(retv) <- 'hensm'
  retv
}
#' @title Friendly interface to softmax regression under Henery model.
#'
#' @description
#'
#' A user friendly interface to the softmax regression under the Henery model.
#'
#' @details
#'
#' Performs a softmax regression by groups, via Maximum Likelihood Estimation.
#' It is assumed that successive sub-races maintain the proportional
#' probability of the softmax, up to some gamma coefficients,
#' \eqn{\gamma_2, \gamma_3, ..., \gamma_n}, which we fit. This model
#' nests the Harville model fit by \code{\link{harsm}}, by fixing all
#' the gammas equal to 1.
#'
#' @inheritParams stats::lm
#' @param na.action  How to deal with missing values in \code{y}, \code{g},
#' \code{X}, \code{wt}, \code{eta0}.
#' @param ngamma  The number of gammas to fit. Should be at least 2.
#' @param fit0    An optional object of class \code{hensm} or of \code{harsm}
#' with the initial fit estimates.
#' These will be used for \sQuote{warm start} of the estimation procedure.
#' A warm start should only speed up estimation, not change the ultimate results.
#' When there is mismatch between the coefficients in \code{fit0} and the model
#' being fit here, the missing coefficients are initialized as zero.
#' If \code{ngamma} is \code{NULL} and \code{fit0} is given,
#' we default to the number of gammas in the initial fit, otherwise
#' we fill any missing gammas with 1.
#' If a \code{harsm} object is given, then \code{ngamma} must be non-null.
#'
#' @template param-weights
#' @template param-regularization
#' @template param-group
#' @template etc
#' @template note-ties
#' @template note-weights
#' @return An object of class \code{hensm}, but also of \code{maxLik} with the
#' fit.
#' @keywords fitting
#' @seealso \code{\link{harsm}}, \code{\link{smlik}}.
#' @template note-normalization
#' @note When regularization is used, the first gamma coefficient is not
#' shrunk, as it always equals one in the Henery model, and is not estimated
#' from the data.
#'
#' @examples
#'
#' nfeat <- 5
#' set.seed(1234)
#' g <- ceiling(seq(0.1,1000,by=0.1))
#' X <- matrix(rnorm(length(g) * nfeat),ncol=nfeat)
#' beta <- rnorm(nfeat)
#' eta <- X %*% beta
#' # 2FIX: do rhenery
#' y <- rsm(eta,g)
#' # now the pretty frontend
#' data <- cbind(data.frame(outcome=y,race=g),as.data.frame(X))
#'
#' fmla <- outcome ~ V1 + V2 + V3 + V4 + V5
#' fitm <- hensm(fmla,data,group=race)
#'
#' # with offset
#' eta0 <- rowMeans(X)
#' data <- cbind(data.frame(outcome=y,race=g,eta0=eta0),as.data.frame(X))
#' fmla <- outcome ~ offset(eta0) + V1 + V2 + V3 + V4 + V5
#' fitm <- hensm(fmla,data,group=race)
#'
#' # on horse race data
#' library(dplyr)
#' data(race_data)
#' df <- race_data %>%
#' 	group_by(EventId) %>%
#' 		mutate(eta0=log(WN_pool / sum(WN_pool))) %>%
#' 	ungroup() %>%
#' 	mutate(weights=ifelse(!is.na(Finish),1,0)) %>%
#' 	mutate(fac_age=cut(Age,c(0,3,5,7,Inf),include.lowest=TRUE))
#'
#' # Henery Model with market efficiency
#' hensm(Finish ~ eta0,data=df,group=EventId,weights=weights,ngamma=3)
#'
#' # look for age effect not captured by consensus odds.
#' fmla <- Finish ~ offset(eta0) + fac_age
#' fit0 <- hensm(fmla,data=df,group=EventId,weights=weights,ngamma=2)
#' # allow warm start.
#' fit1 <- hensm(fmla,data=df,group=EventId,weights=weights,fit0=fit0,ngamma=2)
#' # allow warm start with more gammas.
#' fit2 <- hensm(fmla,data=df,group=EventId,weights=weights,fit0=fit0,ngamma=3)
#' # or a different formula
#' fit3 <- hensm(update(fmla,~ . + PostPosition),data=df,group=EventId,weights=weights,fit0=fit0)
#'
#' # warm start from harsm object
#' fit0_har <- harsm(fmla,data=df,group=EventId,weights=weights)
#' fit4 <- hensm(fmla,data=df,group=EventId,fit0=fit0_har,weights=weights)
#'
#' # regularization examples
#' fmla <- Finish ~ offset(eta0) + fac_age + PostPosition
#' # no regularization
#' fitr0 <- harsm(fmla,data=df,group=year,weights=weight,ngamma=3)
#' # ridge regression on the betas (there are two)
#' fitr2 <- harsm(fmla,data=df,group=year,weights=weight,ngamma=3,
#'   reg_wt=rep(1,2), reg_power=2, reg_zero=0, reg_coef_idx=c(1,2))
#' # l1 regression on the betas (there are two)
#' fitr1 <- harsm(fmla,data=df,group=year,weights=weight,ngamma=3,
#'   reg_wt=rep(1,2), reg_power=1, reg_zero=0, reg_coef_idx=c(1,2))
#' # elasticnet regression on the betas (there are two)
#' fitr12 <- harsm(fmla,data=df,group=year,weights=weight,ngamma=3,
#'   reg_wt=rep(1,4), reg_power=c(1,1,2,2) reg_zero=0, reg_coef_idx=c(1,2,1,2))
#' # l1 regression on the gammas, shrinking to one (there are two estimated)
#' fitrg1 <- harsm(fmla,data=df,group=year,weights=weight,ngamma=3,
#'   reg_wt=rep(1,2), reg_power=1, reg_zero=1, reg_coef_idx=c(3,4))
#' # l2 regression on the betas and gammas, shrinking beta to 0, gamma to 1.
#' fitrg2 <- harsm(fmla,data=df,group=year,weights=weight,ngamma=3,
#'   reg_wt=rep(1,4), reg_power=2, reg_zero=c(0,0,1,1), reg_coef_idx=1:4)
#'
#' @importFrom stats coef formula model.frame model.matrix na.omit
#' @export
#' @rdname hensm
hensm <- function(
  formula,
  data,
  group = NULL,
  weights = NULL,
  ngamma = 4,
  fit0 = NULL,
  reg_wt = NULL,
  reg_zero = NULL,
  reg_power = NULL,
  reg_coef_idx = NULL,
  reg_standardize = FALSE,
  na.action = na.omit
) {
  substitute(formula)

  # I find it highly offensive that this cannot be done reasonably
  # easily in a subfunction because of NSE whatever.

  # https://stackoverflow.com/q/53827563/164611
  mf <- match.call(expand.dots = FALSE)
  #turn weights into symbol if character is passed
  if (is.character(mf$weights)) mf$weights <- as.symbol(mf$weights)
  if (is.character(mf$group)) mf$group <- as.symbol(mf$group)
  m <- match(
    c("formula", "data", "weights", "group", "na.action"),
    names(mf),
    0L
  )
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame()) #evaluate call

  Xs <- model.matrix(formula, mf)
  # remove intercept!
  if (colnames(Xs)[1] == '(Intercept)') {
    Xs <- Xs[, -1, drop = FALSE]
  }
  y <- as.vector(model.response(mf))
  group <- as.vector(model.extract(mf, 'group'))
  eta0 <- model.offset(mf)
  wt <- as.vector(model.weights(mf))

  dat <- list(Xs = Xs, y = y, group = group, eta0 = eta0, wt = wt)

  if (!is.null(fit0)) {
    stopifnot(any(c("hensm", "harsm") %in% class(fit0)))
    feat_names <- colnames(dat$Xs)
    beta0 <- rep(0, ncol(dat$Xs))
    found <- feat_names %in% attr(fit0$beta, 'names')
    beta0[found] <- fit0$beta[feat_names[found]]
    # 2FIX: warn if some names are missing? not needed, I would think.
    # stopifnot(colnames(dat$Xs) == attr(fit0$beta,'names'))
    if (is.null(ngamma)) {
      stopifnot(any(c("hensm") %in% class(fit0)))
      gamma0 <- fit0$gammas
      ngamma <- length(gamma0) + 1
    } else {
      gamma0 <- rep(1, ngamma - 1)
      if (any(c("hensm") %in% class(fit0))) {
        dotake <- min(ngamma - 1, length(fit0$gammas))
        gamma0[1:dotake] <- fit0$gammas[1:dotake]
      }
    }
  } else {
    beta0 <- NULL
    gamma0 <- NULL
  }

  retv <- .hmfit(
    y = dat$y,
    g = dat$group,
    X = dat$Xs,
    wt = dat$wt,
    beta0 = beta0,
    gamma0 = gamma0,
    ngamma = ngamma,
    eta0 = dat$eta0,
    reg_wt = reg_wt,
    reg_zero = reg_zero,
    reg_power = reg_power,
    reg_coef_idx = reg_coef_idx,
    reg_standardize = reg_standardize
  )
  retv <- as.linodds(retv, formula, beta = retv$beta)
  retv
}
#' @export
#' @rdname hensm
#' @importFrom stats vcov
#' @param object  an object of class \code{hensm}.
#' @method vcov hensm
vcov.hensm <- function(object, ...) {
  if (!is.null(object$reg_coef_idx)) {
    warning(
      'Computing vcov on object fit with regularization; statistical properties are dubious.'
    )
  }
  vcov(object$mle)
}

# basically FML
# on print overloading
# https://www.rdocumentation.org/packages/mvbutils/versions/2.7.4.1/topics/print
# and see
# https://stackoverflow.com/questions/8414268/define-a-show-method-for-an-s3-class
# https://stackoverflow.com/questions/23724815/roxygen2-issue-with-exporting-print-method

#' @export
#' @importFrom stats printCoefmat
#' @importFrom methods show
#' @inheritParams base::print
#' @rdname hensm
#' @method print hensm
print.hensm <- function(x, ...) {
  show(summary(x$mle))
  invisible(x)
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
