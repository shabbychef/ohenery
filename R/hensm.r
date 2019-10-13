# /usr/bin/r
#
# Copyright 2018-2019 Steven E. Pav. All Rights Reserved.
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
# Copyright: Steven E. Pav, 2018-2019
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

# apparently necessary to register hensm as an S3 class
setOldClass('hensm')

.hensmlik <- function(theta, group, idx, X, wt, eta0) {
	k <- ncol(X)
	beta <- theta[1:k]
	gamma <- theta[(k+1):length(theta)]
	eta <- X  %*% beta 
	if (!is.null(eta0)) { eta <- eta + eta0 }
	hensmlik(group, idx, eta, gamma=gamma, wt=wt)
}
.hensmgrad <- function(theta, group, idx, X, wt, eta0) {
	k <- ncol(X)
	beta <- theta[1:k]
	gamma <- theta[(k+1):length(theta)]
	eta <- X  %*% beta 
	if (!is.null(eta0)) { eta <- eta + eta0 }
	hval <- hensmlik(group, idx, eta, gamma=gamma, wt=wt, deleta=X)
	c(attr(hval,'gradient'),
		attr(hval,'gradgamma'))
}
#2FIX: why isn't there a experts version of this one?
#  @param ngamma  the number of gammas to model; we model
#        \eqn{\gamma_2} through \eqn{\gamma_n}.
.hmfit <- function(y, g, X, wt=NULL, eta0=NULL, normalize_wt=TRUE,
									 ngamma=4,  method=c('BFGS','NR','CG','NM')) {
	method <- match.arg(method)
	stopifnot(ngamma >= 2)
#2FIX: allow beta0 input.
	k <- ncol(X)
	beta0 <- array(0,k)
	gamma0 <- array(1,ngamma-1)
	theta0 <- c(beta0,gamma0)

	if (!is.null(wt) && normalize_wt) { wt <- wt / abs(mean(wt,na.rm=TRUE)) }  # by having the abs, negative weights still throw an error.
	# turn g into integers?
	if (is.integer(g)) { group <- g } else { group <- match(g,unique(g)) }

	idx <- order(g,y,decreasing=TRUE) - 1
	covadj <- .mean_wt(g=g,idx=idx,wt=wt)
	rv <- maxLik(logLik=.hensmlik,grad=.hensmgrad,hess=NULL,
							 start=theta0,method=method,
							 group=group,idx=idx,X=X,wt=wt,eta0=eta0)

	# adjust it! 
	rv$varcovar <- covadj * vcov(rv)
	retv <- list(mle=rv,
							 beta=rv$estimate[1:k],
							 coefficients=rv,
							 gammas=rv$estimate[(k+1):length(theta0)],
							 gamma2=rv$estimate[k+1],
							 estimate=rv$estimate,  # sigh
							 covadj=covadj,
							 wt=wt,
							 g=g,
							 y=y,
							 formula=NULL,
							 eta0=eta0)

	
	gnames <- paste0('gamma',2:ngamma)
	names(retv$mle$estimate) <- c(colnames(X),gnames)
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
#' @template etc
#' @template note-ties
#' @template note-weights
#' @template param-weights
#' @template param-group
#' @return An object of class \code{hensm}, but also of \code{maxLik} with the
#' fit.
#' @keywords fitting
#' @seealso \code{\link{harsm}}, \code{\link{smlik}}.
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
#' hensm(Finish ~ offset(eta0) + fac_age,data=df,group=EventId,weights=weights,ngamma=2)
#'
#'
#' @importFrom stats coef formula model.frame model.matrix na.omit
#' @export
#' @rdname hensm
hensm <- function(formula,data,group=NULL,weights=NULL,ngamma=4,na.action=na.omit) {
	substitute(formula)

	# I find it highly offensive that this cannot be done reasonably
	# easily in a subfunction because of NSE whatever.

	# https://stackoverflow.com/q/53827563/164611
  mf <- match.call(expand.dots = FALSE)
  #turn weights into symbol if character is passed
  if (is.character(mf$weights)) mf$weights <- as.symbol(mf$weights)
  if (is.character(mf$group)) mf$group <- as.symbol(mf$group)
  m <- match(c("formula", "data", "weights", "group", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE 
  mf[[1L]] <- quote(stats::model.frame) 
  mf <- eval(mf, parent.frame()) #evaluate call

	Xs <- model.matrix(formula,mf)
	# remove intercept!
	if (colnames(Xs)[1] == '(Intercept)') { Xs <- Xs[,-1,drop=FALSE] }
	y <- as.vector(model.response(mf))
  group <- as.vector(model.extract(mf,'group'))
	eta0 <- model.offset(mf)
  wt <- as.vector(model.weights(mf))

	dat <- list(Xs=Xs,y=y,group=group,eta0=eta0,wt=wt)

	retv <- .hmfit(y=dat$y, g=dat$group, X=dat$Xs, wt=dat$wt, ngamma=ngamma, eta0=dat$eta0)
	retv <- as.linodds(retv, formula, beta=retv$beta)
	retv
}
#' @export
#' @rdname hensm
#' @importFrom stats vcov
#' @param object  an object of class \code{hensm}.
#' @method vcov hensm
vcov.hensm <- function(object, ...) {
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
#' @rdname hensm
#' @method print hensm
print.hensm <- function(x, ...) {
	show(summary(x$mle))
	invisible(x)
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
