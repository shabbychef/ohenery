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
# Created: 2018.09.29
# Copyright: Steven E. Pav, 2018-2019
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

# apparently necessary to register harsm as an S3 class
setOldClass('harsm')

.harsmlik <- function(beta, grp, idx, X, wt, eta0) {
	eta <- X  %*% beta 
	if (!is.null(eta0)) { eta <- eta + eta0 }
	harsmlik(grp, idx, eta, wt)
}
.harsmgrad <- function(beta, grp, idx, X, wt, eta0) {
	eta <- X  %*% beta 
	if (!is.null(eta0)) { eta <- eta + eta0 }
	attr(harsmlik(grp, idx, eta, wt, deleta=X),'gradient')
}
.mse <- function(x,y,na.rm=TRUE) { sum((x-y)^2,na.rm=na.rm) }
#' @title Experts only softmax regression under Harville model.
#'
#' @description 
#'
#' An \dQuote{experts only} softmax fitting function for the Harville
#' model.
#'
#' @details
#' 
#' Given a number of events, indexed by group, and a vector \eqn{y} of
#' the ranks of each entry within that group, perform maximum likelihood
#' estimation under the softmax and proportional probability model.
#'
#' The user can optionally supply a vector of \eqn{\eta_0}, which are
#' taken as the fixed, or \sQuote{consensus} odds. The estimation is
#' then conditional on these fixed odds.
#'
#' Weighted estimation is supported.
#'
#' The code relies on the likelihood function of \code{\link{smlik}},
#' and MLE code from \code{\link[maxLik]{maxLik}}.
#'
#' @param y a vector of the ranked outcomes within each group. Only
#' the order within a group matters.
#' @param g a vector giving the group indices. Need not be integers, but
#' that is more efficient. Need not be sorted.
#' Must be the same length as \code{y}.
#' @param X a matrix of the independent variables. Must have as many rows
#' as the length of \code{y}.
#' @param wt  an optional vector of the observation level weights. These must
#' be non-negative, otherwise an error is thrown. Note that the weight of
#' the last ranked outcome within a group is essentially ignored.
#' Must be the same length as \code{y}.
#' @param eta0  an optional vector of the consensus odds. These are added to
#' the fit odds in odds space before the likelihood caclulation. If given,
#' then when the model is used to predict, similar consensus odds must be
#' given.
#' Must be the same length as \code{y}.
#' @param normalize_wt  if \code{TRUE}, we renormalize \code{wt}, if given,
#' to have mean value 1.
#' @inheritParams maxLik::maxLik
#' @return An object of class \code{harsm}, \code{maxLik}, and \code{linodds}.
#' @keywords fitting
#' @seealso the likelihood function, \code{\link{smlik}}, and the
#' expected rank function (the inverse link), \code{\link{erank}}.
#' @template etc
#' @template ref-harville
#' @importFrom maxLik maxLik
#' @importFrom dplyr summarize
#' @importFrom stats sd
#'
#' @examples 
#' nfeat <- 5
#' set.seed(1234)
#' g <- ceiling(seq(0.1,1000,by=0.1))
#' X <- matrix(rnorm(length(g) * nfeat),ncol=nfeat)
#' beta <- rnorm(nfeat)
#' eta <- X %*% beta
#' y <- rsm(eta,g)
#'      
#' mod0 <- harsmfit(y=y,g=g,X=X)
#' summary(mod0)
#' # now upweight finishers 1-5
#' modw <- harsmfit(y=y,g=g,X=X,wt=1 + as.numeric(y < 6))
#' summary(modw)
#' @export
harsmfit <- function(y, g, X, wt=NULL, eta0=NULL, normalize_wt=TRUE,
										 method=c('BFGS','NR','CG','NM')) {
	method <- match.arg(method)
#2FIX: allow beta0 input.
	beta0 <- array(0,ncol(X))
	if (!is.null(wt) && normalize_wt) { wt <- wt / abs(mean(wt,na.rm=TRUE)) }  # by having the abs, negative weights still throw an error.
	# turn g into integers?
	if (is.integer(g)) { grp <- g } else { grp <- match(g,unique(g)) }

	idx <- order(g,y,decreasing=TRUE) - 1
	rv <- maxLik(logLik=.harsmlik,grad=.harsmgrad,hess=NULL,
							 start=beta0,method=method,
							 grp=grp,idx=idx,X=X,wt=wt,eta0=eta0)
	retv <- list(mle=rv,
							 coefficients=rv$estimate,
							 estimate=rv$estimate,  # sigh
							 wt=wt,
							 g=g,
							 y=y,
							 formula=NULL,
							 eta0=eta0)
	# do some summarization
	retv$deviance <- -2 * rv$maximum
	retv$deviance_df <- length(retv$coefficients)
	# now the estimated rank
	deleta <- X %*% retv$coefficients 
	etahat <- deleta + ifelse(!is.null(eta0),eta0,rep(0,length(g)))
	retv$etahat <- etahat
	retv$erank <- harsm_invlink(retv$etahat,g=g)

	SSres <- .mse(retv$erank,y)
	SStot <- data.frame(g=g,y=y) %>%
		group_by(g) %>%
			mutate(dumb_rank=(1 + n())/2) %>%
		ungroup() %>%
		summarize(err=.mse(dumb_rank,y)) %>%
		{ .$err }
	retv$R2 <- 1- SSres / SStot

	if (!is.null(eta0)) {
		erank0 <- harsm_invlink(eta0,g=g)
		SSeta0 <- .mse(erank0,y)
		retv$delta_R2 <- 1- SSres / SSeta0
	} else {
		retv$delta_R2 <- NA_real_
	}
	class(retv) <- 'harsm'
	retv
}

#' @title Friendly interface to softmax regression under Harville model.
#'
#' @description 
#'
#' A user friendly interface to the softmax regression under the Harville
#' model.
#'
#' @details
#'
#' Performs a softmax regression by groups, via Maximum Likelihood Estimation,
#' under the Harville model. 
#' We fit \eqn{\beta} where odds are \eqn{\eta = x^{\top}\beta} for 
#' independent variables \eqn{x}. 
#' The probability of taking first place is then \eqn{\mu=c\exp{\eta}},
#' where the \eqn{c} is chosen so the \eqn{\mu} sum to one.
#' Under the Harville model, conditional on the first place finisher
#' having been observed, the probability model for second
#' (and successive) places with the probabilities of the remaining
#' participants renormalized.
#'
#' The \code{print} method of the \code{harsm} object includes
#' a display of the R-squared. This measures the improvement
#' in squared errors of the expected rank from the model
#' over the null model which posits that all odds are equal.
#' When the formula includes an offset, a \sQuote{delta R-squared}
#' is also output. This is the improvement in predicted
#' ranks over the model based on the offset term.
#' Note that the expected ranks are only easy to produce
#' under the Harville model; under the Henery model, 
#' the summary R-squared statistics are not produced.
#'
#' @inheritParams stats::lm
#' @template param-group
#' @param na.action  How to deal with missing values in the outcomes,
#' groups, weights, etc.
#' @template etc
#' @return An object of class \code{harsm}, but also of \code{maxLik} with the
#' fit.
#' @keywords fitting
#' @seealso \code{\link{harsmfit}}, \code{\link{smlik}}.
#'
#' @examples 
#'
#' nfeat <- 5
#' set.seed(1234)
#' g <- ceiling(seq(0.1,1000,by=0.1))
#' X <- matrix(rnorm(length(g) * nfeat),ncol=nfeat)
#' beta <- rnorm(nfeat)
#' eta <- X %*% beta
#' y <- rsm(eta,g)
#' # now the pretty frontend
#' data <- cbind(data.frame(outcome=y,race=g),as.data.frame(X))
#' 
#' fmla <- outcome ~ V1 + V2 + V3 + V4 + V5
#' fitm <- harsm(fmla,data,group=race)
#' 
#' eta0 <- rowMeans(X)
#' data <- cbind(data.frame(outcome=y,race=g,eta0=eta0),as.data.frame(X))
#' fmla <- outcome ~ offset(eta0) + V1 + V2 + V3 + V4 + V5
#' fitm <- harsm(fmla,data,group=race)
#'
#' # with weights
#' data <- cbind(data.frame(outcome=y,race=g,eta0=eta0,wts=runif(length(y),min=1,max=2)),as.data.frame(X))
#' fmla <- outcome ~ offset(eta0) + V1 + V2 + V3 + V4 + V5
#' fitm <- harsm(fmla,data,group=race,weights=wts)
#'
#' # softmax on the Best Picture data
#' data(best_picture)
#' df <- best_picture
#' df$place <- ifelse(df$winner,1,2)
#' df$weight <- ifelse(df$winner,1,0)
#'
#' fmla <- place ~ nominated_for_BestDirector + nominated_for_BestActor + Drama 
#'
#' harsm(fmla,data=df,group=year,weights=weight) 
#'
#' @importFrom stats coef formula model.frame model.matrix na.omit
#' @template note-ties
#' @template note-weights
#' @template param-weights
#' @export
#' @rdname sm
harsm <- function(formula,data,group=NULL,weights=NULL,na.action=na.omit) {
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

	# call the fit function
	retv <- harsmfit(y=dat$y, g=dat$group, X=dat$Xs, wt=dat$wt, eta0=dat$eta0)
	names(retv$mle$estimate) <- colnames(dat$Xs)
	names(retv$coefficients) <- colnames(dat$Xs)
	retv <- as.linodds(retv, formula, beta=retv$coefficients)
	retv
}
#' @export
#' @rdname harsm
#' @importFrom stats vcov
#' @param object  an object of class \code{harsm}.
#' @method vcov harsm
vcov.harsm <- function(object, ...) {
	vcov(object$mle)
}

# on print overloading 
# https://www.rdocumentation.org/packages/mvbutils/versions/2.7.4.1/topics/print
# and see
# https://stackoverflow.com/questions/8414268/define-a-show-method-for-an-s3-class
# https://stackoverflow.com/questions/23724815/roxygen2-issue-with-exporting-print-method

#' @export
#' @importFrom stats printCoefmat
#' @importFrom methods show
#' @rdname harsm
#' @method print harsm
print.harsm <- function(x, ...) {
	show(summary(x$mle))

	cat('   R2:',x$R2,'\n')
	if (!is.na(x$delta_R2)) {
		cat('delR2:',x$delta_R2,'\n')
	}
	cat(paste0(rep('-',44),collapse=''),'\n')

	invisible(x)
}


#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
