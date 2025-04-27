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
# Created: 2018.09.29
# Copyright: Steven E. Pav, 2018-2024
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav


# apparently necessary to register harsm as an S3 class
setOldClass('harsm')

.harsmlik <- function(beta, grp, idx, X, wt, eta0, ...) {
	eta <- X  %*% beta 
	if (!is.null(eta0)) { eta <- eta + eta0 }
	harsmlik(grp, idx, eta, wt) + .regularization_term(beta, ...)
}
.harsmgrad <- function(beta, grp, idx, X, wt, eta0, ...) {
	eta <- X  %*% beta 
	if (!is.null(eta0)) { eta <- eta + eta0 }
	attr(harsmlik(grp, idx, eta, wt, deleta=X),'gradient') + .regularization_grad(beta, ...)
}
.wmse <- function(x,y,wt=NULL,na.rm=TRUE) { 
	if (!is.null(wt)) {
		retv <- sum(wt*(x-y)^2,na.rm=na.rm) 
	} else {
		retv <- sum((x-y)^2,na.rm=na.rm) 
	}
	retv
}

# to get CRAN checks to not complain about mutated variables
globalVariables(c('dumb_rank','.'))
#. <- NULL
#rm(.)

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
#' The code relies on the likelihood function of \code{\link{harsmlik}},
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
#' @param beta0  an optional vector of the initial estimate of beta for
#' \sQuote{warm start} of the estimation procedure.
#' Must be the same length as number of columns in \code{X}.
#' Should only affect the speed of the computation, not the results.
#' Defaults to all zeroes.
#' @param normalize_wt  if \code{TRUE}, we renormalize \code{wt}, if given,
#' to have mean value 1. Note that the default value has changed
#' since version 0.1.0 of this package. Moreover, non-normalized
#' weights can lead to incorrect inference. Use with caution.
#' @template param-regularization
#'
#' @inheritParams maxLik::maxLik
#' @return An object of class \code{harsm}, \code{maxLik}, and \code{linodds}.
#' @keywords fitting
#' @seealso the likelihood function, \code{\link{harsmlik}}, and the
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
harsmfit <- function(y, g, X, wt=NULL, eta0=NULL, beta0=NULL, normalize_wt=FALSE,
										 reg_wt=NULL, reg_zero=0, reg_power=NULL, reg_coef_idx=NULL,
										 method=c('BFGS','NR','CG','NM')) {
	method <- match.arg(method)
	if (is.null(beta0)) {
		beta0 <- array(0,ncol(X))
	} else {
		stopifnot(length(beta0)==ncol(X))
	}
	reg_zero <- .regularization_default_zero(reg_zero, reg_coef_idx, num_beta=length(beta0))
  .check_regularization(beta0, reg_wt, reg_zero, reg_power, reg_coef_idx) 
	if (!is.null(wt) && normalize_wt) { wt <- wt / abs(mean(wt,na.rm=TRUE)) }  # by having the abs, negative weights still throw an error.
	# turn g into integers?
	if (is.integer(g)) { grp <- g } else { grp <- match(g,unique(g)) }

	idx <- order(g,y,decreasing=TRUE) - 1
	rv <- maxLik(logLik=.harsmlik,grad=.harsmgrad,hess=NULL,
							 start=beta0,method=method,
							 grp=grp,idx=idx,X=X,wt=wt,eta0=eta0,
							 reg_wt=reg_wt, reg_zero=reg_zero, reg_power=reg_power, reg_coef_idx=reg_coef_idx)
	retv <- list(mle=rv,
							 coefficients=rv$estimate,
							 estimate=rv$estimate,  # sigh
							 wt=wt,
							 g=g,
							 y=y,
							 formula=NULL,
							 eta0=eta0,
							 reg_wt=reg_wt, reg_zero=reg_zero, reg_power=reg_power, reg_coef_idx=reg_coef_idx)
	# do some summarization
	retv$deviance <- -2 * rv$maximum
	retv$deviance_df <- length(retv$coefficients)
	# now the estimated rank
	deleta <- X %*% retv$coefficients 
	etahat <- deleta + ifelse(!is.null(eta0),eta0,rep(0,length(g)))
	retv$etahat <- etahat
	retv$erank <- harsm_invlink(retv$etahat,g=g)

	SSres <- .wmse(retv$erank,y,wt=wt)
	ssdf <- data.frame(g=g,y=y)
	if (is.null(wt)) { 
		ssdf$wt <- 1
	} else {
		ssdf$wt <- wt
	}

	SStot <- ssdf %>%
		group_by(g) %>%
		mutate(dumb_rank=(1 + n())/2) %>%
		ungroup() %>%
		summarize(err=.wmse(dumb_rank,y,wt=wt)) %>%
		{ .$err }
	retv$R2 <- 1- SSres / SStot

	if (!is.null(eta0)) {
		erank0 <- harsm_invlink(eta0,g=g)
		SSeta0 <- .wmse(erank0,y,wt=wt)
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
#' Note that this computation uses weighted sums of squares,
#' as the weights contribute to the likelihood term.
#' However, the square sum computation does not take into
#' account the standard error of the rank, and so
#' unlike in linear regression, the softmax regression
#' does not always give positive R-squareds,
#' and the statistic is otherwise hard to interpret.
#'
#' @inheritParams stats::lm
#' @param na.action  How to deal with missing values in the outcomes,
#' groups, weights, etc.
#' @param fit0    An optional object of class \code{harsm} or of \code{hensm} 
#' with the initial fit estimates. 
#' These will be used for \sQuote{warm start} of the estimation procedure. 
#' A warm start should only speed up estimation, not change the ultimate results. 
#' When there is mismatch between the coefficients in \code{fit0} and the model 
#' being fit here, the missing coefficients are initialized as zero. 
#' @template param-weights
#' @template param-regularization
#' @template param-group
#' @template etc
#' @template note-ties
#' @template note-weights
#' @return An object of class \code{harsm}, but also of \code{maxLik} with the
#' fit.
#' @keywords fitting
#' @seealso \code{\link{harsmfit}}, \code{\link{harsmlik}}.
#' @template note-normalization
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
#' data <- cbind(data.frame(outcome=y,race=g,eta0=eta0),as.data.frame(X))
#' data$wts <- runif(nrow(data),min=1,max=2)
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
#' fit0 <- harsm(fmla,data=df,group=year,weights=weight) 
#'
#' # warm start is a thing:
#' sub_fmla <- place ~ nominated_for_BestDirector + nominated_for_BestActor 
#' fit1 <- harsm(sub_fmla,data=df,group=year,weights=weight,fit0=fit0) 
#'
#' # ridge regression.
#' fitr2 <- harsm(sub_fmla,data=df,group=year,weights=weight,
#' reg_wt=rep(1,2), reg_power=2, reg_zero=0, reg_coef_idx=c(1,2))
#'
#' # l1 regularization regression.
#' fitr1 <- harsm(sub_fmla,data=df,group=year,weights=weight,
#' reg_wt=rep(1,2), reg_power=1, reg_zero=0, reg_coef_idx=c(1,2))
#'
#' # elasticnet regularization regression.
#' fitr12 <- harsm(sub_fmla,data=df,group=year,weights=weight,
#' reg_wt=rep(1,4), reg_power=c(1,1,2,2), reg_zero=0, reg_coef_idx=c(1,2,1,2))
#'
#' \donttest{
#' # test against logistic regression
#' if (require(dplyr)) {
#' nevent <- 10000
#' set.seed(1234)
#' adf <- data_frame(eventnum=floor(seq(1,nevent + 0.7,by=0.5))) %>%
#'   mutate(x=rnorm(n()),
#'          program_num=rep(c(1,2),nevent),
#'          intercept=as.numeric(program_num==1),
#'          eta=1.5 * x + 0.3 * intercept,
#'          place=ohenery::rsm(eta,g=eventnum))
#' 
#' # Harville model
#' modh <- harsm(place ~ intercept + x,data=adf,group=eventnum)
#' 
#' # the collapsed data.frame for glm
#' ddf <- adf %>%
#'   arrange(eventnum,program_num) %>%
#'   group_by(eventnum) %>%
#'     summarize(resu=as.numeric(first(place)==1),
#'               delx=first(x) - last(x),
#'               deli=first(intercept) - last(intercept)) %>%
#'   ungroup()
#' 
#' # glm logistic fit
#' modg <- glm(resu ~ delx + 1,data=ddf,family=binomial(link='logit'))
#' 
#' all.equal(as.numeric(coef(modh)),as.numeric(coef(modg)),tolerance=1e-4)
#' all.equal(as.numeric(vcov(modh)),as.numeric(vcov(modg)),tolerance=1e-4)
#'
#' }
#'
#' }
#'
#'
#' @importFrom stats coef formula model.frame model.matrix na.omit model.response model.weights
#' @export
#' @rdname harsm
harsm <- function(formula,data,group=NULL,weights=NULL,fit0=NULL,
									reg_wt=NULL, reg_zero=0, reg_power=NULL, reg_coef_idx=NULL,
									na.action=na.omit) {
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

	if (!is.null(fit0)) {
		stopifnot(any(c("hensm","harsm") %in% class(fit0)))
		feat_names <- colnames(dat$Xs)
		beta0 <- rep(0,ncol(dat$Xs))
		found <- feat_names %in% attr(fit0$beta,'names')
		beta0[found] <- fit0$beta[feat_names[found]]
	} else {
		beta0 <- NULL
	}
	# call the fit function
	retv <- harsmfit(y=dat$y, g=dat$group, X=dat$Xs, wt=dat$wt, eta0=dat$eta0, beta0=beta0,
									 reg_wt=reg_wt, reg_zero=reg_zero, reg_power=reg_power, reg_coef_idx=reg_coef_idx)
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
	if (!is.null(object$reg_coef_idx)) {
		warning('Computing vcov on object fit with regularization; statistical properties are dubious.')
	}
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
#' @inheritParams base::print
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
