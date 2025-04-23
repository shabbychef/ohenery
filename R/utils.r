# /usr/bin/r
#
# Copyright 2018-2024 Steven E. Pav. All Rights Reserved.
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
# Created: 2018.09.17
# Copyright: Steven E. Pav, 2018-2024
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

#' @title Normalize a vector to sum to one.
#'
#' @description 
#'
#' Divide a vector by its sum, resulting in a vector with
#' sum equal to one.
#'
#' @param x vector of input data.
#' @return the input divided by its sum. 
#' For the row-wise version, each row is divided by its sum.
#' @template etc
#' @note This function will return \code{NA} when any elements of the input
#' are \code{NA}. May return \code{Inf} if the elements sum
#' to zero.
#' @export
#' @rdname normalize
normalize <- function(x) { x / sum(x) }

globalVariables(c('rowid'))

# this is a big more complicated than I would have liked, but here we go:
# assume the input is already centered:
.smax_one <- function(xc) {
	if (any(xc > 709.78,na.rm=TRUE)) {  # where overflow occurs
		fsrt <- sort.int(xc,decreasing=TRUE,index.return=TRUE) 
		vals <- normalize(c(1,exp(cumsum(diff(fsrt$x)))))
		retv <- rep(0,length(xc))
		retv[fsrt$ix] <- vals
		return(retv)
	}
	# else
	normalize(exp(xc))
}
# 2FIX: check for infinity?
.smax_two <- function(x) { .smax_one(x - mean(x,na.rm=TRUE)) }

#' @title The softmax function.
#'
#' @description 
#'
#' The softmax function: exponentiate a vector and then
#' normalize.
#'
#' @details
#'
#' Given vector \eqn{\eta} for a single group, essentially
#' computes vector \eqn{\mu} defined by 
#' \deqn{\mu_i = \frac{\exp{\eta_i}}{\sum_j \exp{\eta_j}}.}
#'
#' Note that this computation should be invariant with respect
#' to level shifts of the \eqn{\eta}, and thus we de-mean
#' the odds first.
#'
#' @inheritParams rsm
#' @param eta  numeric array of the odds.
#' The odds are de-meaned within each group.
#' @return the exponentiated data normalized.
#' For the row-wise version, each row is soft maxed.
#' @note This function can deal with overflow in a semi-coherent way.
#' @importFrom dplyr group_by mutate ungroup arrange n
#' @template etc
#' @export
#' @seealso \code{\link{normalize}}, \code{\link{inv_smax}}.
#' @examples
#' # we can deal with large values:
#' set.seed(2345)
#' eta <- rnorm(12,sd=1000)
#' smax(eta)
#' @rdname smax
smax <- function(eta,g=NULL) { 
	if (is.null(g) || (all(g==g[1]))) { return(.smax_two(eta)) }
	stopifnot(length(g)==length(eta))
	rv <- data.frame(g=g,eta=as.numeric(eta)) %>%
		mutate(rowid=seq_len(n())) %>%
		group_by(g) %>%
			mutate(retv=.smax_two(eta)) %>%
		ungroup() %>%
		arrange(rowid)
	rv$retv
}

.inv_smax_two <- function(x) { 
	eta <- log(x)
	iie <- is.infinite(eta)
	if (any(iie)) {
		eta[!iie] <- eta[!iie] - mean(eta[!iie],na.rm=TRUE)
		return(eta)
	} 
	eta - mean(eta,na.rm=TRUE)
}

#' @title The inverse softmax function.
#'
#' @description 
#'
#' The inverse softmax function: take a logarithm and center.
#'
#' @details
#'
#' This is the inverse of the softmax function. Given
#' vector \eqn{\mu} for a single group, finds vector 
#' \eqn{\eta} such that
#' \deqn{\eta_i = \log{\mu_i} + c,}
#' where \eqn{c} is chosen such that the \eqn{\eta} sum 
#' to zero:
#' \deqn{c = \frac{-1}{n} \sum_i \log{\mu_i}.}
#'
#' @inheritParams rsm
#' @return the centered log probabilities.
#' @note This function can deal with overflow in a semi-coherent way.
#' @importFrom dplyr group_by mutate ungroup arrange n
#' @template etc
#' @export
#' @seealso \code{\link{smax}}
#' @examples
#' # we can deal with large values:
#' set.seed(2345)
#' eta <- rnorm(12,sd=1000)
#' mu <- smax(eta)
#' eta0 <- inv_smax(mu)
#' @rdname inv_smax
inv_smax <- function(mu,g=NULL) { 
	stopifnot(all(0 <= mu & mu <= 1))
	if (is.null(g) || (all(g==g[1]))) { return(.inv_smax_two(mu)) }
	stopifnot(length(g)==length(mu))
	rv <- data.frame(g=g,mu=as.numeric(mu)) %>%
		mutate(rowid=seq_len(n())) %>%
		group_by(g) %>%
			mutate(retv=.inv_smax_two(mu)) %>%
		ungroup() %>%
		arrange(rowid)
	rv$retv
}

#' @title The inverse link for the softmax.
#'
#' @description 
#'
#' The inverse link function for the softmax. This function
#' takes the group-wise probabilities, \eqn{\mu}, and computes
#' the expected ranks within each group under the Harville
#' model. That is, it is a groupwise computation of 
#' the \code{\link{erank}} function.
#'
#' @param mu  a vector of the probabilities. Should sum to one,
#' at least per group. 
#' Should be the same size as \code{g} if given. 
#' If both \code{mu} and \code{eta} are given, a warning
#' is issued, and the \code{mu} is used.
#' @inheritParams rsm
#' @inheritParams erank
#' @template etc 
#' @return a vector of the ranks. 
#' @seealso the ungrouped version of this, \code{\link{erank}}.
#' @importFrom dplyr group_by mutate ungroup arrange
#' @examples
#'
#' mus <- runif(12)
#' mus <- mus / sum(mus)
#' harsm_invlink(mus)
#'
#' harsm_invlink(mus,c(rep(1,6),rep(2,6)))
#' @export
harsm_invlink <- function(eta, mu=smax(eta,g), g=NULL) {
	if (!missing(eta) && !missing(mu)) { warning('both mu and eta given; taking mu') }
	if (is.null(g)) { return(erank(mu)) }
	if (all(g==g[1])) { return(erank(mu)) }
	rv <- data.frame(g=g,mu=as.numeric(mu)) %>%
		mutate(rowid=seq_len(n())) %>%
		group_by(g) %>%
			mutate(retv=erank(mu)) %>%
		ungroup() %>%
		arrange(rowid)
	rv$retv
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
