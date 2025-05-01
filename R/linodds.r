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
# Created: 2018.10.18
# Copyright: Steven E. Pav, 2018-2024
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

# apparently necessary to register linodds as an S3 class
setOldClass('linodds')

#' @title An object for modeling linear odds.
#'
#' @description
#'
#' A model for odds linear in some feature.
#'
#' @details
#'
#' An object which holds a formula, some fit coefficients
#' \eqn{\beta} which fit in that formula to generate odds
#' in odds space.
#' The odds can then be converted, via \code{predict.linodds}
#' to probabilities,
#' or to expected ranks under the Harville model.
#' Both \code{\link{harsm}} and \code{\link{hensm}} return
#' objects of class \code{linodds}.
#'
#' We think of linear odds as
#' \eqn{\eta = x^{\top}\beta},
#' for independent variables \eqn{x}. The odds, \eqn{\eta}
#' are converted to probabilities, \eqn{\mu} via
#' \eqn{\mu = c \exp{\eta},} where the constant \eqn{c}
#' is chosen so the \eqn{\mu} for a given matching
#' sum to one.
#'
#' @param object  some list-like object.
#' @inheritParams stats::lm
#' @param beta  the fit coefficients.
#' @rdname linodds
#' @seealso \code{\link{harsm}}, \code{\link{hensm}}.
#' @template etc
#' @export
as.linodds <- function(object, formula, beta) {
  object$formula <- formula
  object$beta <- beta
  class(object) <- c(class(object), 'linodds')
  object
}

#' @rdname linodds
#' @importFrom stats predict
#' @param newdata  a \code{data.frame} from which we can extract a model
#' frame via the formula of the \code{object}.
#' @template param-group
#' @param ... other arguments.
#' @param type  indicates which prediction should be returned:
#' \describe{
#' \item{\code{eta}}{The odds.}
#' \item{\code{mu}}{The probability.}
#' \item{\code{erank}}{The expected rank.}
#' }
#' @param na.action  How to deal with missing values in \code{y}, \code{g},
#' \code{X}, \code{wt}, \code{eta0}.
#' @seealso \code{\link{smax}}, \code{\link{harsm_invlink}}.
#' @importFrom stats delete.response terms model.offset model.matrix model.extract as.formula na.pass
#' @export
#' @method predict linodds
predict.linodds <- function(
  object,
  newdata,
  type = c('eta', 'mu', 'erank'),
  na.action = na.pass,
  group = NULL,
  ...
) {
  type <- match.arg(type)

  fmla <- object$formula
  tt <- terms(fmla)
  Terms <- delete.response(tt)

  # https://stackoverflow.com/q/53827563/164611
  mf <- match.call(expand.dots = FALSE)
  #turn weights into symbol if character is passed
  if (is.character(mf$group)) mf$group <- as.symbol(mf$group)
  # wheee
  names(mf) <- gsub('^newdata$', 'data', names(mf))

  m <- match(c("data", "group", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, 1L, m)]
  mf$drop.unused.levels <- TRUE
  # need this
  # mf$xlev <-
  mf[[1L]] <- quote(stats::model.frame)
  mf[[2L]] <- Terms
  mf <- eval(mf, parent.frame()) #evaluate call
  Xs <- model.matrix(as.formula(Terms), mf)
  # remove intercept!
  if (colnames(Xs)[1] == '(Intercept)') {
    Xs <- Xs[, -1, drop = FALSE]
  }
  group <- as.vector(model.extract(mf, 'group'))
  eta0 <- model.offset(mf)
  wt <- as.vector(model.weights(mf))

  dat <- list(Xs = Xs, group = group, eta0 = eta0, wt = wt)
  if (!all(colnames(dat$Xs) %in% names(object$beta))) {
    stop("some levels in data unknown to fit model")
  }
  # 2FIX: we subset beta to the colnames, but are these promised to be in the right order?
  eta <- as.numeric(dat$Xs %*% matrix(object$beta[colnames(dat$Xs)], ncol = 1))
  # deal with offset
  if (!is.null(eta0)) {
    eta <- eta + eta0
  }
  retval <- switch(
    type,
    eta = eta,
    mu = smax(eta = eta, g = dat$group),
    erank = harsm_invlink(eta = eta, g = dat$group)
  )
  attr(retval, 'na.action') <- attr(dat, 'na.action')
  retval
}

#' @rdname linodds
#' @importFrom stats coef
#' @export
#' @method coef linodds
coef.linodds <- function(object, ...) {
  object$beta
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
