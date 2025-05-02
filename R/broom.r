# /usr/bin/r
#
# Copyright 2025-2025 Steven E. Pav. All Rights Reserved.
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
# Created: 2025.04.30
# Copyright: Steven E. Pav, 2025
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' @title Tidy
#'
#' @description
#'
#' Tidy a harsm, hensm or linodds model.
#'
#' @details
#'
#' Returns a table with information on the fit coefficients
#' and standard errors of an estimated Harville or Henery model.
#'
#' @usage
#'
#' tidy(x, ...)
#'
#' @param x an object of type \code{harsm} or \code{hensm}
#' @param ... arguments for generic consistency.
#' @return A tidy \code{tibble::tibble()} with fields
#' \describe{
#'  \item{term}{The name of the estimated parameter. Betas typically come before gammas.}
#'  \item{estimate}{The estimated parameter.}
#'  \item{std.error}{The standard error of the estimate.}
#'  \item{statistic}{The z-statistic of the estimate.}
#' }
#' @seealso \code{\link{tidy.maxLik}}.
#'
#' @examples
#'
#' # softmax on the Best Picture data
#' data(best_picture)
#' df <- best_picture
#' df$place <- ifelse(df$winner,1,2)
#' df$weight <- ifelse(df$winner,1,0)
#'
#' fmla <- place ~ nominated_for_BestDirector + nominated_for_BestActor + Drama
#' fit0 <- harsm(fmla,data=df,group=year,weights=weight)
#' print(tidy(fit0))
#'
#' @template etc
#' @rdname tidy
#' @export
tidy.harsm <- function(x, ...) {
  result <- tidy(x$mle, ...)
}
#' @rdname tidy
#' @export
tidy.hensm <- function(x, ...) {
  result <- tidy(x$mle, ...)
}

#' @title Glance
#'
#' @description
#'
#' Glance at a harsm, hensm or linodds model.
#'
#' @details
#'
#' Returns a table with information on the overall fit
#' an estimated Harville or Henery model.
#'
#' @usage
#'
#' glance(x, ...)
#'
#' @param x an object of type \code{harsm} or \code{hensm}
#' @param ... arguments for generic consistency.
#' @return A glanced \code{tibble::tibble()} with fields
#' \describe{
#'  \item{df}{The degrees of freedom of the model.}
#'  \item{logLik}{The log-likelihood of the model.}
#'  \item{AIC}{Akaike's Information Criterion for the model.}
#'  \item{nobs}{The number of observations, if this is available, otherwise ‘NA’.}
#' }
#' @note
#' In the future this may include information about the regularization, if any.
#' @seealso \code{\link{glance.maxLik}}.
#'
#' @examples
#'
#' # softmax on the Best Picture data
#' data(best_picture)
#' df <- best_picture
#' df$place <- ifelse(df$winner,1,2)
#' df$weight <- ifelse(df$winner,1,0)
#'
#' fmla <- place ~ nominated_for_BestDirector + nominated_for_BestActor + Drama
#' fit0 <- harsm(fmla,data=df,group=year,weights=weight)
#' print(glance(fit0))
#'
#' @template etc
#' @rdname glance
#' @export
glance.harsm <- function(x, ...) {
  result <- glance(x$mle, ...)
}
#' @rdname glance
#' @export
glance.hensm <- function(x, ...) {
  result <- glance(x$mle, ...)
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
