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

# Created: 2018-09-17
# Copyright: Steven E. Pav, 2018-2019
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

#' Modeling of ordinal outcomes via the softmax function under
#' the Harville and Henery models.
#'
#' @section Harville and Henery models:
#'
#' The Harville and Henery models describe the probability of 
#' ordered outcomes in terms of some parameters. 
#' Typically the ordered outcomes are things like place 
#' in a race, or winner among a large number of contestants.
#' The Harville model could be described as a softmax probability
#' for the first place finish, with a recursive model on the
#' remaining places.
#' The Henery model generalizes that to adjust the remaining
#' places with another parameter.
#' 
#' These are best illustrated with an example.
#' Suppose you observe a race of 20 contestants.
#' Contestant number 11 takes first place,
#' number 6 takes second place, and 17 takes third place,
#' while the fourth through twentieth places are not
#' recorded or not of interest.
#' Under the Harville model, the probability of this outcome
#' can be expressed as
#' \deqn{\frac{\mu_{11}}{\sum_i \mu_i} \frac{\mu_6}{\sum_{i \ne 11} \mu_i}
#' \frac{\mu_{17}}{\sum_{i \ne 11, i \ne 6} \mu_i},}
#' where \eqn{\mu_i = \exp{\eta_i}}.
#' In a softmax regression under the Harville model, 
#' one expresses the odds as \eqn{\eta_i = x_i^{\top}\beta}, where
#' \eqn{x_i} are independent variables, for some
#' \eqn{\beta} to be fit by the regression.
#'
#' Under the Henery model, one adds gammas, \eqn{\gamma_2, \gamma_3, ...} such
#' that the probability of the outcome above is
#' \deqn{\frac{\mu_{11}}{\sum_i \mu_i} \frac{\mu_6^{\gamma_2}}{\sum_{i \ne 11} \mu_i^{\gamma_2}}
#' \frac{\mu_{17}^{\gamma_3}}{\sum_{i \ne 11, i \ne 6} \mu_i^{\gamma_3}}.}
#' There is no reason to model a \eqn{\gamma_1} as anything but one,
#' since it would be redundant. 
#' The Henery softmax regression estimates the \eqn{\beta} as well as 
#' the \eqn{\gamma_j}. 
#' To simplify the regression, the higher order gammas are assumed to equal
#' the last fit value. That is, we usually model
#' \eqn{\gamma_5=\gamma_4=\gamma_3}.
#'
#' The regression supports weighted estimation as well. The weights are
#' applied to the \emph{places}, not to the participants. The
#' weighted likelihood under the example above, for the Harville model
#' is
#' \deqn{\left(\frac{\mu_{11}}{\sum_i \mu_i}\right)^{w_1} \left(\frac{\mu_6}{\sum_{i \ne 11} \mu_i}\right)^{w_2}
#' \left(\frac{\mu_{17}}{\sum_{i \ne 11, i \ne 6} \mu_i}\right)^{w_3}.}
#' The weighting mechanism is how this package deals with unobserved
#' places.
#' Rather than marking all runners-up as tied for fourth place, in this
#' case one sets the \eqn{w_i=0} for \eqn{i > 3}.
#' The regression is then not asked to make distinctions between the
#' tied runners-up.
#'
#' @section Breaking Changes:
#'
#' This package is a work in progress. Expect breaking changes.
#' Please file any bug reports or issues at
#' \url{https://github.com/shabbychef/ohenery/issues}.
#' 
#' @section Legal Mumbo Jumbo:
#'
#' ohenery is distributed in the hope that it will be useful,
#' but WITHOUT ANY WARRANTY; without even the implied warranty of
#' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#' GNU Lesser General Public License for more details.
#'
#' @template etc
#'
#' @name ohenery
#' @rdname ohenery
#' @docType package
#' @keywords package
#' @title The 'ohenery' package.
#' @useDynLib ohenery
#' @template ref-harville
#' @template ref-henery
#' @importFrom Rcpp evalCpp
#' @note
#' 
#' This package is maintained as a hobby. 
#'
NULL

#' @title News for package 'ohenery':
#'
#' @description 
#'
#' News for package \sQuote{ohenery}
#'
#' \newcommand{\CRANpkg}{\href{https://cran.r-project.org/package=#1}{\pkg{#1}}}
#' \newcommand{\ohenery}{\CRANpkg{ohenery}}
#'
#' @section \ohenery{} Initial Version 0.1.1 (2018-10-14) :
#' \itemize{
#' \item Change default in harsm and hensm to use unnormalized weights,
#' correcting inference when not all finishes are observed.
#' }
#'
#' @section \ohenery{} Initial Version 0.1.0 (2018-10-01) :
#' \itemize{
#' \item first CRAN release.
#' }
#'
#' @name ohenery-NEWS
#' @rdname NEWS
NULL

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
