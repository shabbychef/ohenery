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
# Created: 2025.04.27
# Copyright: Steven E. Pav, 2025
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

.check_regularization <- function(
  coef,
  reg_wt,
  reg_zero,
  reg_power,
  reg_coef_idx
) {
  if (is.null(reg_wt)) {
    return
  }
  # 2FIX: check for compatible sizes?
  stopifnot(all((reg_coef_idx >= 1) & (reg_coef_idx <= length(coef))))
  stopifnot(all(reg_wt >= 0))
  stopifnot(all(reg_power >= 0))
}

.regularization_term <- function(
  coef,
  reg_wt,
  reg_zero,
  reg_power,
  reg_coef_idx
) {
  if (is.null(reg_wt)) {
    return(0)
  }
  -sum(reg_wt * (abs(coef[reg_coef_idx] - reg_zero)^reg_power))
}

.regularization_grad <- function(
  coef,
  reg_wt,
  reg_zero,
  reg_power,
  reg_coef_idx
) {
  if (is.null(reg_wt)) {
    return(0)
  }
  gradbits <- -reg_wt *
    reg_power *
    (abs(coef[reg_coef_idx] - reg_zero)^(reg_power - 1)) *
    sign(coef[reg_coef_idx] - reg_zero)
  grad <- rep(0, length(coef))
  for (idx in seq_along(reg_coef_idx)) {
    grad[reg_coef_idx[idx]] <- grad[reg_coef_idx[idx]] + gradbits[idx]
  }
  grad
}

# the reg_zero should default to zero for beta parts of the coefficient,
# and one for the gamma parts of the coefficient. if the user gives a null
# for reg_zero, then fill in those sensible defaults.
.regularization_default_zero <- function(reg_zero, reg_coef_idx, num_beta) {
  if (is.null(reg_zero)) {
    reg_zero <- as.numeric(reg_coef_idx > num_beta)
  }
  reg_zero
}

.zero_to_one <- function(z) {
  ifelse(z == 0, 1, z)
}

# standardize the reg_wt
.regularization_standardize <- function(
  reg_wt,
  reg_power,
  reg_coef_idx,
  reg_standardize,
  X
) {
  if (reg_standardize) {
    if (length(reg_wt) == 1) {
      reg_wt <- rep(reg_wt, length(reg_coef_idx))
    }
    if (length(reg_power) == 1) {
      reg_power <- rep(reg_power, length(reg_coef_idx))
    }
    stds <- apply(X, FUN = sd, MARGIN = 2)
    check_us <- reg_coef_idx[reg_coef_idx <= ncol(X)]
    if (any(stds[check_us] == 0)) {
      warning(
        "Design matrix has some columns with zero standard deviation; will not standardize these"
      )
    }
    for (idx in seq_along(reg_coef_idx)) {
      if (reg_coef_idx[idx] <= ncol(X)) {
        reg_wt[idx] <- reg_wt[idx] *
          .zero_to_one(stds[reg_coef_idx[idx]]**reg_power[idx])
      }
    }
  }
  reg_wt
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
