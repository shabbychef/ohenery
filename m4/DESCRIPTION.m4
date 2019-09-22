dnl divert here just means the output from basedefs does not appear.
divert(-1)
include(basedefs.m4)
divert(0)dnl
Package: PKG_NAME()
Maintainer: Steven E. Pav <shabbychef@gmail.com>
Authors@R: c(person(c("Steven", "E."), "Pav", 
    role=c("aut","cre"),
    email="shabbychef@gmail.com",
    comment = c(ORCID = "0000-0002-4197-6195")))
Version: VERSION()
Date: DATE()
License: LGPL-3
Title: Modeling of Ordinal Random Variables via Softmax
BugReports: https://github.com/shabbychef/PKG_NAME()/issues
Description: Supports the modeling of ordinal random variables, 
    like races, via the softmax function. 
In particular, fits the Henery <doi:10.1111/j.2517-6161.1981.tb01153.x>,
and Harville <doi:10.1080/01621459.1973.10482425> models.
Depends: 
    R (>= 3.0.2)
Imports:
    Rcpp (>= 0.12.3),
    maxLik,
    magrittr,
dnl methods,
    rlang,
dnl matrixcalc,
    tidyr,
    dplyr
LinkingTo: Rcpp
Suggests: 
    microbenchmark,
    testthat, 
    numDeriv,
    knitr
dnl https://github.com/r-lib/devtools/issues/1776
Encoding: UTF-8
URL: https://github.com/shabbychef/PKG_NAME()
dnl VignetteBuilder: knitr
Collate:
m4_R_FILES()
dnl vim:ts=2:sw=2:tw=79:syn=m4:ft=m4:et
