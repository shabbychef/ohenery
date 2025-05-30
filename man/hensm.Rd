% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hensm.r
\name{hensm}
\alias{hensm}
\alias{vcov.hensm}
\alias{print.hensm}
\title{Friendly interface to softmax regression under Henery model.}
\usage{
hensm(
  formula,
  data,
  group = NULL,
  weights = NULL,
  ngamma = 4,
  fit0 = NULL,
  na.action = na.omit
)

\method{vcov}{hensm}(object, ...)

\method{print}{hensm}(x, ...)
}
\arguments{
\item{formula}{an object of class \code{"\link[stats]{formula}"} (or one that
    can be coerced to that class): a symbolic description of the
    model to be fitted.  The details of model specification are given
    under \sQuote{Details}.}

\item{data}{an optional data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    the variables in the model.  If not found in \code{data}, the
    variables are taken from \code{environment(formula)},
    typically the environment from which \code{lm} is called.}

\item{group}{the string name of the group variable in the data,
or a bare character with the group name. The group indices
need not be integers, but that is more efficient. 
They need not be sorted.}

\item{weights}{an optional vector of weights, or the string or bare name of the 
weights in the \code{data} for use in the fitting process. The weights
are attached to the outcomes, not the participant.  Set to \code{NULL}
for none.}

\item{ngamma}{The number of gammas to fit. Should be at least 2.}

\item{fit0}{An optional object of class \code{hensm} or of \code{harsm} 
with the initial fit estimates. 
These will be used for \sQuote{warm start} of the estimation procedure. 
A warm start should only speed up estimation, not change the ultimate results. 
When there is mismatch between the coefficients in \code{fit0} and the model 
being fit here, the missing coefficients are initialized as zero. 
If \code{ngamma} is \code{NULL} and \code{fit0} is given, 
we default to the number of gammas in the initial fit, otherwise
we fill any missing gammas with 1.
If a \code{harsm} object is given, then \code{ngamma} must be non-null.}

\item{na.action}{How to deal with missing values in \code{y}, \code{g},
\code{X}, \code{wt}, \code{eta0}.}

\item{object}{an object of class \code{hensm}.}

\item{...}{For \code{lm()}: additional arguments to be passed to the low level
    regression fitting functions (see below).}

\item{x}{an object used to select a method.}
}
\value{
An object of class \code{hensm}, but also of \code{maxLik} with the
fit.
}
\description{
A user friendly interface to the softmax regression under the Henery model.
}
\details{
Performs a softmax regression by groups, via Maximum Likelihood Estimation.
It is assumed that successive sub-races maintain the proportional
probability of the softmax, up to some gamma coefficients, 
\eqn{\gamma_2, \gamma_3, ..., \gamma_n}, which we fit. This model
nests the Harville model fit by \code{\link{harsm}}, by fixing all
the gammas equal to 1.
}
\note{
This regression may give odd results when the outcomes are tied,
imposing an arbitrary order on the tied outcomes.
Moreover, no warning may be issued in this case.
In future releases, ties may be dealt with differently,
perhaps in analogy to how ties are treated in the
Cox Proportional Hazards regression, using
the methods of Breslow or Efron.

To avoid incorrect inference when only the top
performers are recorded, and all others are 
effectively tied, one should use weighting.
Set the weights to zero for participants who
are tied non-winners, and one for the rest
So for example, if you observe the Gold, Silver,
and Bronze medal winners of an Olympic event
that had a starting field of 12 participants,
set weights to 1 for the medal winners, and 0
for the others. Note that the weights do not
attach to the participants, they attach to
the place they took.

Since version 0.1.0 of this package, the
normalization of weights used in this function
have changed under the hood. This is to
give correct inference in the case where
zero weights are used to signify finishing
places were not observed. If in doubt,
please confirm inference by simulations,
taking as example the simulations in the README.
}
\examples{

nfeat <- 5
set.seed(1234)
g <- ceiling(seq(0.1,1000,by=0.1))
X <- matrix(rnorm(length(g) * nfeat),ncol=nfeat)
beta <- rnorm(nfeat)
eta <- X \%*\% beta
# 2FIX: do rhenery
y <- rsm(eta,g)
# now the pretty frontend
data <- cbind(data.frame(outcome=y,race=g),as.data.frame(X))

fmla <- outcome ~ V1 + V2 + V3 + V4 + V5
fitm <- hensm(fmla,data,group=race)

# with offset
eta0 <- rowMeans(X)
data <- cbind(data.frame(outcome=y,race=g,eta0=eta0),as.data.frame(X))
fmla <- outcome ~ offset(eta0) + V1 + V2 + V3 + V4 + V5
fitm <- hensm(fmla,data,group=race)

# on horse race data
library(dplyr)
data(race_data)
df <- race_data \%>\%
	group_by(EventId) \%>\%
		mutate(eta0=log(WN_pool / sum(WN_pool))) \%>\%
	ungroup() \%>\%
	mutate(weights=ifelse(!is.na(Finish),1,0)) \%>\%
	mutate(fac_age=cut(Age,c(0,3,5,7,Inf),include.lowest=TRUE)) 

# Henery Model with market efficiency
hensm(Finish ~ eta0,data=df,group=EventId,weights=weights,ngamma=3)

# look for age effect not captured by consensus odds.
fmla <- Finish ~ offset(eta0) + fac_age
fit0 <- hensm(fmla,data=df,group=EventId,weights=weights,ngamma=2)
# allow warm start.
fit1 <- hensm(fmla,data=df,group=EventId,weights=weights,fit0=fit0,ngamma=2)
# allow warm start with more gammas.
fit2 <- hensm(fmla,data=df,group=EventId,weights=weights,fit0=fit0,ngamma=3)
# or a different formula
fit3 <- hensm(update(fmla,~ . + PostPosition),data=df,group=EventId,weights=weights,fit0=fit0)

# warm start from harsm object
fit0_har <- harsm(fmla,data=df,group=EventId,weights=weights)
fit4 <- hensm(fmla,data=df,group=EventId,fit0=fit0_har,weights=weights)

}
\seealso{
\code{\link{harsm}}, \code{\link{smlik}}.
}
\author{
Steven E. Pav \email{shabbychef@gmail.com}
}
\keyword{fitting}
