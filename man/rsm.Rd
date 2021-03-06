% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsm.r
\name{rsm}
\alias{rsm}
\title{Generate variates from a softmax distribution.}
\usage{
rsm(eta, g = NULL, mu = NULL, gamma = NULL)
}
\arguments{
\item{eta}{a vector of the odds.
Must be the same length as \code{g} if \code{g} is given.}

\item{g}{a vector giving the group indices. If \code{NULL},
then we assume only one group is in consideration.}

\item{mu}{a vector of the probablities.
Must be the same length as \code{g} if \code{g} is given.
If \code{mu} and \code{eta} are both given, we ignore
\code{eta} and use \code{mu}.}

\item{gamma}{a vector of the gamma parameters for the
Henery model. Typically the first element should be 1.
Omit or set all values to 1 to recover the Harville model.
The last element will be extended if the gamma is
not long enough for a given group. Note that gamma is expected
to be very short, and is not \sQuote{aligned} with
\code{eta} in any way.}
}
\value{
a vector of integers, each a permutation of one through
the number of elements in each group.
}
\description{
Generate variates from a softmax distribution
under Harville or Henery models.
}
\details{
Given the \eqn{\eta} in odds space, and a grouping
variable, returns values in one to the number of
elements in a group. That is, we have a permutation
of 1 through \eqn{n_g} on each group as output.

For a single group, the probability that the \eqn{i}th
element of the output is a 1 is equal to
\deqn{\pi_{1,i} = \frac{\mu_i^{\gamma_1}}{\sum_j \mu_j^{\gamma_1}}.}
Once an element has been selected to have output
1, remove it from the set and iterate, but with the
next \eqn{\gamma} elements.
}
\note{
The output of this function satisfies a kind
of order invariance that \code{\link{rhenery}}
does not, at the cost of some computational inefficiency.
Namely that for a fixed randseed, and for distinct
\code{eta}, the output is equivariant with respect
to permutation of the vector \code{eta} when there
is only one group.

Regarding the \sQuote{direction}, we associate
higher odds with a smaller outcome. That is, the ith element of
the output encodes the place that the ith participant 
took in the simulated \sQuote{race}; it should be small if the
odds for that participant are very high.
}
\examples{
# simple use
set.seed(1234)
g <- ceiling(seq(1,10,by=0.1))
eta <- rnorm(length(g))
y <- rsm(eta,g=g)

# same same:
set.seed(235)
y1 <- rsm(eta,g=g)
set.seed(235)
y2 <- rsm(g=g,mu=smax(eta,g=g))
y1 - y2

# the default model is Harville
set.seed(1212)
y1 <- rsm(eta,g=g)
set.seed(1212)
y2 <- rsm(eta,g=g,gamma=c(1,1,1,1))
y1 - y2

# repeat several times with the cards stack against
# early runners
set.seed(1234)
colMeans(t(replicate(1000,rsm(sort(rnorm(10,sd=1.5)),g=rep(1,10)))))

# illustrate the invariance
mu <- (1:10) / 55
set.seed(1414)
y1 <- rsm(mu=mu,gamma=c(1,1,1))
set.seed(1414)
y2 <- rev(rsm(mu=rev(mu),gamma=c(1,1,1)))
y1 - y2

\donttest{
nfeat <- 5
set.seed(1234)
g <- ceiling(seq(0.1,1000,by=0.1))
X <- matrix(rnorm(length(g) * nfeat),ncol=nfeat)
beta <- rnorm(nfeat)
eta <- X \%*\% beta
y <- rsm(eta,g=g)

idx <- order(g,y,decreasing=TRUE) - 1
fooey <- harsmlik(g,idx,eta,deleta=X)
set.seed(3493)
dib <- rnorm(length(beta))
xvl <- seq(-0.01,0.01,length.out=301)
rsu <- sapply(xvl,
              function(del) {
                beta1 <- beta + del * dib
                eta1 <- X \%*\% beta1
                fooey2 <- harsmlik(g,idx,eta1,deleta=X)
                as.numeric(fooey2) - as.numeric(fooey)
              })
drv <- sapply(xvl,
              function(del) {
                beta1 <- beta + del * dib
                eta1 <- X \%*\% beta1
                fooey2 <- harsmlik(g,idx,eta1,deleta=X)
                sum(attr(fooey2,'gradient') * dib)
              })

if (require('ggplot2') && require('dplyr')) {
  bestx <- xvl[which.max(rsu)]
  ph <- data.frame(x=xvl,lik=rsu,grd=drv) \%>\%
    ggplot(aes(x=x,y=lik)) + 
    geom_point() + 
    geom_line(aes(y=grd/200)) +
    geom_vline(xintercept=bestx,linetype=2,alpha=0.5) + 
    geom_hline(yintercept=0,linetype=2,alpha=0.5)
  print(ph)
}
}

if (require('dplyr') && require('knitr')) {
  # expect this to be very small, almost always 1
  set.seed(1234)
  simdraw <- replicate(10000,{
    rsm(eta=c(100,rnorm(7)))[1]
  })

  as.data.frame(table(simdraw)) \%>\%
    mutate(prob=Freq / sum(Freq)) \%>\%
    knitr::kable()

  # expect this to be uniform on 2 through 8
  set.seed(1234)
  simdraw <- replicate(10000,{
    rsm(eta=c(100,rnorm(7)))[2]
  })

  as.data.frame(table(simdraw)) \%>\%
    mutate(prob=Freq / sum(Freq)) \%>\%
    knitr::kable()
}

}
\author{
Steven E. Pav \email{shabbychef@gmail.com}
}
\keyword{probability}
