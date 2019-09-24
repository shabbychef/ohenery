# Copyright 2018 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav

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

# env var:
# nb: 
# see also:
# todo:
# changelog: 
#
# Created: 2018.09.25
# Copyright: Steven E. Pav, 2018-2019
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

# helpers#FOLDUP
set.char.seed <- function(str) {
	set.seed(as.integer(charToRaw(str)))
}
#UNFOLD

library(dplyr)
library(numDeriv)
context("code runs at all")#FOLDUP

test_that("harsmlik bits",{#FOLDUP
	# travis only?
	#skip_on_cran()

	nfeat <- 8
	set.seed(321)
	g <- ceiling(seq(0.1,100,by=0.1))
	X <- matrix(rnorm(length(g) * nfeat),ncol=nfeat)
	beta <- rnorm(nfeat)
	eta <- X %*% beta
	y <- rsm(eta,g=g)

	idx <- order(g,y,decreasing=TRUE) - 1
	fastval <- attr(harsmlik(g,idx,eta,deleta=X),'gradient') 
	numap <- grad(function(beta,g,idx,X) { 
									eta <- X %*% beta
									as.numeric(harsmlik(g,idx,eta))
		},x=beta,g=g,idx=idx,X=X)

	expect_equal(fastval,numap,tolerance=0.1)

	# now reweight them
	wt <- runif(length(g))
	wt <- wt / mean(wt)   # mean 1 is recommended
	expect_error(foores <- harsmlik(g,idx,eta,wt=wt),NA)
})#UNFOLD
test_that("hensmlik and harsmlik nest",{#FOLDUP
	# travis only?
	#skip_on_cran()

	nfeat <- 8
	set.seed(321)
	g <- ceiling(seq(0.1,100,by=0.1))
	X <- matrix(rnorm(length(g) * nfeat),ncol=nfeat)
	beta <- rnorm(nfeat)
	eta <- X %*% beta
	y <- rsm(eta,g=g)

	# hensmlik and harsmlik should return the same for gamma=1
	gams <- c(1.0,1.0,1.0)
	idx <- order(g,y,decreasing=TRUE) - 1
	expect_error(hooval <- hensmlik(g=g,idx=idx,eta=eta,gamma=gams,deleta=X),NA)
	expect_true(!is.na(hooval))
	expect_error(sooval <- harsmlik(g=g,idx=idx,eta=eta,deleta=X),NA)
	expect_equal(as.numeric(sooval),as.numeric(hooval))
})#UNFOLD
test_that("hensmlik bits",{#FOLDUP
	# travis only?
	#skip_on_cran()

	nfeat <- 8
	set.seed(321)
	g <- ceiling(seq(0.1,100,by=0.1))
	X <- matrix(rnorm(length(g) * nfeat),ncol=nfeat)
	beta <- rnorm(nfeat)
	eta <- X %*% beta
	y <- rsm(eta,g=g)

	gams <- c(1.0,1.0,1.0)
	idx <- order(g,y,decreasing=TRUE) - 1
	expect_error(hooval <- hensmlik(g=g,idx=idx,eta=eta,gamma=gams,deleta=X),NA)
	expect_true(!is.na(hooval))
	fastval <- attr(hooval,'gradient') 
	numap <- grad(function(beta,g,idx,X,gams) { 
									eta <- X %*% beta
									as.numeric(hensmlik(g,idx,eta,gamma=gams))
		},x=beta,g=g,idx=idx,X=X,gams=gams)
	expect_equal(fastval,numap,tolerance=0.1)

	gams <- c(1.0,0.9,0.8)
	idx <- order(g,y,decreasing=TRUE) - 1
	expect_error(hooval <- hensmlik(g=g,idx=idx,eta=eta,gamma=gams,deleta=X),NA)
	expect_true(!is.na(hooval))
	fastval <- attr(hooval,'gradient') 
	numap <- grad(function(beta,g,idx,X,gams) { 
									eta <- X %*% beta
									as.numeric(hensmlik(g,idx,eta,gamma=gams))
		},x=beta,g=g,idx=idx,X=X,gams=gams)
	expect_equal(fastval,numap,tolerance=0.1)

	# now reweight them
	wt <- runif(length(g))
	wt <- wt / mean(wt)   # mean 1 is recommended
	expect_error(foores <- hensmlik(g,idx,eta,wt=wt,gamma=gams),NA)
})#UNFOLD

test_that("smax bits",{#FOLDUP
	# location invariant
	nfeat <- 10
	set.seed(1001)
	eta <- rnorm(nfeat)
	keta <- eta + 10
	expect_error(mu <- smax(eta),NA)
	expect_error(kmu <- smax(keta),NA)
	expect_equal(mu,kmu,tolerance=1e-7)

	expect_error(smax(c(1,0)),NA)
	expect_error(mubig <- smax(c(10000,0,0,0)),NA)
	expect_equal(mubig,c(1,0,0,0))

	nfeat <- 5
	set.seed(4234)
	g <- ceiling(seq(0.1,100,by=0.1))
	X <- matrix(rnorm(length(g) * nfeat),ncol=nfeat)
	beta <- rnorm(nfeat)
	eta <- X %*% beta
	eta <- data.frame(eta=as.numeric(eta),g=g) %>%
			group_by(g) %>%
				mutate(eta=eta-mean(eta)) %>%
			ungroup() %>%
			{ .$eta }

	expect_error(mu <- smax(eta,g=g),NA)
	expect_error(eta0 <- inv_smax(mu,g=g),NA)
	expect_equal(eta,eta0,tolerance=1e-14)

	nfeat <- 8
	set.seed(808)
	expect_error(mu <- normalize(runif(nfeat)),NA)
	expect_error(eta0 <- inv_smax(mu),NA)
	expect_error(mu1 <- smax(eta0),NA)
	expect_equal(mu,mu1,tolerance=1e-10)

	# fingers crossed.
	mu <- c(1,0,0,0)
	expect_error(eta0 <- inv_smax(mu),NA)
	expect_error(mu1 <- smax(eta0),NA)
	expect_equal(mu,mu1,tolerance=1e-10)


})#UNFOLD
test_that("rsm bits",{#FOLDUP
	# travis only?
	#skip_on_cran()
	nfeat <- 5
	set.seed(1234)
	g <- ceiling(seq(0.1,100,by=0.1))
	X <- matrix(rnorm(length(g) * nfeat),ncol=nfeat)
	beta <- rnorm(nfeat)
	eta <- X %*% beta
	mu <- smax(eta,g=g)
	set.seed(5234)
	expect_error(y1 <- rsm(eta,g=g),NA)
	# first, is it deterministic?
	set.seed(5234)
	expect_error(y2 <- rsm(eta,g=g),NA)
	expect_equal(y1,y2)
	
	# secondly, does it also accept mu?
	set.seed(5234)
	expect_error(y3 <- rsm(g=g,mu=mu),NA)
	expect_equal(y1,y3)
})#UNFOLD
test_that("rhenery bits",{#FOLDUP
	# travis only?
	#skip_on_cran()
	nfeat <- 5
	set.seed(1234)
	g <- ceiling(seq(0.1,100,by=0.1))
	X <- matrix(rnorm(length(g) * nfeat),ncol=nfeat)
	eta <- rnorm(10)
	mu <- smax(eta)
	gamma <- runif(length(eta)-1,min=0.5,max=1.5)
	set.seed(5234)
	expect_error(y1 <- rhenery(mu,gamma),NA)
	# first, is it deterministic?
	set.seed(5234)
	expect_error(y2 <- rhenery(mu,gamma),NA)
	expect_equal(y1,y2)
})#UNFOLD
test_that("erank bits",{#FOLDUP
	set.seed(321)
	for (nfeat in c(1,2,3,4,8)) {
		probs <- runif(nfeat)
		probs <- probs / sum(probs)
		expect_error(erank(probs),NA)
	}

	set.seed(2345)
	eta <- rnorm(7)
	val1 <- rowMeans(replicate(5000,rsm(eta)))
	val2 <- erank(smax(eta))

	expect_equal(val1, val2, tolerance=0.05)
})#UNFOLD
test_that("utils",{#FOLDUP
	set.seed(321)
	x <- rnorm(10)
	expect_error(normalize(x),NA)
	expect_error(smax(x),NA)
})#UNFOLD
test_that("normalize ok",{#FOLDUP
	set.seed(111)
	x <- rnorm(10)
	kx <- 4 * x
	expect_error(y <- normalize(x),NA)
	expect_error(ky <- normalize(kx),NA)
	expect_equal(y, ky, tolerance=1e-8)
})#UNFOLD

#UNFOLD
context("group invariance")# FOLDUP

test_that("smax",{#FOLDUP
	nfeat <- 8
	set.seed(123)
	g <- ceiling(seq(0.1,100,by=0.1))
	X <- matrix(rnorm(length(g) * nfeat),ncol=nfeat)
	beta <- rnorm(nfeat)
	eta <- X %*% beta
	mu <- smax(eta,g=g)

	set.seed(789)
	newidx <- sample.int(length(g),length(g))
	gi <- g[newidx]
	ei <- eta[newidx]
	mui <- smax(ei,g=gi)

	expect_equal(mu[newidx],mui)
})#UNFOLD
test_that("harsm_invlink",{#FOLDUP
	nfeat <- 8
	set.seed(123)
	g <- ceiling(seq(0.1,100,by=0.1))
	X <- matrix(rnorm(length(g) * nfeat),ncol=nfeat)
	beta <- rnorm(nfeat)
	eta <- X %*% beta

	ernk <- harsm_invlink(eta,g=g)

	set.seed(789)
	newidx <- sample.int(length(g),length(g))
	gi <- g[newidx]
	ei <- eta[newidx]
	ernki <- harsm_invlink(ei,g=gi)

	expect_equal(ernk[newidx],ernki)
})#UNFOLD
test_that("rsm invariant wrt reordering",{#FOLDUP
	#nfeat <- 8
	for (nfeat in c(4,8)) { 
		set.seed(123)
		for (g in list(rep(1,20),ceiling(seq(0.1,100,by=0.1)))) {
			X <- matrix(rnorm(length(g) * nfeat),ncol=nfeat)
			beta <- rnorm(nfeat)
			eta <- X %*% beta
			set.seed(115)
			y <- rsm(eta,g=g)

			set.seed(789)
			newidx <- sample.int(length(y),length(y))
			gi <- g[newidx]
			ei <- eta[newidx]
			set.seed(115)
			yo <- rsm(ei,g=gi)
			yi <- y[newidx]
			expect_equal(yo,yi)
		}
	}
})#UNFOLD
test_that("rhenery not invariant wrt reordering",{#FOLDUP
	set.seed(432)
	mu <- normalize(runif(40))
	set.seed(1111)
	expect_error(r1 <- rhenery(mu),NA)
	set.seed(1111)
	expect_error(r2 <- rhenery(rev(mu)),NA)
	expect_true(!all(r1==rev(r2)))
})#UNFOLD

# UNFOLD
context("sm makes sense")# FOLDUP
test_that("sm foo",{#FOLDUP
	nfeat <- 8
	set.seed(123)
	g <- ceiling(seq(0.1,1000,by=0.1))
	X <- matrix(rnorm(length(g) * nfeat),ncol=nfeat)
	set.seed(678)
	beta <- rnorm(nfeat,sd=3)
	eta <- X %*% beta
	y <- rsm(eta,g=g)

	# now the consensus
	beta0 <- beta
	beta0[1] <- 0
	beta0[2] <- 0
	beta0[3] <- 0
	eta0 <- X %*% beta0
	expect_error(foo <- harsmfit(y=y,g=g,X=X,eta0=eta0),NA)
	donotuse <- capture.output(expect_error(print(foo),NA))

	#.mse <- function(x,y) sum((x-y)^2)
	#ppooh <- data.frame(y=foo$y,
											#g=foo$g,
											#er=foo$erank) %>%
		#group_by(g) %>%
			#mutate(qo=rank(er),
						 #qdumb=(n() + 1)/2) %>%
		#ungroup() %>%
		#summarize(ssre=.mse(qo,y),
							#ssto=.mse(qdumb,y)) 
	

})#UNFOLD

# UNFOLD
context("modeling functions")#FOLDUP
test_that("harsmfit bits",{#FOLDUP
	# travis only?
	#skip_on_cran()
	nfeat <- 5
	set.seed(1234)
	g <- ceiling(seq(0.1,100,by=0.1))
	X <- matrix(rnorm(length(g) * nfeat),ncol=nfeat)
	beta <- rnorm(nfeat)
	eta <- X %*% beta
	y <- rsm(eta,g=g)
			 
	expect_error(mod0 <- harsmfit(y=y,g=g,X=X),NA)
	expect_equal(as.numeric(coefficients(mod0)),beta,tolerance=0.1)

	# now the pretty frontend
	data <- cbind(data.frame(outcome=y,race=g),as.data.frame(X))

	fmla <- outcome ~ V1 + V2 + V3 + V4 + V5
	expect_error(fitm <- harsm(fmla,group=race,data=data),NA)
	expect_equal(as.numeric(coefficients(mod0)),as.numeric(coefficients(fitm)),tolerance=0.0001)

	# can deal with a single offset
	fmla <- outcome ~ offset(V1) + V2 
	expect_error(fitm <- harsm(fmla,group=race,data=data),NA)
	fmla <- outcome ~ V1 + offset(V2)
	expect_error(fitm <- harsm(fmla,group=race,data=data),NA)

	# with consensus odds, but there
	eta0 <- rowMeans(X)
	data <- cbind(data.frame(outcome=y,race=g,eta0=eta0),as.data.frame(X))
	fmla <- outcome ~ offset(eta0) + V1 + V2 + V3 + V4 + V5
	expect_error(fitm <- harsm(fmla,group=race,data=data),NA)

	# with consensus odds, but there isn't Z1 and Z2
	fmla <- outcome ~ offset(eta0) + Z1 + Z2
	expect_error(fitm <- harsm(fmla,group=race,data=data))

	# fit with weights
	wt <- runif(length(y))
	data <- cbind(data.frame(outcome=y,race=g,wt=wt),as.data.frame(X))
	fmla <- outcome ~ V1 + V2 + V3 + V4 + V5
	expect_error(fitm <- harsm(fmla,data,group=race,weights='wt'),NA)

	# this should error: negative weights
	data <- cbind(data.frame(outcome=y,race=g,wt=rep(-1,length(y))),as.data.frame(X))
	fmla <- outcome ~ V1 + V2 + V3 + V4 + V5
	expect_error(fitm <- harsm(fmla,data,group=race,weights='wt'))
	
	# check on non-numeric race ids
	data <- cbind(data.frame(outcome=y,race=g),as.data.frame(X))
	data <- data[data$race <= 26,]
	data$letrace <- letters[data$race]
	data$facrace <- factor(data$letrace)

	expect_error(fitnum <- harsm(outcome ~ V1 + V2,data,group=race),NA)
	expect_error(fitlet <- harsm(outcome ~ V1 + V2,data,group=letrace),NA)
	expect_error(fitfac <- harsm(outcome ~ V1 + V2,data,group=facrace),NA)

	expect_equal(as.numeric(coefficients(fitnum)),as.numeric(coefficients(fitlet)),tolerance=1e-7)
	expect_equal(as.numeric(coefficients(fitnum)),as.numeric(coefficients(fitfac)),tolerance=1e-7)

})#UNFOLD
test_that("harsmfit prediction",{#FOLDUP
	# travis only?
	#skip_on_cran()
	nfeat <- 2
	set.seed(1234)
	g <- ceiling(seq(0.1,100,by=0.1))
	X <- matrix(rnorm(length(g) * nfeat),ncol=nfeat)
	beta <- rnorm(nfeat)
	eta <- X %*% beta
	y <- rsm(eta,g=g)
	
	# check on non-numeric race ids
	data <- cbind(data.frame(outcome=y,race=g),as.data.frame(X))
	data <- data[data$race <= 26,]
	data$letrace <- letters[data$race]
	data$facrace <- factor(data$letrace)

	expect_error(fitnum <- harsm(outcome ~ V1 + V2,data,group=race),NA)
	expect_error(fitlet <- harsm(outcome ~ V1 + V2,data,group=letrace),NA)
	expect_error(fitfac <- harsm(outcome ~ V1 + V2,data,group=facrace),NA)

	for (ttype in c('eta','mu','erank')) {
		expect_error(fuh <- predict(fitnum,newdata=data,type=ttype),NA)
		expect_error(fuh <- predict(fitlet,newdata=data,type=ttype),NA)
		expect_error(fuh <- predict(fitnum,newdata=data,type=ttype,group=race),NA)
		expect_error(fuh <- predict(fitlet,newdata=data,type=ttype,group=letrace),NA)
	}
	# deal with na actions
	expect_error(fuh <- as.numeric(predict(fitlet,newdata=data,type='eta',group=letrace,na.action=na.pass)),NA)
	expect_equal(length(fuh),nrow(data))

	# outcome should have no affect on the na action
	badata <- data
	badata$outcome[1] <- NA
	expect_error(fuh <- as.numeric(predict(fitlet,newdata=badata,type='eta',group=letrace,na.action=na.pass)),NA)
	expect_equal(length(fuh),nrow(badata))
	expect_error(fuh <- as.numeric(predict(fitlet,newdata=badata,type='eta',group=letrace,na.action=na.omit)),NA)
	expect_equal(length(fuh),nrow(badata))

	# but for V1
	badata <- data
	badata$V1[1] <- NA
	expect_error(fuh <- as.numeric(predict(fitlet,newdata=badata,type='eta',group=letrace,na.action=na.pass)),NA)
	expect_equal(length(fuh),nrow(badata))
	expect_true(all(is.na(fuh[is.na(badata$V1)])))
	expect_true(all(!is.na(fuh[!is.na(badata$V1)])))

	expect_error(fuh <- as.numeric(predict(fitlet,newdata=badata,type='eta',group=letrace,na.action=na.omit)),NA)
	expect_equal(length(fuh),sum(!is.na(badata$V1)))
	expect_true(all(!is.na(fuh)))




})#UNFOLD
test_that("hensm bits",{#FOLDUP
	# travis only?
	#skip_on_cran()
	nfeat <- 5
	set.seed(1234)
	g <- 1 + ((1:10000) %% 1000)
	X <- matrix(rnorm(length(g) * nfeat),ncol=nfeat)
	beta <- rnorm(nfeat)
	eta <- X %*% beta
	y <- rsm(eta,g=g)

	# now the pretty frontend
	data <- cbind(data.frame(outcome=y,race=g),as.data.frame(X))

	fmla <- outcome ~ V1 + V2 + V3 + V4 + V5
	# deterministic?
	expect_error(fitm <- hensm(fmla,data,group=race),NA)
	expect_error(fitm2 <- hensm(fmla,data,group=race),NA)
	expect_equal(fitm$coefficients,fitm2$coefficients)

	# can deal with a single offset
	fmla <- outcome ~ offset(V1) + V2 
	expect_error(fitm <- hensm(fmla,data,group=race),NA)
	fmla <- outcome ~ V1 + offset(V2)
	expect_error(fitm <- hensm(fmla,data,group=race),NA)

	# with consensus odds, but there
	eta0 <- rowMeans(X)
	data <- cbind(data.frame(outcome=y,race=g,eta0=eta0),as.data.frame(X))
	fmla <- outcome ~ offset(eta0) + V1 + V2 + V3 + V4 + V5
	expect_error(fitm <- hensm(fmla,data,group=race),NA)

	# with consensus odds, but there isn't Z1 and Z2
	fmla <- outcome ~ offset(eta0) + Z1 + Z2
	expect_error(fitm <- hensm(fmla,data,group=race))

	# fit with weights
	wt <- runif(length(y))
	data <- cbind(data.frame(outcome=y,race=g,wt=wt),as.data.frame(X))
	fmla <- outcome ~ V1 + V2 + V3 + V4 + V5
	expect_error(fitm <- hensm(fmla,data,group=race,weights='wt'),NA)

	# this should error: negative weights
	data <- cbind(data.frame(outcome=y,race=g,wt=rep(-1,length(y))),as.data.frame(X))
	fmla <- outcome ~ V1 + V2 + V3 + V4 + V5
	expect_error(fitm <- hensm(fmla,data,group=race,weights='wt'))
	
	# check on non-numeric race ids
	data <- cbind(data.frame(outcome=y,race=g),as.data.frame(X))
	data <- data[data$race <= 26,]
	data$letrace <- letters[data$race]
	data$facrace <- factor(data$letrace)

	expect_error(fitnum <- hensm(outcome ~ V1 + V2,data,group=race),NA)
	expect_error(fitlet <- hensm(outcome ~ V1 + V2,data,group=letrace),NA)
	expect_error(fitfac <- hensm(outcome ~ V1 + V2,data,group=facrace),NA)

	expect_equal(as.numeric(coefficients(fitnum)),as.numeric(coefficients(fitlet)),tolerance=1e-7)
	expect_equal(as.numeric(coefficients(fitnum)),as.numeric(coefficients(fitfac)),tolerance=1e-7)

	#for (ttype in c('eta','mu','erank')) {
		#expect_error(fuh <- predict(fitnum,newdata=data,type=ttype),NA)
		#expect_error(fuh <- predict(fitlet,newdata=data,type=ttype),NA)
		#expect_error(fuh <- predict(fitnum,newdata=data,type=ttype,group=race),NA)
		#expect_error(fuh <- predict(fitlet,newdata=data,type=ttype,group=letrace),NA)
	#}
})#UNFOLD
test_that("hensm consistency",{#FOLDUP
	# travis only?
	skip_on_cran()
	nfeat <- 2
	set.seed(1234)
	g <- 1 + ((1:120000) %% 20000)
	X <- matrix(rnorm(length(g) * nfeat),ncol=nfeat)
	beta <- rnorm(nfeat)
	eta <- X %*% beta
	mu <- smax(eta,g=g)
	gammas <- c(0.9,0.8,0.7,rep(1,6))

	# now the pretty frontend
	data <- cbind(data.frame(mu=mu,race=g),as.data.frame(X)) %>%
		group_by(race) %>%
			mutate(outcome=rhenery(mu,gammas[1:(n()-1)])) %>%
		ungroup()

	fmla <- outcome ~ V1 + V2 
	expect_error(fitm <- hensm(fmla,data,group=race,ngamma=5),NA)

	# close enough
	expect_equal(fitm$gammas[1:3],gammas[1:3],tolerance=0.03)
	expect_equal(as.numeric(fitm$beta),beta,tolerance=0.03)
})#UNFOLD
#UNFOLD


#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
