
/*
 
  This file is part of ohenery.
  
  ohenery is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ohenery is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.
  
  You should have received a copy of the GNU Lesser General Public License
  along with ohenery.  If not, see <http://www.gnu.org/licenses/>.

  softmax likelihood function with gradient attached.

  Created: 2018.09.23
  Copyright: Steven E. Pav, 2018-2019
  Author: Steven E. Pav <shabbychef@gmail.com>
  Comments: Steven E. Pav
*/

#ifndef __DEF_SMLIK__
#define __DEF_SMLIK__

#define MAX(_a_,_b_) ((_a_>_b_)? (_a_):(_b_))
#define MIN(_a_,_b_) ((_a_<_b_)? (_a_):(_b_))

#include <math.h>
#include <string.h>

#endif /* __DEF_SMLIK__ */

#include <Rcpp.h>
using namespace Rcpp;

//' @title
//' Softmax log likelihood under Harville and Henery Models.
//'
//' @description
//'
//' Compute the softmax log likelihood and gradient of the same.
//'
//' @details
//'
//' Given vectors \eqn{g}, \eqn{\eta} and optionally the gradient of \eqn{\eta}
//' with respect to some coefficients, computes the log likelihood under the
//' softmax. The user must provide a reverse ordering as well, which is sorted
//' first by the groups, \eqn{g}, and then, within a group, in increasing
//' quality order. For a race, this means that the index is in order from the
//' last place to the first place in that race, where the group numbers
//' correspond to one race.
//'
//' Under the Harville model, the log likelihood on a given group, where we are indexing
//' in \emph{forward} order, is 
//' \deqn{\left(\eta_1 - \log \sum_{j \ge 1} \mu_j\right) + \left(\eta_2 - \log \sum_{j \ge 2} \mu_j\right) + ...}
//' where \eqn{\mu_i = \exp{\eta_i}}.
//' By \dQuote{forward order}, we mean that \eqn{\eta_1} corresponds to the
//' participant taking first place within that group, \eqn{\eta_2} took second
//' place, and so on.
//' 
//' The Henery model log likelihood takes the form
//' \deqn{\left(\eta_1 - \log \sum_{j \ge 1} \mu_j\right) + \left(\gamma_2 \eta_2 - \log \sum_{j \ge 2} \mu_j^{\gamma_2}\right) + ...}
//' for gamma parameters, \eqn{\gamma}.
//' The Henery model corresponds to the Harville model where all the gammas equal 1.
//'
//' Weights in weighted estimation apply to each summand.
//' The weight for the last place participant in a group is irrelevant.
//' The weighted log likelihood under the Harville model is
//' \deqn{w_1\left(\eta_1 - \log \sum_{j \ge 1} \mu_j\right) + w_2\left(\eta_2 - \log \sum_{j \ge 2} \mu_j\right) + ...}
//' One should think of the weights as applying to the outcome,
//' not the participant.
//' 
//' @param g a vector giving the group indices.
//' @param idx a vector of integers, the same length as \code{g}, which
//' describes the reverse sort order on the observations, first by group,
//' then by place within the group. That is, the element \code{idx[0]} is the
//' index of the last place finisher in the group \code{g[0]}; then
//' \code{idx[1]} is the index of the next-to-last place finisher in
//' \code{g[1]} (assuming it equals \code{g[0]}), and so on.
//' The \code{idx} shall be zero based.
//' @param eta  a vector of the odds.
//' Must be the same length as \code{g}.
//' @param wt   an optional vector of non-negative elements, the same length
//' as \code{g}, giving the observation level weights. We then compute a 
//' weighted log likelihood, where the weights are in each summand.
//' The weights should probably have mean 1, but that's just, like, my
//' opinion, man. Negative weights throw an error. Note that the weights
//' for the last place in each group have no effect on the computation.
//' @param deleta  an optional matrix whose row count equals the number of elements
//' of \code{g}, \code{idx}, and \code{eta}. The rows of \code{deleta} are the derivatives of
//' eta with respect to some \eqn{z}. This is used to then maximize
//' likelihood with respect to \eqn{z}.
//'
//' @note
//' The likelihood function does not yet support ties.
//'
//' @return The log likelihood. If \code{deleta} is given, we add an attribute
//' to the scalar number, called \code{gradient} giving the derivative.
//' For the Henery model we also include a term of \code{gradgamma} which is
//' the gradient of the log likelihood with respect to the gamma vector.
//' @examples
//' # a garbage example
//' set.seed(12345)
//' g <- as.integer(sort(ceiling(20 * runif(100))))
//' idx <- as.integer(rev(1:length(g)) - 1L)
//' eta <- rnorm(length(g))
//' foo <- harsmlik(g=g,idx=idx,eta=eta,deleta=NULL)
//'
//' # an example with a real idx
//' nfeat <- 5
//' set.seed(1234)
//' g <- ceiling(seq(0.1,1000,by=0.1))
//' X <- matrix(rnorm(length(g) * nfeat),ncol=nfeat)
//' beta <- rnorm(nfeat)
//' eta <- X %*% beta
//' y <- rsm(eta,g)
//' 
//' idx <- order(g,y,decreasing=TRUE) - 1
//' foores <- harsmlik(g,idx,eta,deleta=X)
//'
//' # now reweight them
//' wt <- runif(length(g))
//' wt <- wt / mean(wt)   # mean 1 is recommended
//' foores <- harsmlik(g,idx,eta,wt=wt)
//' 
//' # try hensmlik 
//' foores <- hensmlik(g,idx,eta,gamma=c(0.9,0.8,1),wt=wt)
//'
//' # check the value of the gradient by numerical approximation
//' \donttest{
//'  nfeat <- 8
//'  set.seed(321)
//'  g <- ceiling(seq(0.1,1000,by=0.1))
//'  X <- matrix(rnorm(length(g) * nfeat),ncol=nfeat)
//'  beta <- rnorm(nfeat)
//'  eta <- X %*% beta
//'  y <- rsm(eta,g)
//'  
//'  idx <- order(g,y,decreasing=TRUE) - 1
//'  if (require(numDeriv)) {
//'  
//'  	fastval <- attr(harsmlik(g,idx,eta,deleta=X),'gradient')
//'  	numap <- grad(function(beta,g,idx,X) { 
//'  			 eta <- X %*% beta
//'  			 as.numeric(harsmlik(g,idx,eta))
//'  			 },
//'  			 x=beta,g=g,idx=idx,X=X)
//'  	rbind(fastval,numap)
//'  }
//' }
//' @template etc
//' @template ref-harville
//' @template note-weights
//' @rdname smlik 
//' @name smlik 
//' @export
// [[Rcpp::export]]
NumericVector harsmlik(IntegerVector g,
											 IntegerVector idx,
											 NumericVector eta,
											 Rcpp::Nullable<NumericVector> wt = R_NilValue,
											 Rcpp::Nullable<NumericMatrix> deleta = R_NilValue) {
	int gel = g.length();
	if (gel != idx.length()) { stop("length mismatch: g and idx"); }
	if (gel != eta.length()) { stop("length mismatch: g and eta"); }
	if (all( idx > 0 ).is_true()) { stop("idx should be zero indexed, not 1 indexed"); }
	if ((any( idx < 0 ).is_true()) or (any( idx >= gel ).is_true())) { stop("idx out of bounds. should be zero indexed, not 1 indexed"); }

	bool has_deriv = deleta.isNotNull();
	int gradl = 1;
	NumericMatrix del_eta;
	if (has_deriv) {
		del_eta = deleta.get();
		if (gel != del_eta.nrow()) { stop("size mismatch: g and rows of deleta"); }
		gradl = del_eta.ncol();
	}
	NumericVector gradi(gradl);
	NumericVector sumgrad(gradl);

	bool has_wt = wt.isNotNull();
	NumericVector wt_val;
	if (has_wt) {
		wt_val = wt.get();
		if (gel != wt_val.length()) { stop("size mismatch: g and wt"); }
	}

	NumericVector mu = exp(eta);  // Rcpp is pretty sweet.
	int jjj = idx[0];
	int prevgrp = g[jjj];
	int thisgrp;
	double loglik = 0;
	double summu = mu[jjj];
	double thiswt;

	for (int iii=1;iii < idx.length();iii++) {
		jjj = idx[iii];
		if (has_wt) {
			thiswt = wt_val[jjj];
			if (thiswt < 0) { stop("negative weight encountered; try again"); }
		}
		thisgrp = g[jjj];
		if (thisgrp == prevgrp) {
			summu += mu[jjj];
			if (has_wt) {
				loglik += thiswt * (eta[jjj] - log(summu));
			} else {
				loglik += (eta[jjj] - log(summu));
			}
			if (has_deriv) { 
				sumgrad += mu[jjj] * del_eta(jjj, _);
				if (has_wt) {
					gradi += thiswt * (del_eta(jjj,_) - ((1.0/summu) * sumgrad));
				} else {
					gradi += (del_eta(jjj,_) - ((1.0/summu) * sumgrad));
				}
			}
		} else {
			// new group created.
			prevgrp = thisgrp;
			summu = mu[jjj];
			if (has_deriv) { sumgrad = summu * del_eta(jjj, _); }
			// this is odd b/c it seems to ignore the weight on the last place case.
			// which is truly odd.
		}
	}

	NumericVector out = NumericVector::create(loglik);
	if (has_deriv) { out.attr("gradient") = gradi; }

	return out;
}

int grpsize(IntegerVector g,IntegerVector idx,int startfrom,int maxl) {
	int prevgrp = g[idx[startfrom]];
	int thisgrp;
	int gsize=1;
	for (int iii=startfrom+1;iii < maxl;iii++) {
		thisgrp = g[idx[iii]];
		if (thisgrp == prevgrp) { ++gsize; } else { return gsize; }
	}
	return gsize;
}

//' @param gamma a vector of the gamma parameters. It is assumed that the
//' first element is \eqn{\gamma_2}, while the last element is applied
//' to all higher order tie breaks.
//' @template ref-henery
//' @rdname smlik 
//' @export
// [[Rcpp::export]]
NumericVector hensmlik(IntegerVector g,
											 IntegerVector idx,
											 NumericVector eta,
											 NumericVector gamma,
											 Rcpp::Nullable<NumericVector> wt = R_NilValue,
											 Rcpp::Nullable<NumericMatrix> deleta = R_NilValue) {
	int gel = g.length();
	if (gel != idx.length()) { stop("length mismatch: g and idx"); }
	if (gel != eta.length()) { stop("length mismatch: g and eta"); }
	if (all( idx > 0 ).is_true()) { stop("idx should be zero indexed, not 1 indexed"); }
	if ((any( idx < 0 ).is_true()) or (any( idx >= gel ).is_true())) { stop("idx out of bounds. should be zero indexed, not 1 indexed"); }
	// dunno.
	//if (any( gamma < 0 ).is_true()) { stop("gammas should be positive?"); }

	int ngamma = gamma.length();

	bool has_deriv = deleta.isNotNull();
	int gradl = 1;
	NumericMatrix del_eta;
	if (has_deriv) {
		del_eta = deleta.get();
		if (gel != del_eta.nrow()) { stop("size mismatch: g and rows of deleta"); }
		gradl = del_eta.ncol();
	}
	NumericVector gradi(gradl, 0.0);
	NumericVector sumgrad(gradl, 0.0);
	NumericVector dgamma(ngamma, 0.0);

	bool has_wt = wt.isNotNull();
	NumericVector wt_val;
	if (has_wt) {
		wt_val = wt.get();
		if (gel != wt_val.length()) { stop("size mismatch: g and wt"); }
	}

	NumericVector mu = exp(eta);  // Rcpp is pretty sweet.
	int jjj = idx[0];
	int prevgrp = g[jjj];
	int thisgrp;
	double loglik = 0;
	double summu = mu[jjj];
	double thiswt = 1.0;
	// the index, in iii space, of the last place entry in this group
	int botidx = 0;
	// the total # of entries in this group
	int fldsize = grpsize(g,idx,0,gel);
	// the place that the current entry took in this group
	int ford = fldsize;
	// the index, in gammas, for this guy's gamma. will be negative when ford=1
	int gidx;
	int ppp;
	// my gamma; take care here;
	double mygamma;
	double thispow;
	double sumdg;
	//Rcout << "field size is " << fldsize << std::endl;

	for (int iii=1;iii < idx.length();iii++) {
		jjj = idx[iii];
		if (has_wt) {
			thiswt = wt_val[jjj];
			if (thiswt < 0) { stop("negative weight encountered; try again"); }
		}
		thisgrp = g[jjj];
		if (thisgrp == prevgrp) {
			--ford;
			gidx = MIN(ford-2,ngamma-1);
			if (gidx < 0) { mygamma = 1.0; } else { mygamma = gamma[gidx]; }
			
			// init for sum
			summu = 0.0;
			if (has_deriv) { 
				std::fill(sumgrad.begin(), sumgrad.end(), 0.0);
				// for (int lll=0;lll < gradl;lll++) { sumgrad[lll] = 0.0; }
				sumdg = 0.0;
			}
			for (int kkk=botidx;kkk <= iii;kkk++) {
				ppp = idx[kkk];
				thispow = pow(mu[ppp],mygamma);
				summu += thispow;
				if (has_deriv) { 
					sumgrad += thispow * mygamma * del_eta(ppp, _);
					sumdg += thispow * eta[ppp];
				}
			}
			loglik += thiswt * (mygamma * eta[jjj] - log(summu));
			if (has_deriv) { 
				gradi += thiswt * ((mygamma * del_eta(jjj,_)) - ((1.0/summu) * sumgrad));
				if (gidx >= 0) { 
					dgamma[gidx] += thiswt * (eta[jjj] - (sumdg / summu)); 
				}
			}
		} else {
			// new group created.
			prevgrp = thisgrp;
			botidx = iii;
			fldsize = grpsize(g,idx,iii,gel);
			ford = fldsize;
		}
	}

	NumericVector out = NumericVector::create(loglik);
	if (has_deriv) { 
		out.attr("gradient") = gradi; 
		out.attr("gradgamma") = dgamma;
	}
	return out;
}

//for vim modeline: (do not edit)
// vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
