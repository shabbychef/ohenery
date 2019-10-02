
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

 
	Compute the expected rank given probabilities.
	Hard to do as it is a big multivariate integration.

  Created: 2018.09.26
  Copyright: Steven E. Pav, 2018-2019
  Author: Steven E. Pav <shabbychef@gmail.com>
  Comments: Steven E. Pav
*/

#ifndef __DEF_ERANK__
#define __DEF_ERANK__

#include <math.h>
#include <string.h>

#define MAX(_a_,_b_) ((_a_>_b_)? (_a_):(_b_))
#define MIN(_a_,_b_) ((_a_<_b_)? (_a_):(_b_))
#define ODDS(_p_) ((_p_) / (1.0 - _p_))

// this is 1 + p / (1-p) = 1 / (1-p)
#define ODDSPLUS(_p_) ((1.0) / (1.0 - _p_))
#define ERANK_THREE(_p_,_opa_,_opb_) (3.0 - _p_ * (_opa_ + _opb_))

#endif /* __DEF_ERANK__ */

#include <Rcpp.h>
using namespace Rcpp;

//' @title
//' Expected rank under the Harville model.
//'
//' @description
//'
//' Compute the expected rank of a bunch of entries based on their probability
//' of winning under the Harville model.
//'
//' @details
//'
//' Given the vector \eqn{\mu}, we compute the expected rank of each
//' entry, under the Harville model, where tail probabilities of winning
//' remain proportional.  
//'
//' Under the Harville model, the probability that the \eqn{i}th element
//' is assigned value 1 is 
//' \deqn{\pi_{1,i} = \frac{\mu_i}{\sum_j \mu_j}.}
//' Once an element has been assigned a 1, the Harville procedure 
//' removes it from the set and iterates.
//' If there are \eqn{k} elements in \eqn{\mu}, then the \eqn{i}th element
//' can be assigned any place between \eqn{1} and \eqn{k}. This
//' function computes the expected value of that random variable.
//'
//' While a naive implementation of this function would take
//' time factorial in \eqn{k}, this implementation takes time quadratic
//' in \eqn{k}, since it can be shown that the expected rank of the \eqn{i}th
//' element takes value
//' \deqn{e_i = k + \frac{1}{2} - \sum_j \frac{\mu_i}{\mu_i + \mu_j}.}
//'
//' @param mu a vector giving the probabilities. Should sum to one.
//'
//' @note
//' we should have the sum of ranks equal to the sum of \code{1:length(mu)}.
//'
//' @return The expected ranks, a vector.
//' @examples
//' # a garbage example
//' set.seed(12345)
//' mus <- runif(12)
//' mus <- mus / sum(mus)
//' erank(mus)
//'
//' # confirm the expected rank via simulation
//' set.seed(123)
//' mus <- runif(6,min=0,max=2)
//' mus <- mus / sum(mus)
//' set.seed(101)
//' emp <- rowMeans(replicate(200,rhenery(mu=mus,gamma=rep(1,length(mus)-1)))) 
//' (emp - erank(mus)) / emp
//'
//' \donttest{
//' if (require(microbenchmark)) {
//'   p10 <- 1:10 / sum(1:10)
//'   p16 <- 1:16 / sum(1:16)
//'   p24 <- 1:24 / sum(1:24)
//'   microbenchmark(erank(p10), erank(p16), erank(p24))
//' }
//' }
//' @template etc
//' @rdname erank
//' @export
// [[Rcpp::export]]
NumericVector erank(NumericVector mu) {
	int nel = mu.length();
	if (nel < 2) { return mu; }
	if (nel == 2) { return 2 - mu; }
	if (nel == 3) { //FOLDUP
		double O_1 = ODDSPLUS(mu[0]);
		double O_2 = ODDSPLUS(mu[1]);
		double O_3 = ODDSPLUS(mu[2]);
		return NumericVector::create( ERANK_THREE(mu[0],O_2,O_3), 
																	ERANK_THREE(mu[1],O_1,O_3), 
																	ERANK_THREE(mu[2],O_1,O_2) );
	}//UNFOLD

	NumericVector output(nel,double(nel));
	double nextv, mui, muj;
	for (int iii=0;iii < nel;iii++) {
		mui = mu[iii];
		for (int jjj=0;jjj < iii;jjj++) {
			muj = mu[jjj];
			nextv = (mui + muj);
			output[iii] -= mui / nextv;
			output[jjj] -= muj / nextv;
		}
	}
	return output;
}

//for vim modeline: (do not edit)
// vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
