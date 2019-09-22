
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
 
  random generation under the Henery Model.

  Created: 2018.11.23
  Copyright: Steven E. Pav, 2018-2019
  Author: Steven E. Pav <shabbychef@gmail.com>
  Comments: Steven E. Pav
*/

#ifndef __DEF_HENERY__
#define __DEF_HENERY__

#include <math.h>
#include <string.h>

#endif /* __DEF_HENERY__ */

#include <Rcpp.h>
using namespace Rcpp;

// max_el is a hack b/c we will use the same vector over again.
int pick_winner(NumericVector probs,int max_el) {
	double apick = R::runif(0.0,1.0);
	if (probs.length() < max_el) {
		max_el = probs.length();
	} 
	if ((apick < 0.0) or (apick > 1.0)) { stop("runif acting weird."); }
	double cumprob=probs[0];
	int retval = 0;
	while (cumprob < apick) {
		retval++;
		if (retval < max_el) {
			cumprob += probs[retval];
		} else {
			break;
		}
	}
	return retval;
}

//' @title
//' Random generation under the Henery (or Harville) softmax model.
//'
//' @description
//'
//' Given base probabilities, and Henery gamma coefficients,
//' performs random generation, using R's built in rand seed,
//' of the final outcome of a race for each participant.
//'
//' @details
//' 
//' Given vectors \eqn{\mu} and \eqn{\gamma}, first computes
//' \deqn{\pi_{1,i} = \frac{\mu_i^{\gamma_1}}{\sum_j \mu_j^{\gamma_1}},}
//' then assigns a \eqn{1} to participant \eqn{i} with probability
//' \eqn{\pi_{1,i}}. The \sQuote{winning} participant is then removed
//' from consideration, and the process is repeated using the remaining
//' \eqn{\mu} and \eqn{\gamma} vectors.
//'
//' Typically one has that \eqn{\mu_i = \exp{\eta_i}}, for some
//' \sQuote{odds}, \eqn{\eta_i}.
//'
//' When the \eqn{\gamma} are all one, you recover the Harville softmax
//' model.
//'
//' @param mu  a vector of the probabilities of taking first place.
//' @param gamma  a vector of the gamma coefficients. Should have length
//' one less than \code{mu}, but if longer the unused elements are ignored.
//' If shorter, we reserve the right to either throw an error or extend out
//' the last gamma element. If not given, the coefficients are assumed
//' to be all one, which is the Harville model.
//' @return A vector, of the same length as the probabilities, giving
//' the entry of each horse. Note that the expected value of this
//' returned thing makes sense, it is \emph{not} the finished rank ordering
//' of a race.
//'
//' @seealso \code{\link{rsm}}
//' @template etc
//' @keywords probability
//' @rdname rhenery
//' @export
// [[Rcpp::export]]
IntegerVector rhenery(NumericVector mu,Rcpp::Nullable<NumericVector> gamma = R_NilValue) {
	int nel = mu.length();
	bool harville = gamma.isNull();
	
	NumericVector gammavals;
	if (!harville) {
		gammavals = gamma.get();
		// 2FIX: cannot extend out yet.
		if (nel > (1 + gammavals.length())) { stop("length mismatch: mu and gamma"); }
	}
	if (any( mu < 0 ).is_true()) { stop("mu out of bounds. must be non-negative."); }
	double summu = sum(mu);
	if (summu < 0) { stop("negative mu makes no sense, or you have NA probs."); }
	int max_el = nel;
	IntegerVector tmpord(nel);
	NumericVector tmpmu(nel);
	for (int iii=0;iii < nel;iii++) {
		tmpord[iii] = iii;
		tmpmu[iii] = mu[iii] / summu;
	}
	int nextwin;
	double sumprob;
	IntegerVector retval(nel);
	for (int iii=0;iii < (nel-1);iii++) {
		nextwin = pick_winner(tmpmu,max_el);
		//Rcout << "winner was " << nextwin << std::endl;
		//assign in
		retval[tmpord[nextwin]] = iii + 1;

		for (int jjj=nextwin;jjj < max_el - 1;jjj++) {
			tmpord[jjj] = tmpord[jjj+1];
		}
		max_el--;
		sumprob = 0;
		for (int jjj=0;jjj < max_el;jjj++) {
			if (harville) {  // harville
				tmpmu[jjj] = mu[tmpord[jjj]];
			} else {  // henery
				tmpmu[jjj] = pow(mu[tmpord[jjj]],gammavals[iii]);
			}
			sumprob += tmpmu[jjj];
		}
		for (int jjj=0;jjj < max_el;jjj++) {
			tmpmu[jjj] = tmpmu[jjj] / sumprob;
		}
	}
	// now the last guy
	retval[tmpord[0]] = nel;
	return retval;
}

//for vim modeline: (do not edit)
// vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
