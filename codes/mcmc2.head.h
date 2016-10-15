#ifndef MCMC2_HEAD_H_
#define MCMC2_HEAD_H_
#include "mcmc.head.h"
#include <RcppArmadillo.h>
#include<Rcpp.h>
#include<math.h>

//[[Rcpp::depends(RcppArmadillo)]]
#include<R.h>
using namespace Rcpp;
using namespace arma;


arma::mat forbackwardMCMC(arma::vec dt, arma::vec tau, arma::vec c,arma::vec w,int iter, int burnin, int thin);

#endif /* MCMC2_HEAD_H_ */
