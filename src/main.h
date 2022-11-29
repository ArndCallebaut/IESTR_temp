#ifndef MODSTRING_H
#define MODSTRING_H

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <iostream>
#include <vector>
//#include <Rcpp.h>
#include <numeric>    
#include <algorithm>    
#include <bits/stdc++.h> 
#include <bits/stdc++.h>
#include <iostream>
#include <stdlib.h> 
#include <string>
#include <cstdlib>
#include "maths_tools.h";
#include "local_optimising_planting_choice.h";

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace std;
using namespace Rcpp;
using namespace Eigen;
using namespace arma;

Eigen::SparseMatrix<double> sqrt_(Eigen::SparseMatrix<double> X);

Rcpp::List optgenam3 (Eigen::SparseMatrix<double> currentPresenceMatrix,
                      std::list<Eigen::SparseMatrix<double>> consecutiveSuitabilityMatrix,
                      Eigen::SparseMatrix<double> costMatrix,
                      Rcpp::NumericMatrix migrationKernel,
                      int threshold,
                      double confidence,
                      int npop,
                      int nsur,
                      int ngen);

#endif