#ifndef LOPC
#define LOPC

#include <RcppEigen.h>
#include <iostream>
#include <vector>
#include <Rcpp.h>
#include <numeric>      
#include <algorithm>    
#include <bits/stdc++.h> 
#include <bits/stdc++.h>
#include <iostream>
#include <stdlib.h> 
#include <string>
#include <cstdlib>
#include "maths_tools.h";

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(Matrix, RcppArmadillo)]]

using namespace std;
using namespace Rcpp;
using namespace Eigen;

int optimise_planting_choice4(SparseMatrix<double> viablesValues2, int threshold, double confidence, Rcpp::NumericMatrix viablesTriplets, Rcpp::NumericMatrix population, int pop, int nbtoplant,  Eigen::SparseVector<double>& current);

void optimise_planting_choice5 (Eigen::SparseMatrix<double> viablesValues2, int threshold, double confidence, Rcpp::NumericMatrix viablesTriplets, Rcpp::NumericMatrix population, int pop, int nbtoplant,  Eigen::SparseVector<double> current);

void voidistheonlyoption();


#endif