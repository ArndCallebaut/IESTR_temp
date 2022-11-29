#ifndef MODSTRING_H
#define MODSTRING_H

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


// [[Rcpp::depends(RcppEigen)]]

using namespace std;
using namespace Rcpp;
using namespace Eigen;

int index_random_choice_non_uniform (Rcpp::NumericVector ununiform_probabilities);

NumericVector eval_probabilityVector (Eigen::SparseVector<double> probabilityVector, int threshold1);

NumericVector eval_probabilityVector_adding (Eigen::SparseVector<double> probabilityVector, Eigen::SparseVector<double> probabilityVector2,int threshold1);

int randomfunc3(int j);

NumericVector generate_permutation3(int permutation_size, int total_size);

NumericVector generate_permutation4(int permutation_size, Rcpp::NumericVector ununiform_probabilities0);

Eigen::SparseMatrix<double> proba_matrix_mult(Eigen::SparseMatrix<double> A,Eigen::SparseMatrix<double> B);

Eigen::SparseMatrix<double> proba_matrix_mult2(Eigen::SparseMatrix<double> A,Eigen::SparseMatrix<double> B);

Eigen::SparseMatrix<double> proba_matrix_mult3(Eigen::SparseMatrix<double> A,Eigen::SparseMatrix<double> B);


#endif