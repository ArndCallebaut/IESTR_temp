#include <RcppEigen.h>
#include <iostream>
#include <vector>
#include <Rcpp.h>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
#include <bits/stdc++.h> //For 2D vector and sort() function
#include <bits/stdc++.h>
#include <iostream>
#include <stdlib.h> 
#include <string>
#include <cstdlib>
#include "rcppeigen_hello_world.h";

// [[Rcpp::depends(RcppEigen)]]

using namespace std;
using namespace Rcpp;
using namespace Eigen;
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// [[Rcpp::export]]
NumericVector timesTwo2(NumericVector x) {
  return x * 2;
}



// [[Rcpp::export]]
NumericVector timesFour(NumericVector x) {
  return 2 * timesTwo(x);
}

