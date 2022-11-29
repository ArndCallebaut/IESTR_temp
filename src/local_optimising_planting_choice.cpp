


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

int optimise_planting_choice4(SparseMatrix<double> viablesValues2, int threshold, double confidence, Rcpp::NumericMatrix viablesTriplets, Rcpp::NumericMatrix population, int pop, int nbtoplant,  Eigen::SparseVector<double>& current){
  
  //We want to eliminate a site
  
  int n = population.cols();
  if (n==0){
    return(-1);
  }
  
  int elimination_index = -1;
  double gain = 0;
  
  // for each in the set of choice
  for (int s=0 ; s < n ; ++s){
    if(population(pop,s)!=-1){
      
      Eigen::SparseVector<double> cur = current;
      
      for (Eigen::SparseMatrix<double>::InnerIterator it(viablesValues2,population(pop,s)); it; ++it) {
        cur.coeffRef(it.index()-1) = (cur.coeffRef(it.index()-1)-it.value())/(1-it.value());
      }
      
      Rcpp::NumericVector vec = eval_probabilityVector(cur,threshold);
      //std::cout << "HOW THE TURN TABLE !!!"<< vec <<  std::endl;
      if (vec(threshold)>=confidence){
        //std::cout << "WELL WELL WELL................... gain="<<gain<<" val+="<<viablesTriplets(population(pop,s),4) << std::endl;
        if(gain<viablesTriplets(population(pop,s),3)){
          //std::cout << "HOW THE TURN TABLE..................." << std::endl;
          elimination_index = s;
          //Rcout << "H" << std::endl;
          gain = viablesTriplets(population(pop,s),3);
          //Rcout << "O" << std::endl;
        }
      }
    }
  }
  
  //
  
  if (elimination_index ==-1){
    //Rcout << "AND NOW THE TRUE TEST 1" << std::endl;
    return(-1);
  }
  
  if (elimination_index >=0){
    //Rcout << "AND NOW THE TRUE TEST 2" << std::endl;
    
    //std::cout << "HOW THE TURN TABLE..................." << elimination_index << std::endl;
    
    //Eigen::SparseVector<double> cur = current;
    
    for (Eigen::SparseMatrix<double>::InnerIterator it(viablesValues2,population(pop,elimination_index)); it; ++it) {
      current.coeffRef(it.index()-1) = (current.coeffRef(it.index()-1)-it.value())/(1-it.value());
    }
    
    population(pop,elimination_index) = -1;
    
    //std::cout << "HOW THE TURN TABLE TWICE..................." << elimination_index << std::endl;
    
    //optimise_planting_choice4(viablesValues2,threshold,confidence,viablesTriplets,population,pop,nbtoplant,cur);
    //optimise_planting_choice2 (probabilityMatrix,threshold,confidence,suitable,population,pop,nbtoplant);
    return(elimination_index);
  }
  
  return(-1);
}

void voidistheonlyoption(){
  return;
}



void optimise_planting_choice5 (Eigen::SparseMatrix<double> viablesValues2, int threshold, double confidence, Rcpp::NumericMatrix viablesTriplets, Rcpp::NumericMatrix population, int pop, int nbtoplant,  Eigen::SparseVector<double> current){

  // We want to add a site

  int n1 = viablesValues2.rows();
  if (n1==0){
    return;
  }
  int n2 = population.cols();

  int added_index = -1;
  double cost = 1000000;
  int s;
  for (s=0 ; s < n2 ; ++s){
    if(population(pop,s)==-1){

      for(int m=0; m<n1;++m){
        Rcpp::NumericVector vec = eval_probabilityVector_adding(current,viablesValues2.col(m),threshold);
        //Rcout << "HOW THE TURN TABLE !!!"<< vec(threshold) << " vs "<< confidence << std::endl;
        if (vec(threshold)>=confidence){
          if(cost>viablesTriplets(m,4)){
            //Rcout << "HOW THE TURN TABLE..................." << std::endl;
            added_index = m;
            //Rcout << "H" << std::endl;
            cost = viablesTriplets(m,4);
            //Rcout << "O" << std::endl;
          }
        }
      }
      break;
    }
  }

  if (added_index ==-1){
    //Rcout << "AND NOW THE TRUE TEST 3" << std::endl;
    return;
  }

  if (added_index >=0){
    //Rcout << "AND NOW THE TRUE TEST 4" << std::endl;
    population(pop,s) = added_index;
    //optimise_planting_choice3 (probabilityMatrix,threshold,confidence,suitable,population,pop,nbtoplant);
    //optimise_planting_choice5(viablesValues2,threshold,confidence,viablesTriplets,population,pop,nbtoplant,current);
    return;
  }
  return;
}