


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
#include <random>
#include <chrono>
typedef Eigen::SparseVector<double> SpVec;
typedef SpVec::InnerIterator InIterVec;

// [[Rcpp::depends(RcppEigen)]]

using namespace std;
using namespace Rcpp;
using namespace Eigen;

//[[Rcpp::export]]
int index_random_choice_non_uniform (NumericVector ununiform_probabilities0){

  Rcpp::NumericVector ununiform_probabilities=ununiform_probabilities0/sum(ununiform_probabilities0);
  Rcpp::NumericVector cumulative_probabilities(ununiform_probabilities.length());
  cumulative_probabilities(0) = ununiform_probabilities(0);
  
  for (int h=1; h<ununiform_probabilities.length() ; h++){
    cumulative_probabilities(h) = ununiform_probabilities(h)+cumulative_probabilities(h-1);  
  }
  
  //std::cout << "LETS TRY THIS"<<cumulative_probabilities << std::endl;
  //if( cumulative_probabilities((cumulative_probabilities.length()-1)) < 1){
    
    //std::cout << "ALERTE ALERTE !"<<cumulative_probabilities(cumulative_probabilities.length()-1) << std::endl;
  //}
  
  double randy = ((double) rand() / (RAND_MAX)) ;
  
  for (int h = 0 ; h < ununiform_probabilities.length() ; h++){
    if (randy<cumulative_probabilities(h)){
      return(int(h));
    }
  }
  //std::cout << "YOU FUCKING DONKEY" << std::endl;
  return(ununiform_probabilities.length()-1);
}




NumericVector eval_probabilityVector (Eigen::SparseVector<double> probabilityVector, int threshold1){
  
  int const threshold = threshold1;
  
  Rcpp::NumericVector nScores(threshold+1);
  
  int h;
  for (h=0; h<=threshold; ++h){
    nScores(h) = 0;
  }
  nScores(0)=1;
  
  for (Eigen::SparseVector<double>::InnerIterator it(probabilityVector); it; ++it)
  {
    
    nScores(threshold) =  it.value() * nScores(threshold-1) +  nScores(threshold);
    for (h=(threshold-1); h>=1; --h){
      nScores(h) =  it.value() * nScores(h-1) +  (1-it.value()) * nScores(h);
    }
    nScores[0] = ((1-it.value()) * nScores[0]);
  }
  return(nScores) ;
}



NumericVector eval_probabilityVector_adding (Eigen::SparseVector<double> probabilityVector, Eigen::SparseVector<double> probabilityVector2,int threshold1){
  
  Eigen::SparseVector<double> cury = probabilityVector;
  
  for (InIterVec i_(probabilityVector2); i_; ++i_){
    //Rcpp::Rcout << " i=" << i_.index() << " value=" << i_.value() << std::endl;
    cury.coeffRef(i_.index()) = cury.coeffRef(i_.index()) + i_.value() - cury.coeffRef(i_.index()) * i_.value()  ;
  }
  
  return(eval_probabilityVector(cury,threshold1)) ;
}

int randomfunc3(int j){
  return rand() % j;
}

NumericVector generate_permutation3(int permutation_size, int total_size){
  
  Rcpp::NumericVector v(total_size);
  for (int i = 0; i < total_size; i++){
    v(i) = i ;
  }
  //std::cout << "OK here it's done." << std::endl;
  
  Rcpp::NumericVector result(permutation_size);
  int indice_max = total_size;
  int alea;
  int saver;
  for (int i = 0; i < permutation_size; i++){
    alea = randomfunc3(indice_max);
    //std::cout << "OK here it's done..."<<v(alea) << std::endl;
    result(i) = v(alea);
    saver = v(indice_max-1);
    v(indice_max-1) = v(alea);
    v(alea) = saver;
    indice_max --;
  }
  return(result);
}

//[[Rcpp::export]]
NumericVector generate_permutation4(int permutation_size, Rcpp::NumericVector ununiform_probabilities0){
  
  //std::cout << "Checkup 1" << std::endl;
  
  Rcpp::NumericVector ununiform_probabilities = ununiform_probabilities0;
  int total_size = ununiform_probabilities.length();
  
  Rcpp::NumericVector v(total_size);
  
  for (int i = 0; i < total_size; i++){
    v(i) = i ;
  }
  //std::cout << "OK here it's done." << std::endl;
  //std::cout << "Checkup 2" << std::endl;
  Rcpp::NumericVector result(permutation_size);
  int indice_max = total_size;
  int alea;
  int saver;
  for (int i = 0; i < permutation_size; i++){
    
    //std::cout << "Checkup 1" << std::endl;
    //std::cout << ununiform_probabilities.length() << std::endl;
    alea = index_random_choice_non_uniform(ununiform_probabilities);
    //std::cout << "Checkup 3.2" << std::endl;
    //std::cout << "OK here it's done..."<<v(alea) << std::endl;
    //std::cout << "Checkup 1.1 - "<< alea << " - ok on TS - "<< total_size << std::endl;
    result(i) = v(alea);
    //std::cout << "Checkup 1.2" << std::endl;
    saver = v(indice_max-1);
    //std::cout << "Checkup 2" << std::endl;
    v(indice_max-1) = v(alea);
    v(alea) = saver;
    indice_max --;
    
    //std::cout << "Checkup 3" << std::endl;
    
    ununiform_probabilities(v(alea)) = ununiform_probabilities(indice_max-1);
    ununiform_probabilities(indice_max-1) = 0;
    
    //ununiform_probabilities.erase((alea));
  }
  return(result);
}




Eigen::SparseMatrix<double> proba_matrix_mult(Eigen::SparseMatrix<double> A,Eigen::SparseMatrix<double> B){
  // probabilistic matrix multiplication of A*B
  int n = A.rows();
  Eigen::SparseMatrix<double>result_matrix(n,n);
  double val1;
  int ind1;
  int ind2;
  double val2;
  
  //Rcout << "A..."<<A.nonZeros()<<" ------- B..."<<B.nonZeros()<<std::endl;
  
  
  for (int k=0; k<A.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(A,k); it; ++it)
    {
      val1 = it.value();
      ind1 = it.row();
      
      if(val1>0.01){
        
        for (int h=0; h<B.outerSize(); ++h)
          for (Eigen::SparseMatrix<double>::InnerIterator it2(B,h); it2; ++it2)
          {
            if (it.col()==it2.row()){
              
              val2 = it2.value()*val1;
              ind2 = it2.col();
              
              if (val2>0.001){
                result_matrix.coeffRef(ind1,ind2) = result_matrix.coeffRef(ind1,ind2) + val2 - val2* result_matrix.coeffRef(ind1,ind2);
              }
            }
          }
      }
    }
    return (result_matrix);
}

Eigen::SparseMatrix<double> proba_matrix_mult2(Eigen::SparseMatrix<double> A,Eigen::SparseMatrix<double> B){
  // probabilistic matrix multiplication of A*B
  int n = A.rows();
  Eigen::SparseMatrix<double>result_matrix(n,n);
  double val1;
  int ind;
  double val2;
  for (int k=0; k<A.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(A,k); it; ++it)
    {
      val1 = it.value();
      ind = it.col();
      for (int h=0; h<n;++h){
        val2 = B.coeffRef(ind,h);
        if(val2!=0){
          result_matrix.coeffRef(it.row(),h)+= val1 + val2 - val1*val2;
        }
      }
    }
    return (result_matrix);
}

/*
//[[Rcpp::export]]
Eigen::SparseMatrix<double> proba_matrix_mult3(Eigen::SparseMatrix<double> AA,Eigen::SparseMatrix<double> B){
  // probabilistic matrix multiplication of A*B
  Eigen::SparseMatrix<double>A = AA.transpose();
  int n = A.rows();
  Eigen::SparseMatrix<double>result_matrix(n,n);
  double val1;
  int ind1;
  int ind2;
  double val2;
  double value;
  int i;
  std::cout << "A..."<<A.nonZeros()<<" ------- B..."<<B.nonZeros()<<std::endl;

  
  for (int k=0; k<B.outerSize(); ++k){
    Eigen::SparseVector<double> local = B.col(k);
    for (int h=0; h<A.outerSize(); ++h){
      
      
      
      value = 0;
    
      for (Eigen::SparseMatrix<double>::InnerIterator it(A,h); it; ++it){
        val1 = it.value();
        ind1 = it.row();
        if(val1>0.01){
          val2 = val1 * local.coeffRef(ind1);
          if ((val2)>0.01){
            value = value + val2 - val2* value;
          }
        }
      }
      
      if (value !=0){
        result_matrix.insert(h,k) = value;
      }
      
    
    }
  }
    return (result_matrix);
}


 //[[Rcpp::export]]
 Eigen::SparseMatrix<double> proba_matrix_mult3(Eigen::SparseMatrix<double> AA,Eigen::SparseMatrix<double> B){
   // probabilistic matrix multiplication of A*B
   Eigen::SparseMatrix<double>A = AA.transpose();
   int n = A.rows();
    
   Eigen::SparseMatrix<double>result_matrix(n,n);
   double val1;
   int ind1;
   int ind2;
   double val2;
   double val3;
   int i;
   std::cout << "A..."<<A.nonZeros()<<" ------- B..."<<B.nonZeros()<<std::endl;
   
   
   for (int k=0; k<B.outerSize(); ++k)
     for (Eigen::SparseMatrix<double>::InnerIterator it(B,k); it; ++it)
     {
       val1 = it.value();
       ind1 = it.row();
       if(val1>0.01){
         for (Eigen::SparseMatrix<double>::InnerIterator it2(A,k); it2; ++it2)
         {
           val2 = it2.value()*val1;
           ind2 = it2.row();
           
           if ((val2)>0.01){
             val3 = result_matrix.coeffRef(ind1,ind2);
             //std::cout << "x..."<<ind1<<" ------- y..."<<ind2<<" ---- val="<<result_matrix.coeffRef(ind1,ind2) + val2 - val2* result_matrix.coeffRef(ind1,ind2)<<std::endl;
             result_matrix.coeffRef(ind1,ind2) = val3+ val2 - val2*val3;
           }
         }
       }
     }
   return (result_matrix);
 }
 */

//[[Rcpp::export]]
Eigen::SparseMatrix<double> proba_matrix_mult3(Eigen::SparseMatrix<double> A,Eigen::SparseMatrix<double> B){
  
  // probabilistic matrix multiplication of A*B
  //Eigen::SparseMatrix<double>A = AA.transpose();
  int n = A.rows();
  Eigen::SparseVector<double> local2(n);
  Eigen::SparseMatrix<double>result_matrix(n,n);
  double val1;
  int ind1;
  int ind2;
  double val2;
  double val3;
  int i;
  std::cout << "A..."<<A.nonZeros()<<" ------- B..."<<B.nonZeros()<<std::endl;
  
  
  for (int k=0; k<B.outerSize(); ++k){
    Eigen::SparseVector<double> local(n);
        for (Eigen::SparseMatrix<double>::InnerIterator it(B,k); it; ++it)
          {
         
         
            val1 = it.value();
            ind1 = it.row();
            //std::cout << "\nval of B calculated - "<<val1<<" k="<<k<<" h="<<ind1<<std::endl;
            if(val1>0.01){
              local2 = val1 * A.col(ind1);
              //std::cout << "A cols resulting of this - "<<local2<<std::endl;
              for (Eigen::SparseVector<double>::InnerIterator it(local2); it; ++it)
                {
                //std::cout << "Now index="<<it.index()<<std::endl;
                local.coeffRef(it.index()) = local.coeffRef(it.index()) +it.value() - local.coeffRef(it.index()) *it.value() ;
                }
              //std::cout << "local then - "<<local<<std::endl;
              //local = local +local2 - local *local2;
              
              //std::cout << "Now locval="<<local<<std::endl;
              
              
            }
          }
        
        for (InIterVec i_(local); i_; ++i_){
          if (i_.value()>0.01){
            result_matrix.insert(i_.index(),k) = i_.value();
          }
        }
        //std::cout << "result_matrix then - "<<result_matrix<<std::endl;
        //result_matrix.col(k) = local;
    }
  result_matrix.pruned(0.01);
  result_matrix.makeCompressed();
  return (result_matrix);
}
/*
//[[Rcpp::export]]
Eigen::SparseMatrix<double> proba_matrix_mult4(Eigen::SparseMatrix<double> A,Eigen::SparseMatrix<double> B){
  
  
  // probabilistic matrix multiplication of A*B
  int n = A.rows();
  
  Eigen::SparseMatrix<double>AA = A.transpose();
  
  Eigen::SparseVector<double> local2(n);
  Eigen::SparseMatrix<double>result_matrix(n,n);
  double val1;
  int ind1;
  int ind2;
  double val2;
  double val3;
  int i;
  
  Rcpp::NumericMatrix tAA(AA.nonZeros(),3);
  Rcpp::NumericMatrix tB(B.nonZeros(),3);
  
  i = 0
  for (int k=0; k<AA.outerSize(); ++k){
    for (Eigen::SparseMatrix<double>::InnerIterator it(AA,k); it; ++it){
      tAA(i,0) = it.row();
      tAA(i,1) = it.col();
      tAA(i,2) = it.value();
    }
    
  i = 0
  for (int k=0; k<B.outerSize(); ++k){
    for (Eigen::SparseMatrix<double>::InnerIterator it(B,k); it; ++it){
      tB(i,0) = it.row();
      tB(i,1) = it.col();
      tB(i,2) = it.value();
    }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  std::cout << "A..."<<A.nonZeros()<<" ------- B..."<<B.nonZeros()<<std::endl;
  
  for (int k=0; k<B.outerSize(); ++k){
    Eigen::SparseVector<double> local(n);
    for (Eigen::SparseMatrix<double>::InnerIterator it(B,k); it; ++it)
    {
      val1 = it.value();
      ind1 = it.row();
      
      if(val1>0.01){
        local2 = val1 * A.col(ind1);
        local = local +local2 - local *local2;
      }
    }
    
    for (InIterVec i_(local); i_; ++i_){
      if (i_.value()>0.05){
        result_matrix.insert(i_.index(),k) = i_.value();
      }
    }
    //result_matrix.col(k) = local;
  }
  result_matrix.pruned(0.01);
  result_matrix.makeCompressed();
  return (result_matrix);
}







*/











