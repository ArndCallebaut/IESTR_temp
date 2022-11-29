#include <RcppEigen.h>
//#include <RcppArmadillo.h>
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
#include <chrono>
#include <thread>
#include "maths_tools.h";
#include "local_optimising_planting_choice.h";
#include "test.h";
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace std;
using namespace Rcpp;
using namespace Eigen;

//using namespace arma;

// [[Rcpp::export]]
Eigen::SparseMatrix<double> rcpp_global_suitable_sites(std::list<Eigen::SparseMatrix<double>> consecutiveSuitabilityMatrix){
  int nbr = consecutiveSuitabilityMatrix.back().rows();
  int nbc = consecutiveSuitabilityMatrix.back().cols();
  
  Eigen::SparseMatrix<double> globalSuitableSites(nbr,nbc);
  
  Rcpp::NumericVector vec_nbpp(consecutiveSuitabilityMatrix.size());
  int k = 0;
  for (auto const& thisSuitabilityMatrix : consecutiveSuitabilityMatrix){
    vec_nbpp(k) = thisSuitabilityMatrix.nonZeros();
    ++k;
  }
  int tnbpp = sum(vec_nbpp);
  Rcpp::NumericMatrix possibilities(tnbpp,5);
  int j=1;
  for (auto const& thisSuitabilityMatrix : consecutiveSuitabilityMatrix){
    for (k=0; k<thisSuitabilityMatrix.outerSize(); ++k)
      for (Eigen::SparseMatrix<double>::InnerIterator it(thisSuitabilityMatrix,k); it; ++it)
      {
        if(globalSuitableSites.coeffRef(it.row(),it.col())==0){
          globalSuitableSites.coeffRef(it.row(),it.col()) = double(j);
          j++;
        }
      }
    }
  globalSuitableSites.makeCompressed();
  std::cout << "In time... you  will know the tragic extends of my fallings"<<std::endl;
  return(globalSuitableSites);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_global_suitable_coordinates(Eigen::SparseMatrix<double> globalSuitableSites){
  int nbsites = globalSuitableSites.nonZeros();
  Rcpp::NumericMatrix globalSuitableCoordinates(nbsites,3);
  for (int k=0; k<globalSuitableSites.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(globalSuitableSites,k); it; ++it)
    {
      globalSuitableCoordinates(it.value()-1,0)=it.row();
      globalSuitableCoordinates(it.value()-1,1)=it.col();
      globalSuitableCoordinates(it.value()-1,2)=int(it.value())-1;
    }
  return(globalSuitableCoordinates);
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> rcpp_local_transition_matrix(Eigen::SparseMatrix<double> globalSuitableSites,Rcpp::NumericMatrix globalSuitableCoordinates,Rcpp::NumericMatrix migrationKernel){
  int nbr = globalSuitableSites.rows();
  int nbc = globalSuitableSites.cols();
  int i;
  int j;  
  int migrationRange = (migrationKernel.rows()-1)/2;
  int nbsites = globalSuitableSites.nonZeros();
  Eigen::SparseMatrix<double> localTransitionMatrix(nbsites,nbsites);
  for (int index=0;index<nbsites;++index){
    i = globalSuitableCoordinates(index,0);
    j = globalSuitableCoordinates(index,1);
    for (int h=max(0,i-migrationRange); h<min(nbr,i+migrationRange+1); ++h){
      for (int k=max(0,j-migrationRange); k<min(nbc,j+migrationRange+1); ++k){
        if (globalSuitableSites.coeffRef(h,k)!=0 && migrationKernel(h+migrationRange-i,k+migrationRange-j)!=0){
          localTransitionMatrix.insert(globalSuitableSites.coeffRef(h,k)-1,index)=migrationKernel(h+migrationRange-i,k+migrationRange-j);
        }
      }
    }
  }
  localTransitionMatrix.pruned(0.01);
  localTransitionMatrix.makeCompressed();
  return (localTransitionMatrix);
}

// [[Rcpp::export]]
std::list<Eigen::SparseMatrix<double>> rcpp_transition_matrices(std::list<Eigen::SparseMatrix<double>> consecutiveSuitabilityMatrix,Eigen::SparseMatrix<double> localTransitionMatrix,Rcpp::NumericMatrix globalSuitableCoordinates){
  std::list<Eigen::SparseMatrix<double>> transitionMatrices;
  Eigen::SparseMatrix<double> transitionMatrix;
  int i=0;
  for (auto thisSuitabilityMatrix : consecutiveSuitabilityMatrix){
    transitionMatrix = localTransitionMatrix;
    for (int k=0; k<localTransitionMatrix.outerSize(); ++k){
      for (Eigen::SparseMatrix<double>::InnerIterator it(localTransitionMatrix,k); it; ++it)
      {
        transitionMatrix.coeffRef(it.row(),it.col()) *= thisSuitabilityMatrix.coeffRef(globalSuitableCoordinates(it.row(),0),globalSuitableCoordinates(it.row(),1));
      }
    }
    transitionMatrix.prune(0.01);
    ++i;
    transitionMatrix.makeCompressed();
    transitionMatrices.push_back(transitionMatrix);
  }
  transitionMatrices.reverse();
  return(transitionMatrices);
}

// [[Rcpp::export]]
std::vector<Eigen::SparseMatrix<double>> rcpp_colonisation_matrices(std::list<Eigen::SparseMatrix<double>> transitionMatrices){
  std::vector<Eigen::SparseMatrix<double>> colonisationMatrices;
  colonisationMatrices.push_back(transitionMatrices.front());
  int k=1;
  bool first = true;
  for (auto const& elt : transitionMatrices) {
    if (first){first = false;}
    else{
      colonisationMatrices.push_back(proba_matrix_mult3(colonisationMatrices.back(),elt));
      ++k;
    }
  }
  //colonisationMatrices.pruned(0.01);
  reverse(colonisationMatrices.begin(),colonisationMatrices.end());
  return(colonisationMatrices);
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> rcpp_viable_sites(std::vector<Eigen::SparseMatrix<double>> colonisationMatrices){
  int nbsites = colonisationMatrices.front().rows();
  int nbperiod = colonisationMatrices.size();
  Eigen::SparseMatrix<double> viableSites(nbperiod,nbsites);
  bool av1;
  bool av2;
  int h;
  int k;
  double var;
  Eigen::SparseVector<int> suit1;
  Eigen::SparseVector<int> suit2;
  bool nonzero1;
  bool nonzero2;
    
  for (int i=0;i<nbsites;++i){
    Rcpp::NumericVector viability (nbperiod);
    for (int j=0;j<nbperiod;++j){
      Eigen::SparseVector<int> suit2 = colonisationMatrices[j].row(i);
      if(suit2.nonZeros()!=0){
        viability(j)=1;
      }
    }
    h=0;
    for (int ss=0; ss<(nbperiod-1); ++ss){
      h=ss;
      suit1 = colonisationMatrices[ss].row(i);
      if ( viability(h)==0 ){
        break;
      }
      for (int sss=ss+1; sss<nbperiod; ++sss){
        k=sss;
        suit2 = colonisationMatrices[sss].row(i);
        nonzero1 = false;
        nonzero2 = false;
        if (viability(k)==0 || viability(h)==0 ){
          break;
        }
        av1 = false;
        av2 = false;
        if (k>=h){
          ++k;
          break;
        }
        for (int j=0;j<nbsites;++j){
          if(suit1.coeffRef(j,i)!=0){
            nonzero1 = true;
          }
          if(suit2.coeffRef(j,i)!=0){
            nonzero2 = true;
          }
          var = suit1.coeffRef(j)-suit2.coeffRef(j);
          if(var>0){
            av1=true;
          }
          if(var<0){
            av2=true;
          }
          if(av1 && av2){
            break;
          }
        }
        if(!nonzero1){
          // si le premier site ne donne aucun resultat
          viability(h)=0;
        }
        if(!nonzero2){
          // si le second site ne donne aucun resultat
          viability(k)=0;
        }
        if (av1 && !av2){
          viability(k)=0;
        }
        if (av2 && !av1){
          viability(h)=0;
        }
        ++k;
      }
      ++h;
    }
    
    for (int j=0;j<nbperiod;++j){
      if (viability[j]==1){
        if((colonisationMatrices.at(j).col(i)).sum()>0.5){
          viableSites.coeffRef(int(j),int(i))=1;
          //viableSites2.coeffRef(int(i),int(j))=1;
        }
      }
    }
  }
  viableSites.makeCompressed();
  return(viableSites);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_viable_triplets(Eigen::SparseMatrix<double> viableSites,
                                         std::list<Eigen::SparseMatrix<double>> colonisationMatrices,
                                         Rcpp::NumericMatrix globalSuitableCoordinates,
                                         Eigen::SparseMatrix<double> globalSuitableSites,
                                         Eigen::SparseMatrix<double> costMatrix){
  int nboc = viableSites.nonZeros();
  Rcpp::NumericMatrix viablesTriplets(nboc,6);
  
  
  int t = 0;
  int index = 0;
  for (auto colma : colonisationMatrices){
    Eigen::SparseVector<double> loc = viableSites.row(t);
    for (Eigen::SparseVector<double>::InnerIterator it(loc); it; ++it){
      viablesTriplets(index,0) = globalSuitableCoordinates(it.index(),0);
      viablesTriplets(index,1) = globalSuitableCoordinates(it.index(),1);
      viablesTriplets(index,2) = t;
      viablesTriplets(index,3) = costMatrix.coeffRef(viablesTriplets(index,0),viablesTriplets(index,1));
      viablesTriplets(index,4) = globalSuitableSites.coeffRef(viablesTriplets(index,0),viablesTriplets(index,1));
      viablesTriplets(index,5) = colma.col(globalSuitableSites.coeffRef(viablesTriplets(index,0),viablesTriplets(index,1))).sum();
      ++index;
    }
    ++t;
  }
  return (viablesTriplets);
}



// [[Rcpp::export]]
Eigen::SparseMatrix<double> rcpp_viable_values(Rcpp::NumericMatrix viablesTriplets,
                                                Eigen::SparseMatrix<double> viableSites,
                                                Eigen::SparseMatrix<double> globalSuitableSites,
                                                std::vector<Eigen::SparseMatrix<double>> colonisationMatrices){
  int nbsites = globalSuitableSites.nonZeros();
  int nboc = viableSites.nonZeros();
  int t = 0;
  int index = 0;
  Eigen::SparseMatrix<double> viableSites2 = viableSites.transpose();
  Eigen::SparseMatrix<double> viablesValues(nbsites,nboc);
  Eigen::SparseVector<double> colo2;
  
  for (int k=0; k<viableSites2.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(viableSites2,k); it; ++it)
    {
      if(it.value()==0){
        continue;
      }
      int xx = viablesTriplets(index,0);
      int yy = viablesTriplets(index,1);
      double val = globalSuitableSites.coeffRef(xx,yy) ;
      colo2 = colonisationMatrices.at(it.col()).col(it.row());
      
      for (Eigen::SparseVector<double>::InnerIterator it2(colo2); it2; ++it2){
        viablesValues.insert(it2.index(),index) = it2.value();
      }
      ++index;
    }
    viablesValues.makeCompressed();
    return(viablesValues);
}    


/*
Eigen::SparseVector<double> rcpp_get_current_vector(Eigen::SparseMatrix<double> currentPresenceMatrix,
                                                    std::vector<Eigen::SparseMatrix<double>> colonisationMatrices,
                                                    Eigen::SparseMatrix<double> globalSuitableSites){
  int nbsites = (globalSuitableSites).nonZeros();
  int x;
  int y;
  int s;
  int v;
  Eigen::SparseVector<double> currentVector(nbsites);
  std::cout << "..............................................."<<std::endl;
  for (int k=0; k<currentPresenceMatrix.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(currentPresenceMatrix,k); it; ++it)
    {
      x = it.row();
      y = it.row();
      s = globalSuitableSites.coeffRef(it.row(),it.col());
      for (Eigen::SparseMatrix<double>::InnerIterator it2(colonisationMatrices.back(),s-1); it2; ++it2) {
        v = it2.value();
        currentVector.coeffRef(it2.row()) = currentVector.coeffRef(it2.row()) +v-v*currentVector.coeffRef(it2.row());
        }
      }
  return(currentVector);
  }
*/

// [[Rcpp::export]]
NumericVector rcpp_eval_current_prob(int threshold,
                                     Eigen::SparseMatrix<double> currentPresenceMatrix,
                                     std::vector<Eigen::SparseMatrix<double>> colonisationMatrices,
                                     Eigen::SparseMatrix<double> globalSuitableSites){
  

  
  int nbsites = (globalSuitableSites).nonZeros();
  int x;
  int y;
  int s;
  double v;
  Eigen::SparseVector<double> currentVector(nbsites);
  //std::cout << "..............................................."<<std::endl;

  for (int k=0; k<currentPresenceMatrix.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(currentPresenceMatrix,k); it; ++it)
    {
      x = it.row();
      y = it.row();
      s = globalSuitableSites.coeffRef(it.row(),it.col());
      //std::cout << "KILLING MACHINE - "<<s<<std::endl;
      
      for (Eigen::SparseMatrix<double>::InnerIterator it2(colonisationMatrices.back(),s-1); it2; ++it2) {
        
        v = it2.value();
        //std::cout << "..............................................."<<it2.row()<<std::endl;
        currentVector.coeffRef(it2.row()) = currentVector.coeffRef(it2.row()) +v-v*currentVector.coeffRef(it2.row());
         
      }
       
       
    }
    

  NumericVector res = eval_probabilityVector(currentVector,threshold);
  //NumericVector res(5) ;
  return(res);
}

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_pheromons(Rcpp::NumericMatrix viablesTriplets){
  int nboc = viablesTriplets.nrow();
  Rcpp::NumericVector pheromons(nboc);
  for (int i=0 ; i<nboc ; i++){
    pheromons[i] = viablesTriplets(i,5)/viablesTriplets(i,3);
  } 
  //pheromons = 1 / (1+pheromons-min(pheromons));
  pheromons = pheromons / sum(pheromons);
  return(pheromons);
}

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_generate_population(Rcpp::NumericVector pheromons,Eigen::SparseMatrix<double> globalSuitableSites,int npop, int nbtoplant){
  
  Rcpp::NumericMatrix population(npop,nbtoplant);
  Rcpp::NumericVector permutation(nbtoplant);
  for (int j=0; j<npop; ++j) {
    //std::cout << "..............................................."<<j<<std::endl;
    permutation=generate_permutation4(nbtoplant,pheromons);
    population(j,_) = permutation ;
  }
  return(population);
}

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_algorithm_opt(Rcpp::NumericVector pheromons,
                                       Rcpp::NumericMatrix viablesTriplets,
                                       Rcpp::NumericMatrix population0,
                                       Eigen::SparseMatrix<double> costMatrix,
                                       Eigen::SparseMatrix<double> currentPresenceMatrix,
                                       std::vector<Eigen::SparseMatrix<double>> colonisationMatrices,
                                       Eigen::SparseMatrix<double> globalSuitableSites,
                                       Eigen::SparseMatrix<double> viablesValues,
                                       int threshold,
                                       double confidence,
                                       int npop,
                                       int nsur,
                                       int ngen,
                                       int nbtoplant){
  
  //return(Rcpp::NumericVector(5));
  
  
  int index;    
  int g;
  int x;
  int y;
  int s;
  double v;
  double cost;
  int indexmini;
  double cout_mini;
  int papa;
  int mama;
  int f;
  double v1;
  double v2;
  
  Rcpp::NumericMatrix population= population0;
  Rcpp::NumericMatrix survivors(npop,nbtoplant);
  Rcpp::NumericVector ranks(npop);
  
  
  
  
  int nboc = viablesTriplets.nrow();
  int nbsites = (globalSuitableSites).nonZeros();
  Eigen::SparseVector<double> currentVector(nbsites);
  //std::cout << "..............................................."<<std::endl;
  for (int k=0; k<currentPresenceMatrix.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(currentPresenceMatrix,k); it; ++it)
    {
      x = it.row();
      y = it.row();
      s = globalSuitableSites.coeffRef(it.row(),it.col());
      //std::cout << "..............................................."<<s<<std::endl;
      
      for (Eigen::SparseMatrix<double>::InnerIterator it2(colonisationMatrices.back(),s-1); it2; ++it2) {
        v = it2.value();
        currentVector.coeffRef(it2.row()) = currentVector.coeffRef(it2.row()) +v-v*currentVector.coeffRef(it2.row());
      }
      
    }


  //return(Rcpp::NumericVector(5));
  //int nbsites = currentVector.size();
  
  Rcpp::NumericMatrix evaluation(npop,2);
  // std::cout << "FINISH HIM" <<npop<<std::endl;
  //
  for (g=0; g<ngen; ++g) {
    std::cout << "Calculating gen number " << g <<std::endl;
    for (int j=0; j<npop;++j) {
      //std::cout << "pop number " <<j<<std::endl;
      Eigen::SparseVector<double> npm = currentVector;
      Eigen::SparseVector<double> lea(nbsites);
      cost = 0;
      
      for (int h=0;h<nbtoplant;++h){
        index = population(j,h);
        if(index != -1){
          
          x = viablesTriplets(index,0);
          y = viablesTriplets(index,1);
          
          
          cost = cost + costMatrix.coeffRef(x,y);
          lea = viablesValues.col(index);

          for (SparseVector<double>::InnerIterator it(lea); it; ++it)
          {
            it.value(); // == vec[ it.index() ]
            it.index();
            v1 = npm.coeffRef(it.index()-1);
            v2 = it.value();
            npm.coeffRef(it.index()-1) = v1 + v2 - v1 * v2;
          }
          
        } 
         
      }
      
      
      
      Rcpp::NumericVector evaluate = eval_probabilityVector(npm, threshold+1);
      evaluation(j,1)=cost;
      
      if((evaluate(threshold)+evaluate(threshold+1))>=confidence)
      {
        evaluation(j,0)=1;
        
        while(evaluate(threshold+1)>=confidence){
          int index_eliminated = optimise_planting_choice4(viablesValues,threshold,confidence,viablesTriplets,population,j,nbtoplant,npm);
          if (index_eliminated==-1){
            break;
          }
          
          if (index_eliminated!=-1){
            evaluate = eval_probabilityVector(npm, threshold+1);
          }
        }
      }
      
      else{
        if(evaluate(threshold)>=confidence){
          optimise_planting_choice5(viablesValues,threshold,confidence,viablesTriplets,population,j,nbtoplant,npm);
          evaluation(j,0)=1;
        }
        else{
          evaluation(j,0)=0;
          evaluation(j,1)=1000000;
        }
      }
    }
    
    //std::cout << "The great evaluation ! : " << threshold<<std::endl;
    //std::cout << evaluation << std::endl;
    
    int saving;
    for (int h=0;h<npop;++h){
      ranks[h]=h;
    }
    for (int h=0;h<nsur;++h){
      indexmini = (h);
      cout_mini = evaluation(ranks(h),1);
      for (int k=h+1;k<npop;++k){
        //Rcout << "Let it flow" << ranks(k) <<std::endl;
        if(cout_mini>evaluation(ranks(k),1)){
          indexmini = (k);
          cout_mini = evaluation(ranks(k),1);
        }
      }
      saving = ranks(h);
      ranks(h)=ranks(indexmini);
      ranks(indexmini)=saving;
    }
    
    for (int h=0;h<nsur;++h){
      survivors(h,_) = population(int(ranks(h)),_);
    }
    
    if (1==1 || g!=ngen){
      int mid = randomfunc3(nbtoplant);
      pheromons = 1 / (1+evaluation(_,1)-min(evaluation(_,1)));
      pheromons = pheromons / sum(pheromons);
      for (int h=nsur;h<npop;++h){
        papa = index_random_choice_non_uniform(pheromons);
        mama = index_random_choice_non_uniform(pheromons);
        for (int j=0;j<nbtoplant;++j){
          f = randomfunc3(2);
          if (f==1){
            survivors(h,j) = population(papa,j);
          }
          else{
            survivors(h,j) = population(mama,j);
          }
          
          f = randomfunc3(10);
          
          if (f==1){
            survivors(h,j) =  randomfunc3(nboc);
          }
        }
      }
      }
    if (true || g!=ngen){
      population = survivors;
      }
  
    }
    
  return(population);
}



// [[Rcpp::export]]  
Rcpp::NumericMatrix rcpp_result_to_choice(Rcpp::NumericMatrix lastPopulation,Rcpp::NumericMatrix viablesTriplets){
  
  int nbtoplant =(lastPopulation).ncol();
  Rcpp::NumericMatrix  choices_matrix(nbtoplant,6);  
  
  //Rcpp::NumericMatrix result(nbtoplant,6);
  
  if(true){
    
    for (int j=0;j<nbtoplant;++j){
      //survivors(h,j) = population(int(ranks(h%nsur)),j);
      if( lastPopulation(1,j)!=-1){
        choices_matrix(j,0) = viablesTriplets(lastPopulation(1,j),0)+1;
        choices_matrix(j,1) = viablesTriplets(lastPopulation(1,j),1)+1;
        choices_matrix(j,2) = viablesTriplets(lastPopulation(1,j),2)+1;
        choices_matrix(j,3) = viablesTriplets(lastPopulation(1,j),3)+1;
        choices_matrix(j,4) = viablesTriplets(lastPopulation(1,j),4)+1;
        //lastPopulation(j,5) = (viablesValues2.col(lastPopulation(1,j))).sum();
      }
      else {
        lastPopulation(j,0) = -1;
        lastPopulation(j,1) = -1;
        lastPopulation(j,2) = -1;
        lastPopulation(j,3) = -1;
      }
    }
  }
  return (choices_matrix);
}


/*
// [[Rcpp::export]]
arma::sp_mat global_suitable_sites(std::list<arma::sp_mat> consecutiveSuitabilityMatrix){
  int nbr = (consecutiveSuitabilityMatrix.back()).n_rows;
  int nbc = (consecutiveSuitabilityMatrix.back()).n_cols;
  arma::sp_mat globalSuitableSites(nbr,nbc);
  Rcpp::NumericVector vec_nbpp(consecutiveSuitabilityMatrix.size());
  int k = 0;
  for (auto const& thisSuitabilityMatrix : consecutiveSuitabilityMatrix){
    vec_nbpp(k) = thisSuitabilityMatrix.n_nonzero;
    ++k;
  }
  int tnbpp = sum(vec_nbpp);
  int j=1;
  for (auto const& thisSuitabilityMatrix : consecutiveSuitabilityMatrix){
    sp_mat::const_iterator start = thisSuitabilityMatrix.begin();
    sp_mat::const_iterator end   = thisSuitabilityMatrix.end();
    
    for(sp_mat::const_iterator it = start; it != end; ++it)
    {
      globalSuitableSites(it.row(),it.col()) = int(j);
      j++;
    }
  }
  return(globalSuitableSites);
}
*/

/*
Rcpp::List transition(std::list<arma::sp_mat> consecutiveSuitabilityMatrix,Rcpp::NumericMatrix migrationKernel) {
  std::list<arma::sp_mat> transitionMatrices;
  
  t = 0;
  std::list<arma::sp_mat> localTransitionMatrix(nbsites,nbsites);
  
  for (index=0;index<nbsites;++index){
    i = globalSuitableCoordinates(index,0);
    j = globalSuitableCoordinates(index,1);
    for (h=max(0,i-migrationRange); h<min(nbr,i+migrationRange+1); ++h){
      for (k=max(0,j-migrationRange); k<min(nbc,j+migrationRange+1); ++k){
        if (globalSuitableSites.coeffRef(h,k)!=0 && migrationKernel(h+migrationRange-i,k+migrationRange-j)!=0){
          
          //std::cout << ".... HOW THE TURNTABLE x="<<(h+migrationRange-i)<< " y="<<(k+migrationRange-j)<<std::endl;
          
          localTransitionMatrix.insert(globalSuitableSites.coeffRef(h,k)-1,index)=migrationKernel(h+migrationRange-i,k+migrationRange-j);
        }
      }
    }
  }
  
  arma::sp_mat transitionMatrix;
  i=0;
  for (auto thisSuitabilityMatrix : consecutiveSuitabilityMatrix){
    transitionMatrix = localTransitionMatrix;
    for (k=0; k<localTransitionMatrix.outerSize(); ++k){
      for (Eigen::SparseMatrix<double>::InnerIterator it(localTransitionMatrix,k); it; ++it)
      {
        
        transitionMatrix.coeffRef(it.row(),it.col()) *= thisSuitabilityMatrix.coeffRef(globalSuitableCoordinates(it.row(),0),globalSuitableCoordinates(it.row(),1));
      }
    }
    transitionMatrix.prune(0.0);
    ++i;
    transitionMatrices.push_back(transitionMatrix);
  }
  
}
*/