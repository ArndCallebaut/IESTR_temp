
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
using namespace arma;

// [[Rcpp::export]]
Eigen::SparseMatrix<double> sqrt_(Eigen::SparseMatrix<double> X) {
  return (X);
}

//[[Rcpp::export]]
Rcpp::List optgenam3 (Eigen::SparseMatrix<double> currentPresenceMatrix,
                      std::list<Eigen::SparseMatrix<double>> consecutiveSuitabilityMatrix,
                      Eigen::SparseMatrix<double> costMatrix,
                      Rcpp::NumericMatrix migrationKernel,
                      int threshold,
                      double confidence,
                      int npop,
                      int nsur,
                      int ngen){
  
  Rcout << "LANCEMENT OPTGENAM 3 !"<<std::endl;
  
  Rcout << "Balise n?0"<<std::endl;
  //std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  
  // **************************************************************************
  // ************************************************ VARIABLES ***************
  // **************************************************************************
  // all matrix, except migrationKernel, have the same size.
  
  // currenPresenceMatrix = matrix of 0 and 1, 1 for presence and 0 for absence.
  // give the presence/absence at the beginning.
  
  // consecutiveSuitabilityMatrix = list of suitability matrix (SM_i).
  // SM_i is suitability at time i, values between 0 and 1.
  
  // costMatrix = cost of planting for each cell.
  
  // migrationKernel = matrix expliciting how a specie migrate.
  // it's a square matrix, the size must be an odd number.
  
  // threshold, confidence = parameters to satisfied for the optimisation.
  // we have to end with more cells occuped by threshold, with a confidence.
  
  // npop, nsur, ngen = parameters for optimization work.
  
  // **************************************************************************
  // **************************************** DEDUCED VARIABLES ***************
  // **************************************************************************
  
  // number of rows and columns in the map.
  int nbr = currentPresenceMatrix.rows();
  int nbc = currentPresenceMatrix.cols();
  
  // number of periods covered.
  int nbperiod = consecutiveSuitabilityMatrix.size();
  
  // working variables.
  int k;
  int t;
  int i;
  int j;
  int h;
  int index;
  
  // max range of migration matrix
  int migrationRange = (migrationKernel.rows()-1)/2;
  
  // **************************************************************************
  // ****************************** I. ORGANISE SUITABLE MATRIX ***************
  // **************************************************************************
  
  Rcout << "Balise n?1"<<std::endl;
  
  // A matrix of allsuitable site, the "working sites"
  // will contains reference to the first occurence in possibilities.
  Eigen::SparseMatrix<int> globalSuitableSites(nbr,nbc);
  
  // Counter for number of planting possibilities in space (nbpp), in each period.
  Rcpp::NumericVector vec_nbpp(nbperiod);
  
  
  // Counting the npps, in each period
  k = 0;
  for (auto const& thisSuitabilityMatrix : consecutiveSuitabilityMatrix){
    vec_nbpp(k) = thisSuitabilityMatrix.nonZeros();
    ++k;
  }
  
  // Total number of planting possibilities
  int tnbpp = sum(vec_nbpp);
  
  // All the choices you can have to do assisted migration.
  // col1=row ; col2=column ; col3 = time ; col4=suitability ;col5=cost
  Rcpp::NumericMatrix possibilities(tnbpp,5);
  
  t=0;
  i=0;
  j=1;
  // for each of the suitability matrix...
  for (auto const& thisSuitabilityMatrix : consecutiveSuitabilityMatrix){
    // for each of the elements in the suitability matrix
    for (k=0; k<thisSuitabilityMatrix.outerSize(); ++k)
      for (Eigen::SparseMatrix<double>::InnerIterator it(thisSuitabilityMatrix,k); it; ++it)
      {
        // construct global suitable sparse matrix.
        if(globalSuitableSites.coeffRef(it.row(),it.col())==0){
          globalSuitableSites.coeffRef(it.row(),it.col()) = int(j);
          j++;
        }
        
        // construct the possibilities matrix.
        possibilities(i,0)=it.row();
        possibilities(i,1)=it.col();
        possibilities(i,2)=t;
        possibilities(i,3)=it.value();
        possibilities(i,4)=costMatrix.coeffRef(it.row(),it.col());
        ++i;
      }
      ++t;
  }
  
  
  
  // Number of working sites, which are at least one time suitable
  int nbsites = globalSuitableSites.nonZeros();
  //Rcout << "TARGET ACCUIRED "<<nbsites<<std::endl;
  
  // coordinates+value from the sparse globalSuitableSites
  Rcpp::NumericMatrix globalSuitableCoordinates(nbsites,3);
  // col1 = row ; col2 = col
  // the order in this matrix is the one used in the matrices manipulation afterward
  index = 0;
  for (k=0; k<globalSuitableSites.outerSize(); ++k)
    for (Eigen::SparseMatrix<int>::InnerIterator it(globalSuitableSites,k); it; ++it)
    {
      globalSuitableCoordinates(it.value()-1,0)=it.row();
      globalSuitableCoordinates(it.value()-1,1)=it.col();
      globalSuitableCoordinates(it.value()-1,2)=it.value()-1;
    }
    
    // std::cout << "Global suitable sites"<<std::endl;
    // std::cout << Eigen::MatrixXd(globalSuitableSites.cast<double>()) << std::endl;
    
    // **************************************************************************
    // ********************************** II. TRANSITION MATRICES ***************
    // **************************************************************************
    
    std::cout<< "Balise n?2"<<std::endl;
  
  //Rcout << "Checking"<<std::endl;
  
  std::list<Eigen::SparseMatrix<double>> transitionMatrices;
  t = 0;
  Eigen::SparseMatrix<double> localTransitionMatrix(nbsites,nbsites);
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
  
  //std::cout<< "Skeleton... "<<std::endl;
  //std::cout << Eigen::MatrixXd( localTransitionMatrix) << std::endl;
  
  
  
  
  
  //std::cout << "Checking"<<std::endl;
  Eigen::SparseMatrix<double> transitionMatrix;
  
  for (auto thisSuitabilityMatrix : consecutiveSuitabilityMatrix){
    
    //std::cout << "I'm here for the computer.."<< thisSuitabilityMatrix.nonZeros()<<std::endl;
    
  }
  
  
  
  
  i=0;
  for (auto thisSuitabilityMatrix : consecutiveSuitabilityMatrix){
    
    
    transitionMatrix = localTransitionMatrix;
    //std::cout << " OK bon on part de "<<transitionMatrix.nonZeros()<<std::endl;
    // for each row...
    
    for (k=0; k<localTransitionMatrix.outerSize(); ++k){
      for (Eigen::SparseMatrix<double>::InnerIterator it(localTransitionMatrix,k); it; ++it)
      {
        
        transitionMatrix.coeffRef(it.row(),it.col()) *= thisSuitabilityMatrix.coeffRef(globalSuitableCoordinates(it.row(),0),globalSuitableCoordinates(it.row(),1));
      }
    }
    
    transitionMatrix.prune(0.0);
    
    
    //std::cout << "Transition Matrice !!! --- "<<i<<std::endl;
    ++i;
    
    transitionMatrices.push_back(transitionMatrix);
    //std::cout << Eigen::MatrixXd(transitionMatrix) << std::endl;
    
  }
  
  
  
  
  //std::cout << ".... HOW THE TURNTABLE "<<localTransitionMatrix.nonZeros()<<std::endl;
  
  // **************************************************************************
  // ******************************* III. COLONISATION MATRICES ***************
  // **************************************************************************
  
  Rcout << "Balise n?3"<<std::endl;
  
  
  std::vector<Eigen::SparseMatrix<double>> colonisationMatrices;
  colonisationMatrices.push_back(transitionMatrices.back());
  
  
  // 
  // std::cout << "Colonisation Matrice !!! --- 0"<<std::endl;
  // std::cout << Eigen::MatrixXd(transitionMatrices.back()) << std::endl;
  
  transitionMatrices.reverse();
  
  k=1;
  bool first = true;
  for (auto const& elt : transitionMatrices) {
    if (first){
      first = false;
    }
    else{
      //Rcout << "EMPTINESS"<<std::endl;
      //Rcout << "TIME TO SUFFER "<< elt.nonZeros()<<"---kill me  "<<colonisationMatrices.back().nonZeros()<<std::endl;
      colonisationMatrices.push_back(proba_matrix_mult3(colonisationMatrices.back(),elt));
      
      //std::cout << "Colonisation Matrice !!! --- "<<k<<std::endl;
      ++k;
      //std::cout << Eigen::MatrixXd(colonisationMatrices.back()) << std::endl;
    }
    
    
  }
  
  //colonisationMatrices.reverse();
  reverse(colonisationMatrices.begin(),colonisationMatrices.end());
  
  // **************************************************************************
  // ****************** IV. ELIMINATE NONE OPTIMAL TIME CHOICES ***************
  // **************************************************************************
  
  // std::cout<< "Balise n?4"<<std::endl;
  // 
  // Eigen::SparseMatrix<int> viableSites(nbperiod,nbsites);
  // Eigen::SparseMatrix<int> viableSites2(nbsites,nbperiod);
  // //viableSites.coeffRef(0,nbsites-1)=1;
  // bool av1;
  // bool av2;
  // double var;
  // 
  // for (i=0;i<nbsites;++i){
  //   
  // 
  //   //std::cout << "i="<< i << "/"<<nbsites<<""<< ".... well well welloooooo "<<std::endl;
  // 
  //   // viability = 1
  //   Rcpp::NumericVector viability (nbperiod);
  // 
  //   for (j=0;j<nbperiod;++j){
  //     viability(j)=1;
  //   }
  //   
  // 
  //   h=0;
  //   for (auto suit1 : colonisationMatrices){
  //     //std::cout << "Your sins will follow you ! "<<h<<std::endl;
  //     //Rcout << "Your sins will follow you "<<h<<std::endl;
  //     k=0;
  //     for (auto suit2 : colonisationMatrices){
  // 
  //       bool nonzero1 = false;
  //       bool nonzero2 = false;
  // 
  //       if (viability(k)==0 || viability(h)==0 ){
  //         break;
  //       }
  // 
  //       av1 = false;
  //       av2 = false;
  // 
  //       if (k>=h){
  //         ++k;
  //         break;
  //       }
  // 
  //       for (j=0;j<nbsites;++j){
  //         
  //         if(suit1.coeffRef(j,i)!=0){
  //           nonzero1 = true;
  //         }
  // 
  //         if(suit2.coeffRef(j,i)!=0){
  //           nonzero2 = true;
  //         }
  // 
  //         var = suit1.coeffRef(j,i)-suit2.coeffRef(j,i);
  // 
  //         if(var>0){
  //           av1=true;
  //         }
  //         if(var<0){
  //           av2=true;
  //         }
  //         if(av1 && av2){
  //           break;
  //         }
  // 
  //       }
  // 
  //       if(!nonzero1){
  //         // si le premier site ne donne aucun r?sultat
  //         viability(h)=0;
  //       }
  // 
  //       if(!nonzero2){
  //         // si le second site ne donne aucun r?sultat
  //         viability(k)=0;
  //       }
  // 
  // 
  // 
  //       if (av1 && !av2){
  // 
  //         viability(k)=0;
  //       }
  //       if (av2 && !av1){
  //         viability(h)=0;
  //       }
  //       ++k;
  //     }
  //     ++h;
  //   }
  // 
  //   //Rcout << "wHOW THE TURNTABLE "<<std::endl;
  // 
  // 
  //   for (j=0;j<nbperiod;++j){
  //     if (viability[j]==1){
  //       //std::cout<< "I.. cut throw.. i="<<i<<" %% j="<<j<<" %% rows="<<viableSites.rows()<<" %% cols="<<viableSites.cols()<<std::endl;
  //       viableSites.coeffRef(int(j),int(i))=1;
  //       viableSites2.coeffRef(int(i),int(j))=1;
  //     }
  //   }
  //   
  //   
  // }
  
  
  std::cout<< "Balise n°4 BIS"<<std::endl;
  
  Eigen::SparseMatrix<int> viableSites(nbperiod,nbsites);
  Eigen::SparseMatrix<int> viableSites2(nbsites,nbperiod);
  //viableSites.coeffRef(0,nbsites-1)=1;
  bool av1;
  bool av2;
  double var;
  
  
  
  
  
  
  for (i=0;i<nbsites;++i){
    
    
    //std::cout << "i="<< i << "/"<<nbsites<<""<< ".... well well welloooooo "<<std::endl;
    
    // viability = 1
    Rcpp::NumericVector viability (nbperiod);
    
    for (j=0;j<nbperiod;++j){
      
      Eigen::SparseVector<int> suit2 = colonisationMatrices[j].row(i);
      
      if(suit2.nonZeros()!=0){
        viability(j)=1;
      }
      
      
      
      //viability(j)=1;
    }
    
    
    h=0;
    
    
    for (int ss=0; ss<(nbperiod-1); ++ss){
      h=ss;
      Eigen::SparseVector<int> suit1 = colonisationMatrices[ss].row(i);
      
      if ( viability(h)==0 ){
        break;
      }
      
      
      //std::cout << "Your sins will follow you ! "<<h<<std::endl;
      //Rcout << "Your sins will follow you "<<h<<std::endl;
      //k=0;
      for (int sss=ss+1; sss<nbperiod; ++sss){
        k=sss;
        Eigen::SparseVector<int> suit2 = colonisationMatrices[sss].row(i);
        
        bool nonzero1 = false;
        bool nonzero2 = false;
        
        if (viability(k)==0 || viability(h)==0 ){
          break;
        }
        
        av1 = false;
        av2 = false;
        
        if (k>=h){
          ++k;
          break;
        }
        
        for (j=0;j<nbsites;++j){
          
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
          // si le premier site ne donne aucun r?sultat
          viability(h)=0;
        }
        
        if(!nonzero2){
          // si le second site ne donne aucun r?sultat
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
    
    //Rcout << "wHOW THE TURNTABLE "<<std::endl;
    
    
    for (j=0;j<nbperiod;++j){
      if (viability[j]==1){
        //std::cout<< "I.. cut throw.. i="<<i<<" %% j="<<j<<" %% rows="<<viableSites.rows()<<" %% cols="<<viableSites.cols()<<std::endl;
        
        if((colonisationMatrices.at(j).col(i)).sum()>0.5){
          viableSites.coeffRef(int(j),int(i))=1;
          viableSites2.coeffRef(int(i),int(j))=1;
        }
        

      }
    }
    
    
  }
  
  
  
  
  
  
  
  // **************************************************************************
  // ************************************ V. REWRITE THE RESULT ***************
  // **************************************************************************
  
  std::cout << "Balise n?5 "<<std::endl;
  int nboc = viableSites.nonZeros();
  // number of optimal choices
  //std::cout<< "Seeing my viablesSites"<<std::endl;
  //std::cout << Eigen::MatrixXd(viableSites.cast<double>()) << std::endl;
  
  // final viable choices
  Rcpp::NumericMatrix viablesTriplets(nboc,5);
  Rcpp::NumericMatrix viablesValues(nbsites,nboc);
  
  t = 0;
  index = 0;
  
  // for each period...
  for (auto colma : colonisationMatrices){
    
    Eigen::SparseVector<int> loc = viableSites.row(t);
    
    
    for (Eigen::SparseVector<int>::InnerIterator it(loc); it; ++it){
      //Rcout << "Ok boss. "<<it.value()<<std::endl;
      //if(it.value()!=0){
      
      
      viablesTriplets(index,0) = globalSuitableCoordinates(it.index(),0);
      viablesTriplets(index,1) = globalSuitableCoordinates(it.index(),1);
      viablesTriplets(index,2) = t;
      viablesTriplets(index,3) = costMatrix.coeffRef(viablesTriplets(index,0),viablesTriplets(index,1));
      viablesTriplets(index,4) = globalSuitableSites.coeffRef(viablesTriplets(index,0),viablesTriplets(index,1));
      
      // for (int j=0; j<nbsites;++j){
      //   viablesValues(j,index) = colma.coeffRef(j,globalSuitableSites.coeffRef(it.row(),it.col()));
      // }
      ++index;
      
    }
    
    
    ++t;
  }
  
  // Rcout << "..............................................."<<std::endl;
  // std::cout << "Viables triples"<<std::endl;
  // std::cout << viablesTriplets<< std::endl;
  // 
  // // **************************************************************************
  // // ************************** VI. REWRITE THE RESULT AGAIN... ***************
  // // **************************************************************************
  // 
  std::cout << "Balise n°6 nbsites="<<nbsites<<" nboc="<<nboc<<" et INDEX_MAX="<<index<<std::endl;
  // 
  // // number of optimal choices
  // std::cout<< "Seeing my viablesSites NBOC="<<nboc<<std::endl;
  // std::cout<< "Seeing my viablesSites NSSITES="<<nbsites<<std::endl;
  // //std::cout << Eigen::MatrixXd(viableSites.cast<double>()) << std::endl;
  
  
  // 
  // std::cout << "TEST BALISING"<<std::endl;
  // 
  //   for (k=0; k<viableSites2.outerSize(); ++k)
  //     for (Eigen::SparseMatrix<int>::InnerIterator it(viableSites2,k); it; ++it)
  //     {
  //       std::cout << "i="<<it.row()<<" j="<<it.col()<<std::endl;
  //     }
  
  
  // 
  // final viable choices
  //Eigen::SparseMatrix<double> viablesValues2(nbsites,nboc);
  
  t = 0;
  index = 0;
  
  //for each period...
  // for (auto colma : colonisationMatrices){
  // 
  //   for (k=0; k<viableSites.outerSize(); ++k)
  //     for (Eigen::SparseMatrix<int>::InnerIterator it(viableSites,k); it; ++it)
  //     {
  //       
  //       //std::cout << "Balise n?6 "<<nbsites<<" "<<nboc<<std::endl;
  //       
  // 
  //       if (it.row()==t){
  // 
  //         if(it.value()==0){
  //           // std::cout << "ET MERDE"<<std::endl;
  //           continue;
  //         }
  //         // std::cout << "ET MERDE..."<<std::endl;
  // 
  // 
  //         int xx = viablesTriplets(index,0);
  //         int yy = viablesTriplets(index,1);
  //         double val = globalSuitableSites.coeffRef(xx,yy) ;
  //         
  //         // std::cout << "On regarde à t="<<t<< " x=" << xx << " et y=" << yy <<std::endl;
  //         
  //         // viablesTriplets(index,2) = t;
  //         // viablesTriplets(index,3) = costMatrix.coeffRef(viablesTriplets(index,0),viablesTriplets(index,1));
  // 
  //         
  //         // for (Eigen::SparseVector<int>::InnerIterator it(colo2); it; ++it){
  //         //   viablesValues2.insert(it.index(),index) = it.value();
  //         // }
  //         
  //         
  // 
  //         for (int j=0; j<nbsites;++j){
  // 
  //           //std::cout << "\n\nJ'essaie... j="<<j<< " in " << val << " time..." << t <<std::endl;
  // 
  //           if (colma.coeffRef(j,val-1)!=0){
  // 
  //             //std::cout << "t="<<t<<" numsite="<<j<<" val="<<val<< " index="<<index/nboc<<std::endl;
  //             viablesValues2.insert(j,index) = colma.coeffRef(j,val-1);
  //             //std::cout << "t="<<t<<" numsite="<<j<<" val="<<val<<std::endl;
  // 
  //           }
  //         }
  // 
  //         ++index;
  // 
  //       }
  //     }
  //     ++t;
  // }
  
  
  std::cout << "Balise n°6 BIS nbsites="<<nbsites<<" nboc="<<nboc<<std::endl;
  //
  // number of optimal choices
  // std::cout<< "Seeing my viablesSites NBOC="<<nboc<<std::endl;
  // std::cout<< "Seeing my viablesSites NSSITES="<<nbsites<<std::endl;
  // std::cout << Eigen::
  
  if (nboc==0){
    List z = List::create();
    return(z);
  }
  
  Eigen::SparseMatrix<double> viablesValues2(nbsites,nboc);
  
  t = 0;
  index = 0;
  
  //for each period...
  for (k=0; k<viableSites2.outerSize(); ++k)
    for (Eigen::SparseMatrix<int>::InnerIterator it(viableSites2,k); it; ++it)
    {
      if(it.value()==0){
        // std::cout << "ET MERDE"<<std::endl;
        continue;
      }
      // std::cout << "ET MERDE..."<<std::endl;
      
      
      
      int xx = viablesTriplets(index,0);
      int yy = viablesTriplets(index,1);
      double val = globalSuitableSites.coeffRef(xx,yy) ;
      
      Eigen::SparseVector<double> colo2 = colonisationMatrices.at(it.col()).col(it.row());
      
      //std::cout << "ET MERDE..."<<std::endl;
      
      for (Eigen::SparseVector<double>::InnerIterator it2(colo2); it2; ++it2){
        if(it2.value()<=0){
          std::cout << "I AM THE AVATAR... valin col="<<colonisationMatrices.at(it.col()).coeffRef(it2.index(),it.row())<< "value="<<it2.value()<<std::endl;
        }
      }

      
      
      for (Eigen::SparseVector<double>::InnerIterator it2(colo2); it2; ++it2){
        
        viablesValues2.insert(it2.index(),index) = it2.value();
        
        //if(it2.value()<0){
          //std::cout << "FOUND YOU, BITCH : insert it2index="<<it2.index()<<" index="<<index<< "value="<<it2.value()<<std::endl;
        //}
        
      }
      
      // std::cout << "On regarde à t="<<t<< " x=" << xx << " et y=" << yy <<std::endl;
      
      // viablesTriplets(index,2) = t;
      // viablesTriplets(index,3) = costMatrix.coeffRef(viablesTriplets(index,0),viablesTriplets(index,1));
      
      
      // for (int j=0; j<nbsites;++j){
      // 
      //   //std::cout << "\n\nJ'essaie... j="<<j<< " in " << val << " time..." << t <<std::endl;
      // 
      //   if (colma.coeffRef(j,val-1)!=0){
      // 
      //     //std::cout << "t="<<t<<" numsite="<<j<<" val="<<val<< " index="<<index/nboc<<std::endl;
      //     viablesValues2.insert(j,index) = colma[it.row()].coeffRef(j,val-1);
      //     //std::cout << "t="<<t<<" numsite="<<j<<" val="<<val<<std::endl;
      // 
      //   }
      // }
      ++index;
    }
    
    
    
    std::cout << "..............................................."<<std::endl;
  //Rcout << "Viables Values"<<std::endl;
  //std::cout << Eigen::MatrixXd(viablesValues) << std::endl;
  
  //return(viablesTriplets);
  
  
  // **************************************************************************
  // ****************** VII. TRANSLATE INPUT INTO USABLE THINGS ***************
  // **************************************************************************
  
  int x;
  int y;
  
  int s;
  int v;
  //arma::sp_mat currentVector(nbsites);
  Eigen::SparseVector<double> currentVector(nbsites);
  std::cout << "..............................................."<<std::endl;
  for (k=0; k<currentPresenceMatrix.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(currentPresenceMatrix,k); it; ++it)
    {
      x = it.row();
      y = it.row();
      s = globalSuitableSites.coeffRef(it.row(),it.col());
      
      //colonisationMatrices.back()
      
      for (Eigen::SparseMatrix<double>::InnerIterator it2(colonisationMatrices.back(),s-1); it2; ++it2) {
        
        v = it2.value();
        currentVector.coeffRef(it2.row()) = currentVector.coeffRef(it2.row()) +v-v*currentVector.coeffRef(it2.row());
        
      }
    }
    
    
    
    
    
    
    
  std::cout << "Balise n?7"<<std::endl;
  //std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  
  
  
  NumericVector initial_state = eval_probabilityVector(currentVector,threshold);
  
  v = 0;
  double p = 0;
  for (h=threshold; h>=0; --h){
    p = p + initial_state[h];
    if(p>=confidence){
      v = h;
      break;
    }
  }
  // How many planting do we need ?
  int nbtoplant = (threshold - v);
  if (nbtoplant>nboc){
    std::cout << "ON PEUT DIRE QUE C'EST FOUTU LES BITCHIS TROP A PLANTER. PLANTEZ TOUT"<<std::endl;
    List z = List::create(0,colonisationMatrices);
    return(z);
  }
  
  if (nbtoplant==0){
    std::cout << "ON PEUT DIRE QUE C'EST TERMINE LES BITCHIS RIEN A PLANTER"<<std::endl;
    List z = List::create(0,colonisationMatrices);
    return(z);
  }
  
  nbtoplant = nbtoplant +10;
  
  
  // **************************************************************************
  // 2. Prepare the possible planting sites
  //
  

  
  std::cout << "Balise n7 ; Preparing ant feromons "<<std::endl;  
  //std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  
  Rcpp::NumericVector pheromons(nboc);
  Rcpp::NumericVector sitespheromons(nbsites);
  //std::cout << " pheromons legth="<<pheromons.length()<<std::endl;
  

  
  
  
  for (i=0 ; i<nboc ; i++){
    //std::cout << "Miam miam my pheromons ="<<viablesTriplets(i,3)<<" and...="<<viablesTriplets(i,4)<<std::endl; 
    pheromons[i] = viablesTriplets(i,3);
    sitespheromons[viablesTriplets(i,4)-1] = viablesTriplets(i,3);
    //sitespheromons(globalSuitableSites.coeffRef(viablesTriplets(i,0),viablesTriplets(i,1))) = viablesTriplets(i,3);
  } 
  

  //std::cout << " pheromons legth="<<pheromons.length()<<std::endl;
  pheromons = 1 / (1+pheromons-min(pheromons));
  sitespheromons = 1 / (1+sitespheromons-min(sitespheromons));
  //std::cout << " pheromons legth="<<pheromons.length()<<std::endl;
  pheromons = pheromons / sum(pheromons);
  sitespheromons = sitespheromons / sum(sitespheromons);
  
  std::cout << "Checking the feromons feromons "<<std::endl;  
  //std::cout << pheromons <<std::endl;  
  std::cout << " pheromons legth="<<pheromons.length()<<std::endl;
  

  std::cout << "Balise n?8 ; nbtoplant : "<<nbtoplant<<" nboc="<<nboc<<std::endl;
  //std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  
  
  // **************************************************************************
  // 3. Genetic algorithm
  //
  
  // Initialise population
  Rcpp::NumericMatrix population(npop,nbtoplant);
  
  Rcpp::NumericVector permutation(nbtoplant);
  for (j=0; j<npop; ++j) {
    // 
    
    //std::cout << "creating population nbpop="<<j<<" pheromons legth="<<pheromons.length()<<std::endl;
    
    //Rcpp::NumericVector blork(10);
    //blork =  generate_permutation3(10,10) ;
    permutation=generate_permutation4(nbtoplant,pheromons);
    //permutation=generate_permutation3(nbtoplant,nboc);
    //std::cout << "Checking the feromons feromons "<< permutation<<std::endl; 
    
    population(j,_) = permutation ;
    //int alea = randomfunc3(10);
    
    // 
  }
  // 
  // 
  
  List z = List::create(population);
  
  //return(z);
  

  
  
  
  
  
  // 
  // 
  // 
    //std::cout << "Ma population de depart.. "<<std::endl;
    //std::cout << population<<std::endl;
    //std::cout << "Our Triplets "<<std::endl;
    //std::cout << (viablesTriplets) << std::endl;
    //std::cout << "resultants values... "<<std::endl;
    //std::cout << Eigen::MatrixXd(currentVector) << std::endl;
  
  
  
  
  
  
  
  int g;
  
  double cost;
  int indexmini;
  double cout_mini;
  int papa;
  int mama;
  int f;
  
  
  double v1;
  double v2;
  Rcpp::NumericMatrix survivors(npop,nbtoplant);
  Rcpp::NumericVector ranks(npop);

  
  
  std::cout << "Balise n....9"<<std::endl;
  //std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  
  
  Rcpp::NumericMatrix evaluation(npop,2);
  // std::cout << "FINISH HIM" <<npop<<std::endl;
  //
  for (g=0; g<ngen; ++g) {
    std::cout << "inner gates... OPENED !" << g <<std::endl;
    //std::cout << population <<std::endl;
    
    // evaluation of each set of choices
    for (j=0; j<npop;++j) {
      
      //djdjkho
      
      //std::cout << "Calling my population number..." <<j<<std::endl;
      // construction of new local vector  tested
      Eigen::SparseVector<double> npm = currentVector;
      Eigen::SparseVector<double> lea(nbsites);
      //
      cost = 0;
      
      for (h=0;h<nbtoplant;++h){
        
        index = population(j,h);
        if(index != -1){
          
          x = viablesTriplets(index,0);
          y = viablesTriplets(index,1);
          cost = cost + costMatrix.coeffRef(x,y);
          //
          //npm = npm + viablesValues2.col(index);
          //
          lea = viablesValues2.col(index);
          //std::cout << "OK BOOMER" <<(lea)<<std::endl;
          
          //npm = npm + lea - npm * lea;
          
          
          /*
          for (SparseVector<double>::InnerIterator it(lea); it; ++it)
          {
            
            v = it.value();
            //std::cout<< "+++" <<v<<std::endl;
            if (v!=0){
              std::cout<< "===" <<v<<std::endl;
              //npm.coeffRef(it.index()) = npm.coeffRef(it.index()) +v-v*npm.coeffRef(it.index());
            }
          }*/
          
          //npm = npm + 0*lea;
          
          
          //std::cout << "I fil a bit uneasy - "<<npm.size()<< "yeah" <<lea.size()<<std::endl;
          
          
          for (SparseVector<double>::InnerIterator it(lea); it; ++it)
          {
            it.value(); // == vec[ it.index() ]
            it.index();
            //std::cout << "GIMMI YOUR GOLD " <<it.value()<<" - "<< npm.coeffRef(it.index()) <<std::endl;
            
            v1 = npm.coeffRef(it.index()-1);
            v2 = it.value();
            //std::cout << "GIMMI YOUR GOLD v1=" <<v1<<" - "<<it.index()<<" - v2="<< v2<< " - chech on that vtot="<<(v1+v2-v1*v2)<<std::endl;
            
            npm.coeffRef(it.index()-1) = v1 + v2 - v1 * v2;
          }
          
          /*
          // We do the addition of the effect one by one
          for (Eigen::SparseMatrix<double>::InnerIterator it(viablesValues2,index); it; ++it) {
            
            v = it.value();
            //std::cout << "Valeur="<< <<npop<<std::endl;
            
            if (v!=0){
              std::cout<< "===" <<v<<std::endl;
              npm.coeffRef(it.row()) = npm.coeffRef(it.row()) +v-v*npm.coeffRef(it.row());
            }
            
            
          }
           */
        }
      }
      
      
      //std::cout << "It's just that..." << npm <<std::endl;
      
      Rcpp::NumericVector evaluate = eval_probabilityVector(npm, threshold+1);
      evaluation(j,1)=cost;
      
      //std::cout << "SON "<<j<<", GIMMI YOUR GRADE " <<std::endl;
      //std::cout << population(j,_) <<std::endl;
      //std::cout << evaluate <<std::endl;
      //std::cout << "what have I done..." << threshold<<std::endl;
      //std::cout << npm <<std::endl;
      //std::cout << "...and ?" << threshold<<std::endl;
     //std::cout << npm <<std::endl;

      
      
      
      
      if((evaluate(threshold)+evaluate(threshold+1))>=confidence)
      {
        evaluation(j,0)=1;
        
        //std::cout << "BLEUHBLEI" << std::endl;
        
        while(evaluate(threshold+1)>=confidence){
          //std::cout << "CHEH " << std::endl;
          //ptimise_planting_choice4 (Eigen::SparseMatrix<double> viablesValues2, int threshold, double confidence, Rcpp::NumericMatrix viablesTriplets, Rcpp::NumericMatrix population, int pop, int nbtoplant,  Eigen::SparseVector<double> current)
          //std::cout << "Stand down" <<std::endl;
          int index_eliminated = optimise_planting_choice4(viablesValues2,threshold,confidence,viablesTriplets,population,j,nbtoplant,npm);
          
          if (index_eliminated==-1){
            break;
          }
          
          if (index_eliminated!=-1){
            //std::cout << "Not good.." <<population(j,index_eliminated) <<std::endl;
            // for (Eigen::SparseMatrix<double>::InnerIterator it(viablesValues2,population(j,index_eliminated)); it; ++it) {
            //   std::cout << "I may not..."<<it.value()<< " --- "<<npm <<std::endl;
            //   npm.coeffRef(it.index()-1) = (npm.coeffRef(it.index()-1)-it.value())/(1-it.value());
            // }
            //std::cout << "I may not... survive..." <<std::endl;
            
            evaluate = eval_probabilityVector(npm, threshold+1);
            //std::cout << "Seems like we can eliminate someone in "<<j<<" - "<< std::endl;
          }
        }
      }
      
      else{
        if(evaluate(threshold)>=confidence){
          optimise_planting_choice5(viablesValues2,threshold,confidence,viablesTriplets,population,j,nbtoplant,npm);
          evaluation(j,0)=1;
        }
        else{
          evaluation(j,0)=0;
          evaluation(j,1)=1000000;
        }
      }
      
      //std::cout << "After all, i'm doing "<< eval_probabilityVector(npm, threshold+1) <<std::endl;
      //
      
      //std::cout << "Balise n?9"<<std::endl;
    }
    //std::cout << "OUR PEOPLE WILL PREVAIL" <<std::endl;
    //std::cout << population <<std::endl;
    //std::cout << "The great evaluation ! : " << threshold<<std::endl;
    //std::cout << evaluation << std::endl;
    
    //std::cout << "Lotus flower, recto !" << npop <<std::endl;
    
    //bubble sort. (waiting for the fusionsort.)
    // Ranking population by survability
  
  
    
    int saving;
    
    for (h=0;h<npop;++h){
      ranks[h]=h;
    }
    for (h=0;h<nsur;++h){
      indexmini = (h);
      cout_mini = evaluation(ranks(h),1);
      for (k=h+1;k<npop;++k){
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

    
    
    
    //survivors...
    
    for (h=0;h<nsur;++h){
      survivors(h,_) = population(int(ranks(h)),_);
    }
    
    //and their children !
    
    
    
    if (1==1 || g!=ngen){
    
    //Rcpp::NumericMatrix survivors(npop,nbtoplant);
      
    int mid = randomfunc3(nbtoplant);
    
    //Rcout << "Achiente !" <<std::endl;

    
    pheromons = 1 / (1+evaluation(_,1)-min(evaluation(_,1)));
    pheromons = pheromons / sum(pheromons);
    
    
    //std::cout << "pherommons N°"<<g<<std::endl;
    //std::cout << pheromons<< std::endl;
    
    for (h=nsur;h<npop;++h){
      
      papa = index_random_choice_non_uniform(pheromons);
      mama = index_random_choice_non_uniform(pheromons);
      //std::cout << "papa N°"<<papa<<std::endl;
      // Rcout << "papa=" <<papa<<" mama="<<mama<<std::endl;
      
      /*
      
      
      for (j=0;j<mid;++j){
        //survivors(h,j) = population(int(ranks(h%nsur)),j);
        survivors(h,j) = population(papa,j);
      }
      //Rcout << "iju !" <<std::endl;
      for (j=mid;j<nbtoplant;++j){
        //survivors(h,j) = population(int(ranks((7*h+1)%nsur)),j);
        survivors(h,j) = population(mama,j);
      }
      
      */
      for (j=0;j<nbtoplant;++j){
        
        f = randomfunc3(2);
        //survivors(h,j) = population(int(ranks(h%nsur)),j);
        
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
      
      
      //survivors(h*%nsur,Range(mid+1,nbtoplant)) = population(int(ranks((h*2+1)%nsur)),Range(mid+1,nbtoplant));
      
      
    }
    
    if (true || g!=ngen){
      population = survivors;
    }
  }
  
  Rcpp::NumericMatrix survivors2(npop,nbtoplant);
  Rcpp::NumericVector values(npop);
  Eigen::SparseVector<double> lea2(nbsites);
  //lea2 = viablesValues2.col(index);
  //Rcpp::NumericVector ranks(npop);
  int saving;
  
// 
//   for (h=0;h<nbtoplant;++h){
//     
//     index = population(0,h);
//     if(index != -1){
//       
//       x = viablesTriplets(index,0);
//       y = viablesTriplets(index,1);
//       cost = costMatrix.coeffRef(x,y);
//       std::cout << "Cost ="<<cost<<std::endl;
//     }
//     
//     if(index == -1){
//       
//       std::cout << "Dead cell bab"<<std::endl;
//     }
//     
//   }
  
  /*
  
  for (h=0;h<npop;++h){
    ranks[h]=h;
  }
  for (h=0;h<nsur;++h){
    indexmini = (h);
    cout_mini = evaluation(ranks(h),1);
    for (k=h+1;k<npop;++k){
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
  
  for (h=0;h<npop;++h){
    survivors2(h,_) = population(int(ranks(h)),_);
    values(h) = evaluation(int(ranks(h)),1);
  }
  */
  
  
  std::cout << "underwolrd1"<<std::endl;
  std::cout << evaluation(0,1) <<std::endl;
  
  
  std::cout << "underwolrd2"<<std::endl;
  std::cout << ranks <<std::endl;
  
  
  //std::cout << "resultants values... "<<std::endl;
  //std::cout << Eigen::MatrixXd(currentVector) << std::endl;
  
  
  //List z2 = List::create(1);
  //return(z2);
  

  
  //Rcpp::NumericVector mess =population(ranks(0),_);
  
  std::cout << "Pop' "<<std::endl;
  for (j=0;j<nbtoplant;++j){
    index = population(0,j);
    std::cout << population(0,j) << std::endl;
    if (population(0,j)!= -1){
      x = viablesTriplets(index,0);
      y = viablesTriplets(index,1);
      cost = costMatrix.coeffRef(x,y);
      std::cout << "paying checked ="<<cost<<std::endl;
    }
  }
  
  //std::cout << "underwolrd1"<<std::endl;
  
  Rcpp::NumericMatrix result(nbtoplant,6);
  
  if(true){
  
    for (j=0;j<nbtoplant;++j){
      //survivors(h,j) = population(int(ranks(h%nsur)),j);
      if( population(1,j)!=-1){
        result(j,0) = viablesTriplets(population(1,j),0)+1;
        result(j,1) = viablesTriplets(population(1,j),1)+1;
        result(j,2) = viablesTriplets(population(1,j),2)+1;
        result(j,3) = viablesTriplets(population(1,j),3)+1;
        result(j,4) = viablesTriplets(population(1,j),4)+1;
        result(j,5) = (viablesValues2.col(population(1,j))).sum();
      }
      else {
        result(j,0) = -1;
        result(j,1) = -1;
        result(j,2) = -1;
        result(j,3) = -1;
      }
      
    }
  
  }
  
  
  
  //Rcpp::NumericMatrix result(6,4);
  
  //fixed list return
  List z2 = List::create( result, colonisationMatrices, viablesTriplets, initial_state, globalSuitableCoordinates) ;

  //List z2 = List::create(1);
  return(z2);
  //, currentVector
  //List z3 = List::create(1);
  //return(z3);

}

