
library("usethis")
library("roxygen2")
library("devtools")
library(Rcpp)
library(RcppEigen)
library(methods)



Rcpp::compileAttributes()
devtools::document()
build()
install()


library(Matrix)
library(optAM)
#sourceCpp("src/functions.cpp")





n = 1000
t = 3
d = 30
off = 1
k = 3




trying = function(t,d,off){
  message("Ready when you are t=",t," & d=",d)
  k = 3
  n = d+off*(t-1)
  suitableMatrices = list()
  h = 1
  i <- rep(((h) * off):((h) * off + d - 1), d)
  j <- rep(((h) * off):((h) * off + d - 1), each = d)
  v <- rep(1, d * d)
  
  costMatrix <- sparseMatrix(i, j, x = v, dims = c(n, n))
  
  
  presenceMatrix <- sparseMatrix(i[1:k], j[1:k], x = rep(1, k), dims = c(n, n))
  
  
  for (h in 1:t) {
    i <- rep(((h) * off):((h) * off + d - 1), d)
    j <- rep(((h) * off):((h) * off + d - 1), each = d)
    v <- rep(1, d * d)
    suitableMatrices[[h]]  <- sparseMatrix(i, j, x = v, dims = c(n, n))
    
    
    v <- rep(h, d * d)
    costMatrix = costMatrix + sparseMatrix(i, j, x = v, dims = c(n, n))
    
  }
  
  migr = array(0, c(3, 3))
  kk = 0.7
  migr[2, 2] = 1
  migr[1, 2] = kk
  migr[2, 1] = kk
  migr[3, 2] = kk
  migr[2, 3] = kk
  
  # confidence
  threshold = 8
  confidence = 0.9
  
  #GA indices
  npop = 100
  nsur = 70
  ngen = 3
  
  T1<-Sys.time()
  
  AAA=optgenam3(presenceMatrix,
                suitableMatrices,
                costMatrix,
                migr,
                threshold,
                confidence,
                npop,
                nsur,
                ngen)
  #print(AAA)
  #print(sum(AAA[[1]][,4]))
  
  T2<-Sys.time()
  return(difftime(T2,T1))
  
}

trying2 = function(t,d,off){
  message("Ready when you are t=",t," & d=",d)
  k = 3
  n = d+off*(t-1)
  suitableMatrices = list()
  h = 1
  i <- rep(((h) * off):((h) * off + d - 1), d)
  j <- rep(((h) * off):((h) * off + d - 1), each = d)
  v <- rep(1, d * d)
  
  costMatrix <- sparseMatrix(i, j, x = v, dims = c(n, n))
  
  
  presenceMatrix <- sparseMatrix(i[1:k], j[1:k], x = rep(1, k), dims = c(n, n))
  
  
  for (h in 1:t) {
    i <- rep(((h) * off):((h) * off + d - 1), d)
    j <- rep(((h) * off):((h) * off + d - 1), each = d)
    v <- rep(1, d * d)
    suitableMatrices[[h]]  <- sparseMatrix(i, j, x = v, dims = c(n, n))
    
    
    v <- rep(h, d * d)
    costMatrix = costMatrix + sparseMatrix(i, j, x = v, dims = c(n, n))
    
  }
  
  migr = array(0, c(3, 3))
  kk = 0.7
  migr[2, 2] = 1
  migr[1, 2] = kk
  migr[2, 1] = kk
  migr[3, 2] = kk
  migr[2, 3] = kk
  
  # confidence
  threshold = 8
  confidence = 0.9
  
  #GA indices
  npop = 800
  nsur = 600
  ngen = 5
  
  T1<-Sys.time()
  

  
  AAA=optgenam3(presenceMatrix,
                suitableMatrices,
                costMatrix,
                migr,
                threshold,
                confidence,
                npop,
                nsur,
                ngen)
  #print(AAA)
  #print(sum(AAA[[1]][,4]))
  
  T2<-Sys.time()
  return(AAA)
  
}











AAA=trying2(5,3,1)

d = seq(20,80,by=10)
t = seq(20,70,by=10)

resu = array(0,c(length(d),length(t)))

for (a in 1:length(t)){
  for (b in 1:length(d)){
    resu[b,a] = trying(t[a],d[b],1)
  }
}

library(lattice)

levelplot((resu))



