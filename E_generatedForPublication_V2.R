
library("usethis")
library("roxygen2")
library("devtools")
library(Rcpp)
library(RcppEigen)
library(methods)
library(lattice)
library(ggplot2)
library(dplyr)
library(purrr)
library(Matrix)
library("viridis")
library("optAM")

col = c('green3','green2','blue1','blue3')
col = terrain.colors(20)


#################################################
##
## General aspects
##
#################################################

nrow = 100
ncol = 100

Tmin_alt_0 = 5
Tmax_alt_0 = 10
Tadd_alt_1 = -10
marge = Tmax_alt_0 - (Tmin_alt_0 + Tadd_alt_1)


T = 7

Tgain = 7

# confidence
threshold = 20
confidence = 0.95

#GA indices
npop = 50
nsur = 35
ngen = 3

i_mp = 8
i_tg = 5


#################################################
##
## Preparing species
##
#################################################

# 1 : Tpref
# 2 : Ttol1
# 3 : Ttol2
# 4 : Abondance
# 5 : MigrationType

esps = matrix(0,3,5)
esps[1,] = c(4,0.3,0.6,0.3,1)
esps[2,] = c(5,1,1,0.3,1)
esps[3,] = c(6,1,1,0.3,1)
nb_species = length(esps[,1])

# CREATE MIGRATION KERNELS

# MK 1
migr1 = array(0.5, c(3, 3))
migr1[2, 2] = 1

# MK 2
migr2 = array(0.3, c(5, 5))
migr2[3, 3] = 1

# MK 3 
migr3 = array(0.2, c(7, 7))
migr3[4, 4] = 0.9

get_migr = function(i){
  if (i==1){
    return(migr1)
  }
  if (i==2){
    return(migr2)
  }
  if (i==3){
    return(migr3)
  }
}
#esp = esps[1,]

#################################################
##
## Preparing map
##
#################################################

ff = function(x,y,cx,cy,r){
  de = sqrt((x-cx)^2+(y-cy)^2)
  if (de>r){
    return(-1)
  }
  return(1)
}

f1 = function(x,y){
  return(-(y+47)*sin(sqrt(abs(x/2+(y+47))))-x*sin(sqrt(abs(x-(y+47)))) - 10*min(abs(x),abs(y)))
}
f2 = function(x,y){
  return( abs(abs(x-y+0.25) * abs(y+x-0.5)))
}
f3 = function(x,y){
  return((sin(y*32)))
}
f4 = function(x,y){
  return(sin(y*64))
}
f5 = function(x,y){
  return(sin(x*6)+cos(y*32)+y-x)
}
f6 = function(x,y){
  return( -abs(x+0.5*cos(y)-1)*abs(sin(x*10)-y)-0.2*x+0.2*cos(y*10) )
}
f7 = function(x,y){
  return( 5*x-y*y*y )
}
f8 = function(x,y){
  #return( min(sin(10*y+15),ff(x,y,0.5,0.5,0.25)))
  #return (sin(6*y+17.5)-((1-sin(6*y+17.5))*(sqrt((x-0.5)^2+(y-0.5)^2))>0.0625))
  #return (sin(6*y+17.5)^2  * (sqrt((x-0.5)^2+(y-0.5)^2)>0.3))
  return ((4*sin(3.14*y)^3-6*sin(3.14*y)^2+3*sin(3.14*y))^2 * (sqrt((x-0.5)^2+(y-0.5)^2)<0.45))
  #return (sin(3.14*y)^2 * (sqrt((x-0.5)^2+(y-0.5)^2)<0.45))
}
height_map = function(nrow,ncol,FUN){
  x = seq(0,1,length=nrow)
  y = seq(0,1,length=ncol)
  v = outer(x,y,FUN)
  v = (v-min(v))/(max(v)-min(v))
  v = v*(v>0.1)
  v = (v-min(v))/(max(v)-min(v))
  v[v==0] = NA
  return(v)
}
col = terrain.colors(10)
map = height_map(nrow,ncol,f8)
image.nan.better(map,col=col,zlim=range(map,na.rm=T),outside.below.color='brown',outside.above.color='brown',na.color='navy')
box()

#################################################
##
## Preparing climate adaptability
##
#################################################

nr = seq(0,1,length=nrow)
nc = seq(0,1,length=ncol)

temperatures = function(x,y){
  return(-(Tmax_alt_0-Tmin_alt_0)*y+Tmax_alt_0)
  #return((Tmax-Tmin)*y+Tmin)
}

climat_map_maker = function(height_map,Tgain,t){
  return(outer(nr,nc,FUN="temperatures")*(height_map!=0)+height_map*Tadd_alt_1+Tgain*t) 
}

# Normal evolution of the temperature
ccmm1 = function(height_map,Tgain,T){
  res = list()
  for (i in 0:T){
    t = i/T
    res[[i+1]] = climat_map_maker(height_map,Tgain,t)
  }
  return(res)
}

# Slow then quick evolution of the temperature
ccmm2 = function(height_map,Tgain,T){
  res = list()
  for (i in 0:T){
    res[[i+1]] = climat_map_maker(height_map,(i/T)^2)
  }
  return(res)
}

# Quick then slow evolution of the temperature
ccmm3 = function(height_map,Tgain,T){
  res = list()
  for (i in 0:T){
    res[[i+1]] = climat_map_maker(height_map,(i/T)^(1/2))
  }
  return(res)
}

#################################################
##
## Preparing suitability
##
#################################################

suit_maps = function(esp,height_map,Tgain,T){
  list_suit = list()
  climats = ccmm1(height_map,Tgain,T)
  for (t in 0:T){
    climat = climats[[t+1]]
    notnull = (abs(climat-esp[1])<esp[3])
    H = 1 + esp[2]/esp[3]
    locval = H - abs(climat-esp[1])/(esp[3]+esp[2])*H
    suit = pmin(locval,1)*notnull
    
    suit[is.na(suit)]=0
    
    #print(suit)
    #message("OK BUDDY")
    list_suit[[t+1]] = as(Matrix(suit,sparse=TRUE),"dgCMatrix")
    #list_suit[[t+1]] =as(as(Matrix(suit,sparse=TRUE),"dgTMatrix"),"dgCMatrix")
    #list_suit[[t]] = Matrix(suit,sparse=TRUE)
  }
  return(list_suit)
}

suits_maps = function(esps,height_map,Tgain,T){
  tot_list = list()
  for(n in 1:nb_species){
    suity = suit_maps(esps[n,],height_map,Tgain,T)
    tot_list[[n]] = suity
  }
  return(tot_list)
}

presence = function(esp,suit_maps){
  suit0 = suit_maps[[1]]
  potpres = which(suit0!=0, arr.ind=TRUE)
  nb_possibles = length(potpres[,1])
  #message(nb_possibles)
  nb_pres = round(nb_possibles*esp[4]+0.5)
  #message(nb_pres)
  xy_pres = as.matrix(potpres[sample(1:nb_possibles,nb_pres),])
  #print(xy_pres)
  v <- rep(1, nb_pres)
  return(sparseMatrix(xy_pres[,1], xy_pres[,2], x = v, dims = c(nrow, ncol)))
}

presences = function(esps,suits_maps){
  res = list()
  len = nb_species
  for (i in 1:len){
    res[[i]] = presence(esps[i,],suits_maps[[i]])
  }
  return(res)
}

cost_mat = function(suitabilities_i){
  cost = Matrix(0,nrow,ncol)
  cost = matrix(rep(1+seq(-1,1,length.out =nrow)^2,ncol),nrow = nrow , ncol = ncol)
  # for (j in 1:length(suitabilities_i)){
  #   
  #   for (i in 1:length(suitabilities_i[[j]]))
  #     xy = which(suitabilities_i[[j]][[i]]!=0,arr.ind=TRUE)
  #     cost[xy[,1],xy[,2]] = 10 + 10*runif(length(xy[,1]))
  # }
  return(Matrix(cost,sparse=TRUE))
}


c_maps = ccmm1(map,Tgain,T)
suities = suits_maps(esps,map,Tgain,T)

#i <- rep(128:133, 6)
#j <- rep(85:90, each=6)
#v <- rep(1, 36)

suitnow = suities[[1]][[1]]
maty <- matrix(runif( ncol*nrow),ncol = ncol)                # Specify number of columns

suitnow>0.9
pres = (suitnow>0.9) * (maty>0.8)

#pres = sparseMatrix(i, j, x = v, dims = c(nrow, ncol))


#pres = presence(esps,suities)

#cost = as(cost_mat(suitabilities[[i]]),"dgCMatrix")
cost = (cost_mat(suities))
if(class(cost)=="dtCMatrix"){
  cost = as(as(cost,"dgCMatrix"),"dgCMatrix")
}else{
  cost = as(cost,"dgCMatrix")
}
cost = as(as(cost,"dgTMatrix"),"dgCMatrix")

if(class(pres)=="dtTMatrix"){
  pres = as(as(pres,"dgTMatrix"),"dgCMatrix")
}else{
  pres = as(pres,"dgCMatrix")
}

pres_list = list()
AAA = list()

#Sys.sleep(1)

for (i in 1:nb_species){
  
  message("Let me know")
  
  suitnow = suities[[i]][[1]]
  maty <- matrix(runif( ncol*nrow),ncol = ncol)                # Specify number of columns
  maty <- t(t((maty) * seq(0,1,length=nrow)^4) * (seq(0,1,length=nrow)>0.1))
  
  suitnow>0.9
  pres = (suitnow>0.9) * (maty>0.2)
  
  
  
  #pres = sparseMatrix(i, j, x = v, dims = c(nrow, ncol))
  
  
  #pres = presence(esps,suities)
  
  #cost = as(cost_mat(suitabilities[[i]]),"dgCMatrix")
  cost = (cost_mat(suities))
  if(class(cost)=="dtCMatrix"){
    cost = as(as(cost,"dgCMatrix"),"dgCMatrix")
  }else{
    cost = as(cost,"dgCMatrix")
  }
  cost = as(as(cost,"dgTMatrix"),"dgCMatrix")
  
  if(class(pres)=="dtTMatrix"){
    pres = as(as(pres,"dgTMatrix"),"dgCMatrix")
  }else{
    pres = as(pres,"dgCMatrix")
  }
  
  pres_list[[i]] = pres
  
  
  lis = suities[[i]]
  for (k in 1:length(lis)){
    if(class(lis[[k]])=="dtTMatrix"){
      matt = as(as(lis[[k]],"dgTMatrix"),"dgCMatrix")
    }else{
      matt = as(lis[[k]],"dgCMatrix")
      
    }
    lis[[k]] = matt
  }
  
  
  message("Bring me a cookie...")
  
  if(sum(pres)>1){
    
    
    message('ALIVE')
    
    message("WHat's the class ?")
    print(class(pres))
    #Sys.sleep(1)
    AAA1=optgenam3(pres,
                   lis,
                   cost,
                   get_migr(esps[i,5]),
                   threshold,
                   confidence,
                   npop,
                   nsur,
                   ngen)
    AAA[[i]] = AAA1
    message('Or is it ..?')
  }
  
  else{
    message('DEAD')
    AAA[[i]] = NA
  }
}


i = 1
pres = pres_list[[i]]
AAA1=AAA[[i]]
chx = AAA1[[1]]
col = AAA1[[2]]
trp = AAA1[[3]]

sususu = suities[[i]][[1]]
su = which(sususu!=0,arr.ind=TRUE)
sususu2 = suities[[i]][[T]]
su2 = which(sususu2!=0,arr.ind=TRUE)

colors = terrain.colors(10)
hea = (viridis(10))

par(mfrow=c(2,2))
HM = map
HM[su] = -100
HM[which(pres!=0,arr.ind=TRUE)] = 100
image.nan.better(HM, col =colors,zlim=c(0,1),outside.below.color='red',outside.above.color='cyan',na.color='navy')
title("Carte de présence de l'espèce en 2021")

cc0 = c_maps[[1]]
image.nan.better(cc0, col =hea,zlim=range(-5,15),outside.below.color='brown',outside.above.color='grey',na.color='navy')
title("Carte des températures en 2021")

HM = map
HM[su2] = -100
image.nan.better(HM, col =colors,zlim=c(0,1),outside.below.color='red',outside.above.color='brown',na.color='navy')
title("Carte prédite de présence de l'espèce en 2061")

ccN = c_maps[[T]]
image.nan.better(ccN, col =hea,zlim=range(-5,15),outside.below.color='brown',outside.above.color='grey',na.color='navy')
title("Carte des températures en 2061")












par(mfrow=c(2,2))

chx = AAA1[[1]]
col = AAA1[[2]]
trp = AAA1[[3]]

t0 = trp[trp[,3]!=-1,c(1,2,5)]
colo = AAA1[[2]]
colo1 = colo[[1]]
coloN = colo[[T]]
resu0 = Matrix(NA,NROW,NCOL)
resuN = Matrix(NA,NROW,NCOL)

for (j in (1:length(t0[,1]))){
  message(j)
  resu0[t0[j,1],t0[j,2]] = sum(colo1[,t0[j,3]])
  resuN[t0[j,1],t0[j,2]] = sum(coloN[,t0[j,3]])
}

val0 = resu0[which(resu0>=0,arr.ind=TRUE)]
valN = resuN[which(resuN>=0,arr.ind=TRUE)]


HM1 = map*0
HM1[which(resu0>=0,arr.ind=TRUE)] = val0
HM2 = map*0
HM2[which(resuN>=0,arr.ind=TRUE)] = valN

image.nan.better(as.matrix(HM1),col=colors,zlim=c(-5,20),outside.below.color='red',outside.above.color='brown',na.color="navy")
title("Zones proposées pour la migration assistée (2021)")


image.nan.better(as.matrix(HM2),col=colors,zlim=c(-5,20),outside.below.color='brown',outside.above.color='brown',na.color="navy")
title("Zones proposées pour la migration assistée (2061)")




choix= chx
choix = choix[choix[,1]!=-1,c(1,2)]
KK=is.na(HM)*1.0

KK[su]=KK[su]+0.2
KK[su2]=KK[su2]+0.4

KK[choix]=NA

HM3 = HM2
HM3[choix] = 10000

image.nan.better(as.matrix(HM3),col=colors,zlim=c(-5,20),outside.below.color='red',outside.above.color='brown',na.color="navy")
title("Choix optimaux cibles pour la migration assistés")

