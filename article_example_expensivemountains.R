graphics.off()
### Definitive script
#dev.off()
### Imports
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
library(hrbrthemes)
library(grid)
library(raster)


### Global Details
col = c('green3','green2','blue1','blue3')
col = terrain.colors(20)
hea = (viridis(10))

#################################################
## Initialisation & Variables
#################################################

# Map - size
nrow = 50
ncol = 50

# Map - time caracteristics
N_cycles = 50
Tadd_cyc = 1.5

# Map - space caracteristics
Tmin_alt_0 = 10
Tmax_alt_0 = 15
Tadd_alt_1 = -8
Trange = Tmax_alt_0 - (Tmin_alt_0 + Tadd_alt_1)
Tmarge = c(Tmin_alt_0 + Tadd_alt_1, Tmax_alt_0+Tadd_cyc)

# Algorithm - genetic algo caracteristics
npop = 300
nsur = 50
ngen = 20

# Algorithm - optimum condition values
threshold = 70
confidence = 0.95

# Species - caracteristics
Trange_spe = c(8.8,9.1,11.6,11.9)
migr_spe = array(0.011, c(3, 3))
migr_spe[2,2] = 0.999

xmin = 0.15
xmax = 0.32
ymin = 0.2
ymax = 0.32

#################################################
## Maps construction
#################################################

# Our cool map looking like Jupiter
mapping = function(x,y){
  return (exp(-3*abs(x-y))* (sqrt((x-0.5)^2+(y-0.5)^2)<0.45))
}
height_map = function(nrow,ncol,FUN){
  x = seq(0,1,length=nrow)
  y = seq(0,1,length=ncol)
  v = outer(x,y,FUN)
  v = (v-min(v))/(max(v)-min(v))
  v = v*(v>0.1)
  v = (v-min(v))/(max(v)-min(v))
  v[v==0] = NA
  v[x]
  return(v)
}
col = terrain.colors(10)
map = height_map(nrow,ncol,mapping)
image.nan.better(map,col=col,zlim=range(map,na.rm=T),outside.below.color='brown',outside.above.color='brown',na.color='navy')
levelplot(map,col.regions = c(terrain.colors(29)), main="Terrain simulated")

#box()

# Temperatures on the map
nr = seq(0,1,length=nrow)
nc = seq(0,1,length=ncol)
temperatures = function(x,y){
  return(-(Tmax_alt_0-Tmin_alt_0)*y+Tmax_alt_0)
}
climat_map_maker = function(height_map,Tgain,t){
  return(outer(nr,nc,FUN="temperatures")*(height_map!=0)+height_map*Tadd_alt_1+Tgain*t) 
}
ccmm1 = function(height_map,Tgain,T){
  res = list()
  for (i in 0:T){
    t = i/T
    matry = climat_map_maker(height_map,Tgain,t)
    res[[i+1]] = matry
    #res[[i+1]] = t(matry)[ncol(matry):1,]
  }
  return(res)
}


c_maps = ccmm1(map,Tadd_cyc,N_cycles)
cc0 = c_maps[[1]]
data <- expand.grid(X=nr, Y=nc)
colnames(data) = c("lontitude","latitude")

library(cowplot)

p1 <-levelplot(cc0,col.regions = c(viridis(29)), main="Temperatures simulated (2020)",at=seq(0, 16, length=30))
#ggplot(data, aes(X, Y, fill= cc0[nrow(cc0):1,])) + geom_tile() + scale_fill_gradient(low="blue", high="red",limits=Tmarge) + theme_ipsum()
cc1 = c_maps[[N_cycles]]
#ggplot(data, aes(X, Y, fill= cc1[,ncol(cc1):1])) + geom_tile() + scale_fill_gradient(low="blue", high="red",limits=Tmarge) + theme_ipsum()
p2 <-levelplot(cc1,col.regions = c(viridis(29)), main="Temperatures simulated (2070)",at=seq(0, 16, length=30))
plot_grid(p1,p2)


# Suitabilities on the map
list_suit = list()
for (t in 0:N_cycles){
  climat = c_maps[[t+1]]
  notnull = (Trange_spe[1]<climat)&(Trange_spe[4]>climat)
  loc1 = (climat-Trange_spe[1]) / (Trange_spe[2]-Trange_spe[1])
  loc2 = (climat-Trange_spe[4]) / (Trange_spe[3]-Trange_spe[4])
  suit = pmin(loc1,loc2,1)
  suit[is.na(suit) | suit<0]=0
  list_suit[[t+1]] = as(Matrix(suit,sparse=TRUE),"dgCMatrix")
}

# Plot suitabilities over time
first_suit = 0*map
for (t in 0:N_cycles){
  back_suit = as.matrix(list_suit[[(N_cycles+1)-t]])
  first_suit[back_suit>0.9]=((N_cycles+1)-(t))
}

last_suit = 0*map
for (t in 0:N_cycles){
  forw_suit = as.matrix(list_suit[[t+1]])
  last_suit[forw_suit>0.9]=t
}

R = ((last_suit-first_suit)/N_cycles)[nrow(first_suit):1,]
G = (first_suit/N_cycles)[nrow(first_suit):1,]
B = (1 - last_suit/N_cycles)[nrow(first_suit):1,]
rvb_tensor = array(0.5,c(nrow,ncol,3))
rvb_tensor[,,1]=1*(first_suit==1 & last_suit==N_cycles) # red
rvb_tensor[,,2]=1*((first_suit!=1 & last_suit==N_cycles)|(first_suit==1 & last_suit!=N_cycles)) # green
rvb_tensor[,,3]=1*(first_suit==1 & last_suit!=N_cycles) # blue

rvb_tensor[(rvb_tensor[,,1]+rvb_tensor[,,2]+rvb_tensor[,,3])==0]=0.5

#rvb_tensor[,,3]= rvb_tensor[,,3] - 1 * (R==0 & G==0 & B==1)
rvb_tensor[rvb_tensor==-1 | FALSE]=0.5
rvb_tensor[is.na(map)]=1


# Cost
cost_mapping = function(x,y){
  return (100+y*100)
}
cost_map = function(nrow,ncol,FUN){
  x = seq(0,1,length=nrow)
  y = seq(0,1,length=ncol)
  v = outer(x,y,FUN)
  return(v)
}
cost = (!is.na(map)) * (100+100*map)
ggplot(data, aes(lontitude, latitude, fill= cost)) + geom_tile() + scale_fill_gradient(low="blue", high="red",limits=c(100,200)) + theme_ipsum()
if(class(cost)=="dtCMatrix"){
  cost = as(as(cost,"dgCMatrix"),"dgCMatrix")
}else{
  cost = as(cost,"dgCMatrix")
}
cost = as(as(cost,"dgTMatrix"),"dgCMatrix")


#pres = as(pres,"dgCMatrix")
rotate <- function(x) t(apply(x, 2, rev))
# Presence
suitnow = (first_suit[nrow:1,]==1)*1

###
maty <- (matrix(runif(ncol*nrow),ncol = ncol))                # Specify number of columns
Xmin = round(xmin * ncol)
Xmax = round(xmax * ncol)
Ymin = round(ymin * nrow)
Ymax = round(ymax * nrow)
maty
maty[1:Xmin,]=0
maty[Xmax:ncol,]=0
maty[,1:Ymin]=0
maty[,Ymax:nrow]=0
###

maty <- (matrix(rep(0,ncol*nrow),ncol = ncol)) 

Xmin = round(xmin * ncol)
Xmax = round(xmax * ncol)
Ymin = round(ymin * nrow)
Ymax = round(ymax * nrow)
maty[Xmin:Xmax,Ymin:Ymax]=1

xmin = 0.75
xmax = 0.9
ymin = 0.75
ymax = 0.9
Xmin = round(xmin * ncol)
Xmax = round(xmax * ncol)
Ymin = round(ymin * nrow)
Ymax = round(ymax * nrow)

maty[Xmin:Xmax,Ymin:Ymax]=1

maty = maty * (matrix(runif(ncol*nrow),ncol = ncol))    


pres = (suitnow>0.99) * (maty>0.3)

if(sum(pres,na.rm=TRUE)==0){
  pres =  (matrix(runif(ncol*nrow),ncol = ncol))     
}


XY_pres = which(pres==1,arr.ind = T)

#a = (vt[order(vt[,4])[1:500],c(1,2)]+0.5)/nrow
#points(a[,2],a[,1],pch = 19,col='purple')

pres[is.na(pres)] = 0
pres = Matrix(pres, sparse = T)
#pres = rotate(pres)
pres = as(pres,"dgCMatrix")


XY_pres = which(pres==1,arr.ind = T)
XY_pres = cbind(XY_pres[,1],XY_pres[,2])

for (i in 1:length(XY_pres[,1])){
  rvb_tensor[nrow-XY_pres[i,1]+1,XY_pres[i,2],1]=0
  rvb_tensor[nrow-XY_pres[i,1]+1,XY_pres[i,2],2]=0
  rvb_tensor[nrow-XY_pres[i,1]+1,XY_pres[i,2],3]=0
}

rvb_tensor2 = round(rvb_tensor*255)
raster_RGB = stack(raster(rvb_tensor2[,,1]),raster(rvb_tensor2[,,2]),raster(rvb_tensor2[,,3]))
plotRGB(flip(raster_RGB))
#p3 <-plotRGB(flip(raster_RGB))

legend("topright",c("Area occupied by species","Area becomming unsuitable", "Area always suitable", "Area becomming suitable"), cex=1.0, bty="y",
       fill=c("black","cyan","red","green"),inset=.04)

#threshold = sum(pres)

gss = rcpp_global_suitable_sites(list_suit)
gsc = rcpp_global_suitable_coordinates(gss)
ltm = rcpp_local_transition_matrix(gss,gsc,migr_spe)
tm = rcpp_transition_matrices(list_suit,ltm,gsc)
cm = rcpp_colonisation_matrices(tm)

vs = rcpp_viable_sites(cm)
vt = rcpp_viable_triplets(vs,cm,gsc,gss,cost)
vv = rcpp_viable_values(vt,vs,gss,cm)
#cv = rcpp_get_current_vector(pres,cm,gss)



ph = rcpp_pheromons(vt)
ecp = rcpp_eval_current_prob(threshold,pres[nrow:1,],cm,gss)
ntp = threshold - which(cumsum(ecp)>0.95)[1] 



po = rcpp_generate_population(ph,gss,npop,ntp)
resultat = rcpp_algorithm_opt(ph,vt,po,cost,pres[nrow:1,],cm,gss,vv,threshold,confidence,npop,nsur,ngen,ntp)
choix = rcpp_result_to_choice(resultat,vt)



#####
##### Afficher la situation SAALEE (si aucune action légitime n'est entreprise)
#####

rvb_tensor = array(0.5,c(nrow,ncol,3))
for (i in 1:length(c(colo3))){
  if (colo3[i]>=0.95){
    #0.5 + colo3[i]/2
    rvb_tensor[possibilities[indices_pres[i],1],possibilities[indices_pres[i],2],1]=0.9#R
    rvb_tensor[possibilities[indices_pres[i],1],possibilities[indices_pres[i],2],2]=0.0 #G
    rvb_tensor[possibilities[indices_pres[i],1],possibilities[indices_pres[i],2],3]=0.9#B    
  }
  
}
for (j in 1:length(cm[[1]][,1])){
  rvb_tensor[possibilities[j,1],possibilities[j,2],1]= 0.6
  rvb_tensor[possibilities[j,1],possibilities[j,2],2]= 0.6
  rvb_tensor[possibilities[j,1],possibilities[j,2],3]= 0.6
}
indices_choix_ma2 = indices_choix_ma[indices_choix_ma!=0] 
for (i in 1:length(indices_pres)){
  for (j in 1:length(cm[[1]][,1])){
    #○cm[[1]][j,indices_pres[i]]
    s = cm[[1]][j,indices_pres[i]]
    if (s!=0){
      rvb_tensor[possibilities[j,1],possibilities[j,2],1]= s + rvb_tensor[possibilities[j,1],possibilities[j,2],1] -  rvb_tensor[possibilities[j,1],possibilities[j,2],1] * s
      #rvb_tensor[possibilities[j,1],possibilities[j,2],2]= 0
      rvb_tensor[possibilities[j,1],possibilities[j,2],2]= s + rvb_tensor[possibilities[j,1],possibilities[j,2],1] -  rvb_tensor[possibilities[j,1],possibilities[j,2],1] * s
    }
  }
}
indices_choix_ma2 = indices_choix_ma[indices_choix_ma!=0] 
for (i in 1:length(indices_choix_ma2)){
  for (j in 1:length(cm[[1]][,1])){
    #○cm[[1]][j,indices_pres[i]]
    s = cm[[1]][j,indices_choix_ma2[i]]
    if (s!=0){
      rvb_tensor[possibilities[j,1],possibilities[j,2],1]= s + rvb_tensor[possibilities[j,1],possibilities[j,2],1] -  rvb_tensor[possibilities[j,1],possibilities[j,2],1] * s
      #rvb_tensor[possibilities[j,1],possibilities[j,2],2]= 0
      rvb_tensor[possibilities[j,1],possibilities[j,2],3]= s + rvb_tensor[possibilities[j,1],possibilities[j,2],1] -  rvb_tensor[possibilities[j,1],possibilities[j,2],1] * s
    }
  }
}
indices_choix_ma2 = indices_choix_ma[indices_choix_ma!=0] 
for (i in 1:length(indices_choix_ma2)){
  x = possibilities[indices_choix_ma2[i],1]
  y = possibilities[indices_choix_ma2[i],2]
  rvb_tensor[x,y,1]=0.3 #R
  rvb_tensor[x,y,2]=0.3 #G
  rvb_tensor[x,y,3]=0.9 #B
}
for (i in 1:length(XY_pres[,1])){
  if (colo3[i]>0.8){
    rvb_tensor[nrow+1-XY_pres[i,1],XY_pres[i,2],1]=0.3 #R
    rvb_tensor[nrow+1-XY_pres[i,1],XY_pres[i,2],2]=0.9 #G
    rvb_tensor[nrow+1-XY_pres[i,1],XY_pres[i,2],3]=0.3 #B
  }
  if (colo3[i]<=0.8){
    rvb_tensor[nrow+1-XY_pres[i,1],XY_pres[i,2],1]=0.9 #R
    rvb_tensor[nrow+1-XY_pres[i,1],XY_pres[i,2],2]=0 #G
    rvb_tensor[nrow+1-XY_pres[i,1],XY_pres[i,2],3]=0 #B
  }
}

rvb_tensor[is.na(map)]=1
rvb_tensor2 = round(rvb_tensor*255)
raster_RGB = stack(raster(rvb_tensor2[,,1]),raster(rvb_tensor2[,,2]),raster(rvb_tensor2[,,3]))

plotRGB(flip(raster_RGB))
legend("topright",c("Surviving sites", "Naturally colonised sites", "Not surviving sites", "Planting sites AM", "Colonised sites due to AM"), cex=1.0, bty="y",
       fill=c("green","yellow","red","blue","purple"),inset=.04)