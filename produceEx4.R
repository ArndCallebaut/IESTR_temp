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
N_cycles = 30
Tadd_cyc = 1.5

# Map - space caracteristics
Tmin_alt_0 = 10
Tmax_alt_0 = 15
Tadd_alt_1 = -8
Trange = Tmax_alt_0 - (Tmin_alt_0 + Tadd_alt_1)
Tmarge = c(Tmin_alt_0 + Tadd_alt_1, Tmax_alt_0+Tadd_cyc)

# Algorithm - genetic algo caracteristics
npop = 1000
nsur = 200
ngen = 30

# Algorithm - optimum condition values
threshold = 50
confidence = 0.95

# Species - caracteristics
Trange_spe = c(12.0,12.5,18.6,18.9)
migr_spe = array(0.03, c(3, 3))
migr_spe[2,2] = 0.99

xmin = 0.18
xmax = 0.25
ymin = 0.18
ymax = 0.25

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
#image.nan.better(map,col=col,zlim=range(map,na.rm=T),outside.below.color='brown',outside.above.color='brown',na.color='navy')
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
#ggplot(data, aes(X, Y, fill= cc0[nrow(cc0):1,])) + geom_tile() + scale_fill_gradient(low="blue", high="red",limits=Tmarge) + theme_ipsum()
cc1 = c_maps[[N_cycles]]
#ggplot(data, aes(X, Y, fill= cc1[,ncol(cc1):1])) + geom_tile() + scale_fill_gradient(low="blue", high="red",limits=Tmarge) + theme_ipsum()

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
  return (100+y*50)
}
cost_map = function(nrow,ncol,FUN){
  x = seq(0,1,length=nrow)
  y = seq(0,1,length=ncol)
  v = outer(x,y,FUN)
  return(v)
}
cost = (!is.na(map)) * cost_map(nrow,ncol,cost_mapping)
#ggplot(data, aes(X, Y, fill= cost)) + geom_tile() + scale_fill_gradient(low="blue", high="red",limits=c(100,150)) + theme_ipsum()
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
maty <- (matrix(runif(ncol*nrow),ncol = ncol))                # Specify number of columns

Xmin = round(xmin * ncol)
Xmax = round(xmax * ncol)
Ymin = round(ymin * nrow)
Ymax = round(ymax * nrow)
maty[1:Xmin,]=0
maty[Xmax:ncol,]=0
maty[,1:Ymin]=0
maty[,Ymax:nrow]=0

pres = (suitnow>0.99) * (maty>0.3)

pres = maty>0.1

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

#for (i in 1:1000){
#  generate_permutation4(100,ph)
#}

index_random_choice_non_uniform(ph)



cm0 = cm[[1]]
sp = apply(cm0,FUN = sum,2)



choix_ma = choix
choix_ma = choix_ma[(choix_ma[,1]!=-1),]

XY_ma = (choix[,c(1,2)]-1/2)/nrow
#points(XY_ma[,2], XY_ma[,1],col="white",pch=19)

possibilities = gsc
possibilities = possibilities +1

XY_ma2 = (possibilities[,c(1,2)]-1/2)/nrow
#points(XY_ma2[,2], XY_ma2[,1],col="brown",pch=19)


colonisation_matrices = cm

choix_ma_xy = cbind(choix_ma[,1],choix_ma[,2]) 
indices_choix_ma = matrix(0,length(choix_ma_xy[,1]))

for (i in (1:length(choix_ma_xy[,1]))){
  h = NA
  for (j in (1:length(possibilities[,1]))){
    if ((possibilities[j,1]==(choix_ma_xy[i,1]))&(possibilities[j,2]==(choix_ma_xy[i,2]))){
      h = possibilities[j,3]
      break
    }
  }
  if (!is.na(h)){
    indices_choix_ma[i] = possibilities[h,3]
  }
  
}

print(indices_choix_ma)








XY_pres = which(pres==1,arr.ind = T)
XY_pres = cbind(XY_pres[,1],XY_pres[,2])

indices_pres = matrix(0,length(XY_pres[,1]))
#XY_pres = XY_pres[nrow:1,]
for (i in (1:length(XY_pres[,1]))){
  h = NA
  for (j in (1:length(possibilities[,1]))){
    if ((possibilities[j,1]==(nrow+1-XY_pres[i,1]))&(possibilities[j,2]==(XY_pres[i,2]))){
      h = possibilities[j,3]
      break
    }
  }
  if (!is.na(h)){
    indices_pres[i] = possibilities[h,3]
  }
  
}

print(indices_pres)









#groa = sparseMatrix(i = (a@i[1:length(a@p)-1]+1), j = a@p[1:length(a@p)-1]+1, x = a@x[1:length(a@p)-1], dims = a@Dim)
coloni = as.matrix(cm[[1]])[,indices_pres]
colo2 = 0.99999999 - coloni
colo3 = round(1-apply(colo2, FUN= prod,2),2)
#colo3 = 1 - prod(colo2)

#####
##### Afficher la situation SAALEE (si aucune action légitime n'est entreprise)
#####

rvb_tensor = array(0.5,c(nrow,ncol,3))

for (i in 1:length(c(colo3))){
  if (colo3[i]>=0.95){
    #rvb_tensor[possibilities[indices_pres[i],1],possibilities[indices_pres[i],2],1]=0.5 + colo3[i]/2#R
    rvb_tensor[possibilities[indices_pres[i],1],possibilities[indices_pres[i],2],2]=1 #G
    #rvb_tensor[possibilities[indices_pres[i],1],possibilities[indices_pres[i],2],3]=0#B    
  }
  
}


for (i in 1:length(indices_pres)){
  for (j in 1:length(cm[[1]][,1])){
    rvb_tensor[possibilities[j,1],possibilities[j,2],2]=cm[[1]][j,indices_pres[i]]
  }
}

for (i in 1:length(XY_pres[,1])){
  
  if (cm[[1]][indices_pres[i],indices_pres[i]]>0.5){
    rvb_tensor[nrow+1-XY_pres[i,1],XY_pres[i,2],1]=0 #R
    rvb_tensor[nrow+1-XY_pres[i,1],XY_pres[i,2],2]=0 #G
    rvb_tensor[nrow+1-XY_pres[i,1],XY_pres[i,2],3]=0.9 #B
  }
  
  if (cm[[1]][indices_pres[i],indices_pres[i]]<=0.5){
    rvb_tensor[nrow+1-XY_pres[i,1],XY_pres[i,2],1]=0.9 #R
    rvb_tensor[nrow+1-XY_pres[i,1],XY_pres[i,2],2]=0.9 #G
    rvb_tensor[nrow+1-XY_pres[i,1],XY_pres[i,2],3]=0 #B
  }
}

# for (i in 1:length(XY_pres[,1])){
#   rvb_tensor[nrow-XY_pres[i,1]+1,XY_pres[i,2],1]=0 #R
#   rvb_tensor[nrow-XY_pres[i,1]+1,XY_pres[i,2],2]=1 #G
#   rvb_tensor[nrow-XY_pres[i,1]+1,XY_pres[i,2],3]=0 #B
# }


rvb_tensor[is.na(map)]=1
rvb_tensor2 = round(rvb_tensor*255)
raster_RGB = stack(raster(rvb_tensor2[,,1]),raster(rvb_tensor2[,,2]),raster(rvb_tensor2[,,3]))

#graphics.off()
plotRGB(flip(raster_RGB))

legend("topright",c("Not surviving sites", "Surviving sites", "Colonised sites", "Area not occupied"), cex=1.0, bty="y",
       fill=c("yellow","blue","green","purple"),inset=.04)











indices_pres2 = matrix(0,N_cycles,length(gsc[,1]))
for (i in 1:N_cycles){
  cm_loc = cm[[i]]
  for (j in 1:length(gsc[,1])){
    indices_pres2[i,j] = sum(cm_loc[,j])
  }
}

opt_each_site = apply(indices_pres2,FUN=max,2)

library("ramify")
#arg_opt_site = apply(indices_pres,FUN=argmax,2)


efficienty = opt_each_site / cost[gsc[,c(1,2)]+1]
eff = efficienty / max(efficienty)
#eff = (log(efficienty)-min(log(efficienty)))/(max(log(efficienty))-min(log(efficienty)))

eff_tensor = array(0.5,c(nrow,ncol,3))

for (i in 1:length(efficienty)){
  #rvb_tensor[nrow+1-XY_pres[i,1],XY_pres[i,2],1]=1 #R
  eff_tensor[gsc[i,1]+1,gsc[i,2]+1,2]= eff[i] #G
  eff_tensor[gsc[i,1]+1,gsc[i,2]+1,3]= 1-eff[i] #B
}

# for (i in 1:length(XY_pres[,1])){
#   rvb_tensor[nrow-XY_pres[i,1]+1,XY_pres[i,2],1]=0 #R
#   rvb_tensor[nrow-XY_pres[i,1]+1,XY_pres[i,2],2]=1 #G
#   rvb_tensor[nrow-XY_pres[i,1]+1,XY_pres[i,2],3]=0 #B
# }

eff_tensor [is.na(map)]=1
eff_tensor2 = round(eff_tensor*255)
raster_RGB = stack(raster(eff_tensor2[,,1]),raster(eff_tensor2[,,2]),raster(eff_tensor2[,,3]))

#graphics.off()
plotRGB(flip(raster_RGB))


legend("topright",c("High effect ratio", "Low effect ratio"), cex=1.0, bty="y",
       fill=c("green","purple"),inset=.04)




















diag_cms = matrix(0,N_cycles,length(gsc[,1]))
for (i in 1:N_cycles){
  cm_loc = cm[[i]]
  diag_cms[i,] = diag(cm_loc)
}

maxi_diags = apply(diag_cms,FUN=max,2)

#plot(maxi_diags,opt_each_site)

### ### ###
efficienty3 = maxi_diags / cost[gsc[,c(1,2)]+1]
eff3 = efficienty3 / max(efficienty3)
#eff = (log(efficienty)-min(log(efficienty)))/(max(log(efficienty))-min(log(efficienty)))

eff_tensor3 = array(0.5,c(nrow,ncol,3))

for (i in 1:length(efficienty3)){
  #rvb_tensor[nrow+1-XY_pres[i,1],XY_pres[i,2],1]=1 #R
  eff_tensor3[gsc[i,1]+1,gsc[i,2]+1,2]= eff3[i]  #G
  eff_tensor3[gsc[i,1]+1,gsc[i,2]+1,3]= 1-eff3[i]#B
}

eff_tensor3 [is.na(map)]=1
eff_tensor3 = round(eff_tensor3*255)
raster_RGB = stack(raster(eff_tensor3[,,1]),raster(eff_tensor3[,,2]),raster(eff_tensor3[,,3]))
#graphics.off()
plotRGB(flip(raster_RGB))
legend("topright",c("High effect ratio", "Low effect ratio"), cex=1.0, bty="y",
       fill=c("green","purple"),inset=.04)






eff_tensor4 = array(0.5,c(nrow,ncol,3))
for (i in 1:length(efficienty3)){
  eff_tensor4[gsc[i,1]+1,gsc[i,2]+1,1]=max(eff3[i]-eff[i],0)  #R
  eff_tensor4[gsc[i,1]+1,gsc[i,2]+1,2]=max(0,eff[i]-eff3[i])  #G
  eff_tensor4[gsc[i,1]+1,gsc[i,2]+1,3]=max(eff3[i]-eff[i],0)  #B
}

eff_tensor4 [is.na(map)]=1
eff_tensor4 = round(eff_tensor4*255)
raster_RGB = stack(raster(eff_tensor4[,,1]),raster(eff_tensor4[,,2]),raster(eff_tensor4[,,3]))
#graphics.off()
plotRGB(flip(raster_RGB))


legend("topright",c("Ration higher", "Ratio lower"), cex=1.0, bty="y",
       fill=c("green","purple"),inset=.04)











#####
##### Afficher les endroits où on plante
#####
rvb_tensor3 = array(0.5,c(nrow,ncol,3))

for (i in 1:length(indices_pres)){
  for (j in 1:length(cm[[1]][,1])){
    rvb_tensor3[possibilities[j,1],possibilities[j,2],1]=cm[[1]][j,indices_pres[i]]
    rvb_tensor3[possibilities[j,1],possibilities[j,2],2]=1 #G
  }
}


for (i in 1:length(choix_ma[,1])){
  
  if (choix[i,3]!=0){
    
    loc_colo = as.matrix(cm[[choix[i,3]]])
    ind=indices_choix_ma[i]
    vect_colo = loc_colo[,ind]
    for (j in 1:length(vect_colo)){
      rvb_tensor3[possibilities[j,1],possibilities[j,2],1]=vect_colo[j]+rvb_tensor3[possibilities[j,1],possibilities[j,2],1]-rvb_tensor3[possibilities[j,1],possibilities[j,2],1]*vect_colo[j] #R
      #rvb_tensor3[possibilities[j,1],possibilities[j,2],2]= #G
      #rvb_tensor3[possibilities[j,1],possibilities[j,2],3]=0 #B
    }
    
    rvb_tensor3[choix[i,1],choix[i,2],1]=1
    rvb_tensor3[choix[i,1],choix[i,2],2]=0
    rvb_tensor3[choix[i,1],choix[i,2],3]=0
  }
  
  #rvb_tensor3[possibilities[ind,1],possibilities[ind,2],1]=loc #R
}
# for (i in 1:length(choix_ma[,1])){
#   rvb_tensor3[choix_ma[i,1],choix_ma[i,2],1]=1 #R
#   rvb_tensor3[choix_ma[i,1],choix_ma[i,2],2]=0 #G
#   rvb_tensor3[choix_ma[i,1],choix_ma[i,2],3]=0 #B
# }

rvb_tensor3[is.na(map)]=1
rvb_tensor2 = round(rvb_tensor3*255)
raster_RGB = stack(raster(rvb_tensor2[,,1]),raster(rvb_tensor2[,,2]),raster(rvb_tensor2[,,3]))

plotRGB(flip(raster_RGB))
legend("topright",c("Presence predicted", "Absence predicted","Planting site"), cex=1.0, bty="y",
       fill=c("yellow","green","red"),inset=.04)

#points(XY_ma[,2], XY_ma[,1],col="red",pch=19)





















arg_opt_each_site = apply(indices_pres2,FUN=which.max,2)
moment = (arg_opt_each_site-1)/(N_cycles-1)
eff_tensor5 = array(0.5,c(nrow,ncol,3))

for (i in 1:length(efficienty)){
  eff_tensor5[gsc[i,1]+1,gsc[i,2]+1,1]=0 #R
  eff_tensor5[gsc[i,1]+1,gsc[i,2]+1,2]=1-moment[i]  #G
  eff_tensor5[gsc[i,1]+1,gsc[i,2]+1,3]=moment[i]  #B
}

for(i in 1:length(efficienty)){
  if((opt_each_site)[i]<0.5){
    eff_tensor5[gsc[i,1]+1,gsc[i,2]+1,1]= 0.45#B
    eff_tensor5[gsc[i,1]+1,gsc[i,2]+1,2]= 0.45  #G
    eff_tensor5[gsc[i,1]+1,gsc[i,2]+1,3]= 0.45  #B
  }
}

eff_tensor5 [is.na(map)]=1
eff_tensor2 = round(eff_tensor5*255)
raster_RGB = stack(raster(eff_tensor2[,,1]),raster(eff_tensor2[,,2]),raster(eff_tensor2[,,3]))
#graphics.off()
plotRGB(flip(raster_RGB))

legend("topright",c("Earlier", "In the period","Later"), cex=1.0, bty="y",
       fill=c("green","darkseagreen","blue"),inset=.04)



arg_opt_diag = apply(diag_cms,FUN=which.max,2)
moment2 = (arg_opt_diag-1)/(N_cycles-1)
eff_tensor6 = array(0.5,c(nrow,ncol,3))

for (i in 1:length(efficienty)){
  eff_tensor6[gsc[i,1]+1,gsc[i,2]+1,1]=0 #R
  eff_tensor6[gsc[i,1]+1,gsc[i,2]+1,2]=1-moment2[i]  #G
  eff_tensor6[gsc[i,1]+1,gsc[i,2]+1,3]=moment2[i]  #B
}

for(i in 1:length(efficienty)){
  if(maxi_diags[i]<0.5){
    eff_tensor6[gsc[i,1]+1,gsc[i,2]+1,1]= 0.45#B
    eff_tensor6[gsc[i,1]+1,gsc[i,2]+1,2]= 0.45  #G
    eff_tensor6[gsc[i,1]+1,gsc[i,2]+1,3]= 0.45  #B
  }
}

eff_tensor6 [is.na(map)]=1
eff_tensor2 = round(eff_tensor6*255)
raster_RGB = stack(raster(eff_tensor2[,,1]),raster(eff_tensor2[,,2]),raster(eff_tensor2[,,3]))
#graphics.off()
plotRGB(flip(raster_RGB))

legend("topright",c("Earlier", "In the period","Later"), cex=1.0, bty="y",
       fill=c("green","darkseagreen","blue"),inset=.04)








eff_tensor7 = array(0.5,c(nrow,ncol,3))

delta_moment = moment - moment2
mod_delta = delta_moment -min(delta_moment)
mod_delta = mod_delta / max(mod_delta)

for (i in 1:length(efficienty)){
  eff_tensor7[gsc[i,1]+1,gsc[i,2]+1,1]=0
  eff_tensor7[gsc[i,1]+1,gsc[i,2]+1,2]=1-mod_delta[i]  #G
  eff_tensor7[gsc[i,1]+1,gsc[i,2]+1,3]=0 #B
}

for(i in 1:length(efficienty)){
  if(maxi_diags[i]<0.5){
    eff_tensor7[gsc[i,1]+1,gsc[i,2]+1,1]= 0.45#B
    eff_tensor7[gsc[i,1]+1,gsc[i,2]+1,2]= 0.45  #G
    eff_tensor7[gsc[i,1]+1,gsc[i,2]+1,3]= 0.45  #B
  }
}
eff_tensor7 [is.na(map)]=1
eff_tensor2 = round(eff_tensor7*255)
raster_RGB = stack(raster(eff_tensor2[,,1]),raster(eff_tensor2[,,2]),raster(eff_tensor2[,,3]))
#graphics.off()
plotRGB(flip(raster_RGB))
legend("topright",c("Optimal is sooner", "Optimal is unchanged"), cex=1.0, bty="y",
       fill=c("green","darkgreen"),inset=.04)





ok_site = apply(indices_pres2,FUN=function(x) { sum(x>(max(x)/2))/length(x) },2)
eff_tensor5 = array(0.5,c(nrow,ncol,3))

for (i in 1:length(efficienty)){
  eff_tensor5[gsc[i,1]+1,gsc[i,2]+1,1]=0 #R
  eff_tensor5[gsc[i,1]+1,gsc[i,2]+1,2]= ok_site[i]  #G
  eff_tensor5[gsc[i,1]+1,gsc[i,2]+1,3]= 0 #B
}

eff_tensor5[is.na(map)]=1
eff_tensor2 = round(eff_tensor5*255)
raster_RGB = stack(raster(eff_tensor2[,,1]),raster(eff_tensor2[,,2]),raster(eff_tensor2[,,3]))
#graphics.off()
plotRGB(flip(raster_RGB))
legend("topright",
       c("Low importance","Medium importance","High importance"), 
       cex=1.0, bty="y",
       fill=c("green","darkgreen","black"),inset=.04)
