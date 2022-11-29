library(lattice)
library(ggplot2)
#library(tidyverse)

library(dplyr)
library(purrr)
library("viridis")

nrow = NROW
ncol = NCOL

col = c('green3','green2','blue1','blue3')
col = terrain.colors(20)


###############################################################################
###
###
### ................................... Carte test 1
###
###
###############################################################################

f1 = function(x,y){
  return(-(y+47)*sin(sqrt(abs(x/2+(y+47))))-x*sin(sqrt(abs(x-(y+47)))) - 10*min(abs(x),abs(y)))
}

height_map1 = function(nrow,ncol){
  
  nr = seq(0,1,length=nrow)
  nc = seq(0,1,length=ncol)
  
  x = 1000*nr-500
  
  y = 1000*nc-500
  
  v = outer(x,y,f1)
  v = v*(v>0.1)
  v = (v-min(v))/(max(v)-min(v))
  
  v[v==0] = NA
  
  return(v)
}

map1 = (height_map1(nrow,ncol))
#levelplot(map1)

image.nan.better(map1,col=col,zlim=range(map1,na.rm=T),outside.below.color='red',outside.above.color='brown',na.color='navy')
box()






###############################################################################
###
###
### ................................... Carte test 2
###
###
###############################################################################

f2 = function(x,y){
  return( abs(abs(x-y+0.25) * abs(y+x-0.5)))
}

height_map2 = function(nrow,ncol){
  
  x = seq(0,1,length=nrow)
  y = seq(0,1,length=ncol)
  
  v = outer(x,y,f2)
  v = (v-min(v))/(max(v)-min(v))
  v = v*(v>0.1)
  v = (v-min(v))/(max(v)-min(v))
  
  v[v==0] = NA
  
  return(v)
}

map2 = height_map2(nrow,ncol)
#levelplot(map2)
image.nan.better(map2,col=col,zlim=range(map2,na.rm=T),outside.below.color='red',outside.above.color='brown',na.color='navy')
box()



###############################################################################
###
###
### ................................... Carte test 3
###
###
###############################################################################

f3 = function(x,y){
  return((sin(y*32)))
}

height_map3 = function(nrow,ncol){
  
  x = seq(0,1,length=nrow)
  y = seq(0,1,length=ncol)
  
  v = outer(x,y,f3)
  v = (v-min(v))/(max(v)-min(v))
  v = v*(v>0.1)
  v = (v-min(v))/(max(v)-min(v))
  
  v[v==0] = NA
  
  return(v)
}

map3 = height_map3(nrow,ncol)
#levelplot(map3)
image.nan.better(map3,col=col,zlim=range(map3,na.rm=T),outside.below.color='brown',outside.above.color='brown',na.color='navy')
box()


###############################################################################
###
###
### ................................... Carte test 4
###
###
###############################################################################

f4 = function(x,y){
  return(sin(y*64))
}

height_map4 = function(nrow,ncol){
  
  x = seq(0,1,length=nrow)
  y = seq(0,1,length=ncol)
  
  v = outer(x,y,f4)
  v = (v-min(v))/(max(v)-min(v))
  v = v*(v>0.1)
  v = (v-min(v))/(max(v)-min(v))
  
  v[v==0] = NA
  
  return(v)
}

map4 = height_map4(nrow,ncol)
#levelplot(map4)
image.nan.better(map4,col=col,zlim=range(map4,na.rm=T),outside.below.color='brown',outside.above.color='brown',na.color='navy')
box()



###############################################################################
###
###
### ................................... Carte test 5
###
###
###############################################################################

f5 = function(x,y){
  return(sin(x*6)+cos(y*32)+y-x)
}

height_map5 = function(nrow,ncol){
  
  x = seq(0,1,length=nrow)
  y = seq(0,1,length=ncol)
  
  v = outer(x,y,f5)
  v = (v-min(v))/(max(v)-min(v))
  v = v*(v>0.2)
  v = (v-min(v))/(max(v)-min(v))
  
  v[v==0] = NA
  
  return(v)
}

map5 = height_map5(nrow,ncol)
#levelplot(map5)
image.nan.better(map5,col=col,zlim=range(map5,na.rm=T),outside.below.color='brown',outside.above.color='brown',na.color='navy')
box()



###############################################################################
###
###
### ................................... Carte test 6
###
###
###############################################################################

f6 = function(x,y){
  return( -abs(x+0.5*cos(y)-1)*abs(sin(x*10)-y)-0.2*x+0.2*cos(y*10) )
}

height_map6 = function(nrow,ncol){
  
  x = seq(0,1,length=nrow)
  y = seq(0,1,length=ncol)
  
  v = outer(x,y,f6)
  v = (v-min(v))/(max(v)-min(v))
  v = v*(v>0.51)
  v = (v-min(v))/(max(v)-min(v))
  
  v[v==0] = NA
  
  return(v)
}

map6 = height_map6(nrow,ncol)
#levelplot(map6)
image.nan.better(map6,col=col,zlim=range(map6,na.rm=T),outside.below.color='brown',outside.above.color='brown',na.color='navy')
box()

h_maps = list(map1,map2,map3,map4,map5,map6)


###############################################################################
###
###
### ................................... Carte test 7
###
###
###############################################################################

f6 = function(x,y){
  return( 5*x-y*y*y )
}

height_map7 = function(nrow,ncol){
  
  x = seq(0,1,length=nrow)
  y = seq(0,1,length=ncol)
  
  v = outer(x,y,f6)
  v = (v-min(v))/(max(v)-min(v))
  v = v*(v>0.51)
  v = (v-min(v))/(max(v)-min(v))
  
  v[v==0] = NA
  
  return(v)
}

map7 = height_map6(nrow,ncol)
#levelplot(map6)
image.nan.better(map7,col=col,zlim=range(map7,na.rm=T),outside.below.color='brown',outside.above.color='brown',na.color='navy')
box()

h_maps = list(map1,map2,map3,map4,map5,map6,map7)
















