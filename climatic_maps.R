
library(lattice)
library(ggplot2)

nrow = NROW
ncol = NCOL

nr = seq(0,1,length=nrow)
nc = seq(0,1,length=ncol)

T = PERIOD
Tgain = TGAIN

###############################################################################
###
###
### ................................... Cartes climatiques, généralités
###
###
###############################################################################

# It's colder on north
# Temperatures "on the sea level"
Tmin = TMIN
Tmax = TMAX

# It's colder on top of the mountains
# Temperature lost at a "1" montain
Tlost = TLOST

temperatures = function(x,y){
  
  return(-(Tmax-Tmin)*y+Tmax)
  #return((Tmax-Tmin)*y+Tmin)
}

climat_map_maker = function(height_map,Tgain){
  return(outer(nr,nc,FUN="temperatures")*(height_map!=0)-height_map*Tlost+Tgain) 
}

climat_map_maker2 = function(i_map,i_temp,t){
  return(climat_map_maker(h_maps[[i_map]],Tgain[[i_temp]]*(t)))
}

ccmm = function(i_map,i_temp){
  res = list()
  for (i in 1:T){
    res[[i]] = climat_map_maker2(i_map,i_temp,i)
  }
}



###############################################################################
###
###
### ................................... Climat Senario 1, 2, 3
###
###
###############################################################################

# Consecutive Climat 

# Normal evolution of the temperature
ccmm1 = function(i_map,i_temp){
  res = list()
  for (i in 1:T){
    res[[i]] = climat_map_maker2(i_map,i_temp,i/T)
  }
  return(res)
}

# Slow then quick evolution of the temperature
ccmm2 = function(i_map,i_temp){
  res = list()
  for (i in 1:T){
    res[[i]] = climat_map_maker2(i_map,i_temp,(i/T)^2)
  }
  return(res)
}

# Quick then slow evolution of the temperature
ccmm3 = function(i_map,i_temp){
  res = list()
  for (i in 1:T){
    res[[i]] = climat_map_maker2(i_map,i_temp,(i/T)^(1/2))
  }
  return(res)
}







