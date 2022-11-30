
### ### ###
### Functions for height map construction 
###

local_height = function(x,y){
  # X & Y in [0,1]
  # function of the height of the island. 
  return (exp(-3*abs(x-y))* (sqrt((x-0.5)^2+(y-0.5)^2)<0.45))
}

height_map = function(nrow,ncol,FUN=local_height){
  # building the raster of the island.
  # FUN = function defining the height of the island. Use local_height by default.
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

### ### ###
### Function for temperatures on the map, at each time stap
###

temperatures = function(x,y){
  # Temperature on the map for altitude=0, for x & y in [0,1]
  return(-(Tmax_alt_0-Tmin_alt_0)*y+Tmax_alt_0)
}
climat_map_maker = function(height_map,Tgain,t){
  # Temperature on the map at each altitude, at a time t.
  # t in [0,1], 0 corresponding to the beginning of the period, 1 to the end.
  return(outer(nr,nc,FUN="temperatures")*(height_map!=0)+height_map*Tadd_alt_1+Tgain*t) 
}
climat_maps_maker = function(height_map,Tgain,T){
  # Temperature on the map at each altitude, at each time t.
  res = list()
  for (i in 0:T){
    t = i/T
    matry = climat_map_maker(height_map,Tgain,t)
    res[[i+1]] = matry
  }
  return(res)
}

### ### ###
### Function for suitabilities on the map, at each time stap
###

# OLD VERSION
# suitability_maps = function(climate_maps,Trange_spe){
#   # Suitability on the map, at each time step
#   list_suit = list()
#   for (t in 0:N_cycles){
#     climat = c_maps[[t+1]]
#     notnull = (Trange_spe[1]<climat)&(Trange_spe[4]>climat)
#     loc1 = (climat-Trange_spe[1]) / (Trange_spe[2]-Trange_spe[1])
#     loc2 = (climat-Trange_spe[4]) / (Trange_spe[3]-Trange_spe[4])
#     suit = pmin(loc1,loc2,1)
#     suit[is.na(suit) | suit<0]=0
#     list_suit[[t+1]] = as(Matrix(suit,sparse=TRUE),"dgCMatrix")
#   }
#   return(list_suit)
# }

make_suitability_maps = function(climate_maps,Trange_spe){
  # Suitability on the map, at each time step
  list_suit = list()
  for (t in 0:N_cycles){
    climat = climate_maps[[t+1]]
    notnull = (Trange_spe[1]<climat)&(Trange_spe[4]>climat)
    loc1 = (climat-Trange_spe[1]) / (Trange_spe[2]-Trange_spe[1])
    loc2 = (climat-Trange_spe[4]) / (Trange_spe[3]-Trange_spe[4])
    suit = pmin(loc1,loc2,1)
    suit[is.na(suit) | suit<0]=0
    list_suit[[t+1]] = as(Matrix(suit,sparse=TRUE),"dgCMatrix")
  }
  return(list_suit)
}

### ### ###
### Function for costs on the map
###

cost_mapping1 = function(x,y){
  return (300+y*0)
}
cost_mapping2 = function(x,y){
  return (100+y*400)
}
cost_mapping3 = function(x,y){
  return (100+local_height(x,y)*400)
}
cost_map = function(nrow,ncol,map,FUN=cost_mapping){
  x = seq(0,1,length=nrow)
  y = seq(0,1,length=ncol)
  v = outer(x,y,FUN)
  cost = (!is.na(map)) * v
  cost[is.na(map)] = NA
  if(any (class(cost) %in%"dtCMatrix")){cost = as(as(cost,"dgCMatrix"),"dgCMatrix")}
  else{cost = as(cost,"dgCMatrix")}
  cost = as(as(cost,"dgTMatrix"),"dgCMatrix")
  return(cost)
}
cost_map_modified = function(nrow,ncol,map,FUN=cost_mapping){
  x = seq(0,1,length=nrow)
  y = seq(0,1,length=ncol)
  v = outer(x,y,FUN) + map*400
  cost = (!is.na(map)) * v
  if(any (class(cost) %in% "dtCMatrix")) {cost = as(as(cost,"dgCMatrix"),"dgCMatrix")}
  else{cost = as(cost,"dgCMatrix")}
  cost = as(as(cost,"dgTMatrix"),"dgCMatrix")
  return(cost)
}

### ### ###
### Function for presence on the map
###

make_presence_map = function(nrow,ncol,suit_maps,height_map,xylim,nb_cell){
  suitnow = (suit_maps[[1]][nrow:1,]==1)*1
  maty <- (matrix(runif(ncol*nrow),ncol = ncol))                # Specify number of columns
  Xmin = round(xylim[1] * ncol)
  Xmax = round(xylim[2] * ncol)
  Ymin = round(xylim[3] * nrow)
  Ymax = round(xylim[4] * nrow)
  mat <- (matrix(rep(0,ncol*nrow),ncol = ncol)) 
  mat[Xmin:Xmax,Ymin:Ymax]=1
  mat[1:Xmin,]=0
  mat[Xmax:ncol,]=0
  mat[,1:Ymin]=0
  mat[,Ymax:nrow]=0
  
  viable_coord = which((mat*(0<height_map))==1,arr.ind=TRUE)
  
  nb_possible = length(viable_coord[,1])
  if(nb_possible<nb_cell){
    message("ERROR : too much cells asked. ")
    pres = mat * 0
    pres[is.na(height_map)] = NA
    pres[viable_coord] = 1
    return(pres)
  }
  indices_presence = sample(nb_possible, nb_cell)
  coord_presence = viable_coord[indices_presence,]
  pres = mat * 0
  pres[is.na(height_map)] = NA
  pres[coord_presence] = 1
  return(pres)
}

# presence_map = function(nrow,ncol,suit_maps,height_map,xylim,nb_cell){
#   suitnow = (suit_maps[[1]][nrow:1,]==1)*1
#   maty <- (matrix(runif(ncol*nrow),ncol = ncol))                # Specify number of columns
#   Xmin = round(xylim[1] * ncol)
#   Xmax = round(xylim[2] * ncol)
#   Ymin = round(xylim[3] * nrow)
#   Ymax = round(xylim[4] * nrow)
#   mat <- (matrix(rep(0,ncol*nrow),ncol = ncol)) 
#   mat[Xmin:Xmax,Ymin:Ymax]=1
#   mat[1:Xmin,]=0
#   mat[Xmax:ncol,]=0
#   mat[,1:Ymin]=0
#   mat[,Ymax:nrow]=0
#   
#   viable_coord = which((mat*(0<height_map))==1,arr.ind=TRUE)
#   
#   nb_possible = length(viable_coord[,1])
#   if(nb_possible<nb_cell){
#     message("ERROR : too much cells asked. ")
#     pres = mat * 0
#     pres[is.na(height_map)] = NA
#     pres[viable_coord] = 1
#     return(pres)
#   }
#   indices_presence = sample(nb_possible, nb_cell)
#   coord_presence = viable_coord[indices_presence,]
#   pres = mat * 0
#   pres[is.na(height_map)] = NA
#   pres[coord_presence] = 1
#   return(pres)
# }

### ### ###
### Function for presence on the map
###



















