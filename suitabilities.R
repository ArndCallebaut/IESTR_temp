library(Matrix)

T=PERIOD

# list suitability
# for (t in 1:T){
#   
#   coord = which( (abs(temps[[t]]-esp[i,1])<esp[i,2])&(abs(humis[[t]]-esp[i,3])<esp[i,4])&(abs(heigth_map-esp[i,5])<esp[i,6]) & heigth_map>0,arr.ind=TRUE)
#   
#   if(length(coord)==0){
#     e <- sparseMatrix(dims = c(nrow,ncol), i={}, j={})
#     message("This one is deadperiod..",i," ok... in t=",t)
#     list_suit[[t]]=e
#     next
#   }
#   
#   loc = abs(temps[[t]]-esp[i,1])/esp[i,2]
#   
#   vals = loc[coord]
#   
#   suit = sparseMatrix(coord[,1], coord[,2], x = vals, dims = c(nrow, ncol))
#   
#   
#   list_suit[t] = suit
  
i_map = i_mp
i_temp = i_tg

suit_maps = function(esp,i_m,i_t){
  list_suit = list()
  climats = ccmm1(i_m,i_t)
  for (t in 1:T){
    climat = climats[[t]]
    notnull = (abs(climat-esp[1])<esp[3])
    H = 1 + esp[2]/esp[3]
    locval = H - abs(climat-esp[1])/(esp[3]+esp[2])*H
    suit = pmin(locval,1)*notnull
    
    suit[is.na(suit)]=0
    
    list_suit[[t]] = as(Matrix(suit,sparse=TRUE),"dgCMatrix")
    #list_suit[[t]] =as(as(Matrix(suit,sparse=TRUE),"dgTMatrix"),"dgCMatrix")
  }
  return(list_suit)
}

suits_maps = function(esps,i_m,i_t){
  
  tot_list = list()
  
  for(n in 1:N){
    
    suity = suit_maps(esps[n,],i_m,i_t)
    
    #tot_list[[n]] = as(suity,"dgCMatrix")
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
  len = NN
  
  for (i in 1:len){
    res[[i]] = presence(esps[i,],suits_maps[[i]])
  }
  
  return(res)
}

cost_mat = function(suitabilities_i){
  cost = Matrix(0,nrow,ncol)
  for (j in 1:length(suitabilities_i)){
    xy = which(suitabilities_i[[j]]!=0,arr.ind=TRUE)
    cost[xy[,1],xy[,2]] = 10 + 10*runif(length(xy[,1]))
  }
  return(Matrix(cost,sparse=TRUE))
}






















