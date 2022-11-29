
esp = esps[1,]

T=50

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
  
i_map = 1
i_temp = 2

suit_maps = function(esp){
  
  list_suit = list()
  climats = ccmm1(i_map,i_temp)
  
  for (t in 1:T){
    climat = climats[[t]]
    notnull = (abs(climat-esp[1])<esp[3])
    H = 1 + esp[2]/esp[3]
    locval = H - abs(climat-esp[1])/(esp[3]+esp[2])*H
    suit = pmin(locval)*notnull
    list_suit[[t]] = suit
  }
  
  return(list_suit)
}
