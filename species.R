
N = NN

Tmax = 15
Tmin = 5

Tlost = TLOST
Tgainmax = max(unlist(TGAIN))

Tminmin = 0
Tmaxmax = 20

marge = (Tmaxmax-Tminmin)


esps = array(0,c(N,6))
# 1 : Tpref
# 2 : Ttol1
# 3 : Ttol2
# 4 : Abondance
# 5 : MigrationType

esps[,1]=runif(N)*marge + Tminmin
esps[,2]=runif(N)*marge*0.1
esps[,3]=runif(N)*marge*0.1+esps[,2]
esps[,4]=runif(N)*0.1
esps[,5]= sample(1:3, N ,replace=TRUE)

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
esp = esps[1,]
