f = function(x){
  y = x / 440
  y[y>1] = 1
  return(y)
}

NN = 100:1100

ff = f(NN) + runif(1000)*0.3
plot(ff)

seq(1,1000,10)


SS= ff[seq(1,1000,30)]
SS = SS/max(SS)

SS[SS>0.85] = 0.85

SS= SS/max(SS)

plot(SS)