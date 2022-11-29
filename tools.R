image.nan.better <- function(z,  zlim, col, na.color='gray', outside.below.color='black', outside.above.color='white',...)
{
  zstep <- (zlim[2] - zlim[1]) / length(col); # step in the color palette
  newz.below.outside <- zlim[1] - 2 * zstep # new z for values below zlim
  newz.above.outside <- zlim[2] + zstep # new z for values above zlim
  newz.na <- zlim[2] + 2 * zstep # new z for NA
  
  z[which(z<zlim[1])] <- newz.below.outside # we affect newz.below.outside
  z[which(z>zlim[2])] <- newz.above.outside # we affect newz.above.outside
  z[which(is.na(z>zlim[2]))] <- newz.na # same for newz.na
  
  zlim[1] <- zlim[1] - 2 * zstep # extend lower limit to include below value
  zlim[2] <- zlim[2] + 2 * zstep # extend top limit to include the two new values above and na
  
  col <- c(outside.below.color, col[1], col, outside.above.color, na.color) #correct by including col[1] at bottom of range
  
  image(z=z,  zlim=zlim, col=col, ...) # we finally call image(...)
}


prod_proba <- function(M1,M2) {
  n = length(M1[,1])
  
  result = Matrix(0,n,n)
  
  for (i in 1:n){
    for (j in 1:n){
      
      
      for (k in 1:n){
        result[i,j] = result[i,j] + M1[i,k] * M2[k,j] - result[i,j] * M1[i,k] * M2[k,j] 
      }
      
      
    }
  }
  
  return(result)
}

M1 = rbind( c(0.5,0.25,0.4,0),c(0.5,1,0.4,0),c(0.4,0.2,0.5,0.1),c(0,0,0.2,1))

M2 = rbind( c(1,0.5,0.8,0),c(0.25,0.5,0.2,0),c(0.8,0.4,1,0.2),c(0,0,0.1,0.5))

res = prod_proba(M1,M2)
