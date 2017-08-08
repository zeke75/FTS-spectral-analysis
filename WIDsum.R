WIDsum <- function(x){
  ss = 0
  BT = 0.2
  
  for (i in -500:500){
    ss = ss + 1/BT*hID((x + 2*pi*i)/BT)
  }
  return(ss)
}

hIDsum <- Vectorize(WIDsum)

plot.function(hIDsum, from = -50, to = 50, n = 1000)

