par(mfrow=c(2,2))

c = 0.05
b = 0.25

### lamID: infinite differentiable flat-top function
lamID <- function(x)
{
  ifelse( abs(x) <= c, 
          1, 
          ifelse( abs(x) < 1, 
                  exp(-b*exp(-b/(abs(x)-c)^2)/(abs(x)-1)^2), 
                  0
          ) 
  )
}
s = seq(-2,2,0.01)
plot(s, lamID(s), type='l')


### WID: inverse fourier transform of lamID
WID <- function(w)
{
  t = seq(-1, 1, 0.01)[-1]-0.005
  rs = 1/2/pi*sum(lamID(t)*cos(w*t))*.01
  
  return(rs)
}
hID <- Vectorize(WID)
plot.function(hID, from = -50, to = 50, n = 1000)


### WIDsum: W^(T)
WIDsum <- function(x){
  ss = 0
  BT = 0.1
  
  for (i in -200:200){
    ss = ss + 1/BT*hID((x + 2*pi*i)/BT)
  }
  return(ss)
}
hIDsum <- Vectorize(WIDsum)
plot.function(hIDsum, from = -10, to = 10, n = 1000)




WID_fs <- function(x){
  fs = 0
  B_T = 0.1
  for (i in -20:20){
    fs = fs + 1/2/pi*lamID(i*B_T)*cos(i*x)
  }
  return(fs)
}
  
hIDfs <- Vectorize(WID_fs)
plot.function(hIDfs, from = -10, to = 10, n = 1000)

