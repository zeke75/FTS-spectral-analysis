par(mfrow=c(2,2))

M = 1
m = 1/2

### lam: trapezoidal flap-top function
lam <- function(x)
{
  ifelse( abs(x) <= m, 
          1, 
          ifelse( abs(x) <= M, 
                  1-(abs(x)-m)/(M-m), 
                  0
          ) 
  )
}
s = seq(-2,2,0.1)
plot(s, lam(s), type='l')


### W: fourier transform of lam
W <- function(w)
{
  t = seq(-M, M, 0.01)
  p = w*t
  rs = sum(lam(t)*cos(p))*.01
  
  return(rs)
}
h <- Vectorize(W)
plot.function(h, from = -50, to = 50, n = 1000)


### Wsum: W^(T)
Wsum <- function(x){
  ss = 0
  BT = 0.1
  
  for (i in -500:500){
    ss = ss + 1/BT*h((x + 2*pi*i)/BT)
  }
  return(ss)
}
hsum <- Vectorize(Wsum)
plot.function(hsum, from = -10, to = 10, n = 1000)
