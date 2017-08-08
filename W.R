M=1

W <- function(w)
{
  t = seq(-M, M, 0.01)
  p = w*t
  rs = sum(lam(t)*cos(p))*.01
  
  return(rs)
}

#s = seq(-2,2,0.1)
#plot(s, W(s))
h <- Vectorize(W)

plot.function(h, from = -50, to = 50, n = 1000)