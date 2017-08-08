WID <- function(w)
{
  t = seq(-1, 1, 0.01)
  rs = sum(lamID(t)*cos(w*t))*.01
  
  return(rs)
}

#s = seq(-2,2,0.1)
#plot(s, W(s))
hID <- Vectorize(WID)

plot.function(hID, from = -50, to = 50, n = 1000)