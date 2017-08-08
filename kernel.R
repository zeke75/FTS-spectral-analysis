Kernel <- function (x)
{
  #BT <- T^{-1/3}
  BT <- 1/20
  h <- 0.5
  M <- 1/BT
  m <- 0.5 * M
  
  fx <- (sin(M*x/2)^2 - sin(m*x/2)^2) / (2*pi*(M-m)*sin(x/2)^2)  
  
  #ifelse( x==0, (M+m)/2/pi, ifelse(abs(x)<=pi, fx, 0))
  ifelse( x==0, (M+m)/2/pi, fx)
   
}