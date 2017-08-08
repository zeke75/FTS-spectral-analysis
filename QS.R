QS <- function (x)
{
  b <- 8
  c <- 0.2
  ifelse(x >= 0, 
         ifelse(x <= c, 1, 3/(b^2*(x-c)^2) * (sin(b*(x-c))/(b*(x-c)) - cos(b*(x-c))) ),
         QS(-x))
}