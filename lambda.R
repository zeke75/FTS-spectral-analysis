M = 1
m = 1/2

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