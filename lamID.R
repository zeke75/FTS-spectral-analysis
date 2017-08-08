c=0.05
b=0.25

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