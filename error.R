
B <- 100
err <- rep(NA,B)

for (i in 1:B) {
  source('~/Dropbox/Research/FTS spectral/Simulations/FTS spectral.R')
  err[i] <- sum((f-f_omega)^2)
  print(i)
}

mse <- median(log(Re(err),2))
boxplot(log(Re(err),2))