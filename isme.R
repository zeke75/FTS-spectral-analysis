### ISME ####

B <- 200
isme <- rep(NA, B)

for (i in 1:B) {
  #ev <- rep(NA, freq-1)
  hs <- rep(NA, freq-1)
  source('~/Dropbox/Research/FTS spectral/Simulations/omega.R')
  for (j in 1:(freq-1)){
    #ev[j] = sum((Re(eigen((f-f_omega)[j,,])$value))^2)
    hs[j] = hilbert.schmidt.norm(Re((f-f_omega)[j,,]))^2
  }
  isme[i] <- 2 * sum(hs) * pi/(freq-1)
  print(i)
}

mse <- median(log(isme,2))
boxplot(log(isme,2))