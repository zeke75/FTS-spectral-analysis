library(ftsspec)

T <- 2^9
K <- 1000

A_0 <- matrix(rnorm(50*K,0,(1:50)^(-2)), 50, K)
A_1 <- matrix(rnorm(50*K,0,(1:50)^(-2)), 50, K)

Xi <- matrix(rnorm(K*(T+1)), K, T+1)
lambda <- 1/((1:K-.5)^2*pi^2)

X <- matrix(rep(0, T*50), T, 50)
eps <- matrix(rep(0, (T+1)*K), T+1, K)

#e <- matrix(rep(0, K*K), K, K)

#for (k in 1:K) {
 # e[k, ] <- sqrt(2)*sin((k-.5)*pi*seq(0.5,K-0.5,1)/K)/10
#}


for (t in 1:(T+1)) {
  for (k in 1:K) {
    #eps[t,] <- eps[t,] + Xi[k,t] * sqrt(lambda[k]) * e[k,]
    eps[t, k] <- Xi[k,t] * sqrt(lambda[k])
  }
  if (t>1){
    X[t-1,] <- A_0 %*% eps[t,] + A_1 %*% eps[t-1,]
  }
}



#### 

ans <- Spec(X, W=hID, trace=FALSE, only.diag=FALSE) # Epanechnikov_kernel, hID
freq <- length(ans$omega)
f <- ans$spec[-freq,,]


J <- freq-1
f_omega <- array(NA,c(J,50,50))

for (j in 1:J) {
  #f_omega[j,,] <- 1/(2*pi)*sum(lambda)*(A_0 %*% t(A_0) + A_1 %*% t(A_1) + exp(-((j-1)*pi/J)*1i)*A_0 %*% t(A_1) + exp(((j-1)*pi/J)*1i)* A_1 %*% t(A_0))
  f_omega[j,,] <- 1/(2*pi)*(A_0 %*% diag(lambda) %*% t(A_0) 
                            + A_1 %*% diag(lambda) %*% t(A_1) 
                            + exp(-((j-1)*pi/J)*1i) * A_0 %*% diag(lambda) %*% t(A_1) 
                            + exp(((j-1)*pi/J)*1i) * A_1 %*% diag(lambda) %*% t(A_0))
}



