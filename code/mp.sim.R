#simulate metapop dynamics
mp.sim <- function(npatch = 75, nyears = 10,    #varaiables
                   psi1 = 0.4, sigma = 0.15,     #parameters
                   delta = 1, gamma = 0.65,     #parameters
                   int.e = 3, slope.e = -2.5){  #parameters
  
  alpha <- log(2)/sigma  
  xy <- cbind(x=runif(npatch, 2, 8),
              y=runif(npatch, 2, 8))
  dij <- as.matrix(dist(xy))
  area <- rlnorm(npatch, 0, 0.2)

  zmat <- matrix(NA, nrow=npatch, ncol=nyears)
  zmat[,1] <- rbinom(npatch, 1, psi1)

  for(i in 2:nyears){
    ss <- c(delta * (exp(-alpha * dij^2) %*% (zmat[,i-1]*area^gamma)))
    cc <- 1-exp(-ss)
    ee <- plogis(int.e + slope.e * area)
    zmat[,i] <- rbinom(npatch, 1, zmat[,i-1]*(1-ee) + (1-zmat[,i-1])*cc)
  }
  return(zmat)
}

