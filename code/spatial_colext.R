###################################################################
#  A Col-Ext metapopulation model SPOM
#
#  psi = z_t-1 * (1-E) + (1-z_t-1)*C
#  E = beta_0 + beta_1 * P^gamma[1]
#  C = 1 - exp(-S)
#  S = delta * sum[exp(-1/(2*sigma^2) * d_ij) * Z * P^gamma]
#  
#  P: metric of choice (Area, pop size etc...)
#
#  Parameters to be estimated:
#  - gamma: a power function for P scaling for extinction
#  - beta_0: ext intercept
#  - beta_1: ext slope
#  - delta: per-capita effective dispersal rate


sp.colext <- nimbleCode({

  #~~~~~~~PRIORS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    delta ~ dnorm(0,0.368)
    #gamma ~ dunif(0,10)
    beta0 ~ dnorm(0,0.368)
    beta1 ~ dnorm(0,0.368)
    alpha <- 0.33

  #~~~~~~~Likelihood~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  for(k in 2:nyear){
    for(j in 1:nsite){
      logit(delta0[j,k-1]) <- delta
      delta0_ik[j,k-1] <- delta0[j,k-1] * z[j,k-1] * P[j,k-1]^gamma #P[j,k-1]
    }
    for(i in 1:nsite){
      #connectivity [col]
      for(j in 1:nsite){
        con[i,j,k-1] <- 1 - (delta0_ik[j,k-1] *
                             exp(-alpha / dmat[i,j]) *
                             (1 - equals(i,j)))
      }

      #transition probs
      col[i,k-1] <- 1 - prod(con[i,1:nsite ,k-1]) 
      logit(ext[i,k-1]) <- beta0 + beta1 * P[i,k-1] #area[i]
      ext.re[i,k-1] <- ext[i,k-1] * (1-col[i,k-1])

      #occupancy
      mu.z[i,k-1] <- z[i,k-1] * (1-ext.re[i,k-1]) +
                     (1 - z[i,k-1]) * col[i,k-1]
      z[i,k] ~ dbern(mu.z[i,k-1])
    }
  }
})
