writeLines("
  model{
    for(i in 1:ns) {
      u[i] ~ dnorm(0,.0001)
      w[i,1] <- 0
      delta[i,t[i,1]] <- 0
      for(k in 1:na[i]) {
       r[i,t[i,k]] ~ dbin(p[i,t[i,k]],n[i,t[i,k]])
       logit(p[i,t[i,k]]) <- u[i]+delta[i,t[i,k]]
       rhat[i,t[i,k]] <- p[i,t[i,k]] * n[i,t[i,k]]
       dev[i,t[i,k]] <- 2 * (r[i,t[i,k]] * (log(r[i,t[i,k]])-log(rhat[i,t[i,k]])) + (n[i,t[i,k]]-r[i,t[i,k]]) * (log(n[i,t[i,k]]-r[i,t[i,k]]) - log(n[i,t[i,k]]-rhat[i,t[i,k]])))
      }
      resdev[i] <- sum(dev[i,t[i,1:na[i]]])
      for (k in 2:na[i]) {
        delta[i,t[i,k]]~dnorm(mean[i,t[i,k]],taud[i,t[i,k]])
        mean[i,t[i,k]] <- d[t[i,k]] - d[t[i,1]] + sw[i,k]
        taud[i,t[i,k]] <- PREC * 2 * (k-1) / k
        w[i,k] <- (delta[i,t[i,k]] - d[t[i,k]] + d[t[i,1]])
        sw[i,k] <- sum(w[i,1:(k-1)]) / (k-1)
      }
    }
    totresdev <- sum(resdev[]) 
    d[ref] <- 0
    for (k in 1:(ref-1)) {
      d[k] ~ dnorm(0,.001)
      logOR.ref[k] <- d[k] - d[ref]
      OR.ref[k] <- exp(logOR.ref[k])
      pred.logOR.ref[k] ~ dnorm(logOR.ref[k],PREC)
      pred.OR.ref[k] <- exp(pred.logOR.ref[k])
    }
    for (k in (ref+1):nt) {

      d[k] ~ dnorm(0,.001)
      logOR.ref[k] <- d[k] - d[ref]
      OR.ref[k] <- exp(logOR.ref[k])
      pred.logOR.ref[k] ~ dnorm(logOR.ref[k],PREC)
      pred.OR.ref[k] <- exp(pred.logOR.ref[k])
    } 
    #tau ~ dnorm(0,1)T(0,)
    tau ~ dunif(0,2)
    PREC <- 1 / pow(tau,2)
    for (c in 1:(nt-1)) {
      for (k in (c+1):nt) {
        logOR[c,k] <- d[c] - d[k]
        OR[c,k] <- exp(logOR[c,k])
        pred.logOR[c,k] ~ dnorm(logOR[c,k],PREC)
        pred.OR[c,k] <- exp(pred.logOR[c,k])
      }
    }
    
    ### if outcome=1 when small.values=undesirable and 0 when small.values are desirable
    
  order[1:nt] <- outcome*rank(-d[1:nt])+(1-outcome)*rank(d[1:nt])
    for (k in 1:nt) {
      most.effective[k] <- equals(order[k],1)
      for (j in 1:nt) {
        effectiveness[k,j] <- equals(order[k],j)
      }
    }
    for (k in 1:nt) {
      for (j in 1:nt) {
        cum.effectiveness[k,j] <- sum(effectiveness[k,1:j])
      }
    }
    for (k in 1:nt) {
      SUCRA[k] <- sum(cum.effectiveness[k,1:(nt-1)]) / (nt-1)
    }    
  }
           ", con="NMA.model.txt")
