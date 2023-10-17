##########################################################
### SWAG Covariance Model Gibbs Sampler- SWAG_GS
## SWAG_GS_cpp: CPP ENABLED version
## classic GS parameters
# S: number of iterations for GS
# burnin: number of iterations to burnin
# thin: thinning mechanism
# save_all: logical to determine what parameters to return. If 0, save group covariances, 
#       degrees of freedoms, and weight. If 1, save and return all parameters.
## optional parameters
# mh.delta.star: delta for random walk proposal for weight(lambda) proposal
# K: Sample space parameter for degrees of freedom (upper bound is p+2)
# df_MH: logical to sample degrees of freedom from a Metropolis Hastings step.
#       If 0, sample degrees of freedom from full conditional. If 1, perform MH step. 
#       Suggest to select 1 if p is very large. If 1, K is not used.
##########################################################
SWAG_GS_cpp = function(S,burnin = round(S*.1),thin = 10,
                  save_all = 0,                        
                  mh.delta.star = .1,                 
                  K = 50,                             
                  df_MH = 0){                         
  
  ###########################
  ## set hyper params
  S0 = eye(p); S0.inv = solve(S0) 
  V0 = S0; V0.inv = S0.inv
  
  eta0 = p + 2
  S1 = eye(p1); S1.inv = solve(S1)
  S2 = eye(p2); S2.inv= solve(S2)
  eta1 = p1 + 2; eta2 = p2 + 2
  
  nu.domain = c((p+2):(p+K)); ND = length(nu.domain)
  
  # priors
  N1 = eye(p1)
  N2 = eye(p2)
  xi1 = p1 + 2
  xi2 =  p2 + 2
  ###########################
  
  ###########################
  ## intialize values
  # level 1
  U = lapply(1:g,function(j)matrix(rnorm(ns[j]*p),ncol=p))
  Sig = Sig.inv = Psi = Psi.inv = lapply(1:g,function(j)eye(p))
  # level 2 and 3
  nu0 = p+2
  nu = p+2
  V0 = V0.inv = eye(p)
  V1 = lapply(1:g,function(j)rwish(S1,eta1))
  V1.inv = lapply(V1,solve)
  V2 = lapply(1:g,function(j)rwish(S2,eta2))
  V2.inv = lapply(V2,solve)
  # level 3
  K1 = K1.inv = eye(p1)
  K2 = K2.inv = eye(p2)
  S0 = S0.inv = eye(p)
  eta0 = p+2
  
  pis = .5
  ###########################
  
  ###########################
  ## create storage  
  acc = 0
  if ( df_MH == 1){
    nu.acc = nu0.acc = eta0.acc = 0
  }
  
  if (save_all == 0){
    len.out = length(seq(from = burnin+1, to = S, by=thin))
    # layer 1
    cov.out = cov.inv = array(NA,dim = c(len.out,p*p*g))
    nu.out = array(NA,dim = c(len.out,2))
    eta0.out = array(NA,dim = c(len.out,1))
    pis.out = array(NA,dim=c(len.out,1))
    index = 1
  } else {
    # layer 1
    Psi.out = array(NA,dim = c(S,p*p*g))
    Sig.out = array(NA,dim = c(S,p*p*g))
    cov.out = cov.inv = array(NA,dim = c(S,p*p*g))
    U.out   = array(NA,dim = c(S,sum(ns*p)))
    pis.out = array(NA,dim=c(S,1))
    # layer 2
    nu.out = array(NA,dim = c(S,2))
    V0.out = array(NA,dim = c(S,p*p))
    V1.out = array(NA,dim = c(S,p1*p1*g))
    V2.out = array(NA,dim = c(S,p2*p2*g))
    # layer 3
    K1.out = array(NA,dim = c(S,p1*p1))
    K2.out = array(NA,dim = c(S,p2*p2))
    eta0.out = array(NA,dim = c(S,1))
  }
  ###########################
  
  ###########################
  ## helpers
  
  d.sig = matrix(NA,nrow=g,ncol=ND)
  
  S.list = list()
  for ( j in 1:g ){
    S.list[[j]] = t(Y.list[[j]]) %*% Y.list[[j]]
  }
  
  cov.inv.temp = cov.temp = list()
  
  ###########################
  
  ###########################
  ## GS
  for( s in 1:S ){
    
    if (s < min(1000,burnin - 500)){
      mh.delta = .5
      df.delta = rep(round(p*3),3) # only using if df_MH == 1
    } else {
      mh.delta = mh.delta.star
      df.delta.star = c(3,1,1) #rep(round(p/4),3) # only using if df_MH == 1
    }
    
    
    ## sample nu and Sig
    if (df_MH == 0){
      # sample nu - uniform on nu.domain
      for ( j in 1:g ){
        Vj.inv = kronecker(V2.inv[[j]],V1.inv[[j]])
        Vj = kronecker(V2[[j]],V1[[j]])
        
        helper = (Y.list[[j]] - pis^(1/2)*U[[j]]) %*% Vj.inv %*% 
          t(Y.list[[j]] - pis^(1/2)*U[[j]]) / (1-pis)
        d.sig[j,] = sapply(nu.domain, function(k)
          dmatT.faster(helper,k-p+1,k,p,ns[j]))
      }
      # same nu for each sigj
      temp = apply(d.sig,2,sum)
      temp = exp(temp-max(temp))
      probs = temp/sum(temp)
      nu = sample(nu.domain,size = 1,prob = probs)
    } else if (df_MH == 1){
      nu.star = reflect(sample((nu - df.delta[2]):(nu + df.delta[2]),1),p+2)
      R=0
      for ( j in 1:g ){
        Vj.inv = kronecker(V2.inv[[j]],V1.inv[[j]])
        helper = (Y.list[[j]] - pis^(1/2)*U[[j]]) %*% Vj.inv %*% 
          t(Y.list[[j]] - pis^(1/2)*U[[j]]) / (1-pis)
        
        R = R + dmatT.faster(helper,nu.star-p+1,nu.star,p,ns[j]) - 
          dmatT.faster(helper,nu-p+1,nu,p,ns[j])
      }
      ## if u<r, set nu to be nu.star
      if (log(runif(1))<R){
        nu = nu.star
        nu.acc = nu.acc + 1
      }
    }
    
    
    ## sample xSigj
    for( j in 1:g ){
      V = kronecker(V2[[j]],V1[[j]])
      Y.tilde = Y.list[[j]] - pis^(1/2)*U[[j]]
      M = csolve(V * (nu - p - 1) + t(Y.tilde) %*% Y.tilde / (1-pis))
      Sig.inv[[j]] = rwish.wrapper(M,nu+ns[j]-1)
      Sig[[j]] = csolve(Sig.inv[[j]])
    }
    
    ## sample pi and U
    # sample pi
    pi.hat = pis[1]
    pis.star = MH_sym_proposal_01(pi.hat, mh.delta)
    ## compute acceptance ratio, joint distn' of all variables and pi
    d.star = sapply(1:g,function(j)dmatnorm(Y.list[[j]],0,
                                            eye(ns[j]),
                                            pis.star * Psi[[j]] + (1-pis.star) * Sig[[j]],
                                            if_log = TRUE))

    d.s = sapply(1:g,function(j)dmatnorm(Y.list[[j]],0,
                                         eye(ns[j]),
                                         pi.hat * Psi[[j]] + (1-pi.hat) * Sig[[j]],
                                         if_log = TRUE))

    R = sum(d.star) - sum(d.s) +
      dbeta(pis.star,1/2,1/2,log=T) - dbeta(pi.hat,1/2,1/2,log=T) 
    ## if u<r, set pis to be pis.star
    if (log(runif(1))<R){
      pis = pis.star
      acc = acc + 1
    }
    
    # sample U
    for ( j in 1:g ){
      W = csolve(pis/(1-pis)*Sig.inv[[j]] + Psi.inv[[j]])
      M = pis^(1/2)/(1-pis) * Y.list[[j]] %*% Sig.inv[[j]] %*% W
      U[[j]] = rmatnorm(M,eye(ns[j]),W)
    }
    
    ## sample nu0 and psi
    if (df_MH == 0){
      # sample nu0 - uniform on nu.domain 
      for ( j in 1:g){
        
        helper = U[[j]] %*% V0.inv %*% t(U[[j]])
        d.sig[j,] = sapply(nu.domain, function(k)
          dmatT.faster(helper,k-p+1,k,p,ns[j]))
        ## same thing; maybe above is more stable/def faster
        # d.sig[j,] = sapply(nu.domain,function(k)
        #   dmatT.propto(U[[j]],k-p+1,0,V0 * (k-p-1),k,
        #                Omega.inv = V0.inv / (k-p-1) ) )
        
      }
      temp = apply(d.sig,2,sum)
      temp = exp(temp-max(temp))
      probs = temp/sum(temp)
      nu0 = sample(nu.domain,size = 1,prob = probs)
    } else if (df_MH == 1){
      nu0.star = reflect(sample((nu0 - df.delta[1]):(nu0 + df.delta[1]),1),p+2)
      R = 0
      for ( j in 1:g){
        helper = U[[j]] %*% V0.inv %*% t(U[[j]])
        
        R = R + dmatT.faster(helper,nu0.star-p+1,nu0.star,p,ns[j]) - 
          dmatT.faster(helper,nu0-p+1,nu0,p,ns[j])
      }
      # if u<r, set nu0 to be nu0.star
      if (log(runif(1))<R){
        nu0 = nu0.star
        nu0.acc = nu0.acc + 1
      }
    }

    # sample Psijs
    for ( j in 1:g){
      M = csolve(V0*(nu0 - p - 1) + t(U[[j]]) %*% U[[j]])
      Psi.inv[[j]] = rwish.wrapper(M,nu0 + ns[j] - 1)
      Psi[[j]] = csolve(Psi.inv[[j]])
    }
    
    
    ### wishart means
    
    ## sample V0
    M0 = csolve(Reduce('+',Psi.inv)*(nu0-p-1) + S0.inv*eta0)
    V0 = rwish.wrapper(M0,eta0 + nu0*g)
    V0.inv = csolve(V0)
    
    ## sample g V1s
    Sig.inv.chol = lapply(Sig.inv,function(l)t(chol(l))) ### chol(S) = LL^T
    for ( j in 1:g ){
      L = array(Sig.inv.chol[[j]],dim=c(p1,p2,p))
      helper = lapply(1:p,function(k) L[,,k] %*% t(V2[[j]]) %*% t(L[,,k]))
      M = csolve(Reduce('+',helper)*(nu-p-1) + S1.inv * eta1)
      V1[[j]] = rwish.wrapper(M,eta1 + nu * p2)
      V1.inv[[j]] = csolve(V1[[j]])
    }
    
    ## sample g V2s 
    for ( j in 1:g ){
      L = array(Sig.inv.chol[[j]],dim=c(p1,p2,p))
      helper = lapply(1:p,function(k) t(L[,,k]) %*% t(V1[[j]]) %*% L[,,k])
      M = csolve(Reduce('+',helper)*(nu-p-1) + S2.inv * eta2)
      V2[[j]] = rwish.wrapper(M,eta2 + nu * p1)
      V2.inv[[j]] = csolve(V2[[j]])
    }
    
    ### layer to shrink homo to kron
    
    # sample K1
    V0.chol = t(chol(V0)) ### chol(S) = LL^T
    L = array(V0.chol,dim=c(p1,p2,p)) 
    helper = lapply(1:p,function(k) (L[,,k]) %*% t(K2.inv) %*% t(L[,,k]))
    M = csolve(Reduce('+',helper)*eta0 + N1 * (xi1 - p1 - 1))
    K1.inv = rwish.wrapper(M, p2*eta0 + xi1)
    K1 = csolve(K1.inv)
    
    # sample K2
    helper = lapply(1:p,function(k) t(L[,,k]) %*% t(K1.inv) %*% (L[,,k]))
    M = csolve(Reduce('+',helper)*eta0 + N2 * (xi2 - p2 - 1))
    K2.inv = rwish.wrapper(M, p1*eta0 + xi2)
    K2 = csolve(K2.inv)
    
    # format to be S0 for code
    S0 = kronecker(K2,K1)
    S0.inv = kronecker(K2.inv,K1.inv)
    
    ## sample eta0
    if (df_MH == 0){
      d.v0 = sapply(nu.domain,function(k)
        dwish(V0,S0/k,k,if_log = T) )
      d.v0 = exp(d.v0-max(d.v0))
      probs = d.v0/sum(d.v0)
      eta0 = sample(nu.domain,size = 1,prob = probs)
    } else if (df_MH == 1){
      eta0.star = reflect(sample((eta0 - df.delta[3]):(eta0 + df.delta[3]),1),p+2)
      R = dwish(V0,S0/eta0.star,eta0.star,if_log=T) - 
        dwish(V0,S0/eta0,eta0,if_log=T)
      ## if u<r, set eta0 to be eta0.star
      if (log(runif(1))<R){
        eta0 = eta0.star
        eta0.acc = eta0.acc + 1
      }
    }
    
    
    ## compute for output
    for ( j in 1:g ){
      cov.temp[[j]] = (1-pis)*Sig[[j]]+pis*Psi[[j]]
      cov.inv.temp[[j]] = csolve(cov.temp[[j]])
    }
    
    ## store output
    if (save_all == 0){
      if((s>burnin)&((s %% thin)==0)){
        # layer 1
        cov.out[index,] = unlist(cov.temp)
        cov.inv[index,] = unlist(cov.inv.temp)
        nu.out[index,] = c(nu0,nu)
        eta0.out[index,] = eta0
        pis.out[index,] = pis
        index = index + 1
      }
    } else {
      # layer 1
      Sig.out[s,]  = unlist(Sig)
      Psi.out[s,] = unlist(Psi)
      cov.out[s,] = unlist(cov.temp)
      cov.inv[s,] = unlist(cov.inv.temp)
      U.out[s,] = unlist(U)
      pis.out[s,] = pis
      # layer 2
      V0.out[s,] = c(V0)
      V1.out[s,] = unlist(V1)
      V2.out[s,] = unlist(V2)
      nu.out[s,] = c(nu0,nu)
      # layer 3- shrink homo to kron
      K1.out[s,] = c(K1)
      K2.out[s,] = c(K2)
      eta0.out[s,] = eta0
    }
    
    
  }
  
  if (save_all == 0){
    out = list("cov.out" = cov.out,
               "cov.inv" = cov.inv,"eta0" = eta0.out,
               "nu" = nu.out,"pis" = pis.out,
               "acc" = acc)
  } else {
    out = list("Sig" = Sig.out, "Psi" = Psi.out,
                "U" = U.out,
                "V0" = V0.out, "V1" = V1.out,
                "V2" = V2.out,"nu" = nu.out,
                "K1" = K1.out, "K2" = K2.out,
                "eta0" = eta0.out,"cov.inv" = cov.inv,
                "pis" = pis.out,"acc" = acc,
                "nu.domain" = nu.domain)
  }
  
  if (df_MH == 0){
    return(out)
  } else if (df_MH == 1){
    return(c(out,list("acc.dfs" = c(nu.acc,nu0.acc,eta0.acc))))
  }
  
}

SWAG_GS = function(S,burnin = round(S*.1),thin = 10,
                   save_all = 0,                        
                   mh.delta.star = .1,                 
                   K = 50,                             
                   df_MH = 0,
                   df.delta.star = c(rep(round(p/4),2),1)){ # delta for reflecting random walk proposal for three df                        
  
  ###########################
  ## set hyper params
  
  ## hyper prior for df
  EPS = .1
  MAX.NU = 2/EPS + (p+3)
  
  DF.PROB = 0.2
  DF.SIZE = round((MAX.NU-(p+2))/2) / ((1-DF.PROB)/DF.PROB) # mean is halfway between lower bound and max
  
  ## rest of hyper priors
  S0 = eye(p); S0.inv = solve(S0) 
  V0 = S0; V0.inv = S0.inv
  
  eta0 = p + 2
  S1 = eye(p1); S1.inv = solve(S1)
  S2 = eye(p2); S2.inv= solve(S2)
  eta1 = p1 + 2; eta2 = p2 + 2
  
  nu.domain = c((p+2):(p+K)); ND = length(nu.domain)
  
  # priors
  N1 = eye(p1)
  N2 = eye(p2)
  xi1 = p1 + 2
  xi2 =  p2 + 2
  ###########################
  
  ###########################
  ## intialize values
  # level 1
  U = lapply(1:g,function(j)matrix(rnorm(ns[j]*p),ncol=p))
  Sig = Sig.inv = Psi = Psi.inv = lapply(1:g,function(j)eye(p))
  # level 2 and 3
  nu0 = p+2
  nu = p+2
  V0 = V0.inv = eye(p)
  V1 = lapply(1:g,function(j)rwish(S1,eta1))
  V1.inv = lapply(V1,solve)
  V2 = lapply(1:g,function(j)rwish(S2,eta2))
  V2.inv = lapply(V2,solve)
  # level 3
  K1 = K1.inv = eye(p1)
  K2 = K2.inv = eye(p2)
  S0 = S0.inv = eye(p)
  eta0 = p+2
  
  pis = .5
  ###########################
  
  ###########################
  ## create storage  
  acc = 0
  if ( df_MH == 1){
    nu.acc = nu0.acc = eta0.acc = 0
  }
  
  if (save_all == 0){
    len.out = length(seq(from = burnin+1, to = S, by=thin))
    # layer 1
    cov.out = cov.inv = array(NA,dim = c(len.out,p*p*g))
    nu.out = array(NA,dim = c(len.out,2))
    eta0.out = array(NA,dim = c(len.out,1))
    pis.out = array(NA,dim=c(len.out,1))
    index = 1
  } else {
    # layer 1
    Psi.out = array(NA,dim = c(S,p*p*g))
    Sig.out = array(NA,dim = c(S,p*p*g))
    cov.out = cov.inv = array(NA,dim = c(S,p*p*g))
    U.out   = array(NA,dim = c(S,sum(ns*p)))
    pis.out = array(NA,dim=c(S,1))
    # layer 2
    nu.out = array(NA,dim = c(S,2))
    V0.out = array(NA,dim = c(S,p*p))
    V1.out = array(NA,dim = c(S,p1*p1*g))
    V2.out = array(NA,dim = c(S,p2*p2*g))
    # layer 3
    K1.out = array(NA,dim = c(S,p1*p1))
    K2.out = array(NA,dim = c(S,p2*p2))
    eta0.out = array(NA,dim = c(S,1))
  }
  ###########################
  
  ###########################
  ## helpers
  
  d.sig = matrix(NA,nrow=g,ncol=ND)
  
  S.list = list()
  for ( j in 1:g ){
    S.list[[j]] = t(Y.list[[j]]) %*% Y.list[[j]]
  }
  
  cov.inv.temp = cov.temp = list()
  
  ###########################
  
  ###########################
  ## GS
  for( s in 1:S ){
    
    if (s < min(1000,burnin - 500)){
      mh.delta = .5
      df.delta = rep(round(p*3),3) # only using if df_MH == 1
    } else {
      mh.delta = mh.delta.star
      df.delta = df.delta.star # only using if df_MH == 1
    }
    
    
    ## sample nu and Sig
    if (df_MH == 0){
      # sample nu - uniform on nu.domain
      for ( j in 1:g ){
        Vj.inv = kronecker(V2.inv[[j]],V1.inv[[j]])
        Vj = kronecker(V2[[j]],V1[[j]])
        
        helper = (Y.list[[j]] - pis^(1/2)*U[[j]]) %*% Vj.inv %*% 
          t(Y.list[[j]] - pis^(1/2)*U[[j]]) / (1-pis)
        d.sig[j,] = sapply(nu.domain, function(k)
          dmatT.faster(helper,k-p+1,k,p,ns[j]))
      }
      # same nu for each sigj
      temp = apply(d.sig,2,sum)
      temp = exp(temp-max(temp))
      probs = temp/sum(temp)
      nu = sample(nu.domain,size = 1,prob = probs)
    } else if (df_MH == 1){
      nu.star = reflect(sample((nu - df.delta[2]):(nu + df.delta[2]),1),p+2)
      R=0
      for ( j in 1:g ){
        Vj.inv = kronecker(V2.inv[[j]],V1.inv[[j]])
        helper = (Y.list[[j]] - pis^(1/2)*U[[j]]) %*% Vj.inv %*% 
          t(Y.list[[j]] - pis^(1/2)*U[[j]]) / (1-pis)
        
        R = R + dmatT.faster(helper,nu.star-p+1,nu.star,p,ns[j]) - 
          dmatT.faster(helper,nu-p+1,nu,p,ns[j]) +
          dnbinom(nu.star - (p+2),DF.SIZE,DF.PROB,log=TRUE) -
          dnbinom(nu - (p+2),DF.SIZE,DF.PROB,log=TRUE)
      }
      ## if u<r, set nu to be nu.star
      if (log(runif(1))<R){
        nu = nu.star
        nu.acc = nu.acc + 1
      }
    }
    
    
    ## sample xSigj
    for( j in 1:g ){
      V = kronecker(V2[[j]],V1[[j]])
      Y.tilde = Y.list[[j]] - pis^(1/2)*U[[j]]
      M = solve(V * (nu - p - 1) + t(Y.tilde) %*% Y.tilde / (1-pis))
      Sig.inv[[j]] = rwish(M,nu+ns[j]-1)
      Sig[[j]] = solve(Sig.inv[[j]])
    }
    
    ## sample pi and U
    # sample pi
    pi.hat = pis[1]
    pis.star = MH_sym_proposal_01(pi.hat, mh.delta)
    ## compute acceptance ratio, joint distn' of all variables and pi
    d.star = sapply(1:g,function(j)dmatnorm(Y.list[[j]],0,
                                            eye(ns[j]),
                                            pis.star * Psi[[j]] + (1-pis.star) * Sig[[j]],
                                            if_log = TRUE))
    
    d.s = sapply(1:g,function(j)dmatnorm(Y.list[[j]],0,
                                         eye(ns[j]),
                                         pi.hat * Psi[[j]] + (1-pi.hat) * Sig[[j]],
                                         if_log = TRUE))
    
    R = sum(d.star) - sum(d.s) +
      dbeta(pis.star,1/2,1/2,log=T) - dbeta(pi.hat,1/2,1/2,log=T) 
    ## if u<r, set pis to be pis.star
    if (log(runif(1))<R){
      pis = pis.star
      acc = acc + 1
    }
    
    # sample U
    for ( j in 1:g ){
      W = solve(pis/(1-pis)*Sig.inv[[j]] + Psi.inv[[j]])
      M = pis^(1/2)/(1-pis) * Y.list[[j]] %*% Sig.inv[[j]] %*% W
      U[[j]] = rmatnorm(M,eye(ns[j]),W)
    }
    
    ## sample nu0 and psi
    if (df_MH == 0){
      # sample nu0 - uniform on nu.domain 
      for ( j in 1:g){
        
        helper = U[[j]] %*% V0.inv %*% t(U[[j]])
        d.sig[j,] = sapply(nu.domain, function(k)
          dmatT.faster(helper,k-p+1,k,p,ns[j]))
        ## same thing; maybe above is more stable/def faster
        # d.sig[j,] = sapply(nu.domain,function(k)
        #   dmatT.propto(U[[j]],k-p+1,0,V0 * (k-p-1),k,
        #                Omega.inv = V0.inv / (k-p-1) ) )
        
      }
      temp = apply(d.sig,2,sum)
      temp = exp(temp-max(temp))
      probs = temp/sum(temp)
      nu0 = sample(nu.domain,size = 1,prob = probs)
    } else if (df_MH == 1){
      nu0.star = reflect(sample((nu0 - df.delta[1]):(nu0 + df.delta[1]),1),p+2)
      R = 0
      for ( j in 1:g){
        helper = U[[j]] %*% V0.inv %*% t(U[[j]])
        
        R = R + dmatT.faster(helper,nu0.star-p+1,nu0.star,p,ns[j]) - 
          dmatT.faster(helper,nu0-p+1,nu0,p,ns[j]) +
          dnbinom(nu0.star - (p+2),DF.SIZE,DF.PROB,log=TRUE) -
          dnbinom(nu0 - (p+2),DF.SIZE,DF.PROB,log=TRUE)
      }
      # if u<r, set nu0 to be nu0.star
      if (log(runif(1))<R){
        nu0 = nu0.star
        nu0.acc = nu0.acc + 1
      }
    }
    
    # sample Psijs
    for ( j in 1:g){
      M = solve(V0*(nu0 - p - 1) + t(U[[j]]) %*% U[[j]])
      Psi.inv[[j]] = rwish(M,nu0 + ns[j] - 1)
      Psi[[j]] = solve(Psi.inv[[j]])
    }
    
    
    ### wishart means
    
    ## sample V0
    M0 = solve(Reduce('+',Psi.inv)*(nu0-p-1) + S0.inv*eta0)
    V0 = rwish(M0,eta0 + nu0*g)
    V0.inv = solve(V0)
    
    ## sample g V1s
    Sig.inv.chol = lapply(Sig.inv,function(l)t(chol(l))) ### chol(S) = LL^T
    for ( j in 1:g ){
      L = array(Sig.inv.chol[[j]],dim=c(p1,p2,p))
      helper = lapply(1:p,function(k) L[,,k] %*% t(V2[[j]]) %*% t(L[,,k]))
      M = solve(Reduce('+',helper)*(nu-p-1) + S1.inv * eta1)
      V1[[j]] = rwish(M,eta1 + nu * p2)
      V1.inv[[j]] = solve(V1[[j]])
    }
    
    ## sample g V2s 
    for ( j in 1:g ){
      L = array(Sig.inv.chol[[j]],dim=c(p1,p2,p))
      helper = lapply(1:p,function(k) t(L[,,k]) %*% t(V1[[j]]) %*% L[,,k])
      M = solve(Reduce('+',helper)*(nu-p-1) + S2.inv * eta2)
      V2[[j]] = rwish(M,eta2 + nu * p1)
      V2.inv[[j]] = solve(V2[[j]])
    }
    
    ### layer to shrink homo to kron
    
    # sample K1
    V0.chol = t(chol(V0)) ### chol(S) = LL^T
    L = array(V0.chol,dim=c(p1,p2,p)) 
    helper = lapply(1:p,function(k) (L[,,k]) %*% t(K2.inv) %*% t(L[,,k]))
    M = solve(Reduce('+',helper)*eta0 + N1 * (xi1 - p1 - 1))
    K1.inv = rwish(M, p2*eta0 + xi1)
    K1 = solve(K1.inv)
    
    # sample K2
    helper = lapply(1:p,function(k) t(L[,,k]) %*% t(K1.inv) %*% (L[,,k]))
    M = solve(Reduce('+',helper)*eta0 + N2 * (xi2 - p2 - 1))
    K2.inv = rwish(M, p1*eta0 + xi2)
    K2 = solve(K2.inv)
    
    # format to be S0 for code
    S0 = kronecker(K2,K1)
    S0.inv = kronecker(K2.inv,K1.inv)
    
    ## sample eta0
    if (df_MH == 0){
      d.v0 = sapply(nu.domain,function(k)
        dwish(V0,S0/k,k,if_log = T) )
      d.v0 = exp(d.v0-max(d.v0))
      probs = d.v0/sum(d.v0)
      eta0 = sample(nu.domain,size = 1,prob = probs)
    } else if (df_MH == 1){
      eta0.star = reflect(sample((eta0 - df.delta[3]):(eta0 + df.delta[3]),1),p+2)
      R = dwish(V0,S0/eta0.star,eta0.star,if_log=T) - 
        dwish(V0,S0/eta0,eta0,if_log=T) +
        dnbinom(eta0.star - (p+2),DF.SIZE,DF.PROB,log=TRUE) -
        dnbinom(eta0 - (p+2),DF.SIZE,DF.PROB,log=TRUE)
      ## if u<r, set eta0 to be eta0.star
      if (log(runif(1))<R){
        eta0 = eta0.star
        eta0.acc = eta0.acc + 1
      }
    }
    
    
    ## compute for output
    for ( j in 1:g ){
      cov.temp[[j]] = (1-pis)*Sig[[j]]+pis*Psi[[j]]
      cov.inv.temp[[j]] = solve(cov.temp[[j]])
    }
    
    ## store output
    if (save_all == 0){
      if((s>burnin)&((s %% thin)==0)){
        # layer 1
        cov.out[index,] = unlist(cov.temp)
        cov.inv[index,] = unlist(cov.inv.temp)
        nu.out[index,] = c(nu0,nu)
        eta0.out[index,] = eta0
        pis.out[index,] = pis
        index = index + 1
      }
    } else {
      # layer 1
      Sig.out[s,]  = unlist(Sig)
      Psi.out[s,] = unlist(Psi)
      cov.out[s,] = unlist(cov.temp)
      cov.inv[s,] = unlist(cov.inv.temp)
      U.out[s,] = unlist(U)
      pis.out[s,] = pis
      # layer 2
      V0.out[s,] = c(V0)
      V1.out[s,] = unlist(V1)
      V2.out[s,] = unlist(V2)
      nu.out[s,] = c(nu0,nu)
      # layer 3- shrink homo to kron
      K1.out[s,] = c(K1)
      K2.out[s,] = c(K2)
      eta0.out[s,] = eta0
    }
    
    
  }
  
  if (save_all == 0){
    out = list("cov.out" = cov.out,
               "cov.inv" = cov.inv,"eta0" = eta0.out,
               "nu" = nu.out,"pis" = pis.out,
               "acc" = acc)
  } else {
    out = list("Sig" = Sig.out, "Psi" = Psi.out,
               "U" = U.out,
               "V0" = V0.out, "V1" = V1.out,
               "V2" = V2.out,"nu" = nu.out,
               "K1" = K1.out, "K2" = K2.out,
               "eta0" = eta0.out,"cov.inv" = cov.inv,
               "pis" = pis.out,"acc" = acc,
               "nu.domain" = nu.domain)
  }
  
  if (df_MH == 0){
    return(c(out,list("acc.dfs" = c(nu.acc,nu0.acc,eta0.acc))))
  } else if (df_MH == 1){
    return(c(out,list("acc.dfs" = c(nu.acc,nu0.acc,eta0.acc))))
  }
  
}

##########################################################
### Alternative Covariance Method Functions
##########################################################
##########################################################

##########################################################
### Sample covariance function
## same as base r cov
# X: n x p matrix that we want a pxp column covariance matrix from
##########################################################
cov.func = function(X){
  X = as.matrix(X)
  N = nrow(X)
  
  X.temp = t(apply(X,1,function(KK) KK - colMeans(X)))
  
  t(X.temp) %*% X.temp / (N - 1)
}
##########################################################
### MLE covariance function
# X: n x p matrix that we want a pxp column covariance matrix from
##########################################################
cov.mle = function(X){
  X = as.matrix(X)
  N = nrow(X)
  
  X.temp = t(apply(X,1,function(KK) KK - colMeans(X)))
  
  t(X.temp) %*% X.temp / (N)
}
##########################################################
### Sample pooled covariance matrix
# X: n x p matrix that we want a pxp column sample covariance matrix from
# group: vector of length n identifying group membership of rows in X
##########################################################
cov.pool = function(X,group){
  # sample pooled covariance matrix
  
  X = as.matrix(X)
  group = as.factor(group)
  N = length(group)
  
  sum.cov = tapply(seq_len(N),group,function(KK)
    length(KK) * cov.mle(X[KK,]))
  
  Reduce("+",sum.cov)/(N-length(unique(group)))
  
}
##########################################################
### MLE pooled covariance matrix
# X: n x p matrix that we want a pxp column sample covariance matrix from
# group: vector of length n identifying group membership of rows in X
##########################################################
cov.pool.mle = function(X,group){
  # MLE pooled covariance matrix
  
  X = as.matrix(X)
  group = as.factor(group)
  N = length(group)
  
  sum.cov = tapply(seq_len(N),group,function(KK)
    length(KK) * cov.mle(X[KK,])) ## sum up Y^TY
  
  Reduce("+",sum.cov)/(N) ## MLE divides by total N- see MKB book
  
}
##########################################################
### Posterior mean of covariance from following hierarchical model:
## Model:
# Y ~ N_(nxp)(0,Sig)
# Therefore, S = Y'Y ~ W(Sig,n)
# Sig ~ IW(S0^(-1)/(nu-p-1),nu)
## Inputs:
# S: sample covariance as described above
# n: sample size
# S0: Wishart prior mean
# nu: Wishart prior degrees of freedom
# de.meaned: logical. If de.meaned = true, use n-1 instead of n in degrees of freedom on S
# loss: "frobenus" or "stein". If loss == stein, return Bayes estimator under Stein loss
##########################################################
cov.shrink.pm = function(S,n,S0,nu,
                         de.meaned = F,
                         loss = "frobenus"){
  nu.star = nu + n
  if ( de.meaned == T ){
    nu.star = nu.star - 1
  }
  
  out = (S0 * (nu - p - 1) + S)/(nu.star - p - 1)
  if ( loss == "stein") {
    out = (S0 * (nu - p - 1) + S) / (nu.star)
  }
  
  return(out)
  
}
##########################################################
### Separable covariance MLE
## block coordinate descent algorithm
# X: p1 x p2 x n array, de-meaned
# output: Psi p1 x p1 row covariance
# output: Sig p2 x p2 column covariance
# de.meaned: logical. If de.meaned = true, use n-1 instead of n in degrees of freedom on S
##########################################################
cov.kron.mle = function(X,itmax = 100,eps = 1e-5,
                        de.meaned = F){

  p1 = dim(X)[1]
  p2 = dim(X)[2]
  N = dim(X)[3] # since removing mean
  
  if (de.meaned == T){
    N = N-1
  }
  
  # initialize sig tilde
  Sig.tilde = matrix(rowMeans(apply(X,3,cov.mle)),ncol = p2)
  Sig.tilde.inv = solve(Sig.tilde)
  
  # initialize stopping checks
  Psi.tilde = Psi.tilde.old = matrix(0,ncol = p1,nrow = p1)
  Sig.tilde.old = matrix(0,ncol = p2, nrow = p2)
  check = F
  it = 0
  while((check == FALSE) & (it < itmax)){
    
    Psi.tilde = matrix(rowSums(apply(X,3,function(KK)KK %*% Sig.tilde.inv %*% t(KK))),
                       ncol = p1)/(N*p2)
    Psi.tilde.inv = solve(Psi.tilde)
    
    Sig.tilde = matrix(rowSums(apply(X,3,function(KK)t(KK) %*% Psi.tilde.inv %*% KK)),
                       ncol = p2)/(N*p1)
    Sig.tilde.inv = solve(Sig.tilde)
    
    if (all(abs(Sig.tilde.old-Sig.tilde)<eps) &
        all(abs(Psi.tilde.old-Psi.tilde)<eps)){
      check = TRUE
    }
    
    # update for next iteration in while loop
    it = it+1 
    Psi.tilde.old = Psi.tilde
    Sig.tilde.old = Sig.tilde
    
  }
  
  return(list("Psi" = Psi.tilde,"Sigma" = Sig.tilde,
              "it" = it))
  
}
##########################################################
### Separable covariance MLE for multiple groups
## block coordinate descent algorithm
# X: p1 x p2 x n array, de-meaned
# group: vector of length n identifying group membership of rows in X
# output: Psi p1 x p1 row covariance
# output: Sig p2 x p2 column covariance
# de.meaned: logical. If de.meaned = true, use n-1 instead of n in degrees of freedom on S
##########################################################
cov.kron.pool.mle = function(X,group,itmax = 100,eps = 1e-5,
                             de.meaned = F){

  # params
  group = as.factor(group)
  
  p1 = dim(X)[1]
  p2 = dim(X)[2]
  N = dim(X)[3]; n = N
  
  if (de.meaned == TRUE){
    N = N - length(unique(group))
  }
  
  # initialize sig tilde
  Sig.tilde = matrix(rowMeans(apply(X,3,cov.mle)),ncol = p2)
  Sig.tilde.inv = solve(Sig.tilde)
  
  # initialize stopping checks
  Psi.tilde = matrix(0,ncol = p1,nrow = p1)
  check = F
  it = 0
  Psi.tilde.old = matrix(0,ncol = p1,nrow = p1)
  Sig.tilde.old = matrix(0,ncol = p2, nrow = p2)
  while((check == FALSE) & (it < itmax)){
    
    Psi.tilde = tapply(seq_len(n),group,function(KK)
      matrix(rowSums(apply(X[,,KK],3,function(MM)MM %*% Sig.tilde.inv %*% t(MM))),
             ncol = p1)
      )
    Psi.tilde = Reduce("+",Psi.tilde)/(N*p2)
    Psi.tilde.inv = solve(Psi.tilde)
    
    Sig.tilde = tapply(seq_len(n),group,function(KK)
      matrix(rowSums(apply(X[,,KK],3,function(MM)t(MM) %*% Psi.tilde.inv %*% MM)),
             ncol = p2)
    )
    Sig.tilde = Reduce("+",Sig.tilde)/(N*p1)
    Sig.tilde.inv = solve(Sig.tilde)
    
    if (all(abs(Sig.tilde.old-Sig.tilde)<eps) &
        all(abs(Psi.tilde.old-Psi.tilde)<eps)){
      check = TRUE
    }
    
    # update for next iteration in while loop
    it = it+1 
    Psi.tilde.old = Psi.tilde
    Sig.tilde.old = Sig.tilde
    
  }
  
  return(list("Psi" = Psi.tilde,"Sigma" = Sig.tilde))
  
}
