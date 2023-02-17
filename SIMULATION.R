###########################
## Simulation Study
# loop through parameters and compare distance between truth and estimate
# save loss output
###########################
###########################
## packages
library(foreach)
library(doParallel)
source("./cov_functions.R")
source("./helpers.R")
###########################

###########################
## model parameters
gs = c(4,10)
Ns = c(1)# 
p1s = c(2,4,8)
p2 = 3
data.type = "hetero, not sep" # "homo, sep", "hetero, sep", "homo, not sep"
## file output identifyer
suffix = "_heteroNOTsep"
###########################

###########################
## GS parameters
S = 28000
burnin = 3000
thin = 10
## simulation parameters
sim = 50
###########################


###########################
## run simmy
final.out=array(list(),dim=c(length(gs),length(p1s),length(Ns)))
loss=list()
for (n.ind in 1:length(Ns)){
  for (g.ind in 1:length(gs)){
    for (p1.ind in 1:length(p1s)){
      
      # get changing params
      g = gs[g.ind]
      p1 = p1s[p1.ind]
      N = p1*p2*Ns[n.ind] + 1
      
      # calculate final params
      p = p1*p2
      ns = rep(N,g)
      
      # storage helper
      nopool.collect = array(NA,dim=c(p,p,g))
      
      ### get true covariance
      Sig.true = Sigj.true = list()
      
      set.seed(123)
      if (data.type == "homo, sep"){
        
        RHO = .7
        V1.true = eye(p1)*(1-RHO) + RHO
        RHO = .2
        V2.true = eye(p2)*(1-RHO) + RHO
        V.true = kronecker((V2.true),(V1.true))
        
        for(i in 1:g){
          Sig.true[[i]] = V.true
        }
        
      } else if (data.type == "hetero, sep"){
        
        for(i in 1:g){
          RHO = sample(seq(from = .35,to =.9,length.out = 20),1)
          V1.true = eye(p1)*(1-RHO) + RHO
          RHO = sample(seq(from = .35,to =.9,length.out = 20),1)
          V2.true = eye(p2)*(1-RHO) + RHO
          
          V.true = kronecker((V2.true),(V1.true))
          Sig.true[[i]] = V.true
        }
        
      } else if (data.type == "hetero, not sep") {
        
        for(i in 1:g){
          RHO = sample(seq(from = .35,to =.9,length.out = 20),1)
          Sig.true[[i]] = eye(p)*(1-RHO) + RHO 
        }
        
      } else if (data.type == "homo, not sep") {
        
        RHO = .7
        V0.true = eye(p)*(1-RHO) + RHO
        
        for(i in 1:g){
          Sig.true[[i]] = V0.true
        }
        
      }
      set.seed(Sys.time())
      ###########################
      ## implement GS
      #setup parallel backend to use many processors
      cores=detectCores()
      cl <- makeCluster(cores[1] - 1)  # dont overload your computer
      registerDoParallel(cl)
      
      
      parallel.out <- foreach(sim.ind=1:sim, .combine=cbind) %dopar% {

        output = list()
        
        ###########################
        ## obtain data set
        set.seed(sim.ind)
        Y.list = list()
        for(i in 1:g){  
          Y.list[[i]] = matrix(rnorm(ns[i]*p),nrow=ns[i]) %*% chol(Sig.true[[i]])
        }
        set.seed(Sys.time())
        
        ## remove mean
        de.mean = function(YY){
          nn = nrow(YY)
          one.nxn = matrix(rep(1,nn*nn),ncol=nn)
          YY.o = (eye(nn) - one.nxn/nn) %*% YY
          return(YY.o)
        }
        Y.list = lapply(Y.list,de.mean)
        ###########################
        
        ###########################
        ## run GS for multiple shrinkage
        model = SWAG_GS(S,burnin,thin,save_all = 0)

        ## get estimate under stein's loss
        output$MS.pm = array(colMeans(model$cov.out),dim = c(p,p,g))
        MS.stein.pm.temp = array(colMeans(model$cov.inv),dim = c(p,p,g))
        MS.stein.pm.temp.inv = array(NA,dim=c(p,p,g))
        for ( k in 1:g ){
          MS.stein.pm.temp.inv[,,k] = solve(MS.stein.pm.temp[,,k])
        }
        output$MS.stein.pm = MS.stein.pm.temp.inv
        ###########################
        
        ###########################
        ## reformat data for the cov functions below
        
        ## reformat to matrix with group
        temp = lbind(Y.list)
        group = temp$group
        Y.matrix = temp$mat
        
        ## reformat to p1xp2xN array 
        N = sum(ns)
        Y.array = array(NA,dim=c(p1,p2,N))
        for ( j in 1:N){
          Y.array[,,j] = matrix(Y.matrix[j,],nrow = p1,ncol = p2,byrow = F)
        }
        ###########################
        
        ###########################
        ## no pooling- shrink to decent prior/regularization
        df = p + 4
        for ( j in 1:g){
          sumsq = t(Y.list[[j]]) %*% Y.list[[j]]
          
          nopool.collect[,,j] = cov.shrink.pm(sumsq,ns[j],eye(p),p+2,de.meaned = T)
        }
        
        output$nopool = nopool.collect
        ###########################
          
        ###########################
        ## kron MLE- hetero
        kron.out = lapply(Y.list,function(j)
          cov.kron.mle(vec.inv.array(j,p1,p2), de.meaned = TRUE))
        kron.out = lapply(kron.out,function(j)kronecker(j$Sigma,j$Psi))
        
        output$hetero.kron = list.to.3d.array(kron.out)
        ###########################
        
        ###########################
        ## kron MLE- homo
        temp = cov.kron.pool.mle(Y.array,group, de.meaned = TRUE)
        kron.homo.out = kronecker(temp$Sigma,temp$Psi)
        
        output$homo.kron = rep.array(kron.homo.out,g)
        ###########################
        
        ###########################
        ## homogeneous - pool 
        pool.out = cov.pool(Y.matrix,group)
        
        output$pool = rep.array(pool.out,g)
        ###########################
        
        ###########################
        ## heterogeneous - standard MLE
        standard.mle = tapply(seq_len(sum(ns)),group,function(KK)
          cov.func(Y.matrix[KK,]))
        
        output$mle = list.to.3d.array(standard.mle)
        ###########################
          
        
        ###########################
        ## get distance from each output and the truth
        loss$stein = unlist(lapply(output,function(k)
          mean(sapply(1:g,function(l) loss_stein(k[,,l],Sig.true[[l]] )))))
        loss$sq = unlist(lapply(output,function(k)
          mean(sapply(1:g,function(l) loss_sq(k[,,l],Sig.true[[l]] )))))
        
        ###########################
        

        ###########################
        ## return this output
        output = loss
      
        output
        
      }
      
      #stop cluster
      stopCluster(cl)
      
      # create new structure
      new.parallel.out = list()
      for (j in 1:2){
        temp = matrix(unlist(parallel.out[j,]),ncol = length(parallel.out[j,][[1]]),byrow=T)
        colnames(temp)=names(parallel.out[j,][[1]])
        new.parallel.out$all[[j]] = temp
      }
      names(new.parallel.out$all) = rownames(parallel.out)
      
      # save output
      final.out[g.ind,p1.ind,n.ind] = new.parallel.out

      
    }
  }
}



dimnames(final.out)[[1]]  = paste0("g",gs)
dimnames(final.out)[[2]]  = paste0("p",p1s)
dimnames(final.out)[[3]] = paste0("N",Ns)



## save output
output.filename = paste0("ms_output",suffix,".RDS")
# saveRDS(final.out,file = output.filename)

## summarize output
source("eval_simulation.R")


