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
## what to save
output.save = "all" # "loss", "all"
## run cpp?
run.cpp = FALSE
###########################

###########################
# which indices to save output of from covariance
cov.save.inds.labels = list(c(1,1,1),c(2,4,2),c(3,6,4),c(5,5,3)) # row, column, group

###########################

###########################
## GS parameters
S = 10#28000
burnin = 2#3000
thin = 2#10
## simulation parameters
sim = 50
###########################


###########################
## run simmy
final.out = toc.swag = sig.out = lambda.out = array(list(),dim=c(length(gs),length(p1s),length(Ns)))
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
      
      # get index for cov elements to save
      cov.save.inds = c(lapply(cov.save.inds.labels, function(a)
        matrix(1:p^2,ncol=p,nrow=p)[a[1],a[2]]+(p^2)*(a[3]-1))) %>% 
        unlist()
      
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
      
      if (run.cpp == TRUE) {
      
	cpp.init = function() {
        	 library(Rcpp)
        	 sourceCpp("fast_matrix_ops.cpp")
      		 }
	
	clusterCall(cl,cpp.init)

      }
      
      registerDoParallel(cl)
      
      
      parallel.out <- foreach(sim.ind=1:sim, .combine=cbind) %dopar% {
             ## if CPP:      #, .noexport=c("csolve","crwish")) %dopar% {
      

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
        tic.swag = Sys.time()
        model = SWAG_GS(S,burnin,thin,save_all = 0)
        toc.swag = Sys.time() - tic.swag

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
        tic.nopool = Sys.time()
        df = p + 4
        for ( j in 1:g){
          sumsq = t(Y.list[[j]]) %*% Y.list[[j]]
          
          nopool.collect[,,j] = cov.shrink.pm(sumsq,ns[j],eye(p),p+2,de.meaned = T)
        }
        toc.nopool = Sys.time() - tic.nopool
        output$nopool = nopool.collect
        ###########################
          
        ###########################
        ## kron MLE- hetero
        tic.kron = Sys.time()
        kron.out = lapply(Y.list,function(j)
          cov.kron.mle(vec.inv.array(j,p1,p2), de.meaned = TRUE))
        kron.out = lapply(kron.out,function(j)kronecker(j$Sigma,j$Psi))
        toc.kron = Sys.time() - tic.kron
        output$hetero.kron = list.to.3d.array(kron.out)
        ###########################
        
        ###########################
        ## kron MLE- homo
        tic.kron.homo = Sys.time()
        temp = cov.kron.pool.mle(Y.array,group, de.meaned = TRUE)
        kron.homo.out = kronecker(temp$Sigma,temp$Psi)
        toc.kron.homo = Sys.time() - tic.kron.homo
        output$homo.kron = rep.array(kron.homo.out,g)
        ###########################
        
        ###########################
        ## homogeneous - pool 
        tic.pool = Sys.time()
        pool.out = cov.pool(Y.matrix,group)
        toc.pool = Sys.time() - tic.pool
        output$pool = rep.array(pool.out,g)
        ###########################
        
        ###########################
        ## heterogeneous - standard MLE
        tic.mle = Sys.time()
        standard.mle = tapply(seq_len(sum(ns)),group,function(KK)
          cov.func(Y.matrix[KK,]))
        toc.mle = Sys.time() - tic.mle
        
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
        ## collect run time in a vector
        toc.out = c(toc.swag,
                    toc.nopool, 
                    toc.kron,toc.kron.homo,
                    toc.pool,
                    toc.mle)
        ###########################
        
        
        ###########################
        ## return this output
        if (output.save == "loss"){
          output = loss
        } else if (output.save == "all"){
          output = list("loss" = loss,
                        "toc" = toc.swag, ## just return swag toc
                        "SWAG_sigma" = model$cov.out[,cov.save.inds],
                        "SWAG_lambda" = c(model$pis))
        }
      
        output
        
      }
      
      #stop cluster
      stopCluster(cl)
      
      if (output.save == "loss"){
        # create new structure
        new.parallel.out = list()
        for (j in 1:2){
          temp = matrix(unlist(parallel.out[j,]),
                        ncol = length(parallel.out[j,][[1]]),byrow=T)
          colnames(temp)=names(parallel.out[j,][[1]])
          new.parallel.out$all[[j]] = temp
        }
        names(new.parallel.out$all) = rownames(parallel.out)
        
        # save output
        final.out[g.ind,p1.ind,n.ind] = new.parallel.out
      } else if (output.save == "all") {
        
        ## collect loss- match format above to align with future output
        new.parallel.out = list()
        for (j in 1:2){
          temp = lapply(parallel.out[1,],function(k)k[[j]])
          temp = matrix(unlist(temp),
                        ncol = length(temp[[1]]),
                        byrow=T)
                               
          colnames(temp)=names(parallel.out[1,1][[1]][[1]])
          new.parallel.out$all[[j]] = temp
        }
        names(new.parallel.out$all) = names(parallel.out[1,][[1]])
        # save loss output
        final.out[g.ind,p1.ind,n.ind] = new.parallel.out
        
        ## save time out
        # toc.mean = Reduce("+",parallel.out[2,])/sim
        # names(toc.mean) = c("MS","nopool","hetero.kron","homo.kron","pool","mle")
        all = list(); all[[1]] = unlist(parallel.out[2,])
        toc.swag[g.ind,p1.ind,n.ind] = all
        
        ## save sig
        all = list(); all[[1]] =  parallel.out[3,]
        sig.out[g.ind,p1.ind,n.ind] = all
        
        ## save lambda 
        temp = matrix(unlist(parallel.out[4,]),
                            ncol = ncol(parallel.out),
                            byrow = T)
        colnames(temp) = paste0("result.",1:ncol(parallel.out))
        all = list(); all[[1]] =  temp
        lambda.out[g.ind,p1.ind,n.ind] = all
      }
      
      print(paste0("Finsished running: n: ",Ns[n.ind],
                   ", p: ",p1s[p1.ind],
                   ", g: ",gs[g.ind],
                   "!!!!!!!!"))
      
      
    }
  }
}

print("Saving output now !!!")

dimnames(final.out)[[1]]  = paste0("g",gs)
dimnames(final.out)[[2]]  = paste0("p",p1s)
dimnames(final.out)[[3]] = paste0("N",Ns)



## save output
if(output.save == "loss"){
  output.filename = paste0("ms_output",suffix,".RDS")
  saveRDS(final.out,file = output.filename)
} else if (output.save == "all"){
  output.filename = paste0("ms_output_all",suffix,".Rdata")
  save(final.out,toc.swag,lambda.out,sig.out,
       file = output.filename)
}

## summarize output
# source("eval_simulation.R")


