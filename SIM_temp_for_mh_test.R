###########################
## Simulation Study
# loop through parameters and compare distance between truth and estimate
# save loss output
###########################
###########################
## packages
library(foreach)
library(doParallel)
library(dplyr)
source("./cov_functions.R")
source("./helpers.R")
###########################

###########################
## model parameters
gs = c(4)
Ns = c(1) # a multiplier times P
p1s = c(4)
p2 = 3
data.type = "homo, sep" # "homo, sep", "hetero, sep", "homo, not sep"
## file output identifyer
suffix_additional = ""
## what to save
output.save = "loss" # loss, all
###########################

###########################
# which indices to save output of from covariance
cov.save.inds.labels = list(c(1,1,1),c(2,4,2),c(3,6,4),c(5,5,3)) # row, column, group
###########################

###########################
## GS parameters
S = 48000
burnin = 8000
thin = 10
## simulation parameters
sim = 50
###########################

###########################
## get suffix for filename

if (data.type == "homo, sep"){
  
suffix = "_homosep"
  
} else if (data.type == "hetero, sep"){
  
suffix = "_heterosep"
  
} else if (data.type == "hetero, not sep") {
  
suffix = "_heteroNOTsep"
  
} else if (data.type == "homo, not sep") {
  
suffix = "_homoNOTsep"
  
}
###########################

###########################
print(paste0("Running regime: ",suffix,"!!!!"))
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
      # setup parallel backend to use many processors
      cores=detectCores()
      cl <- makeCluster(cores[1])
      
      registerDoParallel(cl)
      
      parallel.out <- foreach(sim.ind=1:sim, .combine=cbind) %dopar% {

        ###########################
        ## obtain data set
        set.seed(sim.ind)
        Y.list = list()
        for(i in 1:g){  
          Y.list[[i]] = matrix(rnorm(ns[i]*p),nrow=ns[i]) %*% chol(Sig.true[[i]])
        }
        
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
        model = SWAG_GS(S,burnin,thin,save_all = 0,df_MH = 1)
        toc.swag = Sys.time() - tic.swag

        ## get estimate under stein's loss
        MS.stein.pm.temp = array(colMeans(model$cov.inv),dim = c(p,p,g))
        MS.stein.pm.temp.inv = array(NA,dim=c(p,p,g))
        for ( k in 1:g ){
          MS.stein.pm.temp.inv[,,k] = solve(MS.stein.pm.temp[,,k])
        }
        MS.stein.pm = MS.stein.pm.temp.inv
        ###########################
        
        ###########################
        loss = mean(sapply(1:g, function(l) loss_stein(MS.stein.pm[,,l],Sig.true[[l]] )))
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
        new.parallel.out$stein = unlist(parallel.out)
        # save output
        final.out[g.ind,p1.ind,n.ind] = new.parallel.out
      } else if (output.save == "all") {
        
        ## collect loss- match format above to align with future output
        # create new structure
        new.parallel.out = list()
        new.parallel.out$stein = unlist(parallel.out[1,])
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
  output.filename = paste0("./output/check_MH_ms_output",suffix,suffix_additional,".RDS")
  saveRDS(final.out,file = output.filename)
} else if (output.save == "all"){
  output.filename = paste0("./output/ms_output_all",suffix,suffix_additional,".Rdata")
  save(final.out,toc.swag,lambda.out,sig.out,
       file = output.filename)
}

print(paste0("Saved output at file: ",output.filename,"!!!"))

