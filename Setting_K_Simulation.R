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
g = 4
p1 = 8; p2 = 3; p = p1*p2
N = 2*p # p + 1
ns = rep(N,g)

Ks = c(p+2,2*p,3*p,4*p,5*p)

data.type = "hetero, not sep" # "homo, sep", "hetero, sep", "homo, not sep"

## what to save
output.save = "loss" 
###########################

###########################
print(paste0("Results for p = ",p,", N = ",N,"!!!"))
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

if (N == p * 2){
  n.suffix = "_n2P"
} else if (N == p+1){
  n.suffix = "_nPplus1"
} else {
  stop("Specify n.suffix!")
}
###########################


###########################
## run simmy

## prepare files to save output
final.out = array(list(),dim=c(length(Ks)))
loss=list()


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
  ## run GS for multiple shrinkage for each K
  
  for ( df_ss in 1:length(Ks)){
    
    model = SWAG_GS(S,burnin,thin,save_all = 0,
                    K = Ks[df_ss])
    
    ## get estimate under stein's loss
    output$MS.pm = array(colMeans(model$cov.out),dim = c(p,p,g))
    MS.stein.pm.temp = array(colMeans(model$cov.inv),dim = c(p,p,g))
    MS.stein.pm.temp.inv = array(NA,dim=c(p,p,g))
    for ( k in 1:g ){
      MS.stein.pm.temp.inv[,,k] = solve(MS.stein.pm.temp[,,k])
    }
    
    output[[df_ss]] = MS.stein.pm.temp.inv
    
  }

  ###########################

  ###########################
  ## compute loss
  
  loss$stein = unlist(lapply(output,function(k)
    mean(sapply(1:g,function(l) loss_stein(k[,,l],Sig.true[[l]] )))))
  
  ###########################
  
  
  ###########################
  ## return this output
  
  output = loss
  
  output
  
}

#stop cluster
stopCluster(cl)

final.out = matrix(unlist(parallel.out),ncol = length(Ks),
                   byrow = T)


print("Saving output now !!!")

## save output
output.filename = paste0("./output/K_simulation",suffix,"_p",p,n.suffix,".RDS")
saveRDS(final.out,file = output.filename)

print(paste0("Saved output at file: ",output.filename,"!!!"))


