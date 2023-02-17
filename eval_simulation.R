###########################
## Evaluate Output for Simulation Study
###########################

###########################
library(stringr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(xtable)
###########################


# ###########################
# suffix = "_homoNOTsep"
# identifier = suffix
# ###########################
# 
# 
# ###########################
# output.filename = paste0("ms_output",suffix,".RDS")
# file.out = readRDS(output.filename)

dims = lapply(dimnames(file.out),function(k)as.numeric(str_sub(k,2,-1))) 
names(dims) = c("g","p","N")
gs = dims[[1]]
ps = dims[[2]]
Ns = dims[[3]]

error.names = names(file.out[1,1,1][[1]])

make_print_plots = 0
###########################


###########################
## Get loss density plots for each scenario in simulation
# for Stein's loss

if ( make_print_plots == 1){
  
  for (n.ind in 1:length(Ns)){
    for (g.ind in 1:length(gs)){
      for (p.ind in 1:length(ps)){
        temp.out = file.out[g.ind,p.ind,n.ind][[1]]
        
        for(j in 1){ 
          # scale columns
          loss.j = temp.out[[j]] %>%  log()
          loss.melt = melt(loss.j)
          
          # plot
          pdf(paste0("lossplot",identifier,"_N",Ns[n.ind],"g",gs[g.ind],"p",ps[p.ind],".pdf"),
              family="Times",height = 7,width = 6.6)
          print(ggplot(data = loss.melt, aes(x = value,color = Var2))+
            geom_density() +
            ggtitle(paste0(error.names[j]," Error Loss")))
          dev.off()
          
          # get estimate
          quantiles = rbind(apply(loss.j,2,function(jj)quantile(jj,c(.1,.5,.9))),
                apply(loss.j,2,mean))
          quantiles
          colnames(quantiles)[order(quantiles[4,])]
       
        }
        
      }
    }
  }
  
}
###########################

###########################
## get 
avg_loss_output = c()
index = 1
for (n.ind in 1:length(Ns)){
  for (g.ind in 1:length(gs)){
    for (p.ind in 1:length(ps)){
      temp = c()
      
      for ( j in 1:3 ){
        
        
        temp.out = file.out[g.ind,p.ind,n.ind][[1]]
        loss.j = temp.out[[j]]
        
        temp = cbind(temp,colMeans(loss.j))
        
      }
      
      colnames(temp) = error.names
      rownames(temp) = colnames(loss.j)
      
      TEMP = temp[c(2,7,3,6,4,5),1] # no post mean
      
      avg_loss_output = rbind(avg_loss_output,TEMP)
      rownames(avg_loss_output)[index] = paste0("J = ",gs[g.ind],"; p = ",ps[p.ind]*3)
      
      index = index + 1
      
    }
  }
}
###########################


###########################
## print output to screen
round(avg_loss_output,2)
# xtable(round(avg_loss_output,2))
###########################

