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

make_plots = 0
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

## print output to screen
round(avg_loss_output,2)

###########################


###########################
## Get loss density plots in simulation
# for Stein's loss
if ( make_plots == 1){
  suffixes = c("_homosep","_homoNOTsep","_heterosep","_heteroNOTsep")
  titles = c("Homogenenous, Kronecker", "Homogeneous, not Kronecker",
             "Heterogeneous, Kronecker","Heterogeneous, not Kronecker")
  identifier = suffixes
  
  all.out = list()
  for ( j in 1:length(suffixes)){
    temp.out = file.out[2,2,1][[1]]
    loss.stein = temp.out[[1]] %>% log()
    
    all.out[[j]] = loss.stein[,c(2,7,6,4,5)]
  }
  
  names = colnames(all.out[[1]])
  legend.names.edit = c("SWAG","MLE","pool","K","Kpool")
  colnames(all.out[[1]]) = legend.names.edit
  
  # The palette with black:
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  for ( j in 1:4){
    loss.melt = melt(all.out[[j]])
    theme_set(theme_bw())
    theme_update(text=element_text(family = "Times",size=30),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 plot.title = element_text(hjust = 0.5),
                 legend.title=element_blank(),
                 legend.position = c(0.15, 0.8),
                 axis.text = element_text(size=25))
    if(j ==1){
      print(ggplot(data = loss.melt, aes(x = value,color = Var2,fill=Var2))+
              geom_density(alpha=.3) +
              ggtitle(titles[j]) +
              xlab("") +
              ylab("log Stein loss")+
              scale_fill_manual(values=cbPalette)+
              scale_color_manual(values=cbPalette))
    }else{
      print(ggplot(data = loss.melt, aes(x = value,color = Var2,fill=Var2))+
              geom_density(alpha=.3, show.legend = FALSE) +
              ggtitle(titles[j]) +
              xlab("") +
              ylab("log Stein loss")+
              scale_fill_manual(values=cbPalette)+
              scale_color_manual(values=cbPalette))
    }
  }
  
}
###########################
