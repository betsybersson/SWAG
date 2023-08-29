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
# ###########################
# 
# 
# ###########################
# output.filename = paste0("ms_output",suffix,".RDS")
# file.out = readRDS(output.filename)

file.out = final.out

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
      
      for ( j in 1:2 ){
        
        
        temp.out = file.out[g.ind,p.ind,n.ind][[1]]
        loss.j = temp.out[[j]]
        
        temp = cbind(temp,colMeans(loss.j))
        
      }
      
      colnames(temp) = error.names
      rownames(temp) = colnames(loss.j)
      
      TEMP = temp[c(2,7,3,6,4,5),1] # no post mean; re-order estimates; select stein loss
      
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
## set ggplot theme
theme_set(theme_bw())
theme_update(text=element_text(family = "Times",size=30),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank(),
             plot.title = element_text(hjust = 0.5),
             legend.title=element_blank(),
             legend.position = c(0.15, 0.8),
             axis.text = element_text(size=25))
## get better labels for plots
suffixes = c("_homosep","_homoNOTsep","_heterosep","_heteroNOTsep")
titles = c("Homogenenous, Kronecker", "Homogeneous, not Kronecker",
           "Heterogeneous, Kronecker","Heterogeneous, not Kronecker")
j = which(suffix == suffixes)
###########################


###########################
## summarise toc
ggplot(data.frame("toc" = toc.swag/60),aes(x = toc))+
  geom_density(stat="density") +
  xlab("Run Time (minutes)") +
  ylab("Frequency") +
  ggtitle(titles[j])

###########################


###########################
## plot trace plots of Sig
C = ncol(sig.out[[1]])
for ( c in 1:C){
  sig.out.c = matrix(unlist(lapply(sig.out,function(kk)kk[,c])),
                     byrow=F,
                     ncol = sim)
  sig.melt.c = melt(sig.out.c)
  ggplot(filter(sig.melt.c,Var2 %in% sample(1:sim,4)), # randomly select output from 4 simulations
       aes(x = as.numeric(Var1), y = value )) +
    geom_line() +
    facet_wrap(~as.factor(Var2))+
    theme(legend.position = "none",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    xlab("")  +
    ylab(expression("Trace plots of" ~ Sigma[1*","~(1*","~2)]))
}


###########################




###########################
## Get loss density plots in simulation
# for Stein's loss

# grab data
temp.out = file.out[2,2,1][[1]]
loss.stein = temp.out[[1]] %>% log()
all.out = list()
all.out[[j]] = loss.stein[,c(2,7,6,4,5)]
# get labels
legend.names.edit = c("SWAG","MLE","pool","K","Kpool")
colnames(all.out[[1]]) = legend.names.edit
# color palette with black:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# make plot
loss.melt = melt(all.out[[j]])

pp  = ggplot(data = loss.melt, aes(x = value,color = Var2,fill=Var2)) +
  geom_density(alpha=.3) +
  ggtitle(titles[j]) +
  xlab("") +
  ylab("log Stein loss") +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(values=cbPalette)
if ( j == 1) {
  pp = pp
} else {
  pp = pp + theme(legend.position = "none")
}
print(pp)
  