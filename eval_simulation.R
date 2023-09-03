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


###########################
suffixes = c("_homosep","_homoNOTsep","_heterosep","_heteroNOTsep")
titles = c("Homogenenous, Kronecker", "Homogeneous, not Kronecker",
           "Heterogeneous, Kronecker","Heterogeneous, not Kronecker")
###########################

loss.collect.top = loss.collect.bot = c()
runtime.collect = c()
for ( suffix in suffixes ){
  
  ###########################
  output.filename = paste0("./output/ms_output_all",suffix,".Rdata")
  load(output.filename)
  file.out = final.out
  
  
  dims = lapply(dimnames(file.out),function(k)as.numeric(str_sub(k,2,-1))) 
  names(dims) = c("g","p","N")
  gs = dims[[1]]
  ps = dims[[2]]
  Ns = dims[[3]]
  
  error.names = names(file.out[1,1,1][[1]])
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
  title = titles[which(suffix == suffixes)]
  ###########################
  
  ###########################
  ## get 
  avg_loss_output = avg_runtime_output = c()
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
        
        TEMP = temp[c(2,7,6,4,5),1] # no post mean; re-order estimates; select stein loss
        
        avg_loss_output = rbind(avg_loss_output,TEMP)
        avg_runtime_output = c(avg_runtime_output,mean(toc.swag[g.ind,p.ind,n.ind][[1]]))
        
        rownames(avg_loss_output)[index] = names(avg_runtime_output)[index] = paste0("J = ",gs[g.ind],"; p = ",ps[p.ind]*3)
  
        index = index + 1
        
      }
    }
  }
  
  ## collect output
  if (suffix %in% c("_homosep","_heterosep")){
    loss.collect.top = cbind(loss.collect.top,
                     round(avg_loss_output,2))
  } else {
    loss.collect.bot = cbind(loss.collect.bot,
                             round(avg_loss_output,2))
  }
  
  runtime.collect = cbind(runtime.collect,
                          round(avg_runtime_output/60,2) ) # in min
  
  
  ###########################
  

  ###########################
  
  
  temp = sig.out[2,2,1][[1]][[1]] ## for the first simulation
  C = ncol(temp)
  
  sig.melt = melt(temp)
  pp = ggplot(sig.melt, 
         aes(x = as.numeric(Var1), y = value )) +
    geom_line() +
    facet_wrap(~as.factor(Var2))+
    theme(legend.position = "none",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    xlab("")  +
    ylab(expression("Select " ~ Sigma ~ "trace plots")) + 
    ggtitle(paste0(title))
    #ylab(expression("Trace plots of" ~ Sigma[1*","~(1*","~2)]))
  
  pdf(paste0("./plots/sigma_trace",suffix,".pdf"),
      family="Times",height = 7,width = 6.6)
  print(pp)
  dev.off()
  
  ###########################
  
  
  
  
  ###########################
  ## Get loss density plots in simulation
  # for Stein's loss
  
  # grab data
  temp.out = file.out[2,2,1][[1]]
  loss.stein = temp.out[[1]] %>% log()
  all.out = loss.stein[,c(2,7,6,4,5)]
  # get labels
  legend.names.edit = c("SWAG","MLE","pool","K","Kpool")
  colnames(all.out) = legend.names.edit
  # color palette with black:
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  # make plot
  loss.melt = melt(all.out)
  
  
  pp  = ggplot(data = loss.melt, aes(x = value,color = Var2,fill=Var2)) +
    geom_density(alpha=.3) +
    ggtitle(title) +
    xlab("") +
    ylab("log Stein loss") +
    scale_fill_manual(values=cbPalette) +
    scale_color_manual(values=cbPalette)
  if (suffix ==  "_homosep") {
    pp = pp + theme(legend.position = c(.2,.7))
  } else {
    pp = pp + theme(legend.position = "none")
  }
  pdf(paste0("./plots/lossplot",suffix,".pdf"),
      family="Times",height = 7,width = 6.6)
  print(pp)
  dev.off()
  
  ###########################
  
  
  ###########################
  ## plot lambda densities
  
  
  temp = lambda.out[2,3,1][[1]]
  lambda.melt = melt(temp)
  
  pp = ggplot(lambda.melt, 
              aes(x = value,color = as.factor(Var2))) +
    geom_density() +
    theme(legend.position = "none") +
    xlab("")  +
    ylab(expression(lambda)) + 
    ggtitle(paste0(title))
  #ylab(expression("Trace plots of" ~ Sigma[1*","~(1*","~2)]))
  
  # pdf(paste0("./plots/lambda_density",suffix,".pdf"),
  #     family="Times",height = 7,width = 6.6)
  print(pp)
  # dev.off()
        
  
  ###########################
  
}

###########################
## xtable output

xtable(rbind(loss.collect.top,NA,loss.collect.bot))

colnames(runtime.collect) = titles
xtable(runtime.collect)


###########################
  
