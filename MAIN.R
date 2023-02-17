###########-------------------------------
## Run SWAG model on generic data set
###########-------------------------------
###########-------------------------------
## packages
source("cov_functions.R")
source("helpers.R")
###########-------------------------------

###########-------------------------------
## read in data
load("forest_data.Rdata") # upload data saved as Y.list, dims p1, p2
###########-------------------------------

###########-------------------------------
## extract parameters
p = p1*p2
g = length(Y.list)
ns = unlist(lapply(Y.list,nrow))
###########-------------------------------

###########-------------------------------
## Run GS

model = SWAG_GS(33000,3000,20,
                mh.delta.star = .1,
                save_all = 0)

model$acc/S

###########-------------------------------

###########-------------------------------
## get posterior covariance estimates

# get stein bayes estimates
stein.pm.temp = array(colMeans(model$cov.inv),dim = c(p,p,g))
stein.pm = array(NA,dim=c(p,p,g))
for ( k in 1:g ){
  stein.pm[,,k] = (solve(stein.pm.temp[,,k]))
}

# get K1, K2 estimates
K1 = array(colMeans(model$K1),dim =c(p1,p1))
K2 = array(colMeans(model$K2),dim =c(p2,p2))

K1.temp = matrix(0,nrow=p1,ncol=p1)
K2.temp = matrix(0,nrow=p2,ncol=p2)
for ( k in 1:nrow(model$K1)){
  K1.temp = K1.temp + solve(matrix(model$K1[k,],nrow = p1, ncol = p1))
  K2.temp = K2.temp + solve(matrix(model$K2[k,],nrow = p2, ncol = p2))
}
K1 = solve(K1.temp)
K2 = solve(K2.temp)
###########-------------------------------

###########-------------------------------
## Plot and Explore Output

## plot density of weight and mass of degrees of freedom
par(mfrow=c(2,2), mai = c(.4, 0.4, 0.4, 0.4))
w.star = model$pis
S = length(w.star)

plot(density(w.star),
     xlab = "",ylab = "", main=expression(paste("Density of ",lambda)),
     xlim = c(0,1))

plot(table(model$nu[,1])/S,
     xlim = c(p+2,p+50),
     xaxt = "n",
     main = expression(paste("Probability mass of ", nu)),ylab ="")
axis(side = 1,at = c((p+2):(p+50)),labels = c((p+2):(p+50)))

plot(table(model$nu[,2])/S,
     xlim = c(p+2,p+50),
     xaxt = "n",
     main = expression(paste("Probability mass of ", gamma)),ylab ="")
axis(side = 1,at = c((p+2):(p+50)),labels = c((p+2):(p+50)))

plot(table(model$eta0)/S,
     xlim = c(p+2,p+50),
     xaxt = "n",
     main = expression(paste("Probability mass of ", xi)),ylab ="")
axis(side = 1,at = c((p+2):(p+50)),labels = c((p+2):(p+50)))

## plot covariances
# names for plots
row.names = 1:p1
col.names = 1:p2
for ( j in 1:g){
  
  if(is.element(j,1)){
    par(mfrow=c(1,1), mai = c(1, 1, 0.6, 0.05))
  } else {
    par(mfrow=c(1,1), mai = c(1, .05, 0.6, 0.05))
  }
  mainlab = bquote(~ hat(Sigma)[.(group.names[j])])
  
  image((stein.pm[,(p1*p2):1,j]), 
        col = rev(brewer.pal(9,"YlOrRd")), xaxt="n",yaxt="n")
  
  if(is.element(j,1)){
    axis(side=2,at=seq(1,0,length=p),
         labels=rep(row.names,times=3),tick=FALSE,las=1 ,cex.axis=1.25) 
    axis(2, at = seq(1,0,length=p)[which(rep(row.names,times=3)=="MPB")],
         labels=col.names,tick=FALSE,las=3 ,cex.axis=1.25,line=2.5) }
  axis(side=1,at=seq(1,0,length=p),
       labels=rep(row.names,times=3)[p:1],tick=FALSE,las=3 ,cex.axis=1.25) 
  axis(1, at = seq(1,0,length=p)[which(rep(row.names,times=3)=="MPB")],
       labels=col.names[p2:1],tick=FALSE,las=1 ,cex.axis=1.25,line=2.5) 
  
  mtext(side=3,mainlab,line=.2,cex=2 )
  
  sep<-seq(1.04,-.04,length=4) 
  abline(h=sep,col=gray(.25)) 
  abline(v=sep,col=gray(.25))

}

par(mfrow=c(2,1))
image((K1[,(p1):1]), 
      col = rev(brewer.pal(9,"YlOrRd")), xaxt="n",yaxt="n")
axis(side=1,at=seq(0,1,length=p1),
     labels=row.names,tick=FALSE,las=3 ,cex.axis=1.25)
axis(side=2,at=seq(1,0,length=p1),
     labels=row.names,tick=FALSE,las=1 ,cex.axis=1.25) 
mtext(side=3,expression(hat(P)[1]),line=.2 ,cex=2)

image((K2[,(p2):1]), 
      col = rev(brewer.pal(9,"YlOrRd")), xaxt="n",yaxt="n")
axis(side=1,at=seq(0,1,length=p2),
     labels=col.names,tick=FALSE,las=3 ,cex.axis=1.25)
axis(side=2,at=seq(1,0,length=p2),
     labels=col.names,tick=FALSE,las=1 ,cex.axis=1.25) 
mtext(side=3,expression(hat(P)[2]),line=.2 ,cex=2)
