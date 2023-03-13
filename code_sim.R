# supp material to manuscript 
# "A Riemann Manifold Model Framework for Longitudinal Changes in Physical Activity"
# example code for simulations

rm(list = ls())
set.seed(100)
library(fdapace)

#--------------------------------------------------------------------------
# construct PCs; time window = 7am - 21pm (for demonstrative purposes only)
#--------------------------------------------------------------------------

# PC 1: construct momenta in x and y directions
PC1.x <- rep(0, 60 * 14)
# smooth y direction
PC1.y <- c(rep(0, 60*3), rep(0.003, 60*6), rep(0, 60*5))
plot(smooth.spline(PC1.y, df = 8))
PC1.y <- smooth.spline(PC1.y, df = 8)$y
# assemble PC1
PC1 <- c(PC1.x, PC1.y)


#--------------------------------------------------------------------------
# construct PC 2 as local increase (in the morning)
PC2.x <- rep(0, 60 * 14)
PC2.y <- c(rep(0, 60*3), rep(0.003, 60*1), rep(0, 60*10))
plot(smooth.spline(PC2.y, df = 20))
PC2.y <- smooth.spline(PC2.y, df = 20)$y
PC2 <- c(PC2.x, PC2.y)

#--------------------------------------------------------------------------
# construct PC 3 as temporal shift (x direction)
PC3.y <- rep(0, 14 * 60)
PC3.x <- c(rep(0, 60*3), rep(0.004, 60*2), rep(0, 60*9))
PC3.x <- smooth.spline(PC3.x, df = 10)$y
PC3 <- c(PC3.x, PC3.y)


#--------------------------------------------------------------------------
# PC.sim are the RAW PCs before orthogonalization
PC.sim <- cbind(PC1, PC2, PC3)
dim(PC.sim)
PC.sim <- scale(PC.sim, center = F)


#--------------------------------------------------------------------------
# make PCs orthogonal by generating data from them and extract PCs again
#--------------------------------------------------------------------------
nsim <- 500
n.pc <- dim(PC.sim)[2]

# generating coefficients of the PC scores
coefs <- matrix(rnorm(n = n.pc * nsim , mean = 0, sd = 1), nrow = nsim)
coefs <- coefs %*% diag(c(2.5, 1.5, 1))
dim(coefs)

data.sim <- coefs %*% t(PC.sim)

# apply fPCA to data
s <- c(1: dim(data.sim)[2])
L3 <- MakeFPCAInputs(IDs = rep(1:dim(data.sim)[1], each=length(s)), tVec=rep(s,dim(data.sim)[1]), t(data.sim))
FPCAdense <- FPCA(L3$Ly, L3$Lt, optns = list(usergrid = T))

FPCAdense$obsGrid

plot(FPCAdense)

dim(FPCAdense$phi)
PC.new <- scale(FPCAdense$phi, center = F)

#------------------------------------------------
#  construct baseline PA curves from MENU study
baseline_mean_x <- read.csv('/MENU/baseline_mean_x.csv');
baseline_mean_y <- read.csv('/MENU/baseline_mean_y.csv');

baseline_mean_x <- baseline_mean_x[, 2]
baseline_mean_y <- baseline_mean_y[, 2] 



#--------------------------------------------------------
# construct subject-specific deformations from the PCs
#--------------------------------------------------------
nsim <- 100 
n.pc <- dim(PC.new)[2]

coefs <- matrix(rnorm(n = n.pc * nsim , mean = 0, sd = 1), nrow = nsim)
coefs <- coefs %*% diag(c(0.8, 1.2, 1.5))
coefs <- coefs / 2000




#========>>>>>> deform baseline PA curves to get follow-up curves in matlab ================>>>>>>




#---------------------------------------------------------
# fPCA on estimated deformation momenta (proposed method)
#---------------------------------------------------------

# apply fPCA
# PC scores estimated with estimated momenta
s <- c(1: 1680)
L3 <- MakeFPCAInputs(IDs = rep(1:dim(momenta.est)[1], each=length(s)), tVec=rep(s,dim(momenta.est)[1]), t(momenta.est))
FPCAdense <- FPCA(L3$Ly, L3$Lt, optns = list(usergrid = T))

plot(FPCAdense)
FPCAdense$cumFVE
PC.est <- scale(FPCAdense$phi,center = F)
PC.est <- PC.est[, c(3,2,1)] # re-order PCs according to variance


# PC scores estimated with actual momenta
s <- c(1: 1680)
L.a <- MakeFPCAInputs(IDs = rep(1:dim(momenta.actual)[1], each=length(s)), tVec=rep(s,dim(momenta.actual)[1]), t(momenta.actual))
FPCAdense.a <- FPCA(L.a$Ly, L.a$Lt, optns = list(usergrid = T))

plot(FPCAdense.a)
PCscores.a <- FPCAdense.a$xiEst
PC.new.a <- FPCAdense.a$phi[, c(3,2,1)]



#----------------------------------------
# fPCA on vertical difference only
#----------------------------------------
target.all <- matrix(0, nsim, length(baseline_mean))
target.i.all <- matrix(0, nsim, length(baseline_mean))
diff.all <- matrix(0, nsim, length(baseline_mean_y))

# calculate vertical difference between baseline and follow-up PA curves
for (id in 1: nsim) {
  target <- read.csv(file = paste0('results/target_', id, '.csv'), header = F)
  target.all[id, ] <- c(target[,1], target[,2])

  # interpolated target
  target.ynew <- approx(x=target[,1], y=target[,2], xout = baseline_mean_x)
  target.i <- c(target.ynew$x,target.ynew$y)
  target.i.all[id,] <- target.i
  
  # diff
  diff <- target.ynew$y - baseline_mean_y
  diff.all[id,] <- diff
}

# apply PCA to vertical differences
s <- c(1: 840) # only vertical direction
L3 <- MakeFPCAInputs(IDs = rep(1:dim(diff.all)[1], each=length(s)), tVec=rep(s,dim(diff.all)[1]), t(diff.all))
FPCAdense.diff <- FPCA(L3$Ly, L3$Lt, optns = list(usergrid = T))
plot(FPCAdense.diff)




#----------------------------------------
# with Wrobel's time warping method
#----------------------------------------
source("warp_subj_fn.R")

# load target and run the warping alg in Wrobel
res.all <-list()
# registr all subjects
for (i in 1: nsim) {
  target <- read.csv(file = paste0("results/target_", i, ".csv"), header=F)
  test <- warp_subj(baseline, target, format = T, scale = "global", standardize = T)
  res.all[[i]] <- test
}


#----------------------------------------------------------------------
# bind all warping functions from subjects
newtime.all <- NULL
for (i in 1:nsim){
  newtime.all <- rbind(newtime.all, res.all[[i]]$baseline.time.new)
}
dim(newtime.all)


#----------------------------------------------------------------------
# extract fpca from the warping functions
#----------------------------------------------------------------------
s <- test$control.pts
L3 <- MakeFPCAInputs(IDs = rep(1:nsim, each=length(s)), tVec=rep(s,nsim), t(newtime.all))
FPCAdense <- FPCA(L3$Ly, L3$Lt, optns = list(usergrid = T))
plot(FPCAdense)

# save PC score for PC1
score.pc1t <- FPCAdense$xiEst[,1]

# plot PC representing time warping
pc1.sd <- sqrt(FPCAdense$lambda[1])
c <- 1
pc1.c1 <- c* pc1.sd *pc1.w + FPCAdense$mu
pc1.c2 <- -1*c* pc1.sd *pc1.w + FPCAdense$mu

plot(baseline[,1], baseline[,2], type = 'l', col = "blue")
points(pc1.c1, baseline[,2], type = 'l', col = "red")
points(pc1.c2, baseline[,2], type = 'l', col = "green")

# combine new times with baseline
pc1.time <- data.frame(pc1.c1, pc1.c2, baseline)
rownames(pc1.time) <- NULL
colnames(pc1.time) <- NULL


#---------------------------------------------------------------------------
# calculate vertical change wrt WARPED time and extract PCs
#---------------------------------------------------------------------------
diff.y.all <- NULL
for (i in 1:nsim){
  diff.y.all <- rbind(diff.y.all, res.all[[i]]$diff.y)
}
dim(diff.y.all)

# extract fpca from the warping functions
s <- test$control.pts
L3 <- MakeFPCAInputs(IDs = rep(1:nsim, each=length(s)), tVec=rep(s,nsim), t(diff.y.all))
FPCAdense <- FPCA(L3$Ly, L3$Lt, optns = list(usergrid = T))
plot(FPCAdense)

# save PC score for PC1.m (m means magnitude after warping)
score.pc1m <- FPCAdense$xiEst[,1]
score.pc2m <- FPCAdense$xiEst[,2]

pc1.diff <- FPCAdense$phi[,1]
pc2.diff <- FPCAdense$phi[,2]

pc1.diff.sd <- sqrt(FPCAdense$lambda[1])
pc2.diff.sd <- sqrt(FPCAdense$lambda[2])

c <- 1
pc1.c1 <- c* pc1.diff.sd *pc1.diff + FPCAdense$mu
pc1.c2 <- -1*c* pc1.diff.sd *pc1.diff + FPCAdense$mu

pc2.c1 <- c* pc2.diff.sd *pc2.diff + FPCAdense$mu
pc2.c2 <- -1*c* pc2.diff.sd *pc2.diff + FPCAdense$mu

# PC1 (magnitude after warping) plot
plot(baseline[,1], baseline[,2], type = 'l', col = "blue")
points(test$control.pts, baseline[,2]+pc1.c1, type = 'l', col = "red")
points(test$control.pts, baseline[,2]+pc1.c2, type = 'l', col = "green")

# PC2 (magnitude after warping) plot
plot(baseline[,1], baseline[,2], type = 'l', col = "blue")
points(test$control.pts, baseline[,2]+pc2.c1, type = 'l', col = "red")
points(test$control.pts, baseline[,2]+pc2.c2, type = 'l', col = "green")













