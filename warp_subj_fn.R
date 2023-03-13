#----------------------------------------------------------------
# function to warp baseline to target curve using Wrobel's method
#----------------------------------------------------------------
# input:
# baseline: n_time * 2
# target: n_time * 2

# output:
# tstar (original time) vs. t_hat (warped/new time)

#format = F
#scale = "global"
#standardize = T

warp_subj <- function(baseline, target, format = T, scale = "global", standardize = FALSE) {
  rownames(baseline) <- NULL
  colnames(baseline) <- c("index", "value")
  colnames(target) <- c("index", "value")
  baseline <- data.frame(baseline)
  target <- data.frame(target)
  baseline$id <- 1
  target$id <- 2
  
  # save original magnitudes without standardization
  baseline.y <- baseline$value
  target.y <- target$value
    
  # values for scaling
  value.m <- max(baseline$value, target$value)
  value.m.b <- max(baseline$value)
  value.m.t <- max(target$value)
  
  if (format == T) {
    # re-scale time to [0,1]
    baseline$index <- (baseline$index+1)/2
    target$index <- (target$index+1)/2
    
    baseline$index[baseline$index < 0] <- 0
    target$index[target$index < 0] <- 0
    baseline$index[1] <- 0
    target$index[1] <- 0
  }
  
  if (standardize == T) {
    baseline$value <- scale(baseline$value, center = TRUE, scale = TRUE)
    target$value <- scale(target$value, center = TRUE, scale = TRUE)
  } else if (scale == "global") { # apply scaling only
    baseline$value <- baseline$value/value.m 
    target$value <- target$value/value.m
  } else {
    baseline$value <- baseline$value/value.m.b 
    target$value <- target$value/value.m.t
  }
  
  # before register
  plot(baseline$index, baseline$value, type = 'l', col = "blue")
  points(target$index, target$value, type = 'l', col = "red")
  legend("topright", legend=c("baseline", "target"), col = c("blue", "red"), lty = c(1,1))
  
  # map baseline to target
  # ! use npc_criterion instead of npc, the latter does not seem to work
  reg = register_fpca(Y = baseline, Y_template = target, npc_criterion=0.95, family = "gaussian", fpca_type = "variationalEM", max_iterations = 10)
  
  # map target to baseline
  #reg = register_fpca(Y = target, Y_template = baseline, npc = 2, family = "gaussian", fpca_type = "variationalEM", max_iterations = 10)
  
  newdata <- data.frame(reg$Y)
  # after register
  plot(baseline$index, baseline$value, type = 'l', col = "blue", xaxt="n", xlab = "Time", ylab = "PA Magnitude")
  points(target$index, target$value, type = 'l', col = "red")
  points(newdata$t_hat, newdata$value, type = 'l', col= "orange")
  legend("topright", legend=c("baseline", "target", "registered"), col = c("blue", "red", "orange"), lty = c(1,1))
  
  xtick <- seq(from = 1, to = dim(baseline)[1], by = dim(baseline)[1]/14*3)
  xtick <- (xtick - 1)/dim(baseline)[1]
  
  ticks.label <- c("7:00", "10:00", "13:00", "16:00", "19:00")
  axis(side=1, at=xtick, labels = ticks.label)
  

  #------------------------------------------------
  # control points in original scale (-1 to 1)
  # control points are also the old time index of baseline
  control.pts <- baseline$index * 2 - 1

  baseline.time.new <- newdata$t_hat*2 - 1
  target.time.new <- target$index*2 - 1
  #target.new <- data.frame(index = target.time.new, value = target.y)
  
  # interpolate baseline curve at control points 
  baseline.newy.ctrl <- approx(x=baseline.time.new, y=baseline.y, xout = control.pts)$y
  target.newy.ctrl <- approx(x=target.time.new, y=target.y, xout = control.pts)$y
  
  plot(control.pts, baseline.y, type = 'l', col="blue", xaxt="n", xlab = "Time", ylab = "PA Magnitude")
  points(control.pts, target.newy.ctrl, type = 'l', col="red")
  points(control.pts, baseline.newy.ctrl, type = 'l', col="orange")
  legend("topright", legend=c("baseline", "target", "registered"), col = c("blue", "red", "orange"), lty = c(1,1))
  
  xtick <- seq(from = 1, to = dim(baseline)[1], by = dim(baseline)[1]/14*3)
  xtick <- (xtick - (1+dim(baseline)[1])/2)/((dim(baseline)[1]-1)/2)
  
  ticks.label <- c("7:00", "10:00", "13:00", "16:00", "19:00")
  axis(side=1, at=xtick, labels = ticks.label)
  
  
  # difference in vertical
  diff.y <- target.newy.ctrl - baseline.newy.ctrl
  
  # return
  res <- list(control.pts=control.pts, baseline.time.new=baseline.time.new, target.time.new=target.time.new, baseline.newy.ctrl=baseline.newy.ctrl, target.newy.ctrl=target.newy.ctrl, diff.y=diff.y)
  res
}
