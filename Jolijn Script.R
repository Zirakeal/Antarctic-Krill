library(deSolve)

#hoi

two_patch_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Growth
    growth1 <- r1 * B1 * (1 - B1/K1)
    growth2 <- r2 * B2 * (1 - B2/K2)
    growth3 <- r3 * B3 * (1 - B3/K3)
    growth4 <- r4 * B4 * (1 - B4/K4)
    
    # Harvest
    harvest1 <- q1 * E1 * B1
    harvest2 <- q2 * E2 * B2
    harvest3 <- q3 * E3 * B3
    harvest4 <- q4 * E4 * B4
    
    # Migration
    mig1 <- m21 * B2 + m31 * B3 + m41 * B4 - m12 * B1 - m13 * B1 - m14 * B1
    mig2 <- m12 * B1 + m32 * B3 + m42 * B4 - m21 * B2 - m23 * B2 - m24 * B2
    mig3 <- m13 * B1 + m23 * B2 + m43 * B4 - m31 * B3 - m32 * B3 - m34 * B3
    mig4 <- m14 * B1 + m24 * B2 + m34 * B3 - m41 * B4 - m42 * B4 - m43 * B4
    
    # Differential equations (-B1 to ensure B1 is above)
    dB1 <- max(growth1 - harvest1 + mig1, -B1)
    dB2 <- max(growth2 - harvest2 + mig2, -B2)
    dB3 <- max(growth3 - harvest3 + mig3, -B3)
    dB4 <- max(growth4 - harvest4 + mig4, -B4)
    
    
    list(c(dB1, dB2, dB3, dB4))
  })
}

parameters <- c(
  r1 = 0.2,
  r2 = 0.2,
  r3 = 0.2,
  r4 = 0.2,
  K1 = 12000000,
  K2 = 12000000,
  K3 = 11000000,
  K4 = 11000000,
  q1 = 0.001,
  q2 = 0.001,
  q3 = 0.001,
  q4 = 0.001,
  E1 = 0,
  E2 = 0,
  E3 = 0,
  E4 = 200,
  m12 = 0.04,
  m13 = 0,
  m14 = 0,
  m21 = 0.04,
  m23 = 0.075,
  m24 = 0.05,
  m31 = 0,
  m32 = 0.075,
  m34 = 0.05,
  m41 = 0,
  m42 = 0.05,
  m43 = 0.05
  
)

state <- c(B1 = 10000000, B2 = 10000000, B3 = 10000000, B4 = 30000000)

times <- seq(0, 50, by = 0.1)

out <- ode(y = state, times = times,
           func = two_patch_model,
           parms = parameters)

out <- as.data.frame(out)

# Calculate total biomass at each time step
out$total_biomass <- out$B1 + out$B2 + out$B3 + out$B4

plot(out$time, out$B1, type = "l", col = "blue",
     ylim = c(0,40000000), ylab = "Biomass", xlab = "Time")
lines(out$time, out$B2, col = "red")
lines(out$time, out$B3, col = "orange")
lines(out$time, out$B4, col = "green")
lines(out$time, out$total_biomass, col = "black", lwd = 2)
legend("topright", legend=c("Area 1","Area 2", "Area 3", "Area 4", "Total"),
       col=c("blue","red", "orange", "green", "black"), lty=1, lwd=c(1,1,1,1,2))

effort_seq <- seq(0, 300, by = 5)

yield_curves <- vector("list", 4)
total_curves <- vector("list", 4)
results <- vector("list", 4)

for(s in 1:4){
  
  total_yield <- numeric(length(effort_seq))
  B_total_eq  <- numeric(length(effort_seq))
  
  yield_matrix <- matrix(0, nrow=length(effort_seq), ncol=4)
  
  for(i in seq_along(effort_seq)){
    
    E <- effort_seq[i]
    
    # Reset efforts
    parameters[c("E1","E2","E3","E4")] <- E
    params <- parameters
    
    params["E1"] <- E
    params["E2"] <- E
    params["E3"] <- E
    params["E4"] <- E
    
    # Close one patch
    params[paste0("E", s)] <- 0
    
    out <- try(ode(y = state, times = times,
                   func = two_patch_model,
                   parms = params), silent = TRUE)
    if(inherits(out, "try-error")){
      print(paste("ODE failed at effort =", E, "in scenario", s))
      next
    }
    
    out <- as.data.frame(out)
    
    B_eq <- as.numeric(tail(out[,c("B1","B2","B3","B4")],1))
    
    # Yield per patch
    Y1 <- params["q1"] * params["E1"] * B_eq[1]
    Y2 <- params["q2"] * params["E2"] * B_eq[2]
    Y3 <- params["q3"] * params["E3"] * B_eq[3]
    Y4 <- params["q4"] * params["E4"] * B_eq[4]
    
    yield_matrix[i,] <- c(Y1,Y2,Y3,Y4)
    
    total_yield[i] <- sum(c(Y1,Y2,Y3,Y4))
    B_total_eq[i]  <- sum(B_eq)
  }
  
  # Find MSY
  idx <- which.max(total_yield)
  
  results[[s]] <- list(
    MSY = total_yield[idx],
    E_MSY = effort_seq[idx],
    B_MSY = B_total_eq[idx],
    yield_at_MSY = yield_matrix[idx,]
  )
  
  total_curves[[s]] <- total_yield
  yield_curves[[s]] <- yield_matrix
}

for(s in 1:4){
  cat("\nScenario: Patch", s, "closed\n")
  print(results[[s]])
}


par(mfrow=c(2,2), mar=c(4,4,2,1))  # 2x2 layout

colors <- c("blue","red","orange","green")

for(s in 1:4){
  matplot(effort_seq, yield_curves[[s]], type="l", lty=1, col=colors,
          ylim=c(0, max(yield_curves[[s]], total_curves[[s]])),
          xlab="Fishing Effort", ylab="Yield", main=paste("Scenario", s, "- Patch", s, "closed"))
  
  # Total yield
  lines(effort_seq, total_curves[[s]], col="black", lwd=2)
  
  # MSY vertical line
  abline(v = results[[s]]$E_MSY, col="darkgrey", lty=2)
  
  # Horizontal lines for applied effort
  for(p in 1:4){
    if(p == s){
      applied_effort <- rep(0, length(effort_seq))  # closed patch = 0
    } else {
      applied_effort <- effort_seq  # open patches = full effort
    }
    lines(effort_seq, applied_effort, col=colors[p], lty=2)
  }
  
  legend("topright", legend=c(paste("Area",1:4), "Total","Applied Effort"),
         col=c(colors,"black","grey"), lty=c(1,1,1,1,2,2), lwd=c(1,1,1,1,2,1),
         cex=0.7, pt.cex=0.7, bty="n")
}

