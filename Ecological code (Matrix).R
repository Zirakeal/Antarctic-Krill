library(deSolve)

# ----- Parameters -----
params <- list(
  r = c(0.6, 0.5, 0.55, 0.65, 0.58),  # Reproduction rate
  Kmax = c(1000, 800, 900, 1100, 950), # Carrying capacity
  m = c(0.2, 0.25, 0.22, 0.18, 0.21),  # Natural mortality
  p = c(0.1, 0.15, 0.12, 0.09, 0.11) # Predation mortality 
)

# ----- Migration Matrix -----
# Rows = FROM region
# Columns = TO region

M <- matrix(0, nrow = 5, ncol = 5)

# Region 1 connects to 2
M[1,2] <- 0.03

# Region 2 connects to 1 and 3
M[2,1] <- 0.01
M[2,3] <- 0.02

# Region 3 connects to 2 and 4
M[3,2] <- 0.015
M[3,4] <- 0.025

# Region 4 connects to 3
M[4,3] <- 0.02

# Region 5 is isolated (no touching regions)

# ----- Initial Biomass -----
state <- c(K1 = 500, K2 = 400, K3 = 450, K4 = 600, K5 = 480)

# ----- Time -----
times <- seq(0, 50, by = 0.1)

# ----- Model Function -----
krill_model <- function(t, state, parameters) {
  
  K <- state
  r <- parameters$r
  Kmax <- parameters$Kmax
  m <- parameters$m
  p <- parameters$p
  
  dK <- numeric(5)
  
  # Ecological dynamics
  for (i in 1:5) {
    growth <- r[i] * K[i] * (1 - K[i]/Kmax[i])
    mortality <- m[i] * K[i]
    predation <- p[i] * K[i]
    
    dK[i] <- growth - mortality - predation
  }
  
  # Migration dynamics
  for (i in 1:5) {
    for (j in 1:5) {
      if (i != j) {
        dK[i] <- dK[i] + (M[j,i] * K[j]) - (M[i,j] * K[i])
      }
    }
  }
  
  list(dK)
}

# ----- Run Simulation -----
out <- ode(y = state,
           times = times,
           func = krill_model,
           parms = params)

out <- as.data.frame(out)

# ----- Plot -----
matplot(out$time,
        out[,2:6],
        type = "l",
        lwd = 2,
        lty = 1,
        xlab = "Time",
        ylab = "Krill Biomass")

legend("right",
       legend = paste("Region", 1:5),
       col = 1:5,
       lty = 1,
       lwd = 2)