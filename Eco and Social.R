library(deSolve)

# ----------------------------
# ----- Parameters per region
# ----------------------------
params <- list(
  # Reproduction rates
  r1 = 0.6, r2 = 0.5, r3 = 0.55, r4 = 0.65, 
  # Carrying capacities
  Kmax1 = 1000, Kmax2 = 800, Kmax3 = 900, Kmax4 = 1100,
  # Natural mortality
  m1 = 0.2, m2 = 0.25, m3 = 0.22, m4 = 0.18, 
  # Predation mortality
  p1 = 0.1, p2 = 0.15, p3 = 0.12, p4 = 0.09,
  # Maximum fishing effort
  Fmax1 = 0.1, Fmax2 = 0.15, Fmax3 = 0.12, Fmax4 = 0,
  # Migration rates (FROM -> TO)
  mig12 = 0.03,
  mig21 = 0.01, mig23 = 0.02,
  mig32 = 0.015, mig34 = 0.025,
  mig43 = 0.02
)

# ----------------------------
# ----- Initial State
# ----------------------------
state <- c(
  K1 = 500, K2 = 400, K3 = 450, K4 = 600,  # Krill biomass
  C1 = 0, C2 = 0, C3 = 0, C4 = 0           # Cumulative catch
)

# ----------------------------
# ----- Time
# ----------------------------
times <- seq(0, 50, by = 0.1)

# ----------------------------
# ----- Model Function
# ----------------------------
krill_model <- function(t, state, parameters) {
  
  K <- state[1:4]
  C <- state[5:8]
  
  # ------------------------
  # ----- Dynamic fishing proportional to biomass
  F <- numeric(4)
  F[1] <- parameters$Fmax1 * (K[1]/parameters$Kmax1)
  F[2] <- parameters$Fmax2 * (K[2]/parameters$Kmax2)
  F[3] <- parameters$Fmax3 * (K[3]/parameters$Kmax3)
  F[4] <- parameters$Fmax4 * (K[4]/parameters$Kmax4)
  
  # ------------------------
  # ----- Growth, mortality, predation, fishing
  dK <- numeric(4)
  dC <- numeric(4)
  
  dK[1] <- parameters$r1*K[1]*(1 - K[1]/parameters$Kmax1) - parameters$m1*K[1] - parameters$p1*K[1] - F[1]*K[1]
  dK[2] <- parameters$r2*K[2]*(1 - K[2]/parameters$Kmax2) - parameters$m2*K[2] - parameters$p2*K[2] - F[2]*K[2]
  dK[3] <- parameters$r3*K[3]*(1 - K[3]/parameters$Kmax3) - parameters$m3*K[3] - parameters$p3*K[3] - F[3]*K[3]
  dK[4] <- parameters$r4*K[4]*(1 - K[4]/parameters$Kmax4) - parameters$m4*K[4] - parameters$p4*K[4] - F[4]*K[4]
  
  dC <- F * K  # cumulative catch
  
  # ------------------------
  # ----- Explicit Migration
  # Region 1 → 2
  dK[1] <- dK[1] - parameters$mig12 * K[1]
  dK[2] <- dK[2] + parameters$mig12 * K[1]
  
  # Region 2 → 1 and 3
  dK[2] <- dK[2] - parameters$mig21*K[2] - parameters$mig23*K[2]
  dK[1] <- dK[1] + parameters$mig21*K[2]
  dK[3] <- dK[3] + parameters$mig23*K[2]
  
  # Region 3 → 2 and 4
  dK[3] <- dK[3] - parameters$mig32*K[3] - parameters$mig34*K[3]
  dK[2] <- dK[2] + parameters$mig32*K[3]
  dK[4] <- dK[4] + parameters$mig34*K[3]
  
  # Region 4 → 3
  dK[4] <- dK[4] - parameters$mig43*K[4]
  dK[3] <- dK[3] + parameters$mig43*K[4]
  
  list(c(dK, dC))
}

# ----------------------------
# ----- Run Simulation
# ----------------------------
out <- ode(y = state, times = times, func = krill_model, parms = params)
out <- as.data.frame(out)

# ----------------------------
# ----- Plot Krill Biomass
# ----------------------------
matplot(out$time, out[,2:5], type="l", lwd=2, lty=1, xlab="Time", ylab="Krill Biomass")
legend("right", legend=paste("Region",1:4), col=1:4, lty=1, lwd=2)

# ----------------------------
# ----- Plot Cumulative Catch
# ----------------------------
matplot(out$time, out[,6:9], type="l", lwd=2, lty=1, xlab="Time", ylab="Cumulative Catch")
legend("topleft", legend=paste("Region",1:4), col=1:4, lty=1, lwd=2)