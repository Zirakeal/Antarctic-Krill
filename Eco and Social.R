library(deSolve)

# ----------------------------
# ----- Parameters per region
# ----------------------------
params <- list(
  # Reproduction rates
  r1 = 0.6, r2 = 0.5, r3 = 0.55, r4 = 0.65, 
  # Carrying capacities
  Kmax1 = 13000000, Kmax2 = 20000000, Kmax3 = 220000000, Kmax4 = 8000000,
  # Natural mortality
  m1 = 0.2, m2 = 0.25, m3 = 0.22, m4 = 0.18, 
  # Predation mortality
  p1 = 0.1, p2 = 0.15, p3 = 0.12, p4 = 0.09,
  # Maximum fishing effort
  Fmax1 = 0.1, Fmax2 = 0.15, Fmax3 = 0.12, Fmax4 = 0,
  # Migration rates (FROM -> TO)
  m12 = 0.04,
  m21 = 0.04,
  m23 = 0.075,
  m24 = 0.05,
  m32 = 0.075,
  m34 = 0.05,
  m42 = 0.05,
  m43 = 0.05,
  
  # --- Economic parameters
  q1 = 0.00000002, # catchability
  q2 = 0.00000002,
  q3 = 0.00000002,
  q4 = 0.00000002,
  
  price = 2000,      # $ per ton
  cost = 500,        # cost per unit effort
  alpha = 0.0000001  # effort adjustment speed
)
# ----------------------------
# ----- Initial State
# ----------------------------
state <- c(
  K1 = 15750000,
  K2 = 15750000,
  K3 = 15750000,
  K4 = 15750000,
  E1 = 10,
  E2 = 10,
  E3 = 10,
  E4 = 10
)

# ----------------------------
# ----- Time
# ----------------------------
times <- seq(0, 50, by = 0.1)

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

