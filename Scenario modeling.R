# ================================
# KRILL BIOMASS MODEL (4 REGIONS)
# With fixed CCAMLR catch limits
# ================================

library(deSolve)

# --------------------------------
# 1. PARAMETERS
# --------------------------------

params <- list(
  
  # --- Reproduction rates (per year)
  r1 = 0.6,
  r2 = 0.5,
  r3 = 0.55,
  r4 = 0.65,
  
  # --- Carrying capacities (tons)
  Kmax1 = 13000000,
  Kmax2 = 20000000,
  Kmax3 = 220000000,
  Kmax4 = 8000000,
  
  # --- Natural mortality (per year)
  m1 = 0.2,
  m2 = 0.25,
  m3 = 0.22,
  m4 = 0.18,
  
  # --- Predation mortality (per year)
  p1 = 0.1,
  p2 = 0.15,
  p3 = 0.12,
  p4 = 0.09,
  
  # --- CCAMLR Catch Limits (tons per year)
  TAC1 = 155000,
  TAC2 = 279000,
  TAC3 = 279000,
  TAC4 = 93000,
  
  # --- Migration rates
  m21 = 0.04,
  m23 = 0.075,
  m24 = 0.05,
  m32 = 0.075,
  m34 = 0.05,
  m42 = 0.05,
  m43 = 0.05
)

# --------------------------------
# 2. INITIAL BIOMASS (tons)
# --------------------------------

state <- c(
  K1 = 15750000,
  K2 = 15750000,
  K3 = 15750000,
  K4 = 15750000
)

# --------------------------------
# 3. TIME (years)
# --------------------------------

times <- seq(0, 50, by = 0.1)

# --------------------------------
# 4. MODEL FUNCTION
# --------------------------------

krill_model <- function(t, state, parameters) {
  
  # Extract biomass
  K1 <- state[1]
  K2 <- state[2]
  K3 <- state[3]
  K4 <- state[4]
  
  with(as.list(parameters), {
    
    # --------------------------------
    # Fishing harvest (fixed TAC)
    # Prevent harvesting more than biomass
    # --------------------------------
    
    H1 <- min(TAC1, K1)
    H2 <- min(TAC2, K2)
    H3 <- min(TAC3, K3)
    H4 <- min(TAC4, K4)
    
    # --------------------------------
    # Logistic growth - mortality - fishing
    # --------------------------------
    
    dK1 <- r1*K1*(1 - K1/Kmax1) - m1*K1 - p1*K1 - H1
    dK2 <- r2*K2*(1 - K2/Kmax2) - m2*K2 - p2*K2 - H2
    dK3 <- r3*K3*(1 - K3/Kmax3) - m3*K3 - p3*K3 - H3
    dK4 <- r4*K4*(1 - K4/Kmax4) - m4*K4 - p4*K4 - H4
    
    # --------------------------------
    # Migration between regions
    # --------------------------------
    
    # Region 2 → 1,3,4
    dK2 <- dK2 - m21*K2 - m23*K2 - m24*K2
    dK1 <- dK1 + m21*K2
    dK3 <- dK3 + m23*K2
    dK4 <- dK4 + m24*K2
    
    # Region 3 → 2,4
    dK3 <- dK3 - m32*K3 - m34*K3
    dK2 <- dK2 + m32*K3
    dK4 <- dK4 + m34*K3
    
    # Region 4 → 2,3
    dK4 <- dK4 - m42*K4 - m43*K4
    dK2 <- dK2 + m42*K4
    dK3 <- dK3 + m43*K4
    
    # Return results
    list(c(dK1, dK2, dK3, dK4))
  })
}

# --------------------------------
# 5. RUN SIMULATION
# --------------------------------

out <- ode(y = state, times = times, func = krill_model, parms = params)
out <- as.data.frame(out)

# --------------------------------
# 6. PLOT RESULTS
# --------------------------------

matplot(out$time,
        out[,2:5],
        type = "l",
        lwd = 2,
        lty = 1,
        xlab = "Time (years)",
        ylab = "Krill Biomass (tons)")

legend("right",
       legend = c("Region 1","Region 2","Region 3","Region 4"),
       col = 1:4,
       lty = 1,
       lwd = 2)