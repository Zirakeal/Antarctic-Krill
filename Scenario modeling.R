# ==========================================
# ANTARCTIC KRILL – SPATIAL BIOECONOMIC MODEL
# TRUE OPEN ACCESS DYNAMICS
# ==========================================

library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)

# ------------------------------------------
# 1. PARAMETERS
# ------------------------------------------

params <- list(
  
  # --- Biological parameters
  r1 = 0.5, r2 = 0.6, r3 = 0.6, r4 = 0.5,
  
  Kmax1 = 25000000,
  Kmax2 = 35000000,
  Kmax3 = 40000000,
  Kmax4 = 20000000,
  
  m1 = 0.2, m2 = 0.25, m3 = 0.22, m4 = 0.18,
  p1 = 0.1, p2 = 0.15, p3 = 0.12, p4 = 0.09,
  
  # --- Migration
  m21 = 0.04,
  m23 = 0.075,
  m24 = 0.05,
  m32 = 0.075,
  m34 = 0.05,
  m42 = 0.05,
  m43 = 0.05,
  
  # --- Economic parameters
  q = 2e-8,       # catchability
  price = 2000,   # $ per ton
  cost = 500,     # cost per unit effort
  alpha = 1e-7    # effort adjustment speed
)

# ------------------------------------------
# 2. INITIAL STATE
# ------------------------------------------

state <- c(
  K1 = 15750000,
  K2 = 15750000,
  K3 = 15750000,
  K4 = 15750000,
  E1 = 5,
  E2 = 5,
  E3 = 5,
  E4 = 5
)

# ------------------------------------------
# 3. TIME
# ------------------------------------------

times <- seq(0, 120, by = 1)

# ------------------------------------------
# 4. MODEL FUNCTION
# ------------------------------------------

krill_model <- function(t, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    # --------------------------
    # Harvest
    # --------------------------
    
    H1 <- q * E1 * K1
    H2 <- q * E2 * K2
    H3 <- q * E3 * K3
    H4 <- q * E4 * K4
    
    # --------------------------
    # Profit
    # --------------------------
    
    pi1 <- price * H1 - cost * E1
    pi2 <- price * H2 - cost * E2
    pi3 <- price * H3 - cost * E3
    pi4 <- price * H4 - cost * E4
    
    # --------------------------
    # Biomass dynamics
    # --------------------------
    
    dK1 <- r1*K1*(1 - K1/Kmax1) - m1*K1 - p1*K1 - H1
    dK2 <- r2*K2*(1 - K2/Kmax2) - m2*K2 - p2*K2 - H2
    dK3 <- r3*K3*(1 - K3/Kmax3) - m3*K3 - p3*K3 - H3
    dK4 <- r4*K4*(1 - K4/Kmax4) - m4*K4 - p4*K4 - H4
    
    # --------------------------
    # Migration
    # --------------------------
    
    # From 2
    dK2 <- dK2 - m21*K2 - m23*K2 - m24*K2
    dK1 <- dK1 + m21*K2
    dK3 <- dK3 + m23*K2
    dK4 <- dK4 + m24*K2
    
    # From 3
    dK3 <- dK3 - m32*K3 - m34*K3
    dK2 <- dK2 + m32*K3
    dK4 <- dK4 + m34*K3
    
    # From 4
    dK4 <- dK4 - m42*K4 - m43*K4
    dK2 <- dK2 + m42*K4
    dK3 <- dK3 + m43*K4
    
    # --------------------------
    # Open Access Effort Dynamics
    # --------------------------
    
    dE1 <- alpha * pi1
    dE2 <- alpha * pi2
    dE3 <- alpha * pi3
    dE4 <- alpha * pi4
    
    # Prevent negative effort
    if (E1 <= 0 & dE1 < 0) dE1 <- 0
    if (E2 <= 0 & dE2 < 0) dE2 <- 0
    if (E3 <= 0 & dE3 < 0) dE3 <- 0
    if (E4 <= 0 & dE4 < 0) dE4 <- 0
    
    list(c(dK1, dK2, dK3, dK4,
           dE1, dE2, dE3, dE4))
  })
}

# ------------------------------------------
# 5. RUN MODEL
# ------------------------------------------

out <- ode(y = state, times = times,
           func = krill_model, parms = params)

out <- as.data.frame(out)

# ------------------------------------------
# 6. Biomass per region
# ------------------------------------------
out_long <- out %>%
  pivot_longer(cols = c(K1,K2,K3,K4),
               names_to = "Region",
               values_to = "Biomass")

ggplot(out_long, aes(time, Biomass)) +
  geom_line(size = 1) +
  facet_wrap(~Region, scales = "free_y") +
  theme_minimal() +
  labs(title="Krill Biomass Under Open Access",
       x="Time (years)",
       y="Biomass (tons)")

# ------------------------------------------
# Effort per region
#-----------------------------------------
effort_long <- out %>%
  pivot_longer(cols = c(E1,E2,E3,E4),
               names_to = "Region",
               values_to = "Effort")

ggplot(effort_long, aes(time, Effort)) +
  geom_line(size=1) +
  facet_wrap(~Region, scales="free_y") +
  theme_minimal() +
  labs(title="Fishing Effort (Open Access)",
       x="Time (years)",
       y="Effort")

# ------------------------------------------
# Total Biomass
# ------------------------------------------
out$TotalBiomass <- rowSums(out[,c("K1","K2","K3","K4")])

ggplot(out, aes(time, TotalBiomass)) +
  geom_line(size=1.2) +
  theme_minimal() +
  labs(title="Total Antarctic Krill Biomass",
       x="Time (years)",
       y="Total Biomass (tons)")
