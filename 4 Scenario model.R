# ==========================================
# ANTARCTIC KRILL – SPATIAL POLICY MODEL
# 4 SCENARIO COMPARISON
# ==========================================

library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)

# ------------------------------------------
# 1. PARAMETERS
# ------------------------------------------

base_params <- list(
  
  # Growth rates
  r1 = 0.7, r2 = 0.7, r3 = 0.7, r4 = 0.7,
  
  # Carrying capacities
  Kmax1 = 25000000,
  Kmax2 = 35000000,
  Kmax3 = 40000000,
  Kmax4 = 20000000,
  
  # Natural mortality
  m1 = 0.2, m2 = 0.2, m3 = 0.2, m4 = 0.2,
  
  # Predation
  p1 = 0.1, p2 = 0.1, p3 = 0.1, p4 = 0.1,
  
  # Migration
  m21 = 0.04, m23 = 0.075, m24 = 0.05,
  m32 = 0.075, m34 = 0.05,
  m42 = 0.05, m43 = 0.05,
  
  # Catchability
  q = 2e-8,
  
  # Scenario (will be changed in loop)
  scenario = 0
)

# ------------------------------------------
# 2. INITIAL STATE
# ------------------------------------------

state <- c(
  K1 = 15750000,
  K2 = 15750000,
  K3 = 15750000,
  K4 = 15750000,
  
  # Effort (constant)
  E1 = 5,
  E2 = 5,
  E3 = 5,
  E4 = 5
)

# ------------------------------------------
# 3. TIME
# ------------------------------------------

times <- seq(0, 25, by = 1)

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
    
    # Apply closures
    if (scenario == 1) H1 <- 0
    if (scenario == 2) H2 <- 0
    if (scenario == 3) H3 <- 0
    if (scenario == 4) H4 <- 0
    
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
    
    # From region 2
    dK2 <- dK2 - m21*K2 - m23*K2 - m24*K2
    dK1 <- dK1 + m21*K2
    dK3 <- dK3 + m23*K2
    dK4 <- dK4 + m24*K2
    
    # From region 3
    dK3 <- dK3 - m32*K3 - m34*K3
    dK2 <- dK2 + m32*K3
    dK4 <- dK4 + m34*K3
    
    # From region 4
    dK4 <- dK4 - m42*K4 - m43*K4
    dK2 <- dK2 + m42*K4
    dK3 <- dK3 + m43*K4
    
    # Effort constant
    list(c(dK1, dK2, dK3, dK4,
           dE1 = 0,
           dE2 = 0,
           dE3 = 0,
           dE4 = 0))
  })
}

# ------------------------------------------
# 5. RUN ALL SCENARIOS
# ------------------------------------------

results <- list()

for (s in 1:4) {
  
  params <- base_params
  params$scenario <- s
  
  out <- ode(y = state, times = times,
             func = krill_model, parms = params)
  
  out <- as.data.frame(out)
  
  # Total biomass
  out$TotalBiomass <- rowSums(out[,c("K1","K2","K3","K4")])
  
  out$Scenario <- paste("Scenario", s)
  
  results[[s]] <- out
}

all_results <- bind_rows(results)

# ------------------------------------------
# 6. TOTAL BIOMASS COMPARISON
# ------------------------------------------

ggplot(all_results, 
       aes(time, TotalBiomass, color = Scenario)) +
  geom_line(size = 1.2) +
  theme_minimal() +
  labs(title = "Total Krill Biomass Across Policy Scenarios",
       x = "Time (years)",
       y = "Total Biomass (tons)") +
  theme(legend.position = "bottom")

# ------------------------------------------
# 7. REGIONAL BIOMASS COMPOSITION (STACKED)
# ------------------------------------------

out_long <- all_results %>%
  pivot_longer(cols = c(K1,K2,K3,K4),
               names_to = "Region",
               values_to = "Biomass")

ggplot(out_long, 
       aes(x = time, y = Biomass, fill = Region)) +
  geom_area(alpha = 0.8) +
  facet_wrap(~Scenario) +
  theme_minimal() +
  labs(title = "Spatial Biomass Composition Under Each Scenario",
       x = "Time (years)",
       y = "Biomass (tons)") +
  theme(legend.position = "bottom")