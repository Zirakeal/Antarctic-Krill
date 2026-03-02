# ==========================================
# ANTARCTIC KRILL – SPATIAL BIOECONOMIC MODEL
# TRUE OPEN ACCESS DYNAMICS
# ==========================================

library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)

# ------------------------------------------
# 1. PARAMETERS
# ------------------------------------------

params <- list(
  r1 = 0.7, r2 = 0.7, r3 = 0.7, r4 = 0.7, # growth rate
  Kmax1 = 25000000, # Carry capacity
  Kmax2 = 35000000,
  Kmax3 = 40000000,
  Kmax4 = 20000000,
  m1 = 0.2, m2 = 0.2, m3 = 0.2, m4 = 0.2, # Natural mortality
  p1 = 0.1, p2 = 0.1, p3 = 0.1, p4 = 0.1, # Predation mortality
  
  # Migration
  m21 = 0.04, m23 = 0.075, m24 = 0.05, 
  m32 = 0.075, m34 = 0.05,
  m42 = 0.05, m43 = 0.05,
  
  # catchability coefficient
  q = 2e-8
)

# ------------------------------------------
# 2. INITIAL STATE
# ------------------------------------------

state <- c(
  K1 = 15750000, #initial biomass
  K2 = 15750000,
  K3 = 15750000,
  K4 = 15750000
)

# ------------------------------------------
# 3. TIME
# ------------------------------------------

times <- seq(0, 25, by = 1)

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
    # Return derivatives
    # --------------------------
    
    list(c(dK1, dK2, dK3, dK4,
           dE1 = 0,
           dE2 = 0,
           dE3 = 0,
           dE4 = 0))
  })
}

# ------------------------------------------
# 5. RUN MODEL
# ------------------------------------------

out <- ode(y = state, times = times,
           func = krill_model, parms = params)

out <- as.data.frame(out)


# ------------------------------------------
# Stacked biomass plot
# ------------------------------------------

out_long <- out %>%
  pivot_longer(cols = c(K1,K2,K3,K4),
               names_to = "Region",
               values_to = "Biomass")

ggplot(out_long, aes(x = time, y = Biomass, fill = Region)) +
  geom_area(alpha = 0.8) +
  theme_minimal() +
  labs(title = "Spatial Composition of Krill Biomass",
       x = "Time (years)",
       y = "Biomass (tons)") +
  theme(legend.position = "bottom")

# ------------------------------------------
# Biomass Share Over Time (Relative Dynamics)
# ------------------------------------------
out_share <- out_long %>%
  group_by(time) %>%
  mutate(Share = Biomass / sum(Biomass))

ggplot(out_share, aes(time, Share, color = Region)) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(title = "Regional Biomass Share Over Time",
       x = "Time (years)",
       y = "Proportion of Total Biomass") +
  theme(legend.position = "bottom")

# ------------------------------------------
# Total Biomass + Rolling trend
# ------------------------------------------
out$TotalBiomass <- rowSums(out[,c("K1","K2","K3","K4")])

ggplot(out, aes(time, TotalBiomass)) +
  geom_line(size = 1.2) +
  geom_smooth(se = FALSE, linewidth = 1) +
  theme_minimal() +
  labs(title = "Total Krill Biomass with Trend",
       x = "Time (years)",
       y = "Total Biomass (tons)")

# ------------------------------------------
# Phase plot (System stability View for region 1)
# ------------------------------------------
ggplot(out, aes(x = K1, y = TotalBiomass)) +
  geom_path() +
  theme_minimal() +
  labs(title = "Phase Relationship: Region 1 vs Total Biomass",
       x = "K1",
       y = "Total Biomass")

# ------------------------------------------
# Migration Heat map
# ------------------------------------------
out_rel <- out
out_rel[,c("K1","K2","K3","K4")] <-
  sweep(out[,c("K1","K2","K3","K4")], 2,
        out[1,c("K1","K2","K3","K4")], "/")

out_melt <- melt(out_rel,
                 id.vars = "time",
                 measure.vars = c("K1","K2","K3","K4"))

ggplot(out_melt, aes(time, variable, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low="blue", high="red") +
  theme_minimal() +
  labs(title="Relative Biomass Change (Heatmap)",
       x="Time (years)",
       y="Region",
       fill="Relative Level")