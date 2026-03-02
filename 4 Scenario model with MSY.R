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
  r1 = 0.6, r2 = 0.6, r3 = 0.6, r4 = 0.6,
  
  # Carrying capacities
  Kmax1 = 25000000,
  Kmax2 = 35000000,
  Kmax3 = 40000000,
  Kmax4 = 20000000,
  
  # Natural mortality
  m1 = 0.1, m2 = 0.1, m3 = 0.1, m4 = 0.1,
  
  # Predation
  p1 = 0.1, p2 = 0.1, p3 = 0.1, p4 = 0.1,
  
  # Migration
  m21 = 0.04, m23 = 0.075, m24 = 0.05,
  m32 = 0.075, m34 = 0.05,
  m42 = 0.05, m43 = 0.05,
  
  # Catchability
  q = 2e-5,
  
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
  E1 = 10000,
  E2 = 10000,
  E3 = 10000,
  E4 = 10000
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
  
  # --------------------------
  # Compute Harvest
  # --------------------------
  
  out$H1 <- params$q * out$E1 * out$K1
  out$H2 <- params$q * out$E2 * out$K2
  out$H3 <- params$q * out$E3 * out$K3
  out$H4 <- params$q * out$E4 * out$K4
  
  # Apply scenario closures
  if (s == 1) out$H1 <- 0
  if (s == 2) out$H2 <- 0
  if (s == 3) out$H3 <- 0
  if (s == 4) out$H4 <- 0
  
  # --------------------------
  # Total Catch
  # --------------------------
  
  out$TotalCatch <- out$H1 + out$H2 + out$H3 + out$H4
  
  out$CumulativeCatch <- cumsum(out$TotalCatch)
  
  # --------------------------
  # Total Biomass
  # --------------------------
  
  out$TotalBiomass <- rowSums(out[,c("K1","K2","K3","K4")])
  
  out$Scenario <- paste("Scenario", s)
  
  results[[s]] <- out
}

all_results <- bind_rows(results)

final_summary <- all_results %>%
  group_by(Scenario) %>%
  summarise(
    FinalBiomass = last(TotalBiomass),
    TotalCatch = last(CumulativeCatch)
  )

initial_total_biomass <- sum(state[c("K1","K2","K3","K4")])

sustainable <- final_summary %>%
  filter(FinalBiomass >= initial_total_biomass)

best_policy <- sustainable %>%
  arrange(desc(TotalCatch)) %>%
  slice(1)

best_policy

# ==========================================
# 6. MSY ANALYSIS (NUMERICAL)
# ==========================================

effort_values <- seq(0, 20000, by = 100)

msy_results <- data.frame()

for (e in effort_values) {
  
  # set equal effort in all regions
  state_msy <- state
  state_msy["E1"] <- e
  state_msy["E2"] <- e
  state_msy["E3"] <- e
  state_msy["E4"] <- e
  
  out_msy <- ode(y = state_msy, times = times,
                 func = krill_model,
                 parms = base_params)
  
  out_msy <- as.data.frame(out_msy)
  
  # compute harvest
  out_msy$H1 <- base_params$q * e * out_msy$K1
  out_msy$H2 <- base_params$q * e * out_msy$K2
  out_msy$H3 <- base_params$q * e * out_msy$K3
  out_msy$H4 <- base_params$q * e * out_msy$K4
  
  out_msy$TotalCatch <- out_msy$H1 + out_msy$H2 +
    out_msy$H3 + out_msy$H4
  
  mean_catch <- mean(out_msy$TotalCatch[(length(times)/2):length(times)])
  
  final_biomass <- sum(out_msy[nrow(out_msy), c("K1","K2","K3","K4")])
  
  msy_results <- rbind(msy_results,
                       data.frame(Effort = e,
                                  MeanCatch = mean_catch,
                                  FinalBiomass = final_biomass))
}

# Plot MSY 
ggplot(msy_results, aes(Effort, MeanCatch)) +
  geom_line(size = 1.3) +
  theme_minimal() +
  labs(title = "Yield Curve (MSY Analysis)",
       x = "Effort",
       y = "Average Catch") +
  geom_vline(xintercept =
               msy_results$Effort[which.max(msy_results$MeanCatch)],
             linetype="dashed")

# Plot Biomass at different effort levels
ggplot(msy_results, aes(Effort, FinalBiomass)) +
  geom_line(size = 1.3) +
  theme_minimal() +
  labs(title = "Biomass at Different Effort Levels",
       x = "Effort",
       y = "Final Biomass")

# Identify MSY point
msy_point <- msy_results[which.max(msy_results$MeanCatch), ]

msy_point

# Compare scenario catches to MSY
MSY_value <- max(msy_results$MeanCatch)

ggplot(all_results, aes(time, TotalCatch, color = Scenario)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = MSY_value,
             linetype="dashed") +
  theme_minimal() +
  labs(title = "Scenario Catch vs MSY",
       y = "Total Catch")
# ------------------------------------------
# 7. TOTAL BIOMASS COMPARISON
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
# 8. REGIONAL BIOMASS COMPOSITION (STACKED)
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