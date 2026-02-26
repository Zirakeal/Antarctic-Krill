library(deSolve)

# ----- Parameters -----

params <- c(
  r1 = 0.6, r2 = 0.5, r3 = 0.55, r4 = 0.65, r5 = 0.58, # Reproduction rate
  Kmax1 = 1000, Kmax2 = 800, Kmax3 = 900, Kmax4 = 1100, Kmax5 = 950, # Carrying capacity
  m1 = 0.2, m2 = 0.25, m3 = 0.22, m4 = 0.18, m5 = 0.21, # Natural mortality
  p1 = 0.1, p2 = 0.15, p3 = 0.12, p4 = 0.09, p5 = 0.11 # Predation mortality
)

# ----- Initial biomass -----
state <- c(
  K1 = 500,
  K2 = 400,
  K3 = 450,
  K4 = 600,
  K5 = 480
)

# ----- Time -----
times <- seq(0, 120, by = 1) # Simulate for 120 months time units with a step of 1 month

# ----- Model -----
krill_model <- function(t, state, params) {
  
  with(as.list(c(state, params)), {
    
    dK1 <- r1*K1*(1 - K1/Kmax1) - m1*K1 - p1*K1
    dK2 <- r2*K2*(1 - K2/Kmax2) - m2*K2 - p2*K2
    dK3 <- r3*K3*(1 - K3/Kmax3) - m3*K3 - p3*K3
    dK4 <- r4*K4*(1 - K4/Kmax4) - m4*K4 - p4*K4
    dK5 <- r5*K5*(1 - K5/Kmax5) - m5*K5 - p5*K5
    
    list(c(dK1, dK2, dK3, dK4, dK5))
  })
}

# ----- Run model -----
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
        xlab = "months",
        ylab = "Krill Biomass")

legend("right",
       legend = c("Region 1","Region 2","Region 3","Region 4","Region 5"),
       col = 1:5,
       lty = 1,
       lwd = 2)