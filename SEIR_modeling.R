# SEIR model taken from tutorial written by Peter Gleeson (thanks Peter!)

# Background information ------------------------------------------------->

# the model

# S = "Susceptible" – individuals who have not been exposed to the virus
# E = "Exposed" – individuals exposed to the virus, but not yet infectious
# I = "Infectious" – exposed individuals who go on to become infectious
# R = "Recovered" – infectious individuals who recover and become immune to the virus

# # the parameters

# β is the transmission coefficient. Think of this as the average number of infectious contacts an infectious individual in the population makes each time period. A high value of β means the virus has more opportunity to spread.
# σ is the rate at which exposed individuals become infectious. Think of it as the reciprocal of the average time it takes to become infectious. That is, if an individual becomes infectious after 4 days on average, σ will be 1/4 (or 0.25).
# γ is the rate at which infectious individuals recover. As before, think of it as the reciprocal of the average time it takes to recover. That is, if it takes 10 days on average to recover, γ will be 1/10 (or 0.1).
# μ is an optional parameter to describe the mortality rate of infectious individuals. The higher μ is, the more deadly the virus.

# # the equations

# Susceptible = always decreases over time -> probability of any given contact being between an infectious and susceptible individual is (I / N) * (S / N).
#               calculated by multiplying the transmission coefficient β, by the population size N.
#               
# Exposed = The flow into E will be matched by the flow out of S.
#           Formula = dE/dt = βSI/N - σE
#           
# Infectious = There are two ways to leave the "infectious" stage. Recovery or death.
#           Formula = dI/dt = σE - γI - μI
#           σE - γ (recovery rate) - μ (mortality rate)
# 
# Recovered = This is the simplest formula, it's just the recovery rate. For mortality, its just the rate of μI 
#           Recovery formula = dR/dt = γI
#           Mortality formula = dR/dt = μI
#
# Our friend R-naught! R₀ = β/γ so the basic reproductive number = transmission coefficient β / Recovery rate γ
          
# Setup ----------------------------------------------------- >

#require(deSolve)
library(deSolve)
library(ggplot2)
options(scipen=999) # removes scientific notation

# Epidemic with NO control measures ------------------------- >

# Function specifying the SEIR model

SEIR <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    N <- S+E+I+R
    dS <- -(beta*S*I)/N
    dE <- (beta*S*I)/N - sigma*E
    dI <- sigma*E - gamma*I - mu*I
    dR <- gamma*I
    dM <- mu*I
    
    return(list(c(dS, dE, dI, dR, dM)))
  })
}

# Specify parameters
params <- c(beta=0.5, sigma=0.25, gamma=0.2, mu=0.001)
# Specify the starting points of the model
initial_state <- c(S=999999, E=1, I=0, R=0, M=0)
# Vector going from 0 to 365 days
times <- 0:365
# Model with the "Ode" function
# Ode solves the equation with respect to time
# Solves it simultaneously
model <- ode(initial_state, times, SEIR, params)
# Results!
summary(model)
# Interpretation : 
# Out of a million individuals, 108,264 did not become infected. (Min of S)
# At the peak of the epidemic, 126,516 individuals were infectious simultaneously. (Max of I)
# 887,300 individuals recovered by the end of the model. (Max of R)
# A total of 4436 individuals died during the epidemic.(Max of M)

# Visualise the output
matplot(model, type="l", lty=1, main="SEIR model of an epidemic with no control measures", xlab="Time in days")
legend <- colnames(model)[2:6]
legend("right", legend=legend, col=2:6, lty = 1)

# Find time point with maximum number of infections
infections <- as.data.frame(model)$I
peak <- max(infections)
match(peak, infections) # day 112

# Epidemic WITH control measures ------------------------- >

SEIR_lockdown <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    
    beta = ifelse(
      (time <= start_lockdown || time >= end_lockdown),
      0.5, 0.1
    )
    
    N <- S+E+I+R
    dS <- -(beta*S*I)/N
    dE <- (beta*S*I)/N - sigma*E
    dI <- sigma*E - gamma*I - mu*I
    dR <- gamma*I
    dM <- mu*I
    
    return(list(c(dS, dE, dI, dR, dM)))
  })
}

#lockdown begins on day 90, ends on day 150

params <- c(
  sigma=0.25,
  gamma=0.2,
  mu=0.001,
  start_lockdown=90,
  end_lockdown=150
)

initial_state <- c(S=999999, E=1, I=0, R=0, M=0)
times <- 0:365
model <- ode(initial_state, times, SEIR_lockdown, params)
summary(model)

# Out of one million individuals, 156,886 did not become infected.
# At the peak of the epidemic, 72,444 individuals were infectious at the same time.
# 838,917 individuals recovered.
# A total of 4195 individuals died.

# If we plot this, we see two waves, consistent with what we see over
# the coronavirus pandemic
matplot(
  model, 
  type="l",
  lty=1, 
  main="SEIR model (with intervention)", 
  xlab="Time"
)
legend <- colnames(model)[2:6]
legend("right", legend=legend, col=2:6, lty = 1)

# We see that the max infections peak at day 223 instead of 112 with no control measures
infections <- as.data.frame(model)$I
peak <- max(infections)
match(peak, infections)



