# SIR model Gillespie
# see http://homepages.warwick.ac.uk/~masfz/ModelingInfectiousDiseases/Chapter6/Program_6.4/Program_6_4.py 


stoch_eqns <- function(dens) {
    
    # we have 6 different rates to calculate
    rates <- rep(0, length.out=6)
    
    # 6 x 3 matrix reflecting the change in densities of
    # S, I or R in response to the 
    density_change_matrix <- matrix(
        data=0 # fill it with zeros
        ,nrow=6
        ,ncol=3)
    
    # total density at the start    
    Ntotal <- sum(dens)

    # then fill the rate matrix, we follow the order
    # as in the ODE system as provided on
    # http://homepages.warwick.ac.uk/~masfz/ModelingInfectiousDiseases/Chapter6/Program_6.4/index.html
    rates[1] <- mu * Ntotal  # mu N   
    density_change_matrix[1,] <- c(1,0,0)
    
    
    rates[2] <- beta * dens[1] * dens[2] / N  # beta * X * Y  N, infection
    density_change_matrix[2,] <- c(-1,1,0)
    
    rates[3] <- mu * dens[1] # mu X, death of susceptible
    density_change_matrix[3,] <- c(-1,0,0)
    
    rates[4] <- gamma * dens[2] # gamma Y, loss of infection
    density_change_matrix[4,] <- c(0,-1,1)
    
    rates[5] <- mu * dens[2] # mu Y, death of infection
    density_change_matrix[5,] <- c(0,-1,0)
    
    rates[6] <- mu * dens[3] # mu Z, 
    density_change_matrix[6,] <- c(0,0,-1)
    
    sum_rates <- sum(rates)

    # the r1 and r2 random variables
    r1 <- runif(n=1)
    r2 <- runif(n=1)
    
    # calculate time that nothing happens
    ts <- -log(r2) / (sum_rates)

    # now sample which rate 
    which_event <- sample(x = 1:6,  # a number between 1 and 6
                          prob = rates # weights given by rates
                          )
    
    dens[3] = dens[3] + 
        density_change_matrix[which_event,]
    
} # end stoch eqns

SIR_Gillespie <- function(beta=1.0, 
                  gamma=1.0/10.0, 
                  mu=5e-4, 
                  N0=5000.0,
                  max_time=2*365)
{
    # initial number of infecteds
    
} # end SIR_Gillespie