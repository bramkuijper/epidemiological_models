# SIR model Gillespie
# see http://homepages.warwick.ac.uk/~masfz/ModelingInfectiousDiseases/Chapter6/Program_6.4/Program_6_4.py 


# develops the various rates of the SIR models
# and the resulting changes in density of S, I or R
stoch_eqns <- function(input,
                       beta,
                       gamma,
                       mu
                       ) {
    
    dens <- input[2:4]
    ts <- input[1]
    
    # going over eq. (2.13 - 2.15), 
    # we have 6 different rates to calculate
    # pre-allocate vector to contain the different rates
    rates <- rep(0, length.out=6)
    
    # pre-allocate a
    # 6 x 3 matrix reflecting the change in densities of
    # S, I or R during one event
    density_change_matrix <- matrix(
        data=0 # fill it with zeros
        ,nrow=6
        ,ncol=3)
    
    # total density at the start    
    Ntotal <- sum(dens)

    # then fill the rate matrix, we follow the order
    # as in the ODE system as provided on
    # http://homepages.warwick.ac.uk/~masfz/ModelingInfectiousDiseases/Chapter6/Program_6.4/index.html
    
    # birth of susceptible
    rates[1] <- mu * Ntotal  
    density_change_matrix[1,] <- c(1,0,0)
    
    # infection of susceptible
    rates[2] <- beta * dens[1] * dens[2] / Ntotal  # beta * X * Y  N, infection
    density_change_matrix[2,] <- c(-1,1,0)
    
    # death of susceptible
    rates[3] <- mu * dens[1] # mu X
    density_change_matrix[3,] <- c(-1,0,0)
    
    # loss of infection
    rates[4] <- gamma * dens[2] # gamma Y
    density_change_matrix[4,] <- c(0,-1,1)
    
    # death of infected
    rates[5] <- mu * dens[2] # mu Y
    density_change_matrix[5,] <- c(0,-1,0)
    
    # death of resistant
    rates[6] <- mu * dens[3] # mu Z, 
    density_change_matrix[6,] <- c(0,0,-1)
    
    sum_rates <- sum(rates)

    # random variable to sample from inverse distribution
    # calculating the time that nothing happens
    rand1 <- runif(n=1)
    
    # calculate time that nothing happens
    ts <- ts  -log(rand1) / (sum_rates)
    
    # bounds checking
    stopifnot(rates >= 0)    

    # normally we would use rand2 to sample the event that happens
    # but in R we can simply do this by sampling from   
    # 1:6 using weights that are the rates above
    # these rates just have to be positive which they are 
    which_event <- sample(x = 1:6,  # a number between 1 and 6
                          size = 1,
                          prob = rates # weights given by rates
                          )
    
    # materialise the actual change in the densities
    dens = dens + 
        density_change_matrix[which_event,]
    
    # return both the densities and the time that nothing happened
    return(c(ts,dens))
    
} # end stoch eqns

SIR_Gillespie <- function(beta=1.0, 
                  gamma=1.0/10.0, 
                  mu=5e-4, 
                  initial_values=c(5000-1,1,0),
                  max_time=2*365,
                  skip=10
                  )
{
    # pre-allocate data for the numerical output
    variables <- matrix(data=NA, # initially fill it with NAs
                        nrow=max_time,  # max_time rows
                        ncol=4)
    
    # give names to the columns
    colnames(variables) <- c("T","S","I","R")
    
    # set initial time and initial values
    variables[1,] <- c(0, initial_values)
    
    # now iterate the system for max_time steps
    for (time_idx in 2:max_time)
    {
        # evaluate the stochastic equations
        result <- stoch_eqns(input=variables[time_idx - 1,],
                             beta=beta,
                             gamma=gamma,
                             mu=mu
                             )
        
        # update the time variable
        variables[time_idx,] <- result

        if (time_idx %% skip == 0)        
        {
            print(paste0("time step ",time_idx," out of ",max_time))
        }
    } # end for
    
    return(as.data.frame(variables))
} # end SIR_Gillespie


output <- SIR_Gillespie(max_time = 10*365)

output_l <- pivot_longer(data=output
                         ,cols = c(S,I,R)
                         ,names_to="type"
                         ,values_to="density")

ggplot(data = output_l
       ,mapping=aes(x=T, y=density)) +
    geom_line(mapping=aes(colour=type)) +
    scale_colour_brewer(palette="Set1") +
    theme_classic()
