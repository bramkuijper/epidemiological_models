# SIR model Gillespie
# see http://homepages.warwick.ac.uk/~masfz/ModelingInfectiousDiseases/Chapter6/Program_6.4/Program_6_4.py 


stoch_eqns <- function(input_densities){
    
    # we have 6 different rates to calculate
    rates = rep(0, length.out=6)
    
    # 6 x 3 matrix reflecting the change in densities of
    # S, I or R in response to the 
    
    
    
    # total density at the start    
    Ntotal = sum(input_densities)
    
}



SIR_Gillespie <- function(beta=1.0, 
                  gamma=1.0/10.0, 
                  mu=5e-4, 
                  N0=5000.0,
                  max_time=2*365)
{
    # initial number of infecteds
    
} # end SIR_Gillespie