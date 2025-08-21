# Simulate an epidemic through a network of 30
library(epinet)
source("simulator.R")
set.seed(3)
N <- 30
# Build dyadic covariate matrix (X)
# Have a single covariate for overall edge density; this is the Erdos-Renyi model
examplecoords <- lapply(1:N, function(i) c(runif(1), runif(1)))


spatialkernel <- function (d, kappa=1) {
  return(exp(-kappa * d))
}
set.seed(3)
source("simulator.R")
exampleepidemic <- simulate_epidemic(examplecoords, N = 30, alpha = 1.5,
                                     beta = 0.3, ki = 0.75, thetai = 1, 
                                     latencyperiod = 1, ps=0.5, K=spatialkernel)
