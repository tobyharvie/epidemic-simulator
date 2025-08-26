# Simulate an epidemic through a network of 30
n_pop <- 30
examplecoords <- lapply(1:n_pop, function(i) c(runif(1), runif(1)))
spatialkernel <- function (d, kappa=1) {
  return(exp(-kappa * d))
}
#set.seed(101)
library(epinet)
source("simulator.R")
epi <- simulate_epidemic(
  C = examplecoords, 
  N = n_pop, 
  beta = 1,
  ki = 3, 
  thetai = 7, 
  ke = 1,
  thetae = 7,
  ps = 0.5,
  thetas = 5, 
  ks = 4,
  K = spatialkernel
)

# plot using epinet function
epi_plot <- epi[, -6]
plot.epidemic(epi_plot)

