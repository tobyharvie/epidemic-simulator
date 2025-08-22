# Simulate an epidemic through a network of 30
n_pop <- 30
examplecoords <- lapply(1:N, function(i) c(runif(1), runif(1)))
spatialkernel <- function (d, kappa=1) {
  return(exp(-kappa * d))
}
set.seed(3)
library
source("simulator.R")
epi <- simulate_epidemic(
  C = examplecoords, 
  N = n_pop, 
  alpha = 1.5,
  beta = 10, 
  ki = 3, 
  thetai = 7, 
  ke = 3,
  thetae = 1,
  ps = 0.5,
  thetas = 2, 
  ks = 0.75,
  K = spatialkernel
)

# plot using epinet function
epi_plot <- epi[, -6]
plot.epidemic(epi_plot)

