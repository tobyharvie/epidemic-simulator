library(epinet)
library(ape)
library(phyclust)
source("simulator.R")
source('epi2binarynewick.R')
source('plotting.R')

# simulation
n_pop <- 99
coords <- lapply(1:n_pop, function(i) c(runif(1), runif(1)))
ki <- 3; thetai <- 7; ke <- 3; thetae <- 7; ps <- 0.9; kappa <- 17
spatialkernel <- function (d, kappa=17) {
  return(exp(-kappa * d))
}
epi <- simulate_epidemic(
  C = coords, 
  N = n_pop, 
  beta = 10,
  ki = ki, 
  thetai = thetai, 
  ke = ke,
  thetae = thetae,
  ps = ps,
  K = spatialkernel
)
plot.spatial(coords, epi)
newick <- newick_tree(epi)
nbases <- 100
sequences <- seqgen(opts = paste("-mHKY -l",nbases," -on -s0.1",sep=''), newick.tree = newick)

# plotting 
plot.epidemic(epi[,-6])
plot.spatial(coords, epi)
#plot.phylo(newick)

# epinet (for comparison)
#mycov <- data.frame(id = 1:n_pop, xpos = runif(n_pop), ypos = runif(n_pop))
#dyadCov <- BuildX(mycov, binaryCol = list(c(2, 3)),binaryFunc = "euclidean")
#eta <- c(0, -7)
#net <- SimulateDyadicLinearERGM(N = n_pop, dyadiccovmat = dyadCov, eta = eta)
#epi2 <- SEIR.simulator(M = net, N = n_pop, beta = 1, ki = 3, thetai = 7,ke = 3, latencydist = "gamma")
#plot.epidemic(epi2)

# data pipelines for analysis
source("./data_preparation.R")

# write sampled data as NEXUS format for BEAST and SCOTTI inference

