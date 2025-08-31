library(epinet)
library(ape)
library(phyclust)
source("simulator.R")
source('epi2binarynewick.R')
#set.seed(102)

n_pop <- 30
examplecoords <- lapply(1:n_pop, function(i) c(runif(1), runif(1)))
spatialkernel <- function (d, kappa=1) {
  return(exp(-kappa * d))
}
epi <- simulate_epidemic(
  C = examplecoords, 
  N = n_pop, 
  beta = 1,
  ki = 3, 
  thetai = 7, 
  ke = 1,
  thetae = 7,
  ps = 0.5,
  K = spatialkernel
)

epi <- round(epi, 2)
binarynewick <- gsub('\\[&type = "([A-Z])"\\]','',paste(epi2binarynewick(epi),';',sep=''))

# plotting
#epi_plot <- epi[, -6]
#plot.epidemic(epi_plot[,-6])
#tree <- read.tree(text=binarynewick)
#plot(tree)
#edgelabels(tree$edge.length, bg="black", col="white", font=1)

sequences <- seqgen(opts = "-mHKY -l40 -on -s0.2", newick.tree = binarynewick)

