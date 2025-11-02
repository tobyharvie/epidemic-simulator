library(epinet)
library(ape)
library(phyclust)
library(MASS)
library(fs)
library(TAF)
source('plotting.R')
source('epi2binarynewick.R')
source('simulator.R')
source('data_preparation.R')

dir <- './last/'
if (file.exists(dir)){
  dir_delete(dir)
}
mkdir("last")

# epidemic params
n_pop <- 60
# k shape, theta scale. realistic variables from keeling 2001
ki <- 3; thetai <- 2.17; ke <- 3; thetae <- 1.17
tc <- 2 # transmission constant R0
ka <- 0 # spatial effect. Can reanalyse better later.
spatialkernel <- function (d, kappa=ka) {
  return(exp(-kappa * d))
}
beta <- exp(0.5*ka)*tc/(ki*thetai)

# sequencing params
nbases <- 500

# analyse with differing sampling probabilities

seed_used <- 100
for (ps in seq(from = 0.5, to = 1, by = 0.01)){
  seed_used <- seed_used+100
  # repeat twice for each sampling probability
  for (v in 1:2){
    seed_used <- seed_used+100
    print(paste0("Simulating epidemic ",v," with ps=",ps))
    repeat{
      set.seed(seed_used)
      coords <- lapply(1:n_pop, function(i) c(runif(1), runif(1)))
      epi <- simulate_epidemic(
        C = coords, 
        N = n_pop, 
        beta = beta,
        ki = ki, 
        thetai = thetai, 
        ke = ke,
        thetae = thetae,
        ps = ps,
        K = spatialkernel
      )
      # ensures epidemic spreads to at least 80% of the population
      if (!is.na(epi[round(0.8*n_pop),2])) { 
        #plot.spatial(coords, epi)
        #plot.epidemic(epi[,-6])
        #save
        newick <- newick_tree(epi)
        sequences <- seqgen(opts = paste("-mHKY -l",nbases," -on -s0.0005 -q",sep=''), newick.tree = newick)
        
        dir2 <- paste0(dir,"/ps-",ps*100,'-',v,'/')
        dir.create(dir2)
        save(epi,newick,seed_used,file=paste0(dir2, '/epi-data-',ps,'-',v))
        
        #print(epi)
        
        prepare_data(epi,sequences,ps,dir2,v)
        
        break
      }
      seed_used <- seed_used + 10
    }
  }
}
