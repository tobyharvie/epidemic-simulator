# Simulate an epidemic through a network of 30
library(epinet)
library(ape)
library(phyclust)

n_pop <- 30
examplecoords <- lapply(1:n_pop, function(i) c(runif(1), runif(1)))
spatialkernel <- function (d, kappa=1) {
  return(exp(-kappa * d))
}
#set.seed(102)
source("simulator.R")
source('epi2binarynewick.R')
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

# plot using epinet function
epi_plot <- epi[, -6]

print(epi)
epi <- round(epi, 2)
print(epi)

binarynewick <- gsub('\\[&type = "([A-Z])"\\]','',paste(epi2binarynewick(epi),';',sep=''))
tree <- read.tree(text=binarynewick)
#plot.epidemic(epi_plot[,-6])
#plot(tree)
#edgelabels(tree$edge.length, bg="black", col="white", font=1)

sequences <- seqgen(opts = "-mHKY -l40 -on -s0.2", newick.tree = binarynewick)
#sequences[1]

# save transmission tree in newick format
#epi2newick(epi)
#newick_tree <- paste(gsub('\\[&type = "([A-Z])"\\]','\\1',epi2newick(epi)),';')
#print(newick_tree)
#output_name <- "C:/Users/tobyh/OneDrive - The University of Auckland/_ CS 380/exampletree.trees"
#writeLines(newick_tree, output_name)
#print(epi[,6])
#print(length(epi[,6][epi[,6]!=Inf])/length(epi[,6]))


#library(phangorn)
#base_frequencies <- c(0.1, 0.2, 0.3, 0.4)
#tstv_ratio <- 2
#Q_matrix <- c(1, tstv_ratio, 1, 1, tstv_ratio, 1)
#library(treeio)
#tree <- treeio::read.tree(text=newick_tree)
##tree <- treeio::read.beast(file = "C:/Users/tobyh/OneDrive - The University of Auckland/_ CS 380/exampletree.tree")
#hky_ancestral_alignment <- simSeq(
#  x = tree,
#  l = 10,
#  bf = base_frequencies,
#  Q = Q_matrix,
#  type = "DNA",      # Explicitly specify DNA sequences
#  ancestral = TRUE
#)
#as.character(hky_ancestral_alignment)
#as_tibble(tree)
