plot.spatial <- function(coords, epi) {
point_matrix <- do.call(rbind, coords)
plot(point_matrix, 
     main = "Transmission Tree: Spatial", 
     xlab = "x", 
     ylab = "y", 
     pch = 19, 
     col = "darkgreen")
for (i in 1:nrow(epi)) {
  if (!is.na(epi[i,2])) {
    p <- coords[epi[i,2]]
    c <- coords[epi[i,1]]
    arrows(p[[1]][1], p[[1]][2], c[[1]][1], c[[1]][2], col = "red", length = 0.1, angle = 30,)
  }
}
points(coords[epi[1,1]][[1]][1],coords[epi[1,1]][[1]][2],pch = 19,col='blue')
}

plot.phylo <- function(newick) {
  newick <- gsub(";", "", as.character(newick))
  newick <- paste('(',newick,');')
  tree <- read.tree(text=newick)
  
  plot(tree)
  edgelabels(tree$edge.length, bg="black", col="white", font=1)
}