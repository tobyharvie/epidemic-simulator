table.to.phylo <- function(table) {
  ## Create the edge table
  edge <- table[, c("parent", "node")]
  ## Get the edge lengths
  if(any("length" %in% colnames(table))) {
    edge.length <- table[, "length"]
  } else {
    edge.length <- NULL
  }
  
  ## Get the number of nodes (half the number of edges for a bifurcating tree)
  Nnode <- (nrow(edge)/2)
  
  ## Convert the edge table into integers
  ## Detecting each elements
  tips  <- edge[(!edge[,2] %in% edge[,1]), 2] # tips are in the node column only
  nodes <- edge[(edge[,2] %in% edge[,1]), 2] # nodes are in both columns
  root  <- unique(edge[!(edge[,1] %in% edge[, 2]), 1]) # root is only in the parent column
  ## Get the elements integer value
  root_int <- length(tips)+1
  nodes_int <- paste0("node_temp_", 1:Nnode + root_int)
  tips_int <- 1:length(tips)
  
  ## Update the root as integer
  edge[,1] <- gsub(root, root_int, edge[,1])
  ## Update the nodes as integers
  recursive.gsub <- function(x, y, vector) {
    while(length(x) > 0) {
      vector <- gsub(x[1], y[1], vector)
      x <- x[-1]
      y <- y[-1]
    }
    return(vector)
  }
  edge <- apply(edge, 2, function(x, nodes, nodes_int) recursive.gsub(nodes, nodes_int, x), nodes, nodes_int)
  ## Update the tips as integers
  edge[,2] <- recursive.gsub(tips, tips_int, edge[,2])
  ## Recursively reorder the nodes (root must connect to root+1 and root+2, etc...)
  first_node <- max_node <- root_int
  while(length(grep("node_temp_", edge)) > 0) {
    ## Name the new nodes
    new_node_id <- 1:2 + max_node
    ## Find the descendants of the first node
    to_replace <- edge[grep(first_node, edge[,1]),2]
    ## Check if needs replacement
    needs_replacement <- grep("node_temp_", to_replace)
    if(length(needs_replacement) > 0) {
      to_replace <- to_replace[needs_replacement]
      new_node_id <- new_node_id[1:length(to_replace)]
      ## Replace in the table
      edge <- apply(edge, 2, function(x, to_replace, new_node_id) recursive.gsub(to_replace, new_node_id, x), to_replace, new_node_id)
    }
    ## Update next node to check
    first_node <- new_node_id[1]
    max_node <- max(new_node_id)
  }
  
  ## Get the tip labels
  tip.label <- tips
  node.label <- c(root,nodes)
  
  ## Format it as an integer matrix
  edge <- apply(edge, 2, as.integer)
  colnames(edge) <- NULL
  ## Sort them by root for convenience
  edge <- edge[order(edge[, 1]),]
  
  ## Make the list
  tree <- list(edge = edge, tip.label = tip.label, Nnode = Nnode, edge.length = edge.length, node.label = node.label)
  class(tree) <- "phylo"
  return(tree)
}


# let e be the size 
# node number e indicates time of exposure,
# node number e + len(epi) indicates time of infectiousness
epi = removesusceptibles(epi)
n <- nrow(epi)
# initialize start of infection
epi[1,6] <- 0
edgeList <- list()
edgeLengthList <- list()
tips <- numeric(3*n)
tips[epi[[1,1]]] <- paste(epi[[1,1]],'e',sep='')
for (i in 2:nrow(epi)) {
  # get the row. infected event has already been created with node label 'xi'
  parent <- epi[i, ] 
  # get people infected by this node
  children <- matrix(epi[!is.na(epi[,2]) & epi[,2] == parent[1], ],ncol=6)
  # used for node name when creatingn new nodes
  currentNode <- paste(parent[[1]],'i',sep='')
  # used to ensure correct placement of sampling
  samplingTime <- ifelse(!is.na(parent[[6]]), parent[[6]], Inf)
  
  if (nrow(children)>0) {
    for (j in 1:nrow(children)) {
      child <- children[j,]
      
      # insert an event corresponding to sampling time
      if (samplingTime < child[[3]] ) {
        edgeList[[length(edgeList)+1]] <- c(currentNode, paste(parent[[1]],'s',sep=''))
        edgeLengthList[[length(edgeLengthList)+1]] <- child[[3]]-parent[[4]]
        tips[parent[[1]]+2*n] <- paste(parent[[1]],'s',sep='')
        samplingTime <- Inf
        currentNode <- paste(parent[[1]],'s',sep='')
      }
      
      # create edge from most recent event (exposure of some other child, or sampling) to exposure event of child
      edgeList[[length(edgeList)+1]] <- c(currentNode, paste(child[[1]],'e',sep=''))
      edgeLengthList[[length(edgeLengthList)+1]] <- child[[3]]-parent[[4]]
      tips[child[[1]]] <- paste(child[[1]],'e',sep='')
      # create edge from exposure to infectiousness of child
      edgeList[[length(edgeList)+1]] <- c(paste(child[[1]],'e',sep=''), paste(child[[1]],'i',sep=''))
      edgeLengthList[[length(edgeLengthList)+1]] <- child[[4]]-child[[3]]
      tips[child[[1]]+n] <- paste(child[[1]],'i',sep='')
      
      # iterate
      currentNode <- paste(child[[1]],'e',sep='')
      
    }
    # make sure last node is binary
    edgeList[[length(edgeList)+1]] <- c(currentNode, paste(currentNode,'x',sep=''))
    edgeLengthList[[length(edgeLengthList)+1]] <- 0
    
  }
  else if (nrow(children)==0) {
    edgeList[[length(edgeList)+1]] <- c(paste(parent[[1]],'i',sep=''), paste(parent[[1]],'s',sep=''))
    edgeLengthList[[length(edgeLengthList)+1]] <- ifelse(samplingTime==Inf,parent[[4]],samplingTime)-parent[[4]]
    tips[parent[[1]]+2*n] <- paste(parent[[1]],'s',sep='')
  }
  # else tree is binary where infected time is a leaf
  
  #childtime <- ifelse(!is.na(row[[6]]), row[[6]], row[[3]])
  #parenttime <- ifelse(!is.na(epi[epi[, 1] == row[[2]], 6]), epi[epi[, 1] == row[[2]], 6], epi[epi[, 1] == row[[2]], 3])
  #ededgeLengthList[i] <- (childtime-parenttime)[1]
}

edgeMatrix <- do.call(rbind, edgeList)
edgeLengthMatrix <- do.call(rbind, edgeLengthList)

table <- as.data.frame(cbind(edgeMatrix,edgeLengthMatrix))
colnames(table) <- c("parent", "node", "length")
phy <-table.to.phylo(table)
write.tree(phy)