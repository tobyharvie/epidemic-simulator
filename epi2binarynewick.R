epi2binarynewick.int <- function (epi, node = epi[1, 1], lastrow = 1, infectedatrow = 1, 
          rightchild = FALSE) 
{
  nex = ""
  if (lastrow == 1) {
    epi = removesusceptibles(epi)
    if (nrow(epi) < 2) 
      stop("Need at least two infecteds to convert transmission tree to Newick format")
    return(paste("((", epi2binarynewick.int(epi, node = epi[2, 
                                                      "Parent"], lastrow = 2, infectedatrow = 1), ",", 
                 epi2binarynewick.int(epi, node = epi[2, "Node ID"], lastrow = 2, 
                                infectedatrow = 2, rightchild = TRUE), ",:0.0)[&type = \"I\"]:", 
                 epi[2, "Etime"] - epi[1, "Itime"], ",:0.0)[&type = \"E\"]:", 
                 epi[1, "Itime"] - epi[1, "Etime"], sep = ""))
  }
  atleaf = FALSE
  if (lastrow == nrow(epi)) {
    atleaf = TRUE
  }
  else {
    nextrow = match(node, epi[(lastrow + 1):nrow(epi), "Parent"])
    atleaf = is.na(nextrow)
  }
  if (atleaf) {
    if (rightchild) {
      return(paste("(", node, "[&type = \"I\"]:", epi[infectedatrow, 
                                                      "Rtime"] - epi[infectedatrow, "Itime"], ",:0.0)[&type = \"E\"]:", 
                   epi[infectedatrow, "Itime"] - epi[lastrow, "Etime"],
                   sep = ""))
    }
    else {
      return(paste(node, "[&type = \"I\"]:", epi[infectedatrow, 
                                                 "Rtime"] - epi[lastrow, "Etime"], sep = ""))
    }
  }
  else {
    nextrow = nextrow + lastrow
    if (rightchild) {
      return(paste("((", epi2binarynewick.int(epi, node = node, 
                                        lastrow = nextrow, infectedatrow = infectedatrow), 
                   ",", epi2binarynewick.int(epi, node = epi[nextrow, 
                                                       "Node ID"], lastrow = nextrow, infectedatrow = nextrow, 
                                       rightchild = TRUE), ")[&type = \"I\"]:", epi[lastrow, 
                                                                                    "Itime"] - epi[lastrow, "Etime"], ",:0.0)[&type = \"E\"]:", 
                   epi[nextrow, "Etime"] - epi[lastrow, "Itime"], 
                   sep = ""))
    }
    else {
      return(paste("(", epi2binarynewick.int(epi, node = node, 
                                       lastrow = nextrow, infectedatrow = infectedatrow), 
                   ",", epi2binarynewick.int(epi, node = epi[nextrow, 
                                                       "Node ID"], lastrow = nextrow, infectedatrow = nextrow, 
                                       rightchild = TRUE), ")[&type = \"I\"]:", epi[nextrow, 
                                                                                    "Etime"] - epi[lastrow, "Etime"], sep = ""))
    }
  }
}



epi2binarynewick <- function (epi, node = epi[1, 1], lastrow = 1, infectedatrow = 1, 
          rightchild = FALSE) 
{
  epi = preprocessepi(epi)[,-6]
  nex = ""
  if (lastrow == 1) {
    epi = removesusceptibles(epi)
    if (nrow(epi) < 2) 
      stop("Need at least two infecteds to convert transmission tree to Newick format")
    return(paste("((", epi2binarynewick.int(epi, node = epi[2, 
                                                      "Parent"], lastrow = 2, infectedatrow = 1), ",", 
                 epi2binarynewick.int(epi, node = epi[2, "Node ID"], lastrow = 2, 
                                infectedatrow = 2, rightchild = TRUE), ")[&type = \"I\"]:", 
                 epi[2, "Etime"] - epi[1, "Itime"], ",:0.0)[&type = \"E\"]:", 
                 epi[1, "Itime"] - epi[1, "Etime"], sep = ""))
  }
  atleaf = FALSE
  if (lastrow == nrow(epi)) {
    atleaf = TRUE
  }
  else {
    nextrow = match(node, epi[(lastrow + 1):nrow(epi), "Parent"])
    atleaf = is.na(nextrow)
  }
  if (atleaf) {
    if (rightchild) {
      return(paste("(", node, "[&type = \"I\"]:", epi[infectedatrow, 
                                                      "Rtime"] - epi[infectedatrow, "Itime"], ")[&type = \"E\"]:", 
                   epi[infectedatrow, "Itime"] - epi[lastrow, "Etime"], 
                   sep = ""))
    }
    else {
      return(paste(node, "[&type = \"I\"]:", epi[infectedatrow, 
                                                 "Rtime"] - epi[lastrow, "Etime"], sep = ""))
    }
  }
  else {
    nextrow = nextrow + lastrow
    if (rightchild) {
      return(paste("((", epi2binarynewick.int(epi, node = node, 
                                        lastrow = nextrow, infectedatrow = infectedatrow), 
                   ",", epi2binarynewick.int(epi, node = epi[nextrow, 
                                                       "Node ID"], lastrow = nextrow, infectedatrow = nextrow, 
                                       rightchild = TRUE), ")[&type = \"I\"]:", epi[lastrow, 
                                                                                    "Itime"] - epi[lastrow, "Etime"], ")[&type = \"E\"]:", 
                   epi[nextrow, "Etime"] - epi[lastrow, "Itime"], 
                   sep = ""))
    }
    else {
      return(paste("(", epi2binarynewick.int(epi, node = node, 
                                       lastrow = nextrow, infectedatrow = infectedatrow), 
                   ",", epi2binarynewick.int(epi, node = epi[nextrow, 
                                                       "Node ID"], lastrow = nextrow, infectedatrow = nextrow, 
                                       rightchild = TRUE), ")[&type = \"I\"]:", epi[nextrow, 
                                                                                    "Etime"] - epi[lastrow, "Etime"], sep = ""))
    }
  }
}

preprocessepi <- function (epi) {
  # hacky function, makes sampling times new infectious individuals so that we can use same epi2newick function
  # we must have numeric node indices so these are 100*node label - adjust for larger populations
  # thus sequence for 100*node corresponds to nth node label
  for (node in 1:nrow(epi)) {
    s = epi[node,'Stime']
    if (!is.na(s) && s != Inf) {
      epi <- rbind(epi, c((100*node),node, s, s, s, NA))
    }
  }
  return(epi)
}
