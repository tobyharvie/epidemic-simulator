function (M, N, beta, ki, thetai, ke = ki, thetae = thetai, latencydist = "fixed", 
          latencyperiod = 0) 
{
  if (!is.null(M)) {
    if (!is.matrix(M)) 
      stop("Input error: Network M must be an edgelist matrix.")
    if ((length(dim(M)) != 2) || (dim(M)[2] != 2)) 
      stop("Input error: Network M must an a 2-dimensional edgelist matrix.")
  }
  net = network::network.initialize(N, directed = FALSE)
  network::add.edges(net, head = M[, 1], tail = M[, 2])
  t_net = array(FALSE, network::network.edgecount(net))
  all_edge = unlist(net$mel)
  v1 = all_edge[seq(1, length(all_edge), by = 3)]
  v2 = all_edge[seq(2, length(all_edge), by = 3)]
  is_infectious = is_exposed = is_removed = array(FALSE, N)
  is_susceptible = array(TRUE, N)
  init <- sample(1:N, 1)
  is_infectious[init] = TRUE
  is_susceptible[init] = FALSE
  neighbour_id = network::get.edgeIDs(net, init)
  susc_neighbour_id = neighbour_id[is_susceptible[v1[neighbour_id]] | 
                                     is_susceptible[v2[neighbour_id]]]
  t_net[susc_neighbour_id] = TRUE
  t.time = array(dim = N)
  r.time = array(dim = N)
  t.time[init] <- ifelse(latencydist == "fixed", latencyperiod, 
                         rgamma(1, ke, scale = thetae))
  r.time[init] = rgamma(1, ki, scale = thetai) + t.time[init]
  nextrec = init
  inf.list <- matrix(c(init, NA, 0, t.time[init], NA), nrow = 1)
  time <- cm.time <- t.time[init]
  nexttrans = init
  t.time[init] <- Inf
  count_i = 1
  count_e = 0
  for (i in 2:(N * 3)) {
    n.si = sum(t_net)
    dwt <- ifelse(count_i > 0, r.time[nextrec] - cm.time, 
                  Inf)
    bwt <- ifelse(n.si != 0, rexp(1, beta * n.si), Inf)
    twt <- t.time[nexttrans] - cm.time
    ewt <- min(bwt, dwt, twt, na.rm = TRUE)
    time <- c(time, ewt)
    cm.time <- cm.time + ewt
    if (ewt == bwt) {
      test <- "Infect"
    }
    else if (ewt == dwt) {
      test <- "removal"
    }
    else {
      test <- "transition"
    }
    if (test == "Infect") {
      trans_edge_id = ifelse(n.si > 1, sample(which(t_net), 
                                              1), which(t_net))
      if (is_susceptible[v1[trans_edge_id]]) {
        new.inf = v1[trans_edge_id]
        parent = v2[trans_edge_id]
      }
      else {
        new.inf = v2[trans_edge_id]
        parent = v1[trans_edge_id]
      }
      lat <- ifelse(latencydist == "fixed", latencyperiod, 
                    rgamma(1, ke, scale = thetae))
      t.time[new.inf] <- cm.time + lat
      is_exposed[new.inf] = TRUE
      is_susceptible[new.inf] = FALSE
      count_e = count_e + 1
      t_net[network::get.edgeIDs(net, new.inf)] = FALSE
      inf.list <- rbind(inf.list, c(new.inf, parent, cm.time, 
                                    NA, NA))
      nexttrans <- which(t.time == min(t.time, na.rm = TRUE))
    }
    else if (test == "removal") {
      if (i == 2) {
        inf.list[1, 5] <- cm.time
        break
      }
      new.rec <- nextrec
      is_infectious[new.rec] = FALSE
      is_removed[new.rec] = TRUE
      count_i = count_i - 1
      t_net[network::get.edgeIDs(net, new.rec)] = FALSE
      inf.list[which(inf.list[, 1] == new.rec), 5] <- cm.time
      r.time[nextrec] <- NA
      if (count_i > 0) 
        nextrec <- which(r.time == min(r.time, na.rm = TRUE))
      else if (count_e > 0) {
        nextrec <- which(t.time == min(t.time, na.rm = TRUE))
        r.time[nextrec] <- Inf
      }
    }
    else {
      new.trans <- nexttrans
      is_exposed[new.trans] = FALSE
      is_infectious[new.trans] = TRUE
      count_i = count_i + 1
      count_e = count_e - 1
      neighbour_id = network::get.edgeIDs(net, new.trans)
      susc_neighbour_id = neighbour_id[is_susceptible[v1[neighbour_id]] | 
                                         is_susceptible[v2[neighbour_id]]]
      t_net[susc_neighbour_id] = TRUE
      inf.list[which(inf.list[, 1] == new.trans), 4] <- cm.time
      t.time[nexttrans] <- NA
      nexttrans <- which(t.time == min(t.time, na.rm = TRUE))
      r.time[new.trans] <- cm.time + rgamma(1, ki, scale = thetai)
      if (r.time[new.trans] < r.time[nextrec]) 
        nextrec <- new.trans
    }
    if (count_e + count_i == 0) {
      break
    }
  }
  inf.list[, 3:5] <- inf.list[, 3:5] - min(inf.list[, 5])
  if (any(is_susceptible)) 
    inf.list <- rbind(inf.list, cbind(which(is_susceptible), 
                                      NA, NA, NA, NA))
  colnames(inf.list) <- c("Node ID", "Parent", "Etime", "Itime", 
                          "Rtime")
  class(inf.list) <- c("epidemic", "matrix")
  return(inf.list)
}