simulate_epidemic <- function (C, N, alpha, beta, ki, thetai, ke = ki, thetae = thetai, ps, thetas= thetai, ks= ki, latencyperiod = 0, K) 
{
  # TO DO: choose realisitic parameters
  
  if (!is.null(C)) {
    if (!is.list(C)) 
      stop("Input error: C must be a list.")
    if ((length(C) != N) && all(sapply(C, function(x) is.numeric(x) && length(x) == 2)))
      stop("Input error: Vector C.")
  }
  
  coords_mat <- do.call(rbind, examplecoords)
  # Compute Euclidean distance matrix
  dist_matr <- as.matrix(dist(coords_mat))
  
  # SEIR compartment sets
  infectious <- exposed <- removed <- c()
  susceptible <- array(1:N)
  
  # infectious thresholds
  thresholds <- rexp(N, rate = 1)
  
  # records times of all events
  exposure.time <- rep(Inf, N)
  infectious.time <- rep(Inf, N)
  sampling.time <- rep(Inf, N)
  removal.time <- rep(Inf, N)
  
  # initialize first infection
  init <- sample(1:N, 1)
  infectious <- append(infectious, init)
  susceptible <- susceptible[!(susceptible == init)]
  
  # times for first infected
  exposure.time[init] <- rgamma(1, ke, scale = thetae)
  infectious.time[init] <- rgamma(1, ki, scale = thetai) + exposure.time[init]
  sampling <- rbinom(n = 1, size = 1, prob = ps)
  sampling.time[init] <- ifelse(sampling, rgamma(1, ks, scale = thetas) + exposure.time[init], Inf)
  removal.time[init] <- min(rgamma(1, ki, scale = thetai) + infectious.time[init], sampling.time)
  
  # create infection list (Node index, donor index, exposed time, infectious time, removed time, sampled time)
  infections.list <- matrix(c(init, NA, 0, infectious.time[init], NA, NA), nrow = 1)
  
  # set current time
  t <- infectious.time[init]
  
  # simulate rest of outbreak
  for (x in 2:N+1) {
    
    # waiting time for infection s -> e 
    min.event.time = Inf
    for (j in susceptible){
      # calculate infectious pressure
      pressure <- thresholds[j] - alpha*t -
      beta*sum(sapply(infectious, function(i) K(dist_matr[i, j])*(min(t, removal.time[i])-infectious.time[i])))
      # calculate infection time
      event.time <- pressure / (alpha + beta*sum(sapply(infectious, function(i) K(dist_matr[i, j])))) + t
      # choose soonest infection event
      if (event.time < min.event.time) {
        min.event.time <- event.time
        next.infected <- j
      }
    }
    
    # choose donor
    if (length(infectious) > 1) {
      dists <- dist_matr[next.infected, infectious]
      # inverse distance
      probs <- 1 / dists       
      probs <- probs / sum(probs)
      donor <- sample(infectious, size=1, prob=probs)
    }
    else { donor <- infectious[1]}
    
    t <- min.event.time
    
    # update times for newly infected individual
    exposure.time[next.infected] <- t
    infectious.time[next.infected] <- rgamma(1, ke, scale = thetae) + exposure.time[next.infected]
    sampling <- rbinom(n = 1, size = 1, prob = ps)
    sampling.time[next.infected] <- ifelse(sampling, rgamma(1, ks, scale = thetas) + exposure.time[next.infected], Inf)
    removal.time[next.infected] <- min(rgamma(1, ki, scale = thetai) + infectious.time[next.infected], sampling.time)
    
    infections.list <- rbind(infections.list, c(next.infected, donor, exposure.time[next.infected], infectious.time[next.infected], removal.time[next.infected], sampling.time[next.infected]))
    
    # update compartments
    susceptible <- susceptible[!(susceptible == next.infected)]
    exposed <- append(exposed, next.infected)
    for (ind in exposed) {
      if (infectious.time[ind] <= t) { 
        infectious <- append(infectious, ind)
        exposed <- exposed[!(exposed == ind)]
      }
      if (removal.time[ind] <= t) {
        removed <- append(removed, ind)
        infectious <- infectious[!(infectious == ind)]
      }
    }
    
  }
    
  return(infections.list)
}