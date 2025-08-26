

simulate_epidemic <- function (C, N, beta, ki, thetai, ke = ki, thetae = thetai, ps, thetas= thetai, ks= ki, K) 
{
  
  # Compute Euclidean distance matrix
  coords_mat <- do.call(rbind, examplecoords)
  dist_matr <- as.matrix(dist(coords_mat))
  
  # SEIR compartment sets
  infectious <- exposed <- removed <- c()
  susceptible <- array(1:N)
  
  # infectious thresholds
  thresholds <- rexp(N, rate = 1)
  
  # records times of all events
  exposure_times <- rep(Inf, N)
  infectious_times <- rep(Inf, N)
  sampling_times <- rep(Inf, N)
  removal_times <- rep(Inf, N)

  # initialize first infection
  init <- sample(1:N, 1)
  infectious <- append(infectious, init)
  susceptible <- susceptible[!(susceptible == init)]
  
  # times for first infected
  exposure_times[init] <- 0
  infectious_times[init] <- rgamma(1, ke, scale = thetae) + exposure_times[init]
  sampling <- rbinom(n = 1, size = 1, prob = ps)
  sampling_times[init] <- Inf #ifelse(sampling, rgamma(1, ks+10, scale = thetas) + exposure_times[init], Inf)
  removal_times[init] <- min(rgamma(1, ki, scale = thetai) + infectious_times[init], sampling_times)
  
  # create infection list (Node index, donor index, exposed time, infectious time, removed time, sampled time)
  infections_list <- matrix(c(init, NA, 0, infectious_times[init], removal_times[init], NA), nrow = 1)
  # set current time
  t <- infectious_times[init]

  # simulate rest of outbreak
  for (x in 2:(3*N)) {
    
    if (length(susceptible) == 0 || (length(exposed)==0 && length(infectious)==0)) {
      break
    }
    
    if (length(infectious) > 0){
      # waiting time for infection s -> e 
      transmission_rate <- 0
      for (j in susceptible){
        transmission_rate <- transmission_rate + sum(sapply(infectious, function(i) K(dist_matr[i, j])))
      }
      transmission_rate <- transmission_rate * beta / N
      transmission_time <- t + rexp(1, rate=transmission_rate)
    }
    else { transmission_time <- Inf }
    
    # check event time is before any other event
    if ((length(exposed)==0 || all(transmission_time < infectious_times[exposed])) &&
      (length(infectious)==0 || all(transmission_time < removal_times[infectious]))) {
      # choose donor-recipient pair
      dists <- dist_matr[susceptible, infectious, drop=FALSE]
      probs <- 1 / dists # shorter dist = higher prob
      probs <- probs / sum(probs)
      picked_index <- sample(length(probs), size = 1, prob = as.vector(probs))
      pair <- arrayInd(picked_index, dim(probs))
      next_infected <- susceptible[pair[1]]
      donor <- infectious[pair[2]]
    
      # update times and compartments
      exposure_times[next_infected] <- t
      infectious_times[next_infected] <- rgamma(1, ke, scale = thetae) + exposure_times[next_infected]
      sampling <- rbinom(n = 1, size = 1, prob = ps)
      sampling_times[next_infected] <- ifelse(sampling, rgamma(1, ks, scale = thetas) + exposure_times[next_infected], Inf)
      removal_times[next_infected] <- min(rgamma(1, ki, scale = thetai) + infectious_times[next_infected], sampling_times)
      
      susceptible <- susceptible[susceptible != next_infected, drop=TRUE]
      exposed <- append(exposed, next_infected)
    
      infections_list <- rbind(infections_list, c(next_infected, donor, exposure_times[next_infected], infectious_times[next_infected], removal_times[next_infected], sampling_times[next_infected]))
      t <- transmission_time
    
    }
    # else another event occurs first (e -> i or i -> r)
    else{
      # future fix: refactor to be a bit cleaner
      min_exposed <- min(exposure_times[exposure_times > t], Inf)
      min_infectious <- min(infectious_times[infectious_times > t], Inf)
      min_removal <-min(removal_times[removal_times > t], Inf)
      if (min_exposed < min_infectious && min_exposed < min_removal) {
        ind <- match(min_exposed, exposure_times)
        susceptible <- susceptible[susceptible != ind, drop=TRUE]
        exposed <- append(exposed, ind)
        t <- min_exposed
      }
      else if (min_infectious < min_removal) {
        ind <- match(min_infectious, infectious_times)
        exposed <- exposed[exposed != ind, drop=TRUE]
        infectious <- append(infectious, ind)
        t <- min_infectious
      }
      else {
        ind <- match(min_removal, removal_times)
        infectious <- infectious[infectious != ind, drop=TRUE]
        t <- min_removal
      }
    }
  }
  
  if (length(susceptible)>0) {
    infections_list <- rbind(infections_list, cbind(susceptible, NA, NA, NA, NA, NA))
  }
  colnames(infections_list) <- c("Node ID", "Parent", "Etime", "Itime", "Rtime", "Stime")

  return(infections_list)
}
