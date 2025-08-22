simulate_epidemic <- function (C, N, alpha, beta, ki, thetai, ke = ki, thetae = thetai, ps, thetas= thetai, ks= ki, K) 
{
  # TO DO: choose realisitic parameters
  
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
  exposure_times[init] <- rgamma(1, ke, scale = thetae)
  infectious_times[init] <- rgamma(1, ki, scale = thetai) + exposure_times[init]
  sampling <- rbinom(n = 1, size = 1, prob = ps)
  sampling_times[init] <- ifelse(sampling, rgamma(1, ks, scale = thetas) + exposure_times[init], Inf)
  removal_times[init] <- min(rgamma(1, ki, scale = thetai) + infectious_times[init], sampling_times)
  
  # create infection list (Node index, donor index, exposed time, infectious time, removed time, sampled time)
  infections_list <- matrix(c(init, NA, 0, infectious_times[init], removal_times[init], NA), nrow = 1)
  
  # set current time
  t <- infectious_times[init]

  # simulate rest of outbreak
  for (x in 2:N+1) {
    
    # waiting time for infection s -> e 
    min_event_time = Inf
    for (j in susceptible){
      # calculate infectious pressure
      pressure <- thresholds[j] - alpha*t -
      beta*sum(sapply(infectious, function(i) K(dist_matr[i, j])*(min(t, removal_times[i])-infectious_times[i])))
      # calculate infection time
      event_time <- pressure / (alpha + beta*sum(sapply(infectious, function(i) K(dist_matr[i, j])))) + t
      # choose soonest infection event
      
      if (event_time < min_event_time) {
        min_event_time <- event_time
        next_infected <- j
      }
    }
    
    # choose donor
    if (length(infectious) > 1) {
      dists <- dist_matr[next_infected, infectious]
      # inverse distance
      probs <- 1 / dists       
      probs <- probs / sum(probs)
      donor <- sample(infectious, size=1, prob=probs)
    }
    else { donor <- infectious[1]}
    
    t <- min_event_time
    
    # update times for newly infected individual
    exposure_times[next_infected] <- t
    infectious_times[next_infected] <- rgamma(1, ke, scale = thetae) + exposure_times[next_infected]
    sampling <- rbinom(n = 1, size = 1, prob = ps)
    sampling_times[next_infected] <- ifelse(sampling, rgamma(1, ks, scale = thetas) + exposure_times[next_infected], Inf)
    removal_times[next_infected] <- min(rgamma(1, ki, scale = thetai) + infectious_times[next_infected], sampling_times)
    
    infections_list <- rbind(infections_list, c(next_infected, donor, exposure_times[next_infected], infectious_times[next_infected], removal_times[next_infected], sampling_times[next_infected]))
    
    # update compartments
    susceptible <- susceptible[!(susceptible == next_infected)]
    exposed <- append(exposed, next_infected)
    for (ind in exposed) {
      if (infectious_times[ind] <= t) { 
        infectious <- append(infectious, ind)
        exposed <- exposed[!(exposed == ind)]
      }
      if (removal_times[ind] <= t) {
        removed <- append(removed, ind)
        infectious <- infectious[!(infectious == ind)]
      }
    }
    
  }
    
  return(infections_list)
}