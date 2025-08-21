simulate_epidemic <- function (C, N, alpha, beta, ki, thetai, ke = ki, thetae = thetai, latencyperiod = 0, K) 
{
  # alpha = infection rate
  # beta = contact parameter
  
  if (!is.null(C)) {
    if (!is.list(C)) 
      stop("Input error: C must be a list.")
    if ((length(C) != N) && all(sapply(C, function(x) is.numeric(x) && length(x) == 2)))
      stop("Input error: Vector C.")
  }
  
  coords_mat <- do.call(rbind, examplecoords)
  # Compute Euclidean distance matrix
  dist_matr <- as.matrix(dist(coords_mat))
  
  # infectious thresholds
  thresholds <- rexp(N, rate = 1)
  
  # initialize arrays
  infectious <- exposed <- removed <- c()
  susceptible <- array(1:N)
  
  # initialize first infection
  init <- sample(1:N, 1)
  infectious <- append(infectious, init)
  susceptible <- susceptible[!(susceptible == init)]
  
  # time of exposures
  exposure.time <- rep(Inf, N)
  # time of infections
  infectious.time <- rep(Inf, N)
  # time of recoveries
  removal.time <- rep(Inf, N)
  
  # exposure time of init
  exposure.time[init] <- rgamma(1, ke, scale = thetae)
  # infectious time of init
  infectious.time[init] <- rgamma(1, ki, scale = thetai) + exposure.time[init]
  
  # create infection list
  infections.list <- matrix(c(init, NA, 0, infectious.time[init], NA), nrow = 1)
  
  # get time
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
      #print(event.time)
      # update soonest infection event
      if (event.time < min.event.time) {
        min.event.time <- event.time
        next.infected <- j
      }
    }
    
    t <- min.event.time
    
    # update times for newly infected individual
    exposure.time[next.infected] <- t
    infectious.time[next.infected] = rgamma(1, ke, scale = thetae) + exposure.time[next.infected]
    removal.time[next.infected] = rgamma(1, ki, scale = thetai) + infectious.time[next.infected]
    
    susceptible <- susceptible[!(susceptible == next.infected)]
    exposed <- append(exposed, next.infected)
    
    # update compartments
    for (ind in exposed) {
      if (infectious.time[ind] <= t) { 
        infectious <- append(infectious, ind)
        susceptible <- susceptible[!(susceptible == ind)]
      }
      if (removal.time[ind] <= t) {
        removed <- append(removed, ind)
        infectious <- infectious[!(infectious == ind)]
      }
    }
    
    # update infectious list
    #infections.list[which(infections.list[, 1] == new.trans), 4] <- cm.time
    print(t)
  }
    
  print(exposure.time)
  print(infectious.time)  
  print(removal.time)
  return(infections.list)
}