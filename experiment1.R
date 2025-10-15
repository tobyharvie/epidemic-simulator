library(epinet)
library(ape)
library(phyclust)
source('plotting.R')
source('epi2binarynewick.R')
simulate_epidemic <- function (C, N, beta, ki, thetai, ke = ki, thetae = thetai, ps, K) 
{
  
  # Compute Euclidean distance matrix
  coords_mat <- do.call(rbind, C)
  dist_matr <- as.matrix(dist(coords_mat))
  
  # SEIR compartment sets
  infectious <- exposed <- removed <- c()
  susceptible <- array(1:N)
  
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
  time_infected <- rgamma(1, ke, scale = thetae)
  infectious_times[init] <- time_infected + exposure_times[init]
  #removal_times[init] <- min(rgamma(1, ki, scale = thetai) + infectious_times[init], sampling_times)
  removal_times[init] <- rgamma(1, ki, scale = thetai) + infectious_times[init]
  sampling <- rbinom(n = 1, size = 1, prob = ps)
  #sampling_times[init] <- ifelse(sampling, rgamma(1, ks+10, scale = thetas) + exposure_times[init], NA)
  sampling_rate <- -log(1-ps)/time_infected
  sampling_time <- rexp(1, sampling_rate) + infectious_times[init]
  sampling_times[init] <- ifelse(sampling_time<removal_times[init], sampling_time, NA)
  
  
  # create infection list (Node index, donor index, exposed time, infectious time, removed time, sampled time)
  infections_list <- matrix(c(init, NA, 0, infectious_times[init], removal_times[init], sampling_times[init]), nrow = 1)
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
      # choose donor-recipient pair based on spatial kernel function
      dists <- 1/K(dist_matr[susceptible, infectious, drop=FALSE])
      probs <- 1 / dists # shorter dist = higher prob
      probs <- probs / sum(probs)
      #print(probs)
      picked_index <- sample(length(probs), size = 1, prob = as.vector(probs))
      pair <- arrayInd(picked_index, dim(probs))
      next_infected <- susceptible[pair[1]]
      donor <- infectious[pair[2]]
      
      # update times and compartments
      exposure_times[next_infected] <- transmission_time
      infectious_times[next_infected] <- rgamma(1, ke, scale = thetae) + exposure_times[next_infected]
      time_infected <- rgamma(1, ki, scale = thetai)
      removal_times[next_infected] <- time_infected + infectious_times[next_infected]
      # chosen so that (approx) ps samples occur before end of infectious period
      # assume sampling happens during infectious period (not in exposed period)
      #sampling_rate <- -log(1-ps)/time_infected
      #sampling_time <- rexp(1, sampling_rate) + infectious_times[next_infected]
      #sampling_times[next_infected] <- ifelse(sampling_time<removal_times[next_infected], sampling_time, NA)
      sampling_times[next_infected] <- ifelse(runif(1)<ps, removal_times[next_infected], NA)
      #if (!is.na(sampling_times[next_infected])) {removal_times[next_infected] <- min(removal_times[next_infected], sampling_times[next_infected])}
      
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
  
  return(round(infections_list,3))
}

prepare_data <- function (epi,sequences,kappa){
  dir.create(paste(dir,'BREATH',sep=''))
  dir.create(paste(dir,'SCOTTI',sep=''))
  dir.create(paste(dir,'ScITree',sep=''))
  dir.create(paste(dir,'ScITree/gen_inputs',sep=''))
  dir.create(paste(dir,'ScITree/inputs',sep=''))
  
  
  
  nexus_str<-sequences
  # PREPROCESSING
  # Example string
  # Split into lines
  lines <- unlist(strsplit(nexus_str, "\n"))
  # Keep only labelled rows (start with digits or letters)
  labelled <- grep("^[[:alnum:]_]+", lines, value = TRUE)
  # Extract numeric labels
  labels <- suppressWarnings(as.numeric(sub("^([0-9]+).*", "\\1", labelled)))
  # Filter multiples of (n_pop*100)
  keep <- labelled[!is.na(labels) & labels %% (100) == 0]
  # Replace the numeric labels with divided-by-(n_pop*100) labels
  keep_div <- sapply(keep, function(line) {
    num <- as.numeric(sub("^([0-9]+).*", "\\1", line))
    new_num <- num / (100)
    sub("^[0-9]+", new_num, line)
  })
  # Collapse back into a string
  cleaned_str <- paste(keep_div, collapse = "\n")
  #convert back to NEXUS
  matrix_lines <- unlist(strsplit(cleaned_str, "\n"))
  matrix_lines <- matrix_lines[nchar(matrix_lines) > 0]
  
  # BREATH requires NEXUS
  # switch unit to days
  time_range <- 1000 # in days
  max_time <- max(epi[is.finite(epi)], na.rm = TRUE, finite = TRUE)
  start_date <- as.Date("2020-01-01") # arbitrary starting date
  # include time data
  matrix_lines_new <- sapply(matrix_lines, function(line) {
    node <- as.numeric(sub("^([0-9]+).*", "\\1", line))
    seq <- sub("^[0-9]+\\s+", "", line)
    # Calculate timestamp as a Date
    days <- (epi[epi[, 1] == node, 6] / max_time) * time_range
    ts <- start_date + round(days)
    # Build new label
    paste0(node, ":", ts, "    ", seq)
  })
  nexus_out <- c(
    "#NEXUS","","BEGIN DATA;",
    paste0("    DIMENSIONS NTAX=", length(matrix_lines_new), " NCHAR=", max(nchar(gsub(".*?\\s+", "", matrix_lines_new))), ";"),
    "    FORMAT DATATYPE=DNA MISSING=? GAP=-;",
    "MATRIX",matrix_lines_new,";","END;"
  )
  writeLines(nexus_out, paste(dir,"BREATH/seqs.nex",sep=''))
  
  # calculate sampling shape and rate for sampling hazard
  stimes <- (epi[,6]-epi[,3])/max_time*1000
  sk <-mean(stimes,na.rm=TRUE)
  sv <-var(stimes,na.rm=TRUE)
  sa <-sk^2/sv #shape
  sb <-sk/sv #rate
  print(paste('Sampling prob:', length(epi[!is.na(epi[, 6]), 6])/length(epi[,1])))
  print(paste('Sampling rate:', sb))
  print(paste('Sampling shape:', sa))
  
  ttimes <- (epi[,4]-epi[,3])/max_time*1000
  tk <-mean(ttimes,na.rm=TRUE)
  tv <-var(ttimes,na.rm=TRUE)
  ta <-tk^2/tv #shape
  tb <-tk/tv #rate
  print(paste('Transmission rate:', tb))
  print(paste('Transmission shape:', ta))
  # average number of people that a single node infects
  tc <- mean(table(factor(epi[,2],levels=1:length(epi[,1]))))
  print(paste("transmission constant: ", tc))
  
  # sanity checks
  curve(dgamma(x, shape = sa, scale = 1/sb),from = 0, to = max_time)
  curve(dgamma(x, shape = ta, scale = 1/tb),from = 0, to = max_time)
  
  # ScITree requires FASTA
  matrix_lines_new <- sapply(matrix_lines, function(line) {
    node <- as.numeric(sub("^([0-9]+).*", "\\1", line))
    seq <- sub("^[0-9]+\\s+", "", line)
    # Calculate timestamp as a Date
    days <- (epi[epi[, 1] == node, 6] / max_time) * time_range
    # Build new label
    if (!is.na(epi[epi[, 1] == node, 6])){
      paste0('>', node-1, "_", round(epi[epi[, 1] == node, 6]), "\n", seq)
    }
    else{ 
      paste0('')
    }
  })
  fasta_out <- c(
    sort(matrix_lines_new) # sorting ensures correct ordering for scitree
  )
  writeLines(fasta_out, paste(dir,"ScITree/gen_inputs/seqs.fasta",sep=''))
  
  # -----------------
  # SCITREE inputs
  #-------------------
  
  # covariates
  covariates <- cbind(epi[,1]-1, NA, NA, epi[,3],epi[,4],epi[,5], 1, 1, NA)
  colnames(covariates) <- c("k","coor_x","coor_y","t_e","t_i","t_r","ftype0","herdn","initial_source")
  unobserved <- c()
  for (i in 1:nrow(epi)){
    
    if(is.na(epi[i,6])) { unobserved <- append(unobserved, i)}
    if(is.na(epi[i,3])) { next }
    covariates[i,2] <- coords[epi[i,1]][[1]][1] # x
    covariates[i,3] <- coords[epi[i,1]][[1]][2] # y
    print(epi[epi[,1]==epi[i,2],6])
    if (i==1){ covariates[i,9]<-9999}
    else{
    covariates[i,9] <- ifelse(!is.na(epi[epi[,1]==epi[i,2],6]),epi[i,2]-1,-99)
    }
    # arbitrary extreme value needed for comparison
    if (is.na(covariates[i,4])) { covariates[i,4]=9e+10 }
    if (is.na(covariates[i,5])) { covariates[i,5]=9e+10 }
    if (is.na(covariates[i,6])) { covariates[i,6]=9e+10 }
  }
  if (length(unobserved)>0){
    covariates <- covariates[-unobserved,]
  }
  
  covariates <- as.data.frame(covariates)
  covariates <- covariates[order(covariates$'k'), ]
  print(covariates)
  moves.inputs <- as.data.frame(matrix(1, nrow = 100, ncol = 3))
  pars.aux <- data.frame('n' = nrow(covariates),
                         'kernel_type' = 'exponential',
                         'coord_type' = 'cartesian',
                         't_max' = max_time,
                         'unassigned_time' =  9e+10,
                         'processes' = 1,
                         'n_seq' = 1,
                         'n_base' = nbases,
                         'n_iterations' = 10,
                         'n_frequ' = 1,
                         'n_output_source' = 1,
                         'n_output_gm' = 1,
                         'n_cout' = 1,
                         'opt_latgamma' = 1, # assume gamma dist for latent period
                         'opt_k80' = 0,
                         'opt_betaij' = 0,
                         'opt_ti_update' = 0,
                         'opt_mov' = 0,
                         stringsAsFactors = F)
  # inference parameters
  para.key.inits <- data.frame('alpha' = 0.00025, #IMPORTANT
                               'beta' = tb/2, # this is under development, doesnt matter
                               'lat_mu' = ke/50,
                               'lat_var' = thetae/50,
                               'c' = ki/50,
                               'd' = thetai/50,
                               'k_1' = ka/50,
                               'mu_1' = 3e-05, # next 3 realted to substitution, example
                               'mu_2' = 1e-06,
                               'p_ber' = 0.2,
                               'phi_inf1' = 1, # next related to farms.
                               'phi_inf2' = 1,
                               'rho_susc1' = 1,
                               'rho_susc2' = 1,
                               'nu_inf' = 0.2,
                               'tau_susc'= 0.1,
                               'beta_m'= 1)
  para.priors <- data.frame('t_range' = ke/max_time*time_range,
                            't_back' = max(epi[,4]-epi[,3], na.rm = TRUE)/max_time*time_range, #max assuemd latent period
                            't_bound_hi' = max(epi[,5]-epi[,3], na.rm = TRUE)/max_time*time_range,
                            'rate_exp_prior' = 0.01,
                            'ind_n_base_part' = 0, 
                            'n_base_part' = nbases,
                            'alpha_hi' = 2,
                            'beta_hi' = 30,
                            'mu_lat_hi' = 300,
                            'var_lat_lo' = 0.1,
                            'var_lat_hi' = 300,
                            'c_hi' = 100,
                            'd_hi' = 100,
                            'k_1_lo' = 0,
                            'k_1_hi' = 10,
                            'mu_1_hi' = 44,
                            'mu_2_hi' = 0.1,
                            'p_ber_hi' = 1.0,
                            'phi_inf1_hi' = 500,
                            'phi_inf2_hi' = 500,
                            'rho_susc1_hi' = 500,
                            'rho_susc2_hi' = 500,
                            'nu_inf_lo' = 0,
                            'nu_inf_hi' = 1,
                            'tau_susc_lo' = 0,
                            'tau_susc_hi' = 1,
                            'beta_m_hi' = 5,
                            'trace_window' = 20)
  para.sf <- data.frame('alpha_sf' = 0.00015, # do some more work on this
                        'beta_sf' = tb*0.55,
                        'mu_lat_sf'   = 1,
                        'var_lat_sf' = 1,
                        'c_sf' = 1.25,
                        'd_sf' = 0.75,
                        'k_1_sf' = 0.01,
                        'mu_1_sf' = 0.2,
                        'mu_2_sf' = 2.5e-6,
                        'p_ber_sf' = 0.02,
                        'phi_inf1_sf' = 1.75,
                        'phi_inf2_sf' = 1.5,
                        'rho_susc1_sf' = 1,
                        'rho_susc2_sf' = 1.25,
                        'nu_inf_sf' = 0.25,
                        'tau_susc_sf' = 0.25,
                        'beta_m_sf' = 1)
  outdir <- paste(dir,"ScITree/inputs/",sep='')
  write.csv(covariates,     paste(outdir,"covariates.csv",sep=''),     row.names = FALSE)
  write.csv(moves.inputs,   paste(outdir,"moves.inputs.csv",sep=''),   row.names = FALSE)
  write.csv(pars.aux,       paste(outdir,"pars.aux.csv",sep=''),       row.names = FALSE)
  write.csv(para.key.inits, paste(outdir,"para.key.inits.csv",sep=''), row.names = FALSE)
  write.csv(para.priors,    paste(outdir,"para.priors.csv",sep=''),    row.names = FALSE)
  write.csv(para.sf,        paste(outdir,"para.sf.csv",sep=''),        row.names = FALSE)
  
  # ---------------
  # SCOTTI Pipeline
  # ---------------
  # FASTA
  matrix_lines_new <- sapply(matrix_lines, function(line) {
    node <- as.numeric(sub("^([0-9]+).*", "\\1", line))
    seq <- sub("^[0-9]+\\s+", "", line)
    # Build new label
    if (!is.na(epi[epi[, 1] == node, 6])){
      paste0('>', node, "\n", seq)
    }
    else{ paste0('')}
  })
  fasta_out <- c(
    matrix_lines_new
  )
  writeLines(fasta_out, "SCOTTI/seqs.fasta")
  # sampling dates
  write.table(cbind(epi[complete.cases(epi),][,1],epi[complete.cases(epi),][,6]), file = paste(dir,"SCOTTI/dates.csv",sep=''),sep = ",",row.names = FALSE, col.names = FALSE)
  # hosts .. unneeded?
  write.table(cbind(epi[complete.cases(epi),][,1],epi[complete.cases(epi),][,1]), file = paste(dir,"SCOTTI/hosts.csv",sep=''),sep = ",",row.names = FALSE, col.names = FALSE)
  # earliest and latest times host in infection
  # assume sampled during infectious period
  # thus earliest possible time infectious period + exposed period at 95th percentile
  lower <- qgamma(p = 0.99, shape=ke, scale = thetae) + qgamma(p = 0.99, shape=ki, scale = thetai)
  # for upper, we assume that we remove after sampling
  upper <- 1
  write.table(cbind(epi[complete.cases(epi),][,1],pmax(0,epi[complete.cases(epi),][,6]-lower),epi[complete.cases(epi),][,6]+upper), file = paste(dir,"SCOTTI/hostTimes.csv",sep=''),sep = ",",row.names = FALSE, col.names = FALSE)
  # use py script to generate xml
}

dir <- './spatial-experiment/'

# epidemic params
n_pop <- 60
# k shape, theta scale. realistic variables from keeling 2001
ki <- 3; thetai <- 2.17; ke <- 3; thetae <- 1.17; ps <- 0.9;
tc <- 2 # transmission constant R0

# sequencing params
nbases <- 500

# grid search of kappa
for (ka in seq(from = 0, to = 0, by = 0)){
  # in order to recreate.
  spatialkernel <- function (d, kappa=ka) {
    return(exp(-kappa * d))
  }
  beta <- exp(0.5*ka)*tc/(ki*thetai)
  print(paste('Simulating epidemic with kappa=',ka,'beta=',beta))
  seed_used <- 1
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
      # manual inspection
      plot.spatial(coords, epi)
      plot.epidemic(epi[,-6])
      #save
      newick <- newick_tree(epi)
      sequences <- seqgen(opts = paste("-mHKY -l",nbases," -on -s0.0005",sep=''), newick.tree = newick)
      
      dir <- paste("spatial-experiment/kappa-",ka*10,'/',sep='')
      dir.create(dir)
      save(epi,newick,seed_used,file=paste(dir, '/epi-data',sep=''))
      
      prepare_data(epi,sequences,ka)
      
      break
    }
    else{
      print(epi)
    }
    seed_used <- seed_used + 1
  }
}


#newick <- newick_tree(epi)
#nbases <- 100
#sequences <- seqgen(opts = paste("-mHKY -l",nbases," -on -s0.001",sep=''), newick.tree = newick)

# plotting 
#plot.epidemic(epi[,-6])
#plot.spatial(coords, epi)
#plot.phylo(newick)

# epinet (for comparison)
#mycov <- data.frame(id = 1:n_pop, xpos = runif(n_pop), ypos = runif(n_pop))
#dyadCov <- BuildX(mycov, binaryCol = list(c(2, 3)),binaryFunc = "euclidean")
#eta <- c(0, -7)
#net <- SimulateDyadicLinearERGM(N = n_pop, dyadiccovmat = dyadCov, eta = eta)
#epi2 <- SEIR.simulator(M = net, N = n_pop, beta = 1, ki = 3, thetai = 7,ke = 3, latencydist = "gamma")
#plot.epidemic(epi2)

# data pipelines for analysis
#source("./data_preparation.R")

# write sampled data as NEXUS format for BEAST and SCOTTI inference

