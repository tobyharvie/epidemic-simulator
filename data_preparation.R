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
writeLines(nexus_out, "BREATH/seqs.nex")

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
    paste0('>', node, "_", round(days), "\n", seq)
  }
  else{ paste0('')}
})
fasta_out <- c(
  matrix_lines_new
)
writeLines(fasta_out, "ScITree/gen_inputs/seqs.fasta")

# -----------------
# SCITREE inputs
#-------------------

# covariates
covariates <- cbind(epi[,1], NA, NA, epi[,3],epi[,4],epi[,5], 1, 1, NA)
colnames(covariates) <- c("k","coor_x","coor_y","t_e","t_i","t_r","ftype0","herdn","initial_source")
unobserved <- c()
for (i in 1:nrow(epi)){
  if(is.na(epi[i,6])) { unobserved <- append(unobserved, i)}
  covariates[i,2] <- coords[epi[i,1]][[1]][1] # x
  covariates[i,3] <- coords[epi[i,1]][[1]][2] # y
  covariates[i,9] <- 9999
  # arbitrary extreme value needed for comparison
  if (is.na(covariates[i,4])) { covariates[i,4]=9e+10 }
  if (is.na(covariates[i,5])) { covariates[i,5]=9e+10 }
  if (is.na(covariates[i,6])) { covariates[i,6]=9e+10 }
}
if (length(unobserved)>0){
  covariates <- covariates[-unobserved,]
}
covariates <- as.data.frame(covariates)
moves.inputs <- as.data.frame(matrix(NA, nrow = 100, ncol = 3))
pars.aux <- data.frame('n' = nrow(covariates),
                       'kernel_type' = 'exponential',
                       'coord_type' = 'cartesian',
                       't_max' = max_time,
                       'unassigned_time' =  9e+10,
                       'processes' = 1,
                       'n_seq' = 1,
                       'n_base' = nbases,
                       'n_iterations' = 1e5,
                       'n_frequ' = 10,
                       'n_output_source' = 1000,
                       'n_output_gm' = 2000,
                       'n_cout' = 1000,
                       'opt_latgamma' = 1, # assume gamma dist for latent period
                       'opt_k80' = 0,
                       'opt_betaij' = 0,
                       'opt_ti_update' = 0,
                       'opt_mov' = 0,
                       stringsAsFactors = F)
# inference parameters
para.key.inits <- data.frame('alpha' = tb, #IMPORTANT
                             'beta' = tb/2, # this is under development, doesnt matter
                             'lat_mu' = ke,
                             'lat_var' = thetae,
                             'c' = ki,
                             'd' = thetai,
                             'k_1' = ka,
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
para.priors <- data.frame('t_range' = ke,
                          't_back' = max(epi[,4]-epi[,3], na.rm = TRUE), #max assuemd latent period
                          't_bound_hi' = max(epi[,5]-epi[,3], na.rm = TRUE),
                          'rate_exp_prior' = 0.001,
                          'ind_n_base_part' = 0, 
                          'n_base_part' = nbases,
                          'alpha_hi' = 0.1,
                          'beta_hi' = 30,
                          'mu_lat_hi' = 50,
                          'var_lat_lo' = 0.1,
                          'var_lat_hi' = 50,
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
para.sf <- data.frame('alpha_sf' = tb*1.5, # do some more work on this
                      'beta_sf' = tb*0.55,
                      'mu_lat_sf'   = 2,
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
outdir <- "./ScITree/inputs/"
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
write.table(cbind(epi[complete.cases(epi),][,1],epi[complete.cases(epi),][,6]), file = "SCOTTI/dates.csv",sep = ",",row.names = FALSE, col.names = FALSE)
# hosts .. unneeded?
write.table(cbind(epi[complete.cases(epi),][,1],epi[complete.cases(epi),][,1]), file = "SCOTTI/hosts.csv",sep = ",",row.names = FALSE, col.names = FALSE)
# earliest and latest times host in infection
write.table(cbind(epi[complete.cases(epi),][,1],epi[complete.cases(epi),][,3],epi[complete.cases(epi),][,5]), file = "SCOTTI/hostTimes.csv",sep = ",",row.names = FALSE, col.names = FALSE)
# use py script to generate xml