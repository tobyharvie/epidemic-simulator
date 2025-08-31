data_preparation <- function(epi,nexus_str) {
  
  # Example string
  # Split into lines
  lines <- unlist(strsplit(nexus_str, "\n"))
  # Keep only labelled rows (start with digits or letters)
  labelled <- grep("^[[:alnum:]_]+", lines, value = TRUE)
  # Extract numeric labels
  labels <- suppressWarnings(as.numeric(sub("^([0-9]+).*", "\\1", labelled)))
  # Filter multiples of 100
  keep <- labelled[!is.na(labels) & labels %% 100 == 0]
  # Replace the numeric labels with divided-by-100 labels
  keep_div <- sapply(keep, function(line) {
    num <- as.numeric(sub("^([0-9]+).*", "\\1", line))
    new_num <- num / 100
    sub("^[0-9]+", new_num, line)
  })
  # Collapse back into a string
  cleaned_str <- paste(keep_div, collapse = "\n")
  #convert back to NEXUS
  matrix_lines <- unlist(strsplit(cleaned_str, "\n"))
  matrix_lines <- matrix_lines[nchar(matrix_lines) > 0]
  
  # BREATH counts in days
  time_range <- 1000 # in days
  max_time <- max(epi[is.finite(epi)], na.rm = TRUE, finite = TRUE)
  start_date <- as.Date("2020-01-01") # arbitrary starting date
  # include time data
  matrix_lines_new <- sapply(matrix_lines, function(line) {
    node <- as.numeric(sub("^([0-9]+).*", "\\1", line))
    seq <- sub("^[0-9]+\\s+", "", line)
    # Calculate timestamp as a Date
    days <- (epi[[node,6]] / max_time) * time_range
    ts <- start_date + round(days)
    # Build new label
    paste0(node, ":", ts, "    ", seq)
  })
  
  nexus_out <- c(
    "#NEXUS",
    "",
    "BEGIN DATA;",
    paste0("    DIMENSIONS NTAX=", length(matrix_lines_new), " NCHAR=", max(nchar(gsub(".*?\\s+", "", matrix_lines_new))), ";"),
    "    FORMAT DATATYPE=DNA MISSING=? GAP=-;",
    "MATRIX",
    matrix_lines_new,
    ";",
    "END;"
  )
  writeLines(nexus_out, "filtered_sequences.nex")
  
  # calculate sampling shape and rate for sampling hazard
  stimes <- (epi[,6]-epi[,3])/max_time*1000
  sk <-mean(stimes,na.rm=TRUE)
  sv <-var(stimes,na.rm=TRUE)
  sa <-sk^2/sv #shape
  sb <-sk/sv #variance
  print(paste('Sampling rate:', sa))
  print(paste('Sampling shape:', sb))
  
  itimes <- (epi[,4]-epi[,3])/max_time*1000
  ik <-mean(itimes,na.rm=TRUE)
  iv <-var(itimes,na.rm=TRUE)
  ia <-ik^2/iv #shape
  ib <-ik/iv #variance
  print(paste('Transmission rate:', ia))
  print(paste('Transmission shape:', ib))
  print(paste("transmission constant: ", mean(table(factor(epi[,2],levels=1:30)))))
  
  curve(dgamma(x, shape = sa, scale = 1/sb),
        from = 0, to = 100)
  curve(dgamma(x, shape = ia, scale = 1/ib),
        from = 0, to = 100)
  
  # calculate transmission shape and rate for transmission hazard

}