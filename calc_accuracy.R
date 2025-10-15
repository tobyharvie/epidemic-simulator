handlescottiformat <- function (filename) {
  # Read the file
  lines <- readLines(filename)
  
  # Find section for "Probabilities direct transmission"
  start_idx <- grep("^Probabilities direct transmission", lines)
  end_idx <- grep("^Probabilities", lines)
  end_idx <- end_idx[end_idx > start_idx][1] - 1
  
  # Extract only that section
  direct_lines <- lines[start_idx:end_idx]
  direct_lines <- direct_lines[direct_lines != ""]
  
  # Get all hosts listed at the top
  hosts_line <- grep("^Hosts:", lines, value = TRUE)
  hosts <- unlist(strsplit(gsub("Hosts:|,", "", hosts_line), "\\s+"))
  hosts <- hosts[hosts != ""]
  
  # Initialize an empty matrix
  mat <- matrix(0, nrow = length(hosts), ncol = length(hosts),
                dimnames = list(hosts, hosts))
  
  # Parse the direct transmission section
  current_from <- NULL
  for (line in direct_lines) {
    # Detect "From host X to :"
    if (grepl("^From host", line)) {
      current_from <- sub("^From host (\\d+) to.*", "\\1", line)
    } else if (nchar(line) > 0 && !is.null(current_from)) {
      # Split comma-separated pairs like "52 0.035, 24 0.25"
      entries <- unlist(strsplit(line, ","))
      for (entry in entries) {
        entry <- trimws(entry)
        if (nchar(entry) > 0) {
          parts <- unlist(strsplit(entry, "\\s+"))
          if (length(parts) == 2) {
            to <- parts[1]
            val <- as.numeric(parts[2])
            if (!is.na(val) && to %in% hosts && current_from %in% hosts) {
              mat[current_from, to] <- val
            }
          }
        }
      }
    }
  }
  
  # Convert to data frame for easier viewing
  #direct_transmission_df <- as.data.frame(mat)
  
  # Save table (optional)
  #write.csv(direct_transmission_df, "direct_transmission_matrix.csv", row.names = TRUE)
  return(mat)
}

handlebreathformat <- function (filename) {
  lines <- readLines(filename)
  header <- strsplit(trimws(lines[1]), "\\s+")[[1]]
  header <- header[header != ""]
  col_ids <- sub("_.*", "", header)
  df <- read.table(filename, header = FALSE, skip = 1, stringsAsFactors = FALSE)
  df_ids <- sub("_.*", "", df[[1]])
  mat <- as.matrix(df[, -1])
  mode(mat) <- "numeric"
  rownames(mat) <- df_ids
  colnames(mat) <- col_ids
  mat <- t(mat) #transpose to get parent-child
  transmission_df <- as.data.frame(mat)
  #write.csv(transmission_df, "parsed_matrix.csv", row.names = TRUE)
  #print(head(transmission_df[, 1:10]))
  return(mat)
}

calib <- function (filename, epi, type) {
  
  # parse
  if (type == "scotti") { m <- handlescottiformat(paste0('./',filename)) }
  else if (type == "breath") { m <- handlebreathformat(paste0('./',filename)) }
  
  # Example: epi has columns Node, Parent
  # Ensure they’re numeric and comparable
  epi <- as.data.frame(epi)
  epi$Node <- as.numeric(epi$'Node ID')
  epi$Parent <- as.numeric(epi$Parent)
  hosts <- as.numeric(colnames(m))
  true_mat <- matrix(0, nrow = nrow(m), ncol = ncol(m),
                     dimnames = list(rownames(m), colnames(m)))
  
  for (i in seq_len(nrow(epi))) {
    if (!is.na(epi$Parent[i])) {
      parent <- as.character(epi$Parent[i])
      child  <- as.character(epi$Node[i])
      if (parent %in% rownames(true_mat) && child %in% colnames(true_mat)) {
        true_mat[parent, child] <- 1
      }
    }
  }
  # Flatten both matrices for comparison ---
  pred <- as.numeric(m)
  true <- as.numeric(true_mat)
  # Remove self-infections or NAs if any
  valid <- !is.na(pred) & row(m) != col(m)
  pred <- pred[valid]
  true <- true[valid]
  # Bin probabilities into 10% intervals ---
  bins <- cut(pred, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE, right = FALSE)
  # Compute calibration per bin ---
  library(dplyr)
  calibration <- data.frame(bin = bins, correct = true) |>
    group_by(bin) |>
    summarise(
      count = n(),
      mean_pred = mean(as.numeric(sub("\\[(.*),(.*)\\)", "\\1", bin))) + 0.05, # midpoint of bin
      proportion_correct = mean(correct),
      .groups = "drop"
    )
  # plot
  plot(calibration$mean_pred,calibration$proportion_correct,xlim=c(0,1),ylim=c(0,1),col="red",
       main = paste(filename,type))
  abline(b=1,a=0,color="gray")
  
}

sampling_experiment <- function () {
  for (ps in c(6,7,8,9)){
    load(paste0("sampling-experiment/ps-",ps,"/epi-data"))
    calib(paste0("ps/scotti-ps-0.",ps,"_network.txt"),epi,"scotti")
    calib(paste0("ps/breath-ps-0.",ps,"-matr"),epi,"breath")
  }
  load(paste0("sampling-experiment/ps-10/epi-data"))
  calib("ps/breath-ps-1-matr",epi,"breath")
  
}

sampling_experiment()


#calib("kappa-0-1-matr.txt", "breath")
#calib("ps/scotti-8-10-full_network.txt", "scotti")

