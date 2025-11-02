library(dplyr)
library(tidyr)
library(purrr)

create_prob_mats <- function (folder) {
  # make probability matrices
  files <- list.files(folder, full.names = TRUE)
  for (f in files) {
    # Skip if not a -1.trees file
    if (!grepl("\\-1.seqs.trees$", f)) next
    fname <- basename(f)
    if (grepl("scotti", fname, ignore.case = TRUE)) {
      # Call the Python script located in the parent directory
      base_name <- sub("\\.trees$", "", f)
      args <- c(
        "Make_transmission_tree_alternative.py",
        "--input", f,
        "--outputF", base_name,
        "--burnin", "10"
      )
      message("Running SCOTTI Python command for: ", fname)
      system2("python", args = args, stdout = TRUE, stderr = TRUE)
    } else if (grepl("breath", fname, ignore.case = TRUE)) {
      #if(file.exists(normalizePath(matr_out_path))) return()
      partition <- sub(".*-(.*?)\\.trees$", "\\1", fname)
      out_path <- sub("\\.trees$", ".svg", normalizePath(f))
      matr_out_path <- sub("\\.trees$", ".matr", normalizePath(f))
      exe_path <- normalizePath("C:/Users/tobyh/OneDrive - The University of Auckland/_ CS 380/BEAST.v2.7.7.Windows/BEAST/AppLauncher.exe")
      args <- c("WIWVisualiser",
                "-tree", normalizePath(f),
                "-partition", partition,
                "-out", normalizePath(out_path),
                "-matrix", normalizePath(matr_out_path))
      message("Running: ", exe_path, " ", paste(args, collapse = " "))
      system2(exe_path, args)
    }
  }
}
handlescottiformat <- function (filename) {
  # Read the file
  lines <- readLines(filename)
  
  # Find section for "Probabilities direct transmission"
  start_idx <- grep("^Probabilities of direct transmittor to each sampled host", lines)
  end_idx <- grep("^Probabilities", lines)
  end_idx <- end_idx[end_idx > start_idx][1] - 1
  
  # Extract only that section
  direct_lines <- lines[start_idx:length(lines)]
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
    # Detect "To host X to :"
    if (grepl("^To host", line)) {
      current_from <- sub("^To host (\\d+) from.*", "\\1", line)
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
            if (is.na(val)) val <- parts[2]
            if (!is.na(val) && to %in% hosts && current_from %in% hosts) {
              mat[current_from, to] <- val
            }
          }
        }
      }
    }
  }
  
  mat <- mat[rownames(mat) != 'Unsampled', ]
  
  # Convert to data frame for easier viewing
  #direct_transmission_df <- as.data.frame(mat)
  
  # Save table (optional)
  #write.csv(direct_transmission_df, "direct_transmission_matrix.csv", row.names = TRUE)
  return(t(mat))
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

max_prob_accuracy <- function (epi,m){
  predicted_parent <- apply(m, 2, function(col) {
    if (all(is.na(col))) return(NA)
    if (max(col)<(1-sum(col))) return(NA) # probabiliy unsampled is highest
    rownames(m)[which.max(col)]
  })
  # Align with actual nodes
  node_ids <- colnames(m)
  pred_df <- data.frame(Node = as.numeric(node_ids),
                        PredParent = as.numeric(predicted_parent),
                        stringsAsFactors = FALSE)
  # Merge with true epidemic data
  compare_df <- merge(epi[, c("Node", "Parent")], pred_df, by = "Node", all.x = TRUE)
  # Compute accuracy (exclude cases where Parent is NA)
  accuracy1 <- mean(compare_df$Parent == compare_df$PredParent, na.rm = TRUE)

  ## get unsampled accuracy
  parent_col <- names(epi)[1]
  sampled_col <- names(epi)[6]
  # Filter only rows where PredParent is NA
  na_preds <- compare_df[is.na(compare_df$PredParent), ]
  # Join with epi to get unsampled info for the true Parent
  merged <- merge(
    na_preds,
    epi[, c(parent_col, sampled_col)],
    by.x = "Parent",
    by.y = parent_col,
    all.x = TRUE
  )
  # Compute accuracy: proportion of rows where the parent is unsampled
  accuracy2 <- mean(is.na(merged[[sampled_col]]))
  
  if (all(is.na(accuracy1))) return(accuracy2)
  if (all(is.na(accuracy2))) return(accuracy1)
  
  accuracy <- (accuracy1*nrow(na.omit(compare_df))+accuracy2*(nrow(compare_df)-nrow(na.omit(compare_df))))/nrow(compare_df)
  
  return(accuracy)
}
calibration <- function (epi,m) {
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
  # Remove NAs if any
  valid <- !is.na(pred) & row(m) != col(m)
  pred <- pred[valid]
  true <- true[valid]
  # Bin probabilities into 10% intervals ---
  bins <- cut(pred, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE, right = FALSE)
  # Compute calibration per bin ---

  calibration <- data.frame(bin = bins, correct = true) |>
    group_by(bin) |>
    summarise(
      count = n(),
      proportion_correct = mean(correct),
      .groups = "drop"
    ) |>
    complete(
      bin = unique(bins, na.rm=TRUE),
      fill = list(count = 0,proportion_correct = NA)
    ) |>
    drop_na(bin)
  
  # account for unsampled parent probabilities
  pred_probs <- 1 - colSums(m)
  child_ids <- colnames(m)
  true2 <- numeric(length(child_ids))
  for (i in seq_along(child_ids)) {
    child_id <- child_ids[i]
    child_row <- which(epi[[1]] == child_id)
    parent_id <- epi[child_row[1],2]
    parent_row <- which(epi[[1]] == parent_id)
    true2[i] <- if (is.na(epi[parent_row[1],6])) 1 else 0
  }
  bins2 <- cut(pred_probs, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE, right = FALSE)
  calibration2 <- data.frame(bin = bins2, correct = true2) |>
    group_by(bin) |>
    summarise(
      count = n(),
      proportion_correct = mean(correct),
      .groups = "drop"
    ) |>
    complete(
      bin = unique(bins, na.rm=TRUE),
      fill = list(count = 0,proportion_correct = NA)
    ) |>
    drop_na(bin)
  
  combined <- bind_rows(calibration, calibration2) %>%
    group_by(bin) %>%
    summarise(
      proportion_correct = {
        x <- proportion_correct
        w <- count
        valid <- !is.na(x) & !is.na(w)
        if (!any(valid)) NA_real_ else weighted.mean(x[valid], w[valid])
      },
      count = sum(count, na.rm = TRUE),
      .groups = "drop"
    )
  return(combined)
}
percentile_accuracy <- function(epi,m) {
  # Initialize results vector
  in_top_95 <- rep(NA, ncol(m))
  names(in_top_95) <- colnames(m)
  
  # Loop through each child (column)
  for (child_idx in 1:ncol(m)) {
    child_name <- colnames(m)[child_idx]
    parent_name <- epi[which(epi[, 1] == child_name),][[2]]
    #print(parent_name)
    true_parent_row <- which(epi[,1]==parent_name)
    true_parent_name <- ifelse(!is.na(epi[true_parent_row, 6]),parent_name,'not_sampled')
    
    # Get the column of predictions for this child
    col_probs <- m[, child_idx]
    prob_not_sampled <- 1 - sum(col_probs)
    all_probs <- c(col_probs, "not_sampled" = prob_not_sampled)
    sorted_probs <- sort(all_probs, decreasing = TRUE)
    
    # Calculate cumulative probabilities
    cum_probs <- cumsum(sorted_probs)

    # Find position of true parent in sorted list
    
    true_parent_position <- which(names(sorted_probs) == true_parent_name)
    #print(sorted_probs)
    
    # Check if true parent is in top 95%
    if (length(true_parent_position) > 0) {
      in_top_95[child_idx] <- cum_probs[true_parent_position] <= 0.95
    } else {
      # True parent not in predictions (shouldn't happen if data is consistent)
      in_top_95[child_idx] <- NA
    }
    
    #print(child_name)
    #print(cum_probs)
    #print(true_parent_name)

    
  }
  in_top_95 <- in_top_95[!is.na(in_top_95)]
  return (sum(in_top_95)/length(in_top_95))
}

analyse <- function (filename, epi, type, ps, v) {
  
  if (!file.exists(paste0('./',filename))){
    print(paste0("Skipped ",filename, ": file not found"))
    return()
  }
  
  # parse
  if (type == "scotti") { m <- handlescottiformat(paste0('./',filename)) }
  else if (type == "breath") { m <- handlebreathformat(paste0('./',filename)) }
  
  epi <- as.data.frame(epi)
  epi$Node <- as.numeric(epi$'Node ID')
  epi$Parent <- as.numeric(epi$Parent)

  accuracy <- max_prob_accuracy(epi,m)
  print(paste0("Accuracy (most probable parent correct): ", round(accuracy,3)))
  
  calibration <- calibration(epi,m)
  if (type == "scotti") { scotti_calib[2*(ps-50)+v] <- calibration }
  else if (type == "breath") { breath_calib[2*(ps-50)+v] <- calibration }
  #barplot(names.arg=paste0(calibration$bin,calibration$count),height=calibration$proportion_correct,col="blue",
  #     main = paste(filename,type),ylim=c(0,1),las=2)
  #abline(a=0,b=0.1)
  
  percentile <- percentile_accuracy(epi,m)
  print(paste0("Accuracy (predictioni in 95th percentile): ", round(percentile,3)))
}

combine_calib <- function(calib) {
  # Combine all into one data frame
  combined <- bind_rows(calib) %>%
    group_by(bin) %>%
    summarise(
      total_count = sum(count, na.rm = TRUE),
      # Weighted average of proportion_correct by count
      weighted_proportion_correct = ifelse(
        sum(count, na.rm = TRUE) > 0,
        sum(proportion_correct * count, na.rm = TRUE) / sum(count, na.rm = TRUE),
        NA_real_
      ),
      .groups = "drop"
    )
  return(combined)
}

plot_acc <- function(accs,title="Average Accuracy per ps (50–100)",ylab="Average accuracy") {
  ps_values <- 50:100
  split_acc <- split(accs, rep(ps_values, each = 2))
  print(split_acc)
  avg_acc <- sapply(split_acc, function(x) mean(x, na.rm = TRUE))
  plot(
    ps_values,
    avg_acc,
    type = "b",                     
    pch = 19,                       
    col = "blue",
    xlab = "Sampling probability (%)",
    ylab = ylab,
    ylim=c(0,1),
    main = title
  )
}
plot_accs <- function(accs1, accs2, 
                     title = "Average Accuracy per ps (50–100)", 
                     ylab = "Average accuracy",
                     labels = c("Accuracy Set 1", "Accuracy Set 2")) {
  
  ps_values <- 50:100
  
  # Split each accuracy vector by ps values (2 per ps)
  split_acc1 <- split(accs1, rep(ps_values, each = 2))
  split_acc2 <- split(accs2, rep(ps_values, each = 2))
  
  # Compute average accuracy for each ps, ignoring NAs
  avg_acc1 <- sapply(split_acc1, function(x) mean(x, na.rm = TRUE))
  avg_acc2 <- sapply(split_acc2, function(x) mean(x, na.rm = TRUE))
  
  # Plot first accuracy set
  plot(
    ps_values, avg_acc1,
    type = "b",
    pch = 19,
    col = "blue",
    ylim = c(0, 1),
    xlab = "Sampling probability (%)",
    ylab = ylab,
    main = title
  )
  
  # Add second accuracy set
  lines(ps_values, avg_acc2, type = "b", pch = 17, col = "red")
  
  # Add legend
  legend(
    "bottomright",
    legend = labels,
    col = c("blue", "red"),
    pch = c(19, 17),
    lty = 1,
    cex = 0.9
  )
}
plot_calib <- function(calib, title = "Calibration Plot") {
  # Calculate bin midpoints from labels like [0.1,0.2)
  mid_vals <- sapply(strsplit(gsub("\\[|\\)|\\]", "", calib$bin), ","), 
                     function(x) mean(as.numeric(x)))
  
  # Create bar plot without axis labels first
  bp <- barplot(
    height = calib$weighted_proportion_correct,
    col = "skyblue",
    ylim = c(0, 1),
    main = title,
    ylab = "Proportion of predictions correct",
    axes = TRUE,
    names.arg = rep("", nrow(calib)) # suppress for manual labeling
  )
  
  # Add red points at midpoints
  points(bp, mid_vals, pch = 19, col = "red")
  
  # Add nicer slanted labels manually
  text(
    x = bp,
    y = par("usr")[3] - 0.05, # below x-axis
    labels = paste0(calib$bin, "\n(n=", calib$total_count, ")"),
    srt = 80,                # rotate 45 degrees
    adj = 1,
    xpd = TRUE,
    cex = 0.7                # smaller font
  )
  
  # Add legend
  legend("topleft",
         legend = c("Observed", "Ideal"),
         fill = c("skyblue", NA),
         border = NA,
         pch = c(NA, 19),
         col = c("black", "red"),
         bty = "n")
}

calc_statistics_from_data <- function (){
  #create_prob_mats('breath-trees')
  #create_prob_mats('scotti-trees')
  breath_calib <<- c()
  scotti_calib <<- c()
  breath_acc <<- rep(NA, 102)
  scotti_acc <<- rep(NA, 102)
  breath_perc <<- rep(NA, 102)
  scotti_perc <<- rep(NA, 102)
  for (ps in 50:100){
    for (v in 1:2) {
      print(paste0("Analysing epidemic ",ps,"-",v))
      load(paste0("last2/ps-",ps,"-",v,"/epi-data-",ps/100,"-",v))
      
      epi <- as.data.frame(epi)
      epi$Node <- as.numeric(epi$'Node ID')
      epi$Parent <- as.numeric(epi$Parent)
      
      breath_filename <- paste0("breath-trees3/breath-",ps,"-",v,"-seqs.matr")
      scotti_filename <- paste0("scotti-trees3/scotti-",ps,"-",v,"_network.txt")
      
      if (file.exists(paste0('./',breath_filename))){
        m <- handlebreathformat(paste0('./',breath_filename))
        accuracy <- max_prob_accuracy(epi,m)
        #print(paste0("Accuracy (most probable parent correct): ", round(accuracy,3)))
        breath_acc[2*(ps-50)+v] <- accuracy
        calib <- calibration(epi,m)

        breath_calib <<- append(breath_calib, list(calib))
        
        percentile <- percentile_accuracy(epi,m)
        breath_perc[2*(ps-50)+v] <- percentile
        #print(paste0("Accuracy (predictioni in 95th percentile): ", round(percentile,3)))
      }
      else { print(paste0("Skipped ",breath_filename, ": file not found")) }
      
      if (file.exists(paste0('./',scotti_filename))){
        m <- handlescottiformat(paste0('./',scotti_filename))
        accuracy <- max_prob_accuracy(epi,m)
        #print(paste0("Accuracy (most probable parent correct): ", round(accuracy,3)))
        scotti_acc[2*(ps-50)+v] <- accuracy
        calib <- calibration(epi,m)
        #barplot(names.arg=paste0(calibration$bin,calibration$count),height=calibration$proportion_correct,col="blue",
        #     main = paste(filename,type),ylim=c(0,1),las=2)
        #abline(a=0,b=0.1)
        scotti_calib <<- append(scotti_calib, list(calib))
        
        percentile <- percentile_accuracy(epi,m)
        scotti_perc[2*(ps-50)+v] <- percentile
        #print(paste0("Accuracy (predictioni in 95th percentile): ", round(percentile,3)))
      }
      else { print(paste0("Skipped ",scotti_filename, ": file not found")) }
    }
  }
}

calc_statistics_from_data()

plot_acc(breath_acc,"breath accuracy (max likelihood)","accuracy")
plot_acc(scotti_acc,"scotti accuracy (max likelihood)","accuracy")
plot_acc(breath_perc,"breath accuracy (95 percentile)","accuracy")
plot_acc(scotti_perc,"scotti accuracy (95 percentile)","accuracy")

plot_accs(breath_acc,scotti_acc,title="Average Accuracy (truth is max-likelihood prediction)", labels=c('breath','scotti'))
plot_accs(breath_perc,scotti_perc,title="Average Accuracy (truth in 95th percentile of predictions)", labels=c('breath','scotti'))

breath_combined_calib <- combine_calib(breath_calib)
scotti_combined_calib <- combine_calib(scotti_calib)
plot_calib(breath_combined_calib, "BREATH Calibration")
plot_calib(scotti_combined_calib, "SCOTTI Calibration")