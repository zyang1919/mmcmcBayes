compare_dmrs <- function(rst1, rst2) {
  # Ensure the results have the necessary columns
  required_cols <- c("Chromosome", "Start_CpG", "End_CpG")
  
  if (!all(required_cols %in% colnames(rst1)) || !all(required_cols %in% colnames(rst2))) {
    stop("Both result datasets must contain the columns: Chromosome, Start_CpG, and End_CpG.")
  }
  
  # Convert to data frames
  rst1 <- as.data.frame(rst1)
  rst2 <- as.data.frame(rst2)
  
  # Function to safely convert CpG site positions to numeric while keeping original IDs
  clean_numeric <- function(x) {
    x <- gsub("[^0-9]", "", x)  # Remove non-numeric characters
    as.numeric(x)  # Convert to numeric
  }
  
  # Apply numeric conversion for calculations, while keeping original IDs
  rst1$Start_CpG_numeric <- clean_numeric(rst1$Start_CpG)
  rst1$End_CpG_numeric <- clean_numeric(rst1$End_CpG)
  rst2$Start_CpG_numeric <- clean_numeric(rst2$Start_CpG)
  rst2$End_CpG_numeric <- clean_numeric(rst2$End_CpG)
  
  # Remove rows where numeric conversion resulted in NA
  rst1 <- rst1[complete.cases(rst1[, c("Start_CpG_numeric", "End_CpG_numeric")]), ]
  rst2 <- rst2[complete.cases(rst2[, c("Start_CpG_numeric", "End_CpG_numeric")]), ]
  
  # Initialize a list to store overlaps
  overlap_results <- list()
  
  # Loop through each detected region in rst1 and check for overlap with rst2
  for (i in 1:nrow(rst1)) {
    chr1 <- rst1$Chromosome[i]
    start1_num <- rst1$Start_CpG_numeric[i]
    end1_num <- rst1$End_CpG_numeric[i]
    start1_orig <- rst1$Start_CpG[i]
    end1_orig <- rst1$End_CpG[i]
    
    # Find overlaps in rst2 where Chromosome matches
    overlaps <- rst2[rst2$Chromosome == chr1 & 
                       ((rst2$Start_CpG_numeric <= end1_num & rst2$End_CpG_numeric >= start1_num) |  # Partial overlap
                          (rst2$Start_CpG_numeric >= start1_num & rst2$End_CpG_numeric <= end1_num)), ]  # Complete overlap
    
    if (nrow(overlaps) > 0) {
      for (j in 1:nrow(overlaps)) {
        start2_num <- overlaps$Start_CpG_numeric[j]
        end2_num <- overlaps$End_CpG_numeric[j]
        start2_orig <- overlaps$Start_CpG[j]
        end2_orig <- overlaps$End_CpG[j]
        
        # Compute overlap size
        overlap_start <- max(start1_num, start2_num)
        overlap_end <- min(end1_num, end2_num)
        overlap_size <- max(0, overlap_end - overlap_start + 1)
        
        # Compute total region size for normalization
        total_region_size <- max(end1_num - start1_num + 1, end2_num - start2_num + 1)
        overlap_percentage <- (overlap_size / total_region_size) * 100
        
        # Store results without overlap_size
        overlap_results <- rbind(overlap_results, data.frame(
          Chromosome = chr1,
          Start_CpG_Method1 = start1_orig,
          End_CpG_Method1 = end1_orig,
          Start_CpG_Method2 = start2_orig,
          End_CpG_Method2 = end2_orig,
          Overlap_Percentage = round(overlap_percentage, 2)
        ))
      }
    }
  }
  
  # Convert result to data frame
  overlap_results <- as.data.frame(overlap_results)
  
  # Print results
  if (nrow(overlap_results) > 0) {
    return(overlap_results)
  } else {
    print("No overlapping DMRs detected between the two methods.")
  }
  
}
