library(readr)

rewrite_spruce_input <- function(filename){
  mtx <- read_lines(filename)
  mtx <- mtx[-c(1:2)]
  
  # new_filename <- strsplit(filename, '/')[[1]]
  # new_filename <- new_filename[length(new_filename)]
  new_filename <- strsplit(filename, '\\.')[[1]][1]
  new_filename <- paste(new_filename, "cleaned.csv", sep = "_")
  # print(new_filename)
  write_lines(mtx, new_filename)
  return(new_filename)
}

read_spruce_input <- function(filename, reads){
  # read input
  input <- read_delim(filename, 
                             " ", escape_double = FALSE, col_names = FALSE, 
                             trim_ws = TRUE)
  # sest up column names & get # of states
  n_col <- ncol(input)
  fixed_names <- c('sample_idx', 'sample_name', 'snv_idx', 'snv_name', 'prob_lb', 'prob_avg', 'prob_ub')
  if ((n_col-length(fixed_names))%%3 != 0){
    print('error')
    break
  }
  new_names <- c()
  nState <- (n_col-length(fixed_names))/3
  for (i in 1:nState){
    new_names <- c(new_names, paste(c('m_', 'p_', 'mu_'), i, sep = ""))
  }
  colnames(input) <- c(fixed_names, new_names)  
  
  # get # of SNVs
  samples <- unique(input$sample_idx)
  SNVs <- unique(input$snv_idx)
  nSample <- length(samples)
  nSNV <- length(SNVs)
  
  # analysis setup (arrays)
  # total_reads <- matrix(nrow = nSample, ncol = nSNV)
  input_VAF <- matrix(nrow = nSample, ncol = nSNV)
  # 
  # # calculation: vaf = var/total
  for (i in samples){
    sample_rows <- which(input$sample_idx == i)
    # sanity check
    if (length(sample_rows) != nSNV){
      print('sample/SNV# error')
      break
    }

    for (j in SNVs){
      # print(paste('j:',j))
      row <- which(input$snv_idx[sample_rows] == j)
      input_VAF[i+1,j+1] <- input$prob_avg[row]*reads
    }
  }
  # final output
  # return(input$prob_avg)
  return(input_VAF)
}

read_spruce_output <- function(filename){
  output <- read_csv(filename)
  output <- output[,-c(ncol(output))]
  output <- data.matrix(output)
  return(output)
}

read_spruce_readcount <- function(filename){
  reads <- read_csv(filename)
  reads <- as.numeric(reads)
  return(reads)
}

total_solution_likelihood <- function(reads, input_VAF, output_VAF){
  if (nrow(input_VAF) != nrow(output_VAF) | ncol(input_VAF) != ncol(output_VAF)){
    print('input & output incompatible')
    break
  }
  
  a <- reads*input_VAF
  
  probability <- 1
  
  for (i in 1:nrow(input_VAF)){
    for (j in 1:ncol(input_VAF)){
      print('new')
      print(input_VAF[i,j])
      print(output_VAF[i,j])
      
      prob <- pbinom(input_VAF[i,j], size = reads, prob = output_VAF[i,j])
      
      probability <- probability*prob # need to check this
    }
  }
  probability <- log(probability) # should i take the exponent instead? or negative log?
  return(probability)
}

calculate_solution_likelihood <- function(input_filename, output_filename, reads_filename){
  new_filename <- rewrite_spruce_input(input_filename)
  
  output <- read_spruce_output(output_filename)
  reads <- read_spruce_readcount(reads_filename)
  input <- read_spruce_input(new_filename, reads)
  
  return(total_solution_likelihood(reads, input, output))
}

iterate_over_solutions <- function(input_filename, base_sol_file, idxs){
  reads_filename <- paste(base_sol_file, "reads.csv", sep = "_")
  sol_likelihood <- matrix(nrow = length(idxs), ncol = 2)
  i <- 1
  for (idx in idxs){
    output_filename <- paste(base_sol_file, "sol", idx, "parsed.csv", sep = "_")
    print(output_filename)
    likelihood <- calculate_solution_likelihood(input_filename, output_filename, reads_filename)
    sol_likelihood[i,1] <- idx
    sol_likelihood[i,2] <- likelihood
    i <- i+1
  }
  return(sol_likelihood)
}

# input_filename <- "~/Spruce/data/sims_r0_m2_n5_c50.data"
# output_filename <- "~/Spruce/data/sims_r0_m2_n5_c50_parsed.csv"
# reads_filename <- "~/Spruce/data/sims_r0_m2_n5_c50_reads.csv"
# print(calculate_solution_likelihood(input_filename, output_filename, reads_filename))
nsol <- 205267
likelihoods <- iterate_over_solutions("~/Spruce/data/sims_r0_m2_n5_c50.data", "~/GitHub/spruce/build/data/sims_r0_m2_n5_c50", 0:(nsol-1))