library(readr)

get_filenames <- function(m, n, c){
  core <- paste('sims_r0_m', m, '_n', n, '_c', c, sep = '')
  base <- paste("~/GitHub/spruce/build/data/", core, sep = '')
  input <- paste("~/GitHub/spruce/data/sims/n_5_noisy/", core, '.data', sep = '')
  # "~/GitHub/spruce/data/sims/n_5_noisy/sims_r0_m5_n5_c50.data", "~/GitHub/spruce/build/data/sims_r0_m5_n5_c50"
  list(input = input, base = base)
}

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
  input <- read_delim(filename, " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  # print(input)
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
  input_vars <- matrix(nrow = nSample, ncol = nSNV)
  # 
  # # calculation: vaf = var/total
  for (i in samples){
    # print(i)
    sample_rows <- which(input$sample_idx == i)
    # print(sample_rows)
    # sanity check
    if (length(sample_rows) != nSNV){
      print('sample/SNV# error')
      break
    }

    for (j in SNVs){
      # print(paste('j:',j))
      row <- which(input$snv_idx[sample_rows] == j)+(i*nSNV)
      # print(row)
      input_vars[i+1,j+1] <- input$prob_avg[row]#*reads
    }
  }
  # final output
  # return(input$prob_avg)
  return(input_vars)
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
  
  # a <- reads*input_VAF
  # print(input_VAF)
  # print(reads)
  # print(output_VAF)
  
  probability <- 1
  # ones <- 0
  # zeros <- 0
  for (i in 1:nrow(input_VAF)){
    for (j in 1:ncol(input_VAF)){
      # if (output_VAF[i,j] == 0){next}
      # trial <- pbinom(input_VAF[i,j], size = reads, prob = output_VAF[i,j])
      successes <- floor(reads*input_VAF[i,j])
      # print(paste(successes, reads, output_VAF[i,j]))
      # prob <- sum(dbinom(floor(successes-(0.1*reads)):floor(successes+(0.1*reads)), size = reads, prob = output_VAF[i,j]))
      # for (s in (successes-(0.1*reads)):(successes+(0.1*reads))){
      #   prob <- prob + dbinom(s, size = reads, prob = output_VAF[i,j])
      # }
      prob <- pbinom(successes, size = reads, prob = output_VAF[i,j])#, lower.tail = TRUE, log.p = FALSE)
      # print(prob)
      # if (prob == 0){print('zero'); zeros <- zeros + 1}
      # if (prob == 1){print('one'); ones <- ones + 1}
      # a = input_VAF[i,j]
      # B = 1 - a
      # s = a*B/((a+B)*(a+B)*(a+B+1))
      probability <- probability*prob # need to check this
    }
  }
  probability <- probability # should i take the exponent instead? or negative log?
  return(probability)
}

test_solution_likelihood <- function(reads, input_VAF, output_VAF){
  if (nrow(input_VAF) != nrow(output_VAF) | ncol(input_VAF) != ncol(output_VAF)){
    print('input & output incompatible')
    break
  }
  
  prob_diff <- input_VAF - output_VAF
  # print(prob_diff)
  # print(nrow(prob_diff))
  # print(dist(prob_diff[1,]))
  prob <- sum(prob_diff) #sum(vapply((1:nrow(prob_diff)), function(x){sqrt(sum(prob_diff[x,]^2))}, c(1)))
  return(prob)
}

full_likelihood_output <- function(reads, input_VAF, output_VAF){
  if (nrow(input_VAF) != nrow(output_VAF) | ncol(input_VAF) != ncol(output_VAF)){
    print('input & output incompatible')
    break
  }
  
  # a <- reads*input_VAF
  # print(input_VAF)
  # print(reads)
  # print(output_VAF)
  
  # probability <- 0
  likelihood_vars <- matrix(ncol = length(input_VAF), nrow = 1)
  prob_vars <- matrix(ncol = length(input_VAF), nrow = 1)
  
  idx <- 1
  
  for (i in 1:nrow(input_VAF)){
    for (j in 1:ncol(input_VAF)){
      successes <- floor(reads*input_VAF[i,j])
      # prob <- sum(dbinom(floor(successes-(0.05*reads)):floor(successes+(0.05*reads)), size = reads, prob = output_VAF[i,j]))
      prob <- pbinom(floor(successes), size = reads, prob = output_VAF[i,j])
      
      likelihood_vars[idx] <- prob
      prob_vars[idx] <- input_VAF[i,j]-output_VAF[i,j]
      idx <- idx + 1
    }
  }
  
  soln_vars <- matrix(ncol = length(input_VAF), nrow = 1)
  soln_vars[1:length(input_VAF)] <- likelihood_vars
  soln_vars[length(input_VAF):(2*length(input_VAF))] <- prob_vars
  
  return(soln_vars)
}

calculate_solution_likelihood <- function(input_filename, output_filename, reads_filename){
  new_filename <- rewrite_spruce_input(input_filename)
  
  output <- read_spruce_output(output_filename)
  reads <- read_spruce_readcount(reads_filename)
  input <- read_spruce_input(new_filename, reads)
  
  return(total_solution_likelihood(reads, input, output)) #
}

calculate_solution_distance <- function(input_filename, output_filename, reads_filename){
  new_filename <- rewrite_spruce_input(input_filename)

  output <- read_spruce_output(output_filename)
  reads <- read_spruce_readcount(reads_filename)
  input <- read_spruce_input(new_filename, reads)

  return(test_solution_likelihood(reads, input, output)) #
}

calculate_full_solution_likelihood <- function(input_filename, output_filename, reads_filename){
  new_filename <- rewrite_spruce_input(input_filename)
  
  output <- read_spruce_output(output_filename)
  reads <- read_spruce_readcount(reads_filename)
  input <- read_spruce_input(new_filename, reads)
  
  return(full_likelihood_output(reads, input, output))
}

iterate_over_solutions <- function(input_filename, base_sol_file, idxs){
  reads_filename <- paste(base_sol_file, "reads.csv", sep = "_")
  sol_likelihood <- matrix(nrow = length(idxs), ncol = 3)
  i <- 1
  for (idx in idxs){
    output_filename <- paste(base_sol_file, "sol", idx, "parsed.csv", sep = "_")
    # print(output_filename)
    likelihood <- calculate_solution_likelihood(input_filename, output_filename, reads_filename)
    distance <- calculate_solution_distance(input_filename, output_filename, reads_filename)
    sol_likelihood[i,1] <- idx
    sol_likelihood[i,2] <- likelihood
    sol_likelihood[i,3] <- distance
    i <- i+1
  }
  return(sol_likelihood)
}

iterate_soln_vars <- function(input_filename, base_sol_file, idxs){
  reads_filename <- paste(base_sol_file, "reads.csv", sep = "_")
  sol_likelihood <- matrix(nrow = length(idxs), ncol = 21)
  i <- 1
  for (idx in idxs){
    output_filename <- paste(base_sol_file, "sol", idx, "parsed.csv", sep = "_")
    # print(output_filename)
    vars <- calculate_full_solution_likelihood(input_filename, output_filename, reads_filename)
    sol_likelihood[i,1] <- idx
    sol_likelihood[i,(2:21)] <- vars
    i <- i+1
  }
  return(sol_likelihood)
}

nsol <- matrix(data = c(205267, 54044, 9170, 6227, 791, 2832, 1572, 292, 92, 12, 40, 33, 2, 2, 2), nrow = 5, ncol = 3)

recall_m2 <- c()
recall_m5 <- c()
recall_m10 <- c()
recall_m2[[1]] <- sims_r0_m2_n5_c50_recall$X1
recall_m2[[2]] <- sims_r0_m2_n5_c100_recall$X1
recall_m2[[3]] <- sims_r0_m2_n5_c500_recall$X1
recall_m2[[4]] <- sims_r0_m2_n5_c1000_recall$X1
recall_m2[[5]] <- sims_r0_m2_n5_c10000_recall$X1
recall_m5[[1]] <- sims_r0_m5_n5_c50_recall$X1
recall_m5[[2]] <- sims_r0_m5_n5_c100_recall$X1
recall_m5[[3]] <- sims_r0_m5_n5_c500_recall$X1
recall_m5[[4]] <- sims_r0_m5_n5_c1000_recall$X1
recall_m5[[5]] <- sims_r0_m5_n5_c10000_recall$X1
recall_m10[[1]] <- sims_r0_m10_n5_c50_recall$X1
recall_m10[[2]] <- sims_r0_m10_n5_c100_recall$X1
recall_m10[[3]] <- sims_r0_m10_n5_c500_recall$X1
recall_m10[[4]] <- sims_r0_m10_n5_c1000_recall$X1
recall_m10[[5]] <- sims_r0_m10_n5_c10000_recall$X1

recall <- cbind(recall_m2, recall_m5, recall_m10)

m_list <- c(2, 5,10)
c_list <- c(50, 100, 500, 1000, 10000)

likelihood_thresh <- matrix(data = -1, nrow = 3, ncol = 5)
percent_captured <- matrix(data = -1, nrow = 3, ncol = 5)

for (m in 3){
  for (c in 1:5){
    m_val <- m_list[m]
    c_val <- c_list[c]
    nsol_val <- nsol[c,m]
    filenames <- get_filenames(m_val,5,c_val)
    likelihood <- iterate_over_solutions(filenames$input, filenames$base, 0:(nsol_val-1))
    recall_vals <- recall[[c,m]]
    max <- which(recall_vals == 1)
    log_likelihoods <- log(likelihood[,2])
    thresh <- log_likelihoods[2]
    capture <- length(which(log_likelihoods >= log_likelihoods[max]))/length(log_likelihoods)
    
    likelihood_thresh[m,c] <- thresh
    percent_captured[m,c] <- capture
  }
}

# plot(log(m5_n5_likelihood_1000[,2]), sims_r0_m5_n5_c1000_recall$X1)
m <- 2
c <- 2
m_val <- m_list[m]
c_val <- c_list[c]
nsol_val <- nsol[c,m]
filenames <- get_filenames(m_val,5,c_val)
likelihood <- iterate_over_solutions(filenames$input, filenames$base, 0:(nsol_val-1))
recall_vals <- recall[[c,m]]
max <- which(recall_vals == 1)
log_likelihoods <- log(likelihood[,2])
thresh <- log_likelihoods[2]
capture <- length(which(log_likelihoods >= log_likelihoods[max]))/length(log_likelihoods)
