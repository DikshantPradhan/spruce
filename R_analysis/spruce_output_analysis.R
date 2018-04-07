# WORKSPACE SETUP
library(readr)
sample_input <- read_delim("~/Spruce/sample_input.csv", 
                           " ", escape_double = FALSE, col_names = FALSE, 
                           trim_ws = TRUE)
sample_ouput_mix <- read_delim("~/Spruce/sample_ouput_mix.csv", 
                               " ", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE)
sample_output_mix <- data.matrix(sample_ouput_mix)
sample_output_tree <- read_delim("~/Spruce/sample_output_tree.csv", 
                                 " ", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE)
input <- data.matrix(sample_input)
input <- input[,c(1,3:19)]

# SETUP MIXTURE
mixture <- sample_output_mix

nSample <- 2
nSNV <- ncol(mixture)
nState <- nrow(mixture)/nSample

mixture_array <- array(data = 0, dim = c(nState, nSNV, nSample))

for (i in 1:nState){
  for (j in 1:nSNV){
    for (k in 1:nSample){
      mixture_array[i,j,k] <- mixture[2*i-(2-k),j]
    }
  }
}

# SETUP TREE
tree <- sample_output_tree

state_array <- array(data = NA, dim = c(nState, nSNV, 3))
for (i in 1:nState){
  for (j in 1:nSNV){
    # strsplit
    nums <- strsplit(tree[[2*j,i]], split = c('[()]'))[[1]][2]
    nums <- as.numeric(strsplit(nums, split = c(','))[[1]])
    # print(nums)
    for (k in 1:3){
      # print(paste(i, j , k))
      # print(nums[k])
      state_array[i,j,k] <- nums[k]
    }
  }
}

snv_mtx <- matrix(data = NA, nrow = nState, ncol = nSNV)
for (i in 1:nState){
  for (j in 1:nSNV){
    snv_mtx[i,j] <- state_array[i,j,3]/(state_array[i,j,1] + state_array[i,j,2])
  }
}

# SETUP INPUT

nCopyProp <- 4
mCopy <- array(dim = c(nSNV, nCopyProp, nSample))
pCopy <- array(dim = c(nSNV, nCopyProp, nSample))

start_col <- 7
for (i in 1:nSNV){
  for (j in seq(from = start_col, to = ncol(input), by = 3)){
    mCopy[i, (j-(nCopyProp-1))/3, 1] <- input[i,j]
    pCopy[i, (j-(nCopyProp-1))/3, 1] <- input[i,j+1]
    #copyNumProp[i, (j-(nCopyProp-1))/3, 1] <- mCopy+pCopy
  }
}
for (i in (1+nSNV):(2*nSNV)){
  for (j in seq(from = start_col, to = ncol(input), by = 3)){
    mCopy[i-nSNV, (j-(nCopyProp-1))/3, 2] <- input[i,j]
    pCopy[i-nSNV, (j-(nCopyProp-1))/3, 2] <- input[i,j+1]
    # copyNumProp[i-nSNV, (j-(nCopyProp-1))/3, 2] <- mCopy+pCopy
  }
}

copyPropFreq <- array(dim = c(nSNV, nCopyProp, nSample))
for (i in 1:nSNV){
  for (j in seq(from = (start_col+2), to = ncol(input), by = 3)){
    copyPropFreq[i, (j-(nCopyProp+1))/3, 1] <- input[i,j]
  }
}
for (i in (1+nSNV):(2*nSNV)){
  for (j in seq(from = (start_col+2), to = ncol(input), by = 3)){
    copyPropFreq[i-nSNV, (j-(nCopyProp+1))/3, 2] <- input[i,j]
  }
}

# GET SNV OUTPUT FREQUENCY

VAF_1 <- mixture_array[,,1]*snv_mtx
output_VAF_1 <- colSums(VAF_1)
VAF_2 <- mixture_array[,,2]*snv_mtx
output_VAF_2 <- colSums(VAF_2)

# GET SNV INPUT FREQUENCY

estVarCount_1 <- input[1:5,5]
estVarCount_2 <- input[6:10,5]

copyNumProp <- array(dim = c(nSNV, nCopyProp, nSample)) # need to finish using this

refProportions <- copyNumProp*copyPropFreq
input_refCount_1 <- rowSums(varProportions[,,1])
input_refCount_2 <- rowSums(varProportions[,,2])

betaVAF_1 <- estVarCount_1/(estVarCount_1+input_refCount_1)
betaVAF_2 <- estVarCount_2/(estVarCount_2+input_refCount_2)
