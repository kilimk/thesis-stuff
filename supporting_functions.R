# 21.07.2015
# Author: David Toubiana
# some important supporting functions for network analysis


get.igraph <- function(data){
  # function taking in the correlation adjacency matrix and returning the network as an igraph object for subsequent analysis
  # load library
  library(igraph)
  sig.netw.zero <- data
  # sig.netw.zero[which(sig.netw.zero!=0)]  <- 1
  # convert diagonal to 0's
  for(i in 1:nrow(sig.netw.zero)){sig.netw.zero[i,i] <- 0 }
  # make network graph
  gr <- graph.adjacency(sig.netw.zero,mode="lower",diag=F,add.colnames=T,weighted=T)
  gr$names <- row.names(sig.netw.zero)
  # find nonconnected nodes
  zero <- which(degree(gr)==0)
  # delete nonconnected nodes from graph
  gr <- delete.vertices(gr,zero)
  if(length(zero)!=0) gr$names <- gr$names[-zero]
  return(gr) # undirected
} # end function
################################################################################################################################################################
################################################################################################################################################################
triangle <- function(data){
  # get only the lower triangle of a symmetric matrix and store in a vector
    # if it is single celled
  if(dim(data)[1]<2){
      return(data[1,1])
  }
    ret.vec <- NULL
    for(i in 1:(ncol(data)-1)){
       # print(i)
        ret.vec <- c(ret.vec,data[(i+1):nrow(data),i])
    }
    nas <- which(is.na(ret.vec))
    if(length(nas)!=0)
        ret.vec <- ret.vec[-nas]
    return(ret.vec)  
#   first.col <- 2
#   
#   for(k in 1:nrow(data)){
#     for(i in first.col:ncol(data)){
#       ret.vec <- rbind(ret.vec,data[k,i])
#     }# end nested loop
#     first.col <- first.col+1
#     if(first.col==ncol(data)+1){break}
#   } # end for loop
#   
} # end function
################################################################################################################################################################
################################################################################################################################################################
# FDR function showing p values, the frequency, and the actual number of false positives
# within the significantly identified instances

fdr.fun <- function(p.value){
  # load library
  library(qvalue)
  q.vec <- qvalue(p.value) # determining qvalue
  # preparing return matrix
  fdr.matrix <- matrix(ncol=5,nrow=0)
  colnames(fdr.matrix) <- c("pvalue","qvalue","frequency","False positives","%")
  
  seq <- seq(0.001,0.05,0.001)  # for fdr
  
  for(i in seq){
      #print(i)
    fdr.matrix <- rbind(fdr.matrix,c(i,max(q.vec$qvalues[q.vec$pvalues <= i]),length(which(q.vec$pvalues<=i)),0))
    fdr.matrix[nrow(fdr.matrix),4] <- fdr.matrix[nrow(fdr.matrix),2]*fdr.matrix[nrow(fdr.matrix),3]
    fdr.matrix[nrow(fdr.matrix),5] <- fdr.matrix[nrow(fdr.matrix),4]/fdr.matrix[nrow(fdr.matrix),3]
  }
  
  return(fdr.matrix)
  
} # end function
################################################################################################################################################################
################################################################################################################################################################
# removing all columns with NAs only

remove_zero_col <- function(matrix){
  x  <- NULL
  for(i in 1:ncol(matrix)){
    # print(colnames(matrix)[i])
    for(j in 1:nrow(matrix)){
      if(matrix[j,i]!=0){
        break;    	# break the loop if the entry is not 0
      }
      # this if will only be entered if the last entry was 0
      if(j==nrow(matrix)){
        x <- rbind(x,i)
      }
    } # end nested for loop
  } # end for loop
  # remove all columns with NAs only
  
  if(length(x)>0){
    matrix <- matrix[,-x]
  }
  
  matrix   # return variable
  
} # end function
################################################################################################################################################################
################################################################################################################################################################

# removing all rows with NAs only

remove_zero_row <- function(matrix){
  x  <- NULL
  for(i in 1:nrow(matrix)){
    for(j in 1:ncol(matrix)){
      if(matrix[i,j]!=0){
        break;			# break the loop if the entry is not 0
      }
      # this if will only be entered if the last entry was 0
      if(j==ncol(matrix)){
        x <- rbind(x,i)
      }
    } # end nested for loop
  } # end for loop
  # remove all rows with NAs only
  if(length(x)>0){
    matrix <- matrix[-x,]
  }
  matrix   # return variable
} # end function
################################################################################################################################################################
################################################################################################################################################################
# removing all columns with NAs only

remove_na_col <- function(matrix){
  x  <- NULL
  for(i in 1:ncol(matrix)){
    for(j in 1:nrow(matrix)){
      if(is.na(matrix[j,i])==F){
        break;			# break the loop if the entry is not NA
      }
      # this if will only be entered if the last entry was NA
      if(j==nrow(matrix)){
        x <- rbind(x,i)
      }
    } # end nested for loop
  } # end for loop
  # remove all columns with NAs only
  
  if(length(x)>0){
    matrix <- matrix[,-x]
  }
  
  matrix   # return variable
  
} # end function
################################################################################################################################################################
################################################################################################################################################################

# removing all rows with NAs only

remove_na_row <- function(matrix){
  x  <- NULL
  for(i in 1:nrow(matrix)){
    for(j in 1:ncol(matrix)){
      if(is.na(matrix[i,j])==F){
        break;			# break the loop if the entry is not NA
      }
      # this if will only be entered if the last entry was NA
      if(j==ncol(matrix)){
        x <- rbind(x,i)
      }
    } # end nested for loop
  } # end for loop
  # remove all rows with NAs only
  if(length(x)>0){
    matrix <- matrix[-x,]
  }
  matrix   # return variable
} # end function
################################################################################################################################################################
################################################################################################################################################################
sync.rows <- function(dataset.1,dataset.2){
  # function that takes in two datasets and synchronizes
  # their rows for the same name, so that at the end both
  # datasets contain the same rownames
  # prerequisite: a rowname can be present only once for each dataset
  
  # create output variables
  output.dataset.1 <- matrix(ncol=ncol(dataset.1),nrow=0)
  output.dataset.2 <- matrix(ncol=ncol(dataset.2),nrow=0)
  
  # set the colnames to input parameters
  colnames(output.dataset.1) <- colnames(dataset.1) 
  colnames(output.dataset.2) <- colnames(dataset.2) 
  
  # now go over all rows of dataset 1 and compare to dataset 2
  for(i in 1:nrow(dataset.1)){
    # read entries
    x <- which(row.names(dataset.2)==row.names(dataset.1)[i])
    
    # see if entry exists
    if(length(x)>0) { # it exists in second dataset
      output.dataset.1 <- rbind(output.dataset.1,dataset.1[i,])
      output.dataset.2 <- rbind(output.dataset.2,dataset.2[x,])
    } # end if
  } # end for loop
  
  # return variable
  y <- cbind(output.dataset.1,output.dataset.2)
  y
} # end function
################################################################################################################################################################
################################################################################################################################################################
sync.cols <- function(dataset.1,dataset.2){
  # function that takes in two datasets and synchronizes
  # their cols for the same name, so that at the end both
  # datasets contain the same rownames
  # prerequisite: a colname can be present only once for each dataset
  
  # create output variables
  output.dataset.1 <- matrix(nrow=nrow(dataset.1),ncol=0)
  output.dataset.2 <- matrix(nrow=nrow(dataset.2),ncol=0)
  
  # set the rownames to input parameters
  rownames(output.dataset.1) <- rownames(dataset.1) 
  rownames(output.dataset.2) <- rownames(dataset.2) 
  
  # now go over all cols of dataset 1 and compare to dataset 2
  for(i in 1:ncol(dataset.1)){
    # read entries
    x <- which(colnames(dataset.2)==colnames(dataset.1)[i])
    
    # see if entry exists
    if(length(x)>0) { # it exists in second dataset
      output.dataset.1 <- cbind(output.dataset.1,dataset.1[,i])
      colnames(output.dataset.1)[ncol(output.dataset.1)] <- colnames(dataset.1)[i]
      output.dataset.2 <- cbind(output.dataset.2,dataset.2[,x])
    } # end if
  } # end for loop
  
  # return variable
  y <- rbind(output.dataset.1,output.dataset.2)
  y
} # end function
################################################################################################################################################################
################################################################################################################################################################
# function converting all entries in a matrix to numeric instead of character

# matrix parameter is the matrix to convert
convert_to_num <- function(matrix){
  
  # empty matrix
  conv_mat <- matrix(nrow=nrow(matrix),ncol=ncol(matrix),dimnames=list(rownames(matrix),colnames(matrix)))
  
  # fill matrix con_mat with values from matrix
  for(i in 1:nrow(matrix)){
    for(j in 1:ncol(matrix)){
      conv_mat[i,j] <- as.numeric(as.character(matrix[i,j]))
    } # end nested for loop
    
  } # end for loop
  
  # return variable
  conv_mat
} # end function
################################################################################################################################################################
################################################################################################################################################################

# matrix parameter is the matrix to convert
convert.to.character<- function(matrix){
  
  # empty matrix
  conv_mat <- matrix(nrow=nrow(matrix),ncol=ncol(matrix),dimnames=list(rownames(matrix),colnames(matrix)))
  
  # fill matrix con_mat with values from matrix
  for(i in 1:nrow(matrix)){
    for(j in 1:ncol(matrix)){
      conv_mat[i,j] <- as.character(matrix[i,j])
    } # end nested for loop
    
  } # end for loop
  
  # return variable
  conv_mat
} # end function
################################################################################################################################################################
################################################################################################################################################################

# collection of descriptive functions: mean, sd, se, var
# e.g.
# calc_mean receives a data matrix with rows as the samples and columns as variables
# the number of replicates can vary. the function of this function to detect the different number of
# of rows for samples and calculate automatically the mean of each sample for each variable

# the replicates must be grouped together

# it works the same for sd, se, and var
calc_mean <- function(input_matrix){
  # function variables
  nrow <- nrow(input_matrix); ncol <- ncol(input_matrix)
  r_names <- rownames(input_matrix)
  
  sample_name <- r_names[1]  # variable to contain the sample name - initialize with the first row
  row_mean_vector <- NULL
  sample_name_rows <- NULL
  matrix_output_row_names <- sample_name # initialize with first row name
  result_matrix <- NULL # containing all the averages
  # index for matrix
  index  <- 0
  # find all samples with the same name and calculate the averages for each variable
  # in an iterative mode
  while(TRUE){
    # find all identical sample name rows
    sample_name_rows <- which(r_names==sample_name)
    
    # now calculate the means for each variable
    for(i in 1:ncol){
      row_mean_vector <- cbind(row_mean_vector,mean(input_matrix[sample_name_rows,i],na.rm=T))
    }# end for loop 
    
    # append row vector with averages to result_matrix
    result_matrix  <- rbind(result_matrix,row_mean_vector)
    index  <- index + length(sample_name_rows)
    #print(index)
    #print(nrow)
    if(index >= nrow){
      break;   # break the loop if reached the end of input matrix
    } # end if
    else{
      sample_name <- r_names[index+1]  # reinitialize sample_name variable with next sample
      # append sample_name
      matrix_output_row_names <- rbind(matrix_output_row_names,sample_name)
      row_mean_vector <- NULL
      
    } # end else
    
  } # end while loop
  # return variable
  output <- matrix(result_matrix,nrow=nrow(result_matrix),ncol=ncol(result_matrix),byrow=F,
                   dimnames <- list(matrix_output_row_names,colnames(input_matrix)))
  output
  
}# end function

################################################################################################################################################################
################################################################################################################################################################

# the same functionality for standard deviation

calc_sd <- function(input_matrix){
  # function variables
  nrow <- nrow(input_matrix); ncol <- ncol(input_matrix)
  r_names <- rownames(input_matrix)
  
  sample_name <- r_names[1]  # variable to contain the sample name - initialize with the first row
  row_mean_vector <- NULL
  sample_name_rows <- NULL
  matrix_output_row_names <- sample_name # initialize with first row name
  result_matrix <- NULL # containing all the averages
  # index for matrix
  index  <- 0
  # find all samples with the same name and calculate the averages for each variable
  # in an iterative mode
  while(TRUE){
    # find all identical sample name rows
    sample_name_rows <- which(r_names==sample_name)
    
    # now calculate the means for each variable
    for(i in 1:ncol){
      row_mean_vector <- cbind(row_mean_vector,sd(input_matrix[sample_name_rows,i],na.rm=T))
    }# end for loop 
    
    # append row vector with averages to result_matrix
    result_matrix  <- rbind(result_matrix,row_mean_vector)
    index  <- index + length(sample_name_rows)
    #print(index)
    #print(nrow)
    if(index >= nrow){
      break;   # break the loop if reached the end of input matrix
    } # end if
    else{
      sample_name <- r_names[index+1]  # reinitialize sample_name variable with next sample
      # append sample_name
      matrix_output_row_names <- rbind(matrix_output_row_names,sample_name)
      row_mean_vector <- NULL
      
    } # end else
    
  } # end while loop
  # return variable
  output <- matrix(result_matrix,nrow=nrow(result_matrix),ncol=ncol(result_matrix),byrow=F,
                   dimnames <- list(matrix_output_row_names,colnames(input_matrix)))
  output
  
}# end function
################################################################################################################################################################
################################################################################################################################################################

# the same functionality for variance
calc_var <- function(input_matrix){
  # function variables
  nrow <- nrow(input_matrix); ncol <- ncol(input_matrix)
  r_names <- rownames(input_matrix)
  
  sample_name <- r_names[1]  # variable to contain the sample name - initialize with the first row
  row_mean_vector <- NULL
  sample_name_rows <- NULL
  matrix_output_row_names <- sample_name # initialize with first row name
  result_matrix <- NULL # containing all the averages
  # index for matrix
  index  <- 0
  # find all samples with the same name and calculate the averages for each variable
  # in an iterative mode
  while(TRUE){
    # find all identical sample name rows
    sample_name_rows <- which(r_names==sample_name)
    
    # now calculate the means for each variable
    for(i in 1:ncol){
      row_mean_vector <- cbind(row_mean_vector,var(input_matrix[sample_name_rows,i],na.rm=T))
    }# end for loop 
    
    # append row vector with averages to result_matrix
    result_matrix  <- rbind(result_matrix,row_mean_vector)
    index  <- index + length(sample_name_rows)
    #print(index)
    #print(nrow)
    if(index >= nrow){
      break;   # break the loop if reached the end of input matrix
    } # end if
    else{
      sample_name <- r_names[index+1]  # reinitialize sample_name variable with next sample
      # append sample_name
      matrix_output_row_names <- rbind(matrix_output_row_names,sample_name)
      row_mean_vector <- NULL
      
    } # end else
    
  } # end while loop
  # return variable
  output <- matrix(result_matrix,nrow=nrow(result_matrix),ncol=ncol(result_matrix),byrow=F,
                   dimnames <- list(matrix_output_row_names,colnames(input_matrix)))
  output
  
}# end function
################################################################################################################################################################
################################################################################################################################################################

# the same functionality for median

calc_median <- function(input_matrix){
  # function variables
  nrow <- nrow(input_matrix); ncol <- ncol(input_matrix)
  r_names <- rownames(input_matrix)
  
  sample_name <- r_names[1]  # variable to contain the sample name - initialize with the first row
  row_mean_vector <- NULL
  sample_name_rows <- NULL
  matrix_output_row_names <- sample_name # initialize with first row name
  result_matrix <- NULL # containing all the averages
  # index for matrix
  index  <- 0
  # find all samples with the same name and calculate the averages for each variable
  # in an iterative mode
  while(TRUE){
    # find all identical sample name rows
    sample_name_rows <- which(r_names==sample_name)
    
    # now calculate the means for each variable
    for(i in 1:ncol){
      row_mean_vector <- cbind(row_mean_vector,median(input_matrix[sample_name_rows,i],na.rm=T))
    }# end for loop 
    
    # append row vector with averages to result_matrix
    result_matrix  <- rbind(result_matrix,row_mean_vector)
    index  <- index + length(sample_name_rows)
    #print(index)
    #print(nrow)
    if(index >= nrow){
      break;   # break the loop if reached the end of input matrix
    } # end if
    else{
      sample_name <- r_names[index+1]  # reinitialize sample_name variable with next sample
      # append sample_name
      matrix_output_row_names <- rbind(matrix_output_row_names,sample_name)
      row_mean_vector <- NULL
      
    } # end else
    
  } # end while loop
  # return variable
  output <- matrix(result_matrix,nrow=nrow(result_matrix),ncol=ncol(result_matrix),byrow=F,
                   dimnames <- list(matrix_output_row_names,colnames(input_matrix)))
  output
  
}# end function
################################################################################################################################################################
################################################################################################################################################################

# the same functionality for standard error
calc_se <- function(input_matrix){
  # function variables
  nrow <- nrow(input_matrix); ncol <- ncol(input_matrix)
  r_names <- rownames(input_matrix)
  
  sample_name <- r_names[1]  # variable to contain the sample name - initialize with the first row
  row_mean_vector <- NULL
  sample_name_rows <- NULL
  matrix_output_row_names <- sample_name # initialize with first row name
  result_matrix <- NULL # containing all the averages
  # index for matrix
  index  <- 0
  # find all samples with the same name and calculate the averages for each variable
  # in an iterative mode
  while(TRUE){
    # find all identical sample name rows
    sample_name_rows <- which(r_names==sample_name)
    
    # now calculate the means for each variable
    for(i in 1:ncol){
      row_mean_vector <- cbind(row_mean_vector,sd(input_matrix[sample_name_rows,i],na.rm=T)/sqrt(length(sample_name_rows)))
    }# end for loop 
    
    # append row vector with averages to result_matrix
    result_matrix  <- rbind(result_matrix,row_mean_vector)
    index  <- index + length(sample_name_rows)
    #print(index)
    #print(nrow)
    if(index >= nrow){
      break;   # break the loop if reached the end of input matrix
    } # end if
    else{
      sample_name <- r_names[index+1]  # reinitialize sample_name variable with next sample
      # append sample_name
      matrix_output_row_names <- rbind(matrix_output_row_names,sample_name)
      row_mean_vector <- NULL
      
    } # end else
    
  } # end while loop
  # return variable
  output <- matrix(result_matrix,nrow=nrow(result_matrix),ncol=ncol(result_matrix),byrow=F,
                   dimnames <- list(matrix_output_row_names,colnames(input_matrix)))
  output
  
}# end function
################################################################################################################################################################
################################################################################################################################################################
#==========================================================================================================
test.cor <- function(data){
  # calculating the p-values related to the correlation
  # remember that you have to edit the files in excel first - delete first column(numbers)
  # by default the method uses pearson correlation, change the code inside the function
  
  mat <- matrix(nrow=ncol(data),ncol=ncol(data),dimnames=list(colnames(data),colnames(data)))
  
  for(j in 1:ncol(data)){
    for(k in 1:ncol(data)){
      c <- cor.test(data[,j],data[,k],method="pearson",use="pairwise.complete.obs")
      #c <- cor.test(data[,j],data[,k],method="spearman")
      mat[j,k] <- c$p.value
      mat[k,j] <- c$p.value
    } # end nested nested for loop
  }# end nested for loop
  # return variable
  return(mat)
} # end function
#==========================================================================================================


