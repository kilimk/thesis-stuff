correlation.thresholds <- function(data, min.r,max.r,method){
#   function that takes in a normalized dataset param =data and the minimum and maximum r values
#   constructs a correlation matrix and corresponding p-value matrix.
#   networks are constructed from varying thresholds and then subjected to network generation
#   the networks are being analyzed for network measures
#   the return value contains a dataframe that contains the different constellation of r
#   and p-values with corresponding network measures
#   parameter method holds which correlation to run: pearson, spearman etc.
  
  # first create correlation matrix and corresponding p-value matrix by calling function correlation.network
  print("function: correlation thresholds")
  cor.network <- correlation.network(data,method)
  cor <- cor.network@cor.matrix  # by default pearson correlation
  # for the next call, the test_cor function must be preloaded
  p.cor <- cor.network@p.value.mat
  
  # prepare a matrix that contains return values
  network.measures <- matrix(nrow=55,ncol=8)  # this is hard coded
  colnames(network.measures) <- c("r-value","p-value","nnode","nedges",
                                  "degree","density","clust.coef","diameter")
  
  index.seq <- seq(1,51,5)
  r.seq <- seq(min.r,max.r,(max.r-min.r)/10) # divide into 11 bins, from r.min to r.max
  network.measures[index.seq,1] <-  r.seq
  network.measures[,2] <- rep(seq(0.01,0.05,0.01),11)
  #network.measures[,2] <- rep(seq(0.006,0.01,0.001),11)
  # now run nested for loops to go over thresholds r
  count  <- 1 # count variable for the row index of network measure table
  for (r in r.seq){
    # for(p in seq(0.006,0.01,0.001)){
    for(p in seq(0.01,0.05,0.01)){  # hard coded to always go from p-values 0.01 to 0.05
      # do the networks here
     
      # create r coefficient network for specific r and p values
      # call function correlation.network.thresh
      sig.cor <- correlation.network.thresh(cor,p.cor,r,p)
    
      # now that we have the matrix only with "significant" correlations
      # we can construct the correlation networks
      # replace all the significant correlations with value 1
      sig.cor.ones <- sig.cor
      x <- which(sig.cor!=0)
      sig.cor.ones[x] <- 1
      # construct network using igraph method -> load library igraph
      library(igraph)
      
      # convert to igraph network
      for(m in 1:nrow(sig.cor.ones)){sig.cor.ones[m,m] <- 0}
      gr1 <- graph.adjacency(sig.cor.ones,mode="lower", diag=F)
      
      # calculate network parameters
      num.nodes <- vcount(gr1)
      num.edges <- ecount(gr1)
      degree <- mean(degree(gr1))
      density <- graph.density(simplify(gr1),loops=F)
      clus.coef <- transitivity(gr1,type="global")
      diameter <- diameter(gr1,unconnected=T)
      
      #fill values into return table: network.measures
      network.measures[count,3:8] <- c(num.nodes,num.edges,degree,
                                       density,clus.coef,diameter)
      # increment count
      count  <- count + 1
    }# end nested loop
  } # end loop
  
  # prepare return variable - S4 class
  setClass("network.properties",representation(prop="matrix",
                                              cor.matrix="matrix",p.value.mat="matrix"))
  # prepare return object 
  ret.object  <- new("network.properties",prop=network.measures,cor.matrix=cor,p.value.mat=p.cor)
  # return
  return(ret.object)
} # end function

#==========================================================================================================
correlation.network <- function(data,method='pearson'){
  # create correlation matrix and corresponding p-value matrix
  # parameter method holds type of correlation to be run
  print("function: correlation network")
  cor <- cor(data,method=method)  # by default pearson correlation
  # for the next make call to function test.cor
  p.cor <- test.cor(data,method)
  
  # define return variable - S4 class
  setClass("cor.network",representation(cor.matrix="matrix",p.value.mat="matrix"))
  # prepare return object 
  ret.object  <- new("cor.network",cor.matrix=cor,p.value.mat=p.cor)
  return(ret.object)
} # end function

#==========================================================================================================
correlation.network.thresh <- function(cor,p.cor,r,p){
  # same as function correlation.network, but returns an r coefficient matrix for specific 
  # r coefficient and p-value thresholds
  # vector containing indeces for specific r threshold
    print("function: correlation network thresh")
  x.cor <- which(abs(cor)>=r,arr.ind=T)
  # return variable holding matrix for specific r and p values
  sig.cor <- cor
  sig.cor[,] <- 0
  
  # looping over x.cor to identify p values at same indeces and filling sig.cor
  for(k in 1:nrow(x.cor)){
    if(is.na(p.cor[x.cor[k,1],x.cor[k,2]]))
      sig.cor[x.cor[k,1],x.cor[k,2]] <- 0
    else if(p.cor[x.cor[k,1],x.cor[k,2]]<=p){
      sig.cor[x.cor[k,1],x.cor[k,2]]  <- cor[x.cor[k,1],x.cor[k,2]]
    } #end if
    else if(is.na(p.cor[x.cor[k,1],x.cor[k,2]]))
      sig.cor[x.cor[k,1],x.cor[k,2]] <- 0
    
  } #end loop
  return(sig.cor)
} # end function

#==========================================================================================================
test.cor <- function(data,method){
  # calculating the p-values related to the correlation
  # remember that you have to edit the files in excel first - delete first column(numbers)
  # by default the method uses pearson correlation, change the code inside the function
    print("function: test cor")
  mat <- matrix(nrow=ncol(data),ncol=ncol(data),dimnames=list(colnames(data),colnames(data)))
  
  for(j in 1:ncol(data)){
      print(j)
    for(k in 1:ncol(data)){
      c <- cor.test(data[,j],data[,k],method=method,use="pairwise.complete.obs")
      #c <- cor.test(data[,j],data[,k],method="spearman")
      mat[j,k] <- c$p.value
      mat[k,j] <- c$p.value
    } # end nested nested for loop
  }# end nested for loop
  # return variable
  return(mat)
} # end function
#==========================================================================================================

normalize.data <- function(data){
  # function that normalizes a data matrix by chromatogram = rows
  # and then by metabolite = cols
  
  # first complete dataset by imputation
  # load pcaMethods library for imputation 
  #!!!!!!!!!!! this library can be downloaded from Bioconductor only!!!!!!!!!!!!!!
  library(pcaMethods)
  # first remove all NA cols and rows 
  data <- remove.na.col(data); data <- remove.na.row(data)
  # use pca method for imputation
  pca <- pca(data,method="ppca")
  data.complete <- pca@completeObs
  
  # average rows with trimmed mean 0.05
  row.mean <- apply(data.complete,1,mean,trim=0.05)
  # normalize each cell of by corresponding mean value of row.mean
  for(i in 1:length(row.mean)){
      data.complete[i,] <- data.complete[i,]/row.mean[i]
    } # end for loop
  
  # do the same for cols now: average cols with trimmed mean 0.05
  col.mean <- apply(data.complete,2,mean,trim=0.05)
  # normalize each cell of by corresponding mean value of row.mean
  for(i in 1:length(col.mean)){
      data.complete[,i] <- data.complete[,i]/col.mean[i]
    } # end for loop
  # return variable
  return(data.complete)  
} # end function

#==========================================================================================================

remove.na.col <- function(matrix){
  # removing all columns with NAs only
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
  } # end if
  return(matrix)   # return variable
} # end function

#==========================================================================================================

remove.na.row <- function(matrix){
  # function removing all rows with NAs only
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
  return(matrix)   # return variable
} # end function

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
    x <- which(tolower(row.names(dataset.2))==tolower(row.names(dataset.1)[i]))
    
    # see if entry exists
    if(length(x)>0) { # it exists in second dataset
      output.dataset.1 <- rbind(output.dataset.1,dataset.1[i,])
      row.names(output.dataset.1)[nrow(output.dataset.1)] <- row.names(dataset.1)[i]
      output.dataset.2 <- rbind(output.dataset.2,dataset.2[x,])
    } # end if
  } # end for loop
  
  # return variable
  y <- cbind(output.dataset.1,output.dataset.2)
  return(y)
} # end function

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
    x <- which(tolower(colnames(dataset.2))==tolower(colnames(dataset.1)[i]))
    
    # see if entry exists
    if(length(x)>0) { # it exists in second dataset
      output.dataset.1 <- cbind(output.dataset.1,dataset.1[,i])
      colnames(output.dataset.1)[ncol(output.dataset.1)] <- colnames(dataset.1)[i]
      output.dataset.2 <- cbind(output.dataset.2,dataset.2[,x])
    } # end if
  } # end for loop
  
  # return variable
  y <- rbind(output.dataset.1,output.dataset.2)
  return(y)
} # end function

netw.membership <- function(cor,r,p){
  # cor = an object of function correlation.thresholds
  # r = correlation coefficient, p = correlation p value
  # function defining memberships of the network according to the walktrap.community algorithm
  # returns S4 object containing edgelist and membership vector and a matrix with r coefficient
  sig.netw <- correlation.network.thresh(cor@cor.matrix,cor@p.value.mat,r,p)
  gr <- get.igraph(sig.netw)
  # find membership
  wc <- walktrap.community(gr)
  mem <- as.matrix(membership(wc))
  mem <- cbind(mem,degree(igraph))
  row.names(mem) <- gr$names
  for(i in 1:nrow(sig.netw)){sig.netw[i,i] <- 0 }
  setClass("network",representation(cor.matrix="matrix",membership="matrix",degree="matrix",edgelist = "matrix"))
  # prepare return object 
  ret.object  <- new("network",cor.matrix=sig.netw,membership=mem[,1],degree=mem[,2],edgelist = get.edgelist(gr))
  return(ret.object)
} # end function

netw.membership.igraph <- function(igraph){
  if(vcount(igraph)==0 | ecount(igraph)==0){
    print("No exisiting edges to be grouped!")
    return(0)
  }
  # function defining memberships of the network according to the walktrap.community algorithm 
  # and returning the degree of respective node
  # parameter = igraph object
  # returns membership vector as matrix
  # find membership
  wc <- walktrap.community(igraph)
  #lec <- leading.eigenvector.community(igraph,weights=T)
  mem <- as.matrix(membership(wc))
  mem <- cbind(mem,degree(igraph))
  colnames(mem) <- c("membership","degree")
  row.names(mem) <- igraph$names
  return(mem)
} # end function


#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

# this function prepares a matrix for cytoscape analysis, such as a correlation matrix and generates three files: the connection between metabolites, an attribute file for the metabolite, and an attribute file for the edges to be read in into cytoscape

# incoming parameters: data_matrix - the matrix with the correlation significant data
# comp_class: an object of the comp_class library
# membership: a matrix object containing the membership value for each node of the network according to the walktrap.community algorithm
# file_name: the file_name to be given for the output files; type string
cytoscape_preparation <- function(data_matrix,comp_class_file,membership,file_name){
  
  # variable initilization
  
  # open files
  con <- file(paste(file_name, "_connections.txt",sep=""),"a")
  cat(paste("Met.1","\t","Connection","\t","Met.2","\n",sep=""),file=con)
  ea <- file(paste(file_name,"_ea.txt",sep=""),"a")
  cat(paste("Connection","\t","abs_r_value","\t","pos_neg","\n",sep=""),file=ea)
  na <- file(paste(file_name,"_na.txt",sep=""),"a")
  cat(paste("Metabolite","\t","Compound_class","\t","Membership","\t","Type","\t","Degree","\n",sep=""),file=na)
  
  compound_value <- NULL
  degree <- NULL
  mem <- NULL
  type <- NULL
  row_value=2
  #row_value=1
  
  # loop over columns of data_matrix to find connections between metabolites
  for(i in 1:(ncol(data_matrix))){
    #for(i in 1:(ncol(data_matrix))){	
    # find membership of metabolite
    x <- which(row.names(membership)==colnames(data_matrix)[i])
    if(length(x)==0){mem <- 0; degree <- 0} 
    else{mem <- membership[x,1]; degree <- membership[x,2]}
    
    # find metabolite in file comp_class_file and assign number according to compound class
    cc <- which(comp_class_file == colnames(data_matrix)[i],arr.ind=T)
    if(is.na(cc[2])){compound_value <- 0} else{compound_value <- cc[2]}
    # for testing purposes
    if(length(x)==0)
    print(paste(colnames(data_matrix)[i]," ,Memb: ",mem," ,Degree: ", degree, sep=""))
    # specific for Asfaw and Uris data
    if(compound_value==0){type <- 0}
    else if(compound_value<9 & compound_value>0){type <- 1} # primary
    else if(compound_value<13 & compound_value>9){type <- 2} # secondary
    else{type  <- 3} # physiological
    str_n <- paste(colnames(data_matrix)[i],"\t",compound_value,"\t",mem,"\t",type,"\t",degree,"\n",sep="")
    cat(str_n,file=na)
    
    ## break the loop if the i == ncol data_matrix because we don't need to check 
    ## the connections of the last column but we need in the node property file
    if(i==ncol(data_matrix)){break}
    # check for actual value
 
    for(j in row_value:nrow(data_matrix)){
      # first find all connections
     
      if(data_matrix[j,i]==0 && j<nrow(data_matrix)){next}
      
      if(data_matrix[j,i]!=0) {
        #if(j==102){print("yes")}
      # append connections one by one to data.frame containing the connections between metabolites
      connection <- paste(colnames(data_matrix)[i],"\t",
                          colnames(data_matrix)[i],"-",rownames(data_matrix)[j],"\t",
                          rownames(data_matrix)[j],"\n",sep="")
      cat(connection,file=con)
      # for the connection found, prepare output files for edges
     
      p_n <- 0
      # determine whether it is a positive or negative correlation
      if(data_matrix[j,i]>0){p_n=1} else{p_n=-1}
      # format string
      str_e <- paste(colnames(data_matrix)[i],"-",rownames(data_matrix)[j],"\t",
                     abs(data_matrix[i,j]),"\t",p_n,"\n",sep="")
      #str_e <- paste(colnames(data_matrix)[i],"-",rownames(data_matrix)[j],"\t",
      #					abs(data_matrix[j,i]),"\t",p_n,"\n",sep="")					
      cat(str_e,file=ea)		
      } # end if		
     
    } # end nested for loop
    row_value <- row_value+1
    #		row_value <- 1
  }# end for loop
  
  
  # closing files
  close(con); close(ea); close(na)
  
  
} # end function
##############################################################################################
cytoscape_preparation_sym_dif <- function(data_matrices,comp_class_file,file_name){
  # same as function cytoscape_preparation - the difference is that this takes in multiple
  # netw matrices as a list and appends them to each other
  # to run this function you need to run function net.differnce in file grape.R

  
  # variable initilization
  
  # open files
  con <- file(paste(file_name, "_connections.txt",sep=""),"a")
  cat(paste("Met.1","\t","Connection","\t","Met.2","\n",sep=""),file=con)
  ea <- file(paste(file_name,"_ea.txt",sep=""),"a")
  cat(paste("Connection","\t","abs_r_value","\t","pos_neg","\t","edge_affil","\n",sep=""),file=ea)
  na <- file(paste(file_name,"_na.txt",sep=""),"a")
  cat(paste("Metabolite","\t","Compound_class","\t","Type","\n",sep=""),file=na)
  
  metabolite.list <- NULL
  # start here
  for(d in 1:length(data_matrices)){
    if(d==1) metabolite.list <- colnames(data_matrices[[d]])
    else{
      # loop over other matrices
      for(i in 1:ncol(data_matrices[[d]])){
        met.name <- which(tolower(metabolite.list)==tolower(colnames(data_matrices[[d]]))[i])
        # if it doesn't exist append
        if(length(met.name)==0) 
          metabolite.list <- c(metabolite.list,colnames(data_matrices[[d]])[i])
      } # end for loop
    } # end else
    } # end for loop
  
  compound_value <- NULL
  # prepare na file
  for(i in 1:length(metabolite.list)){
    cc <- which(comp_class_file==metabolite.list[i],arr.ind=T)
    if(is.na(cc[2])){compound_value <- 0} else{compound_value <- cc[2]}
    
    # specific for Asfaw and Uris data
    if(compound_value==0){type <- 0}
    else if(compound_value<9 & compound_value>0){type <- 1} # primary
    else if(compound_value<13 & compound_value>3){type <- 2} # secondary
    else{type  <- 3} # physiological
    str_n <- paste(metabolite.list[i],"\t",compound_value,"\t",type,"\n",sep="")
    cat(str_n,file=na)
  } # end for loop
  
  for(d in 1:length(data_matrices)){
    
  type <- NULL
  row_value=2
  #row_value=1
  data_matrix <- data_matrices[[d]]
  
  # loop over columns of data_matrix to find connections between metabolites

    # check for actual value
  for(i in 1:(ncol(data_matrix))){  
      ## break the loop if the i == ncol data_matrix because we don't need to check 
      ## the connections of the last column but we need in the node property file
      if(i==ncol(data_matrix)){break}
      
    for(j in row_value:nrow(data_matrix)){
      # first find all connections
      
      if(data_matrix[j,i]==0 && j<nrow(data_matrix)){next}
      
      if(data_matrix[j,i]!=0) {
        #if(j==102){print("yes")}
        # append connections one by one to data.frame containing the connections between metabolites
        connection <- paste(colnames(data_matrix)[i],"\t",
                            colnames(data_matrix)[i],"-",rownames(data_matrix)[j],"\t",
                            rownames(data_matrix)[j],"\n",sep="")
        cat(connection,file=con)
        # for the connection found, prepare output files for edges
        
        p_n <- 0
        # determine whether it is a positive or negative correlation
        if(data_matrix[j,i]>0){p_n=1} else{p_n=-1}
        # format string
        str_e <- paste(colnames(data_matrix)[i],"-",rownames(data_matrix)[j],"\t",
                       abs(data_matrix[i,j]),"\t",p_n,"\t",d, "\n",sep="")
        #str_e <- paste(colnames(data_matrix)[i],"-",rownames(data_matrix)[j],"\t",
        #					abs(data_matrix[j,i]),"\t",p_n,"\n",sep="")					
        cat(str_e,file=ea)		
      } # end if		
      
    } # end nested for loop
    row_value <- row_value+1
    #		row_value <- 1
  }# end for loop
  } # end for loop for list
  
  
  # closing files
  close(con); close(ea); close(na)
  
  
} # end function
##############################################################################################
##############################################################################################

identical_row_names <- function(data){
  # function to make sure rownames are unique by labelling them 
  for(i in 1:nrow(data)){
    x <- which(rownames(data)==rownames(data)[i])
    if(length(x)>1){
      print(rownames(data)[i])
      counter <- 1
      for(j in x){
        rownames(data)[j] <- paste(rownames(data)[j],".",counter,sep="")
        counter <- counter+1
      }# end if
    }
  }
  return(data)
} # end function

identical_col_names <- function(data){
  # function to make sure rownames are unique by labelling them 
  for(i in 1:nrow(data)){
    x <- which(colnames(data)==colnames(data)[i])
    if(length(x)>1){
      print(colnames(data)[i])
      counter <- 1
      for(j in x){
        colnames(data)[j] <- paste(colnames(data)[j],".",counter,sep="")
        counter <- counter+1
      }# end if
    }
  }
  return(data)
} # end function
################################################################################################################################################################
network.intersect <- function(sig.cor.ls){
    # function determining the network intersect between all supplied networks in paratmeter sig.cor.ls of type list
    
    # first convert all matrices omittiing the actual r values
    for(i in 1:length(sig.cor.ls)){
        sig.cor.ls[[i]][sig.cor.ls[[i]]>0] <- 1 # for positive correlation
        sig.cor.ls[[i]][sig.cor.ls[[i]]<0] <- 2 # for negative correlation
    } # end for loop
    
    # find all connections in the first matrix and use it as originator
    ret.mat <- sig.cor.ls[[1]]
    links <- which(ret.mat!=0,arr.ind=T)
    
    # now look if these links are also present in the other matrices
    # !!!! very important the nodes must have same names
    
    for(i in 1:nrow(links)){
        row.n <- tolower(rownames(ret.mat)[links[i,1]])
        col.n <- tolower(colnames(ret.mat)[links[i,2]])
        
        for(j in 2:length(sig.cor.ls)){
            row.num <- which(tolower(rownames(sig.cor.ls[[j]]))==row.n)
            col.num <- which(tolower(colnames(sig.cor.ls[[j]]))==col.n)
            # check if the connection is existent in the other matrix
            if(length(row.num)>=1 & length(col.num)>=1 ){
                # nodes exist
                # now check if the connection between nodes exist
                con <- sig.cor.ls[[j]][row.num,col.num]
                if(con==0){
                    # it doesn't exist so set connection to 0 in originator matrix and break the loop
                    ret.mat[links[i,1],links[i,2]] <- 0
                    break;
                }
                else{
                    if(con==1)
                        ret.mat[links[i,1],links[i,2]] <- (ret.mat[links[i,1],links[i,2]])+1
                    else
                        ret.mat[links[i,1],links[i,2]] <- (ret.mat[links[i,1],links[i,2]])+2
                }
            } # end if
            else{
                # connection does not exist
                # set connection to 0
                ret.mat[links[i,1],links[i,2]] <- 0
                break;
            } # end else
        } # end nested loop
    } # end for loop
    
    # some explanation for the coming code:
    # the dataset returned contains only those connection that are present in all networks, i.e. positive and negative
    # the positive connections are denoted by value 1 and the negative by value 2
    # if we add up the values and the same connection was positive in all networks, we will get a value of 
    # datasets provided. Therefore, if value/number of datasets = 1 - only under this scenario we have the same positive connections
    # if value/number of datasets = 2 then we have all negative correlation in all datasets
    # any value inbetween 1 and 2 indicates that the connetion exist in all datasets but they are not the same
    # positive and negtive
    # thus
    ret.mat[ret.mat != 0] <- (ret.mat[ret.mat != 0])/(length(sig.cor.ls))
    # set diagonal to zero - just to make sure
    for(i in 1:nrow(ret.mat)) ret.mat[i,i] <- 0
    
    return(ret.mat)
    
} # end function
################################################################################################################################################################
network.difference <- function(sig.cor.ls){
    # function determine the network difference for each network submitted to other networks submitted
    # with parameter of type list
    # prepare a return variable of type list
    ret.list <- sig.cor.ls
    
    for(i in 1:length(sig.cor.ls)){
        # determine sequence for network comparison within list vector
        seq <- NULL
        if(i==1)  seq  <- c(2:length(sig.cor.ls))
        else if(i==length(sig.cor.ls)){
            if(i==2) seq <- 1
            else seq <- c(1:(length(sig.cor.ls)-1))
        }
        else seq <- c(1:(i-1),(i+1):length(sig.cor.ls))  
        
        links <- which(sig.cor.ls[[i]]!=0,arr.ind=T)
        
        # now make the comparison
        for(j in 1:nrow(links)){
            row.n <- tolower(rownames(sig.cor.ls[[i]]))[links[j,1]]
            col.n <- tolower(colnames(sig.cor.ls[[i]]))[links[j,2]]
            for(k in seq){
                row.num <- which(tolower(rownames(sig.cor.ls[[k]]))==row.n)
                col.num <- which(tolower(colnames(sig.cor.ls[[k]]))==col.n)
                # if exists then set corresponding cell value of ret.list to zero
                if(length(row.num)>=1 & length(col.num)>=1){
                    # check for the value
                    if(sig.cor.ls[[k]][row.num,col.num]!=0){
                        # it does exist
                        ret.list[[i]][links[j,1],links[j,2]] <- 0
                        break; # no need to look into the other networks
                    } # end if
                } # end if
            } # end nested nested loop
        } # end nested loop
    } # end for loop
    return(ret.list)
} # end function
################################################################################################################################################################