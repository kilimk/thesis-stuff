
main = function(file,cor.search.bins=9){
  
  # before  first use install following libraries: install.packages(c("qvalue","igraph","Hmisc",
  #"reshape2","plyr","dplyr"))
  #if one of libraries fails to install with error on R version , try to run this setRepositories() and choose bio based repos
  
  library(igraph);library(Hmisc);library(reshape2);
  library(pcaMethods);library(qvalue);library(plyr);library(dplyr)
  
  # Assumption 1 - data is in csv format with components (metabolites/apk function) as columns
  # and samples (genotypes/apk) in rows. first column is string type with entity names
  raw_data <- load_and_validate(file)
  
  # TODO : add normalization function call
  
  col_num = ncol(raw_data)
  # aggregate all columns by sample names ( in case of tomato pericarps there were 6 batches)
  samples_avg=aggregate(raw_data[,2:col_num], by=list(raw_data[,1]), FUN=mean, na.rm=TRUE)
  samples_avg = samples_avg[rowSums(is.na(samples_avg[,-1]))<(col_num-1),]
  # percentag of NAs
  print(paste("Percentage of missing values in given dataset % : " ,100*sum(sapply(samples_avg, function(x) sum(is.na(x)))) / (nrow(samples_avg)*(col_num-1))))  
  
  #aggregate function fills NAs with NaN, so it should be assigned NA for next imputation function
  samples_avg[is.na(samples_avg)] <- NA
  
  #TODO : add another options of imputation 
  # impute data
  pca <- pca(samples_avg,method = "nipals")
  filled_matrix <- pca@completeObs
  
  #calculate pair wise correlation matrix. Pearson or Spearman can be used
  #pearson is more common choice but it assumes normal distribution . Example to check.
  #norm_index = which(unlist(lapply(data.frame(k),function(x) shapiro.test(x)$p.value))>0.05)
  #print(paste("Percentage of normally distributed variables: ",length(norm_index)/col_num))
  
  cor <- correlation.thresholds(filled_matrix,min.r = 0,max.r = 0.9,
                                cor.method = 'pearson',num.bins = cor.search.bins)
  
  net.gold = graph_from_edgelist(as.matrix(cor[,1:2]),directed=F)
  E(net.gold)$weight = cor$Weight
  
  #construct graph based on correlation matrix as an adjacency matrix
  return(list(net=net.gold,cor.matrix=cor))
}

correlation.thresholds <- function(profile, min.r=0,max.r=0.9,cor.method="pearson",num.bins,logger=0){
  #   function that takes in a normalized dataset param =data and the minimum and maximum r values
  #   constructs a correlation matrix and corresponding p-value matrix.
  #   networks are constructed from varying thresholds and then subjected to network generation
  #   the networks are being analyzed for network measures
  #   the return value contains a dataframe that contains the different constellation of r
  #   and p-values with corresponding network measures
  
  # first create correlation matrix and corresponding p-value matrix by calling function rcorr from Hmisc library
  cor.list = rcorr(as.matrix(profile),type=cor.method)
  cor.r = as.matrix(cor.list$r)
  cor.p = as.matrix(cor.list$P)
  cor.r[upper.tri(cor.r,diag=T)]=0
  cor.p[upper.tri(cor.p,diag=T)]=0
  cor.df = data.frame(melt(cor.r))
  
  colnames(cor.df)=c("nodeA","nodeB","Weight")
  cor.df$Sign = melt(cor.p)$value
  cor.df = filter(cor.df,Weight!=0 & Sign<0.05)
  full.num.edges = nrow(cor.df)
  # not must but better for debug and perception
  cor.df.order = arrange(cor.df,abs(Weight),-Sign)
  
  r.values = seq(min.r,max.r,(max.r-min.r)/num.bins)
  p.values = seq(0.01,0.05,0.01)
  net.meta = arrange(expand.grid(r_value=r.values,p_value=p.values),r_value,-p_value)
  net.features = c("nnode","nedges","degree","density","clust.coef","diameter")
  net.meta[,net.features]=NA
  
  for (i in 1:nrow(net.meta)){
    r=net.meta$r_value[i]
    p=net.meta$p_value[i]
    temp.cor = filter(cor.df.order,abs(Weight)>r & Sign<p)
    
    net = graph_from_edgelist(as.matrix(temp.cor[,1:2]))
    #fill values into graph main features df
    net.meta[i,net.features] <- c(vcount(net),ecount(net),mean(degree(net)),
                                  graph.density(simplify(net),loops=F),
                                  transitivity(net,type="global"),
                                  diameter(net,unconnected=T))
  }
  
  print(net.meta)
  
  # find optimal r and adjusted p according to multitest hypothesis (fdr from qvalues package)
  sd_stat = aggregate(net.meta,by=list(net.meta$r_value),FUN=sd,na.rm=T)
  threshold.sensitivity = 0.005
  r.opt = min(sd_stat[which(sd_stat$nedges<threshold.sensitivity*full.num.edges),1])
  
  fdr_mat = fdr.fun(cor.p[lower.tri(cor.p,diag = F)])
  fp.sensitivity = 0.05
  p.adj=min(fdr_mat[,1][which(fdr_mat[,2]>=fp.sensitivity)])
  print(paste("Optimal correlation threshold: ",r.opt))
  print(paste("P value after correction: ",p.adj))
  return (filter(cor.df.order,abs(Weight)>r.opt & Sign<p.adj))
  #return(cor.df.order)
}
  
  

load_and_validate = function(file){
  
  df = read.csv(file)
  col_class = unlist(lapply(df, class))
  class_counts = table(col_class[-1])
  if (col_class[1]!="factor")
    warning("First column should be of factor class")
  if (length(class_counts)>1 || names(class_counts)[1]!="numeric")
    stop("Starting from second column all columns in data set should be numeric. Edit your data and try again")
  return (df)
  
}
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

compound.match  =function (x,y,num=2){
  len=length(x)
  
  match_Cname = character(len)
  score = numeric(len)
  others_top = character(len) 
  y = unique(y)
  
  for (i in 1:len ){
    
    compound = x[i]
    scores = vector()
    for (str2comp in y){
      scores = c(scores,levenshteinSim(as.character(compound),as.character(str2comp)))
    }
    ind = order(scores,decreasing = T)
    match_Cname[i] = as.character(y[ind[1]])
    score[i] = scores[ind[1]]
    if (num>=1)
      others_top[i] = paste(y[ind[2:num+1]],collapse="#")
   
  }
  
  return (data.frame(compound=x,match_Cname,score,others_top))
  
}