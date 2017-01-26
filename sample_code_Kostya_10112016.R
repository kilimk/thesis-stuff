## sample code for Kostya
## to run correlation analysis and extract significant correlations only

dat <- read.csv(file.choose(),check.names=F,row.names = 1)  # this reads in your dataset

## here load the source files containing the functions and function calls:
## files: 
## 1. network_feature_functions.R contains all the functions to run feature analysis on an igraph object
## 2. network_features_07092916 makes a function call to network_feature_functions: pay attention that some functions
## are substituted by NAs you have to remove comment symbols. There will be written comments in th file where
## I didn't use them as they were taking to long to run, especially on the grape network
## 3. correlation.thresholds.R contains all the functions to only extract significant correlations based on given thresholds
## you will also find a function that converts the network into cytoscape files
## 4. supporting_functions.R contains all kind of functions that are usefull
## you will have to find the corresponding libraries to run some functions

## running correlation analysis
## the next step assumes that you have normalized the dat file
## metabolites should be columns and the varieties/genotypes rows

cor <- correlation.thresholds(dat,0,1,'spearman')   # the last parameter can also be 'pearson'
# looking at prep variable to determine thresholds for r
cor@prop
r <- 0.4  # just an example
p <- triangle(cor@p.value.mat)
fdr <- fdr.fun(p.value=p)
p.adj <- 0.012

cor.sig <- correlation.network.thresh(cor@cor.matrix,cor@p.value.mat,r,p.adj)

# making graph object
library(igraph)
gr <- graph_from_adjacency_matrix(cor.sig,mode="lower",weighted=T,diag=F,add.colnames=NULL)

## running network feature calculations  - sample code
## you need to adjust this to your needs
subgraph.node.index <- c(1,2,3)  # can be a single number of a vector of numbers corresponding to the node ids
subgraph <- induced.subgraph(gr,subgraph.node.index)
nonzero <- which(get.edge.ids(gr,c(subgraph.node.index,subgraph.node.index))!=0)
subgraph.edge.index <- get.edge.ids(gr,c(subgraph.node.index,subgraph.node.index))[nonzero]
features <- network_features(gr,subgraph,subgraph.node.index,subgraph.edge.index,community.detection=T)

