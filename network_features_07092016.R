## AUTHOR: DAVID TOUBIANA
## Date original: 23.01.2015
## Date modifications: 07.09.2016
## Modifications are improvements to the code so that it runs faster, especially if 
## many instances must be run
## container file containing all the functions to compute features based on network properties
## for each property one seperate function
## here, especially when computing the 4 mathematical moments, one function will be available for all
## one function to compute the property/feature for the subgraph of the network
## one function to compute the property/feature of the entire network
network_features  <- function(graph, subgraph, subgraph.node.index=1,subgraph.edge.index=1,
                              community.detection=F){
    # this is the function that will be called, which will make calls to all other functions
    # using parameters graph and subgraph to compute features
    # library igraph has to be imported for the graph property functions to run
    #library(igraph)
    # library e1071 has to be imported for skewness and kurtosis measures to run
    library(e1071)
    #source('~/PostDoc April 2016/CNA-ML_R/network_feature_functions.R', echo=TRUE)
    # variable to be returned
    #------------------------------------------------------------------------------------------------
    # precalculating community structures if parameter community.detection=T - by default it is FALSE
    # only want to do it once for the graph
    if(community.detection){
        gr <- graph
        E(gr)$weight <- abs(E(gr)$weight)
        # fast greedy, non-weighted and weighted
        fgc <- fastgreedy.community(gr,weights=NULL)
        fgc.w <- fastgreedy.community(gr,weights=E(gr)$weight)
        # walktrap, non-weighted and weighted
        wtc <- walktrap.community(gr,weights=NULL)
        wtc.w <- walktrap.community(gr,weights=E(gr)$weight)
        # edge betweenness community
        ebc <- walktrap.community(gr,weights=NULL)
        ebc.w <- walktrap.community(gr,weights=E(gr)$weight)
        # label propagation community
        lpc <- walktrap.community(gr,weights=NULL)
        lpc.w <- walktrap.community(gr,weights=E(gr)$weight)
        # leading eigenvector community
        lec <- walktrap.community(gr,weights=NULL)
        lec.w <- walktrap.community(gr,weights=E(gr)$weight)
        # multi level community
        mlc <- walktrap.community(gr,weights=NULL)
        mlc.w <- walktrap.community(gr,weights=E(gr)$weight)
    }
    #------------------------------------------------------------------------------------------------
    return.variable = NULL
    
    return.variable <- cbind(return.variable,subgraph.num.edges(subgraph))                    # no. 1   ok
    return.variable <- cbind(return.variable,t(subgraph.degree(subgraph)))                    # no. 2-5   ok
    return.variable <- cbind(return.variable,t(subgraph.weighted.degree(subgraph)))           # no. 6-10   ok
    return.variable <- cbind(return.variable,t(subgraph.abs.weighted.degree(subgraph)))       # no. 11-15   ok
    print(15)
    return.variable <- cbind(return.variable,geodesic.distance.subgraph(subgraph))            # no. 16 ok
    return.variable <- cbind(return.variable,t(geodesic.distance.subgraph.moments(subgraph)))    # no 17-20 
    return.variable <- cbind(return.variable,t(geodesic.distance.weighted.subgraph.moments(subgraph)))    # no 21-25
    #return.variable <- cbind(return.variable,group.betweenness.subgraph(graph,subgraph.node.index)) # no 26 - uses lists - slower
    # only for testing purposes - uncomment when done
    # is too slow put NAs instead
    return.variable <- cbind(return.variable,NA)
    return.variable <- cbind(return.variable,NA)
   # return.variable <- cbind(return.variable,group.betweenness.subgraph(graph,subgraph.node.index)) # no 26 - uses loops - faster
   # return.variable <- cbind(return.variable,weighted.group.betweenness.subgraph(graph,subgraph.node.index)) # no 27 - uses lists - faster
    # return.variable <- cbind(return.variable,weighted.group.betweenness.subgraph.2(graph,subgraph.node.index)) # no 27 - uses loops - slower
    return.variable <- cbind(return.variable,t(closeness.centrality.subgraph(subgraph)))    # no 28-31
    (print(31))
    return.variable <- cbind(return.variable,t(weighted.closeness.centrality.subgraph(subgraph)))    # no 32-35
    return.variable <- cbind(return.variable,t(closeness.centrality.graph.subgraph(graph,subgraph.node.index)))   # no 36-39  
    return.variable <- cbind(return.variable,t(weighted.closeness.centrality.graph.subgraph(graph,subgraph.node.index)))  # no 40-43
    return.variable <- cbind(return.variable,diameter.subgraph(subgraph))  # no 44
    return.variable <- cbind(return.variable,weighted.diameter.subgraph(subgraph)) # no 45
    print(45)
    # only for testing purposes - uncomment when done
    # is too slow put NAs instead
    return.variable <- cbind(return.variable,NA)
    return.variable <- cbind(return.variable,NA)
    # return.variable <- cbind(return.variable,diameter.centrality(graph,subgraph.node.index))  # no 46
    # return.variable <- cbind(return.variable,weighted.diameter.centrality(graph,subgraph.node.index)) # no 47
    return.variable <- cbind(return.variable,t(clustering.coefficient.subgraph(subgraph)))  # no 48-52
    return.variable <- cbind(return.variable,t(weighted.clustering.coefficient.subgraph(subgraph)))  # no 53-57
    return.variable <- cbind(return.variable,
                             t(clustering.coefficient.subgraph.graph(graph,subgraph.node.index))) # no 58-61
    print(61)
    return.variable <- cbind(return.variable,
                             t(weighted.clustering.coefficient.subgraph.graph(graph,subgraph.node.index))) # no 62-65
    return.variable <- cbind(return.variable,t(node.betweenness.subgraph(graph,subgraph.node.index))) # no 66-69
    return.variable <- cbind(return.variable,t(weighted.node.betweenness.subgraph(graph,subgraph.node.index))) # no 70-73
    return.variable <- cbind(return.variable,t(edge.betweenness.subgraph(graph,subgraph.edge.index))) # no 74-77
    print(77)
    return.variable <- cbind(return.variable,t(weighted.edge.betweenness.subgraph(graph,subgraph.edge.index))) # no 78-81
    return.variable <- cbind(return.variable,assortativity.subgraph(subgraph)) # no 82
    return.variable <- cbind(return.variable,density.subgraph(subgraph)) # no 83
    return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,fgc)) # no 84
    return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,fgc.w)) # no 85
    return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,wtc)) # no 86
    return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,wtc.w)) # no 87
    return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,ebc)) # no 88
    return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,ebc.w)) # no 89
    return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,lpc)) # no 90
    print(90)
    return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,lpc.w)) # no 91
    return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,lec)) # no 92
    return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,lec.w)) # no 93
    return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,mlc)) # no 94
    return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,mlc.w)) # no 95
    # friends measures: common friends, total friend, distinct friends, mixed friends, and
    # preferential attachment score
    return.variable <- cbind(return.variable,t(friend.measures(graph,subgraph.node.index))) # no 96 - 100
    print(100)
    #return(return.variable)
    # the jaccard coefficient is the commonfriends measure / totalfriends measure
    # features 101-105
    if(length(subgraph.node.index)<2){
        sni <- unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))
        return.variable <- cbind(return.variable,t(jaccard.coefficient(graph,sni,
                                                    return.variable[length(return.variable)-4],
                                                    return.variable[length(return.variable)-3])))
    } else{
        return.variable <- cbind(return.variable,t(jaccard.coefficient(graph,subgraph.node.index,
                                                    return.variable[length(return.variable)-4],
                                                    return.variable[length(return.variable)-3])))
    }
    print(105)
    return.variable <- cbind(return.variable,t(friends.measure.subgraph(subgraph)))  # no 106-110
    return.variable <- cbind(return.variable,t(friends.measure.graph(graph,subgraph.node.index)))  # no 111-115
    print(115)
    # whole network features
    return.variable <- cbind(return.variable,num.edges.graph(graph,subgraph.node.index)) # no 116
    #print(116)
    return.variable <- cbind(return.variable,t(degree.graph(graph,subgraph.node.index))) # no 117-120
    #print(120)
    return.variable <- cbind(return.variable,t(weighted.degree.graph(graph,subgraph.node.index))) # no 121-125
    print(125)
    return.variable <- cbind(return.variable,t(abs.weighted.degree.graph(graph,subgraph.node.index))) # no 126-130
    #print(130)
    return.variable <- cbind(return.variable,t(geodesic.distances.graph(graph,subgraph.node.index))) # no 131-135
    return.variable <- cbind(return.variable,t(geodesic.distances.weighted.graph(graph,
                                                subgraph.node.index))) # no 136-140
    # is too slow put NAs instead
    return.variable <- cbind(return.variable,t(rep(NA,8)))
    #return.variable <- cbind(return.variable,NA)
    # return.variable <- cbind(return.variable,t(stress.centrality(graph,subgraph.node.index))) # no 141-144
    # return.variable <- cbind(return.variable,t(weighted.stress.centrality(graph,subgraph.node.index))) # no 145-148
    print(148)
    return(return.variable)
} # end function

