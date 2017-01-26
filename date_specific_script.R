# script: 06.09.2015
# edit all pathways
# create a list for all compounds
all.species.pathways <- pathways(species = "all")
# in order to get function compounds going, define a variable of type string, adding the different 
# species/pathways you are interested in
pathways <- "all.species.pathways"
all.species.compounds <- compounds(pathways = pathways)
write.csv(unlist(all.species.compounds),"all_species_compound_list_060915.csv")
# now for tomato
tomato.pathways <- pathways(species = "Solanum lycopersicum")
pathways <- "tomato.pathways"
tomato.compounds <- compounds(pathways = pathways)
write.csv(unlist(tomato.compounds),"tomato.compounds_060915.csv")
# 20.07.2015  ended working here, begin with next line

# for 2001 dataset
# reading in datafiles tom.01 = fruit pericarp metabolites of 2001 etc.
tom.01 <- read.delim(file.choose(),check.names=F)
temp <- tom.01[order(tom.01[,2]),]
temp <- as.matrix(temp)
rownames(temp) <- temp[,2] 
temp <- temp[,-c(1,2)]
tom.sorted.01 <- temp
tom.sorted.01 <- convert_to_num(tom.sorted.01)
tom.avg.01 <- calc_mean(tom.sorted.01)
# percentag of NAs
length(which(is.na(tom.avg.01))) / length(tom.avg.01)  # = 10.04%
nas <- (which(is.na(tom.avg.01)))
tom.avg.01[nas] <- NA
library(pcaMethods)
# impute data
pca <- pca(tom.avg.01,method = "nipals")
tom.avg.01 <- pca@completeObs

# checking compound names with tomato compound list
match(tolower(colnames(tom.avg.01)),tolower(unlist(tomato.compounds))) # poor results
# use grep instead
for(i in 1:dim(tom.avg.01)[2]){
  comp <- grep(tolower(colnames(tom.avg.01)[i]),unlist(tomato.compounds))
  print(paste(i," compound: ",tolower(colnames(tom.avg.01)[i])))
  print(comp)
  print(unlist(tomato.compounds)[comp])
} # end for loop
cnames.01 <- c("cis-aconitate","L-alanine","arabinose","L-arginine","L-asparagine","L-aspartate","2-oxoglutarate",  
            "alpha;-tocopherol","benzoate","beta;-alanine","citramalate","citrate","L-cysteine","dehydroascorbate (bicyclic form)",
            "erythritol","fa16_0","fa18_0","fa18_2","fructose","fructose 1,6-bisphosphate","GDP-beta;-L-fucose","fumarate",
            "gaba","galactose","UDP-alpha;-D-galacturonate","D-gluconate","glucose","glucose_6_phosphate", "L-glutamate",
            "L-glutamine","D-glycerate","glycerate_3_Phosphate","glycerol","glycerol_1p",           
            "glycine","L-homoserine","inositol","inositol_1p","D-threo-isocitrate","L-isoleucine","isomaltose","L-leucine","L-lysine",          
            "l-ascorbate","(S)-malate","maleate","maltitol","maltose","D-mannitol","GDP-alpha;-D-mannose","melezitose","L-methionine",            
            "5-oxoproline","L-phenylalanine","phosphate","L-proline","putrescine","L-quinate","rhamnose","D-ribofuranose", 
            "saccharate","L-serine","shikimate","sorbitol","succinate","sucrose","Se-methyl-L-selenocysteine","t4_ho_proline",
            "threonate","L-threonine","alpha;,alpha;-trehalose","L-tryptophan","L-tyrosine","uracil","L-valine")

special_compounds <- c("GDP-beta;-L-galactose","UDP-alpha;-D-galactose","alpha;-D-galactose","ADP-alpha;-D-glucose",
                       "GDP-alpha;-D-glucose","UDP-alpha;-D-glucose","alpha;-D-glucose","beta;-D-glucose","UDP-L-rhamnose",
                       "dTDP-beta;-L-rhamnose" )
omit_compounds <- read.delim(file.choose(),check.names=F)

# checking again
length(which(is.na(match(cnames.01,unlist(tomato.compounds)))))  # = 20 NAs - sugars = 17 NAs out of 75 compounds
cnames.01[which(is.na(match(cnames.01,unlist(tomato.compounds))))]

#=========================================================================================================================
# 06.09.2015 - double instance of same pathways are no longer omitted
IL.2001.pathways <- find.pathways(tomato.pathways,special.compounds = special_compounds,omit.compounds = omit_compounds,cnames = cnames.01)


#=========================================================================================================================
# 06.09.2016
# for 2003 dataset
# reading in datafiles tom.01 = fruit pericarp metabolites of 2001 etc.
tom.03 <- read.delim(file.choose(),check.names=F)
temp <- tom.03[order(tom.03[,1]),]
temp <- as.matrix(temp)
rownames(temp) <- temp[,1] 
temp <- temp[,-1]
tom.sorted.03 <- temp
tom.sorted.03 <- convert_to_num(tom.sorted.03)
tom.avg.03 <- calc_mean(tom.sorted.03)
# percentag of NAs
length(which(is.na(tom.avg.03))) / length(tom.avg.03)  # = 10.11%
nas <- (which(is.na(tom.avg.03)))
tom.avg.03[nas] <- NA

library(pcaMethods)
# impute data
pca <- pca(tom.avg.03,method = "nipals")
tom.avg.03 <- pca@completeObs

metabolites <- strsplit(colnames(tom.avg.03),".",T)
mets <- NULL
for(i in 1:length(metabolites)){mets <- rbind(mets,metabolites[[i]][2])}
# checking compound names with tomato compound list
match(tolower(mets),tolower(unlist(tomato.compounds))) # poor results
# use grep instead
for(i in 1:length(mets)){
  comp <- grep(tolower(mets[i]),unlist(tomato.compounds))
  print(paste(i," compound: ",tolower(mets[i])))
  print(comp)
  print(unlist(tomato.compounds)[comp])
} # end for loop
cnames.03 <- c("cis-aconitate","L-alanine","arabinose","L-arginine","L-asparagine","2-oxoglutarate",
               "alpha;-tocopherol","benzoate","beta;-alanine","citramalate","L-cysteine","dehydroascorbate (bicyclic form)",
               "erythritol","fa16_0","fa18_0","fa18_2","fructose_6_phosphate","GDP-beta;-L-fucose","fumarate",
               "gaba","UDP-alpha;-D-galacturonate","D-gluconate","glucose_6_phosphate","D-glycerate",  "Glycerate_3_Phosphate",
               "glycerol","Glycerol_1_Phosphate","glycine","L-homoserine","inositol","inositol_1_phosphate","L-isoleucine",
               "isomaltose","L-leucine","L-lysine","L_Ascorbate",
               "(S)-malate","2-isopropylmaleate","maltitol","maltose","D-mannitol","melezitose","L-methionine",
               "L-phenylalanine","L-proline","putrescine","L-quinate","rhamnose","D-ribofuranose","saccharate", "L-serine","shikimate",
               "sorbitol","succinate","sucrose","Se-methyl-L-selenocysteine","t4_ho_proline","threonate","L-threonine", 
               "alpha;,alpha;-trehalose","L-tryptophan","L-tyrosine","uracil","L-valine","L-aspartate", "citrate","fructose","galactose",
               "glucose","L-glutamate","L-glutamine","D-threo-isocitrate","GDP-alpha;-D-mannose","pyroglutamate","phosphate")

special_compounds <- c("GDP-beta;-L-galactose","UDP-alpha;-D-galactose","alpha;-D-galactose","ADP-alpha;-D-glucose",
                       "GDP-alpha;-D-glucose","UDP-alpha;-D-glucose","alpha;-D-glucose","beta;-D-glucose","UDP-L-rhamnose",
                       "dTDP-beta;-L-rhamnose" )
omit_compounds <- read.delim(file.choose(),check.names=F)

# checking again
length(which(is.na(match(cnames.03,unlist(tomato.compounds)))))  # = 26 NAs - sugars = 23 NAs out of 75 compounds
cnames.03[which(is.na(match(cnames.03,unlist(tomato.compounds))))]

#=========================================================================================================================

IL.2003.pathways <- find.pathways(tomato.pathways,special.compounds = special_compounds,omit.compounds = omit_compounds,cnames = cnames.03)

#=========================================================================================================================



# for 2004 dataset
# reading in datafiles tom.01 = fruit pericarp metabolites of 2001 etc.
tom.04 <- read.delim(file.choose(),check.names=F)
temp <- tom.04[order(tom.04[,1]),]
temp <- as.matrix(temp)
rownames(temp) <- temp[,1] 
temp <- temp[,-1]
tom.sorted.04 <- temp
tom.sorted.04 <- convert_to_num(tom.sorted.04)
tom.avg.04 <- calc_mean(tom.sorted.04)
# percentag of NAs
length(which(is.na(tom.avg.04))) / length(tom.avg.04)  # = 30.31%
nas <- (which(is.na(tom.avg.04)))
tom.avg.04[nas] <- NA

tom.avg.04 <- remove_na_col(tom.avg.04)
tom.avg.04 <- remove_na_row(tom.avg.04)
length(which(is.na(tom.avg.04))) / length(tom.avg.04)  # = 0 %
library(pcaMethods)
# impute data
pca <- pca(tom.avg.04,method = "nipals")
tom.avg.04 <- pca@completeObs

metabolites <- strsplit(colnames(tom.avg.04),".",T)
mets <- NULL
for(i in 1:length(metabolites)){mets <- rbind(mets,metabolites[[i]][2])}
# checking compound names with tomato compound list
match(tolower(mets),tolower(unlist(tomato.compounds))) # poor results
# use grep instead
for(i in 1:length(mets)){
  comp <- grep(tolower(mets[i]),unlist(tomato.compounds))
  print(paste(i," compound: ",tolower(mets[i])))
  print(comp)
  print(unlist(tomato.compounds)[comp])
} # end for loop
# cnames.04 <- c("cis-aconitate","L-alanine","arabinose","L-arginine","L-asparagine","2-oxoglutarate",
#                "alpha;-tocopherol","benzoate","beta;-alanine","citramalate","L-cysteine","dehydroascorbate (bicyclic form)",
#                "erythritol","fa16_0","fa18_0","fa18_2","fructose_6_phosphate","GDP-beta;-L-fucose","fumarate",
#                "gaba","UDP-alpha;-D-galacturonate","D-gluconate","glucose_6_phosphate","D-glycerate",  "Glycerate_3_Phosphate",
#                "glycerol","Glycerol_1_Phosphate","glycine","L-homoserine","inositol","inositol_1_phosphate","L-isoleucine",
#                "isomaltose","L-leucine","L-lysine","L_Ascorbate",
#                "(S)-malate","2-isopropylmaleate","maltitol","maltose","D-mannitol","melezitose","L-methionine",
#                "L-phenylalanine","L-proline","putrescine","L-quinate","rhamnose","D-ribofuranose","saccharate", "L-serine","shikimate",
#                "sorbitol","succinate","sucrose","Se-methyl-L-selenocysteine","t4_ho_proline","threonate","L-threonine", 
#                "alpha;,alpha;-trehalose","L-tryptophan","L-tyrosine","uracil","L-valine","L-aspartate", "citrate","fructose","galactose",
#                "glucose","L-glutamate","L-glutamine","D-threo-isocitrate","GDP-alpha;-D-mannose","pyroglutamate","phosphate")



cnames.04 <- c("adenosine_5_phsophate" , "L-alanine" ,  "L-arginine" , "L-asparagine" , "L-aspartate" , "alpha;-tocopherol",
               "benzoate" , "gaba" , "beta;-alanine" ,"calystegine A3" , "calystegine B2" , "citramalte" , "citrate",
               "L-cysteine" , "dehydroascorbate (bicyclic form)" , "erythritol" , "fructose" , "frucoste_6_phosphate" , 
               "GDP-beta;-L-fucose" , "fumarate" , "galactinol" , "1" , "galactose" , "UDP-alpha;-D-galacturonate" ,
               "D-gluconate" , "glucose" , "glucose_6_phosphate" , "L-glutamate" , "L-glutamine" , "2-oxoglutarate" , 
               "D-glycerate" , "glycerate_3_phosphate" , "glycerol" , "glycerol_1_phosphate" , "glycine" , "guanine" , 
               "L-histidine" , "L-homoeserine" , "inositol" , "D-threo-isocitrate" , "L-isoleucine" , "L-leucine" ,
               "L-lysine" , "L-ascrobate" , "(S)-malate" , "maltitol" , "maltose" , "GDP-alpha;-D-mannose" , 
               "L-methionine" ,  "nicotinate" , "nicotinamid" , "octadeconate" , "L-phenylalanine" , "phenylamin" , 
               "phosphate" , "L-proine" , 
               "putrescine" , "L-pyruvate" , "L-quinate" , "raffinose" , "riobse" , "saccharate" , "salyicylate" ,
               "L-serine" , "shikimate" , "succinate" , "sucrose" , "suqalene" , "threonate" , "L-threonine" , 
               "alpha;,alpha;-trehalose" , "L-tryptophan" , "L-tyrosine" , "uracil" , "urea" , "L-valine" , "xylital",
               "UDP-alpha;-D-xylose")


special_compounds <- c("GDP-beta;-L-galactose","UDP-alpha;-D-galactose","alpha;-D-galactose","ADP-alpha;-D-glucose",
                       "GDP-alpha;-D-glucose","UDP-alpha;-D-glucose","alpha;-D-glucose","beta;-D-glucose","UDP-L-rhamnose",
                       "dTDP-beta;-L-rhamnose" )
omit_compounds <- read.delim(file.choose(),check.names=F)

# checking again
length(which(is.na(match(cnames.04,unlist(tomato.compounds)))))  # = 31 NAs - sugars = 28 NAs out of 75 compounds
cnames.04[which(is.na(match(cnames.04,unlist(tomato.compounds))))]

#=========================================================================================================================

IL.2004.pathways <- find.pathways(tomato.pathways,special.compounds = special_compounds,omit.compounds = omit_compounds,cnames = cnames.04)

#=========================================================================================================================

#06.09.2015
all.season.pathways <- vector("list",3)
all.season.pathways[[1]] <- IL.2001.pathways
all.season.pathways[[2]] <- IL.2003.pathways
all.season.pathways[[3]] <- IL.2004.pathways

overlap.indices <- sync.seasons(pathway.seasons = all.season.pathways)  # this the index variable to be used for network analysis
overlap.ind <- sync.seasons(pathway.seasons = all.season.pathways)

two.season.pathways <- vector("list",2)
two.season.pathways[[1]] <- IL.2001.pathways
two.season.pathways[[2]] <- IL.2004.pathways

overlap.indices.03.04 <- sync.seasons(pathway.seasons = two.season.pathways)
na <- which(is.na(match(1:575,overlap.indices)))
# test
for(i in na){
  for(j in 1:3){
    print(paste(j,": ",i))
    print(all.season.pathways[[j]][[i]]$overlap)
  }
}



#==============================================================================================================
# 06.09.2015
# finding pathway features for pathways, such as enzymes, genes, reactants, reactions etc.
all.pathways.feature.list <- pathway_features() 
counter <- 0
for(i in 1: length(all.pathways.feature.list)){
  x <- which(all.pathways.feature.list[[i]]$species=="Solanum lycopersicum")
  if(length(x)>0){
    counter <- counter+1
  }
}
print(counter)
counter <- 0
for(i in 1: length(tomato.pathways)){
  x <- which(tomato.pathways[[i]]$species=="Solanum lycopersicum")
  if(length(x)>0){
    counter <- counter+1
  }
}
print(counter)
all.pathway.feature.list <- pathway_features()
counter <- 0
for(i in 1:length(all.pathways.feature.list)){
  if(length(which(all.pathways.feature.list[[i]]$species=="Solanum lycopersicum"))>0){
    print(counter)
    counter <- counter+1
  } 
}

print(counter); print(i)

#06.09.2015
# tom 2001
tomato.pathway.feature.list <- species.feature.list(all.pathways.feature.list,"Solanum lycopersicum",overlap.indices)

tomato.output.variable <- species.feature.table(tomato.pathway.feature.list)

write.csv(tomato.output.variable,"tom_pathway_features_060915.csv")

# network construction for tom.01, tom.03, tom.04 fruit metabolite set
colnames(tom.avg.01) <- cnames.01
netw.cor.01 <- correlation.thresholds(tom.avg.01,0.1,1)
# p.value = 0.01 and r.value 0.3
net.thresh.01 <- correlation.network.thresh(netw.cor.01@cor.matrix,netw.cor.01@p.value.mat,0.3,0.01)
tom.01.graph <- create.igraph(net.thresh.01)
vcount(tom.01.graph); ecount(tom.01.graph)

# tom 2003

colnames(tom.avg.03) <- cnames.03
netw.cor.03 <- correlation.thresholds(tom.avg.03,0.1,1)
# p.value = 0.01 and r.value 0.3
net.thresh.03 <- correlation.network.thresh(netw.cor.03@cor.matrix,netw.cor.03@p.value.mat,0.3,0.01)
tom.03.graph <- create.igraph(net.thresh.03)
vcount(tom.03.graph); ecount(tom.03.graph)

# tom 2004

colnames(tom.avg.04) <- cnames.04
netw.cor.04 <- correlation.thresholds(tom.avg.04,0.1,1)
# p.value = 0.01 and r.value 0.3
net.thresh.04 <- correlation.network.thresh(netw.cor.04@cor.matrix,c,0.3,0.01)
tom.04.graph <- create.igraph(net.thresh.04)
vcount(tom.04.graph); ecount(tom.04.graph)


########################################################################################################
########################################################################################################
# find negative controls, which are not part of the tomato metabolom
# full collection = all.pathways.feature.list

count <-0 
negative.control <- vector("list",0)
for(i in 1:length(all.pathways.feature.list)){
 
  if(length(which(all.pathways.feature.list[[i]]$species=="Solanum lycopersicum")>0)){
        next}
      
  count <- count+1
  negative.control[[count]] <- all.pathways.feature.list[[i]]
  names(negative.control[[count]]) <- c("pathway","compounds_of_pathway", "enzymes_of_pathway","genes_of_pathway",
                                        "products_of_pathway","super_pathway","sub_pathways","species",
                                        "enzymes_not_used","reactants_of_pathway", "reaction_of_pathway")
     
} # end for loop

# combine cnames from tomato for nontomato
cnames <- c(cnames.01, cnames.03, cnames.04)
# a compound needs to exist in all three datasets
# delete the ones who don't
cnames.species <- NULL
while(length(cnames>=3)){
  x <- which(tolower(cnames)==tolower(cnames[[1]]))
  if(length(x)>=3){
    cnames.species <- c(cnames.species,cnames[[1]])
    cnames <- cnames[-x]
  }
  else
    cnames <- cnames[-x]
}
# delete redundants
# cnames.species <- unique(cnames.species)
non.tomato.pathways <- find.pathways(negative.control,special_compounds,omit_compounds,cnames.species)
sync.non.tomato.pathways.index <- sync.pathways.intern(non.tomato.pathways) # returns index variable which ones to delete
synced.non.tomato.list <- vector("list",0)
counter <- 0
for(i in sync.non.tomato.pathways.index){
  counter <- counter+1
  synced.non.tomato.list[[counter]] <- non.tomato.pathways[i] 
}
sync.tomato.non.tomato.pathways <-cross.reference.pathways.sync(all.season.pathways,overlap.indices,synced.non.tomato.list)
# make a proper table for the negative control with function pathway.features and find the best negative controls
# this contains only part of the features
synced.tomato.non.tomato.list <- vector("list",0)
counter <- 0
for(i in sync.tomato.non.tomato.pathways){
  counter <- counter+1
  synced.tomato.non.tomato.list[[counter]] <- synced.non.tomato.list[i] 
}
# therefore create a list containing all features
synced.tomato.non.tomato.full.feature.list <- vector("list",0)
counter <- 0
for(i in 1:length(synced.tomato.non.tomato.list)){
  for(j in 1:length(all.pathways.feature.list)){
    if(as.character(synced.tomato.non.tomato.list[[i]][[1]][[1]]$pathway) == as.character(all.pathway.feature.list[[j]]$pathway)){
      counter <- counter+1
      synced.tomato.non.tomato.full.feature.list[[counter]] <- all.pathway.feature.list[[j]]
      break
    }
  }
}


negative.control.feature.table <- species.feature.table(synced.tomato.non.tomato.full.feature.list)


########################################################################################################
########################################################################################################
# 06.09.2015
# combine here variables: tomato.output.variable (169x15) and negative.control.feature.table (33x15)
# on top of each other for positive and negative controls of pathway features = total = 202x15
all.species.pathway.features <- rbind(tomato.output.variable,negative.control.feature.table)
########################################################################################################
########################################################################################################

# 06.09.2015
library(igraph)
# generating network feature tables
V(tom.01.graph)$names <- rownames(net.thresh.01) 
V(tom.03.graph)$names <- rownames(net.thresh.03) 
V(tom.04.graph)$names <- rownames(net.thresh.04) 
# use those variables: IL.2001.pathways..., overlap.indices
network.feature.colnames <- read.delim(file.choose(),check.names=F,header=F)
gr.features.01 <- matrix(ncol = length(network.feature.colnames [[1]]))
colnames(gr.features.01) <- network.feature.colnames[[1]]
# preparing graph objects
for(i in 1:length(overlap.indices)){
  print(paste (IL.2001.pathways[[overlap.indices[i]]]$pathway," index: ",i,sep=""))
  
  subgraph.node.index <- match(IL.2001.pathways[[overlap.indices[i]]]$overlap,V(tom.01.graph)$names)
  # find appropriate subgraph edge index
  sni <- NULL
  for(j in 1:(length(subgraph.node.index)-1)){
    for(k in (j+1):length(subgraph.node.index)){
      sni <- c(sni,subgraph.node.index[j],subgraph.node.index[k])
    } # end for
  } # end for
  subgraph <- induced.subgraph(tom.01.graph,vids = subgraph.node.index)
  nonzero <- which(get.edge.ids(tom.01.graph,sni)!=0)
  subgraph.edge.index <- get.edge.ids(tom.01.graph,sni)[nonzero]
  if(i==1){
    gr.features.01[i,1:115] <- network_features(tom.01.graph,subgraph,subgraph.node.index,subgraph.edge.index)
    gr.features.01[i,116:148] <- whole_network_features(tom.01.graph,subgraph.node.index)
    rownames(gr.features.01)[i] <- as.character(IL.2001.pathways[[overlap.indices[i]]]$pathway)
  }
  else{
    gr.feat.a <- network_features(tom.01.graph,subgraph,subgraph.node.index,subgraph.edge.index)
    gr.feat.b <-  whole_network_features(tom.01.graph,subgraph.node.index)
    gr.feat <- cbind(gr.feat.a,gr.feat.b)
    gr.features.01 <- rbind(gr.features.01,gr.feat)
                          #  network_features(tom.01.graph,subgraph,subgraph.node.index,subgraph.edge.index))
    rownames(gr.features.01)[i] <- as.character(IL.2001.pathways[[overlap.indices[i]]]$pathway)
  }
} # end for loop
###############################################################################################################
# for 2003
gr.features.03 <- matrix(ncol = length(network.feature.colnames [[1]]))
colnames(gr.features.03) <- network.feature.colnames[[1]]
# preparing graph objects

# for 2003
for(i in 1:length(overlap.indices)){
  print(paste (IL.2003.pathways[[overlap.indices[i]]]$pathway," index: ",i,sep=""))
  
  subgraph.node.index <- match(IL.2003.pathways[[overlap.indices[i]]]$overlap,V(tom.03.graph)$names)
  # find appropriate subgraph edge index
  sni <- NULL
  for(j in 1:(length(subgraph.node.index)-1)){
    for(k in (j+1):length(subgraph.node.index)){
      sni <- c(sni,subgraph.node.index[j],subgraph.node.index[k])
    } # end for
  } # end for
  subgraph <- induced.subgraph(tom.03.graph,vids = subgraph.node.index)
  nonzero <- which(get.edge.ids(tom.03.graph,sni)!=0)
  subgraph.edge.index <- get.edge.ids(tom.03.graph,sni)[nonzero]
  if(i==1){
    gr.features.03[i,1:115] <- network_features(tom.03.graph,subgraph,subgraph.node.index,subgraph.edge.index)
    gr.features.03[i,116:148] <- whole_network_features(tom.03.graph,subgraph.node.index)
    rownames(gr.features.03)[i] <- as.character(IL.2003.pathways[[overlap.indices[i]]]$pathway)
  }
  else{
    gr.feat.a <- network_features(tom.03.graph,subgraph,subgraph.node.index,subgraph.edge.index)
    gr.feat.b <-  whole_network_features(tom.03.graph,subgraph.node.index)
    gr.feat <- cbind(gr.feat.a,gr.feat.b)
    gr.features.03 <- rbind(gr.features.03,gr.feat) 
                            #network_features(tom.03.graph,subgraph,subgraph.node.index,subgraph.edge.index))
    rownames(gr.features.03)[i] <- as.character(IL.2003.pathways[[overlap.indices[i]]]$pathway)
  }
} # end for loop
###############################################################################################################
# for 2004
gr.features.04 <- matrix(ncol = length(network.feature.colnames [[1]]))
colnames(gr.features.04) <- network.feature.colnames[[1]]
# preparing graph objects
for(i in 1:length(overlap.indices)){
  print(paste (IL.2004.pathways[[overlap.indices[i]]]$pathway," index: ",i,sep=""))
  
  subgraph.node.index <- match(IL.2004.pathways[[overlap.indices[i]]]$overlap,V(tom.04.graph)$names)
  # find appropriate subgraph edge index
  sni <- NULL
  for(j in 1:(length(subgraph.node.index)-1)){
    for(k in (j+1):length(subgraph.node.index)){
      sni <- c(sni,subgraph.node.index[j],subgraph.node.index[k])
    } # end for
  } # end for
  subgraph <- induced.subgraph(tom.04.graph,vids = subgraph.node.index)
  nonzero <- which(get.edge.ids(tom.04.graph,sni)!=0)
  subgraph.edge.index <- get.edge.ids(tom.04.graph,sni)[nonzero]
  if(i==1){
    gr.features.04[i,1:115] <- network_features(tom.04.graph,subgraph,subgraph.node.index,subgraph.edge.index)
    gr.features.04[i,116:148] <- whole_network_features(tom.04.graph,subgraph.node.index)
    rownames(gr.features.04)[i] <- as.character(IL.2004.pathways[[overlap.indices[i]]]$pathway)
  }
  else{
    gr.feat.a <- network_features(tom.04.graph,subgraph,subgraph.node.index,subgraph.edge.index)
    gr.feat.b <-  whole_network_features(tom.04.graph,subgraph.node.index)
    gr.feat <- cbind(gr.feat.a,gr.feat.b)
    gr.features.04 <- rbind(gr.features.04,gr.feat) 
                           # network_features(tom.04.graph,subgraph,subgraph.node.index,subgraph.edge.index))
    rownames(gr.features.04)[i] <- as.character(IL.2004.pathways[[overlap.indices[i]]]$pathway)
  }
} # end for loop
###############################################################################################################
###############################################################################################################
all.season.network.features <- cbind(gr.features.01,gr.features.03,gr.features.04)
###############################################################################################################
###############################################################################################################
###############################################################################################################
# now add the features calculated for the negative controls in background of the correlation matrices using
# variable synced.tomato.non.tomato.list[[1]][[1]][[1]]$overlap

# something wrong in the overlap 02.08.15 pay attention - has been fixed

# for 2001

neg.ctrl.features.01 <- matrix(ncol = length(network.feature.colnames [[1]]))
colnames(neg.ctrl.features.01) <- network.feature.colnames[[1]]
# preparing graph objects
for(i in 1:length(synced.tomato.non.tomato.list)){
  print(paste (synced.tomato.non.tomato.list[[i]][[1]][[1]]$pathway," index: ",i,sep=""))
  
  subgraph.node.index <- match(synced.tomato.non.tomato.list[[i]][[1]][[1]]$overlap,V(tom.01.graph)$names)
#   nas <- which(is.na(subgraph.node.index))
#   if(length(nas>0))
#     subgraph.node.index <- subgraph.node.index[-nas]
  # find appropriate subgraph edge index
  sni <- NULL
  for(j in 1:(length(subgraph.node.index)-1)){
    for(k in (j+1):length(subgraph.node.index)){
      sni <- c(sni,subgraph.node.index[j],subgraph.node.index[k])
    } # end for
  } # end for
  subgraph <- induced.subgraph(tom.01.graph,vids = subgraph.node.index)
  nonzero <- which(get.edge.ids(tom.01.graph,sni)!=0)
  subgraph.edge.index <- get.edge.ids(tom.01.graph,sni)[nonzero]
  if(i==1){
    neg.ctrl.features.01[i,1:115] <- network_features(tom.01.graph,subgraph,subgraph.node.index,subgraph.edge.index)
    neg.ctrl.features.01[i,116:148] <- whole_network_features(tom.01.graph,subgraph.node.index)
    #neg.ctrl.features.01[i,] <- network_features(tom.01.graph,subgraph,subgraph.node.index,subgraph.edge.index)
    rownames(neg.ctrl.features.01)[i] <- as.character(synced.tomato.non.tomato.list[[i]][[1]][[1]]$pathway)
  }
  else{
    gr.feat.a <- network_features(tom.01.graph,subgraph,subgraph.node.index,subgraph.edge.index)
    gr.feat.b <-  whole_network_features(tom.01.graph,subgraph.node.index)
    gr.feat <- cbind(gr.feat.a,gr.feat.b)
    neg.ctrl.features.01  <- rbind(neg.ctrl.features.01 ,gr.feat) 
    #neg.ctrl.features.01 <- rbind(neg.ctrl.features.01, 
    #        network_features(tom.01.graph,subgraph,subgraph.node.index,subgraph.edge.index))
    rownames(neg.ctrl.features.01)[i] <- as.character(synced.tomato.non.tomato.list[[i]][[1]][[1]]$pathway)
  }
} # end for loop
###############################################################################################################

# for 2003

neg.ctrl.features.03 <- matrix(ncol = length(network.feature.colnames [[1]]))
colnames(neg.ctrl.features.03) <- network.feature.colnames[[1]]
# preparing graph objects
for(i in 1:length(synced.tomato.non.tomato.list)){
  print(paste (synced.tomato.non.tomato.list[[i]][[1]][[1]]$pathway," index: ",i,sep=""))
  
  subgraph.node.index <- match(synced.tomato.non.tomato.list[[i]][[1]][[1]]$overlap,V(tom.03.graph)$names)
  #   nas <- which(is.na(subgraph.node.index))
  #   if(length(nas>0))
  #     subgraph.node.index <- subgraph.node.index[-nas]
  # find appropriate subgraph edge index
  sni <- NULL
  for(j in 1:(length(subgraph.node.index)-1)){
    for(k in (j+1):length(subgraph.node.index)){
      sni <- c(sni,subgraph.node.index[j],subgraph.node.index[k])
    } # end for
  } # end for
  subgraph <- induced.subgraph(tom.03.graph,vids = subgraph.node.index)
  nonzero <- which(get.edge.ids(tom.03.graph,sni)!=0)
  subgraph.edge.index <- get.edge.ids(tom.03.graph,sni)[nonzero]
  if(i==1){
    neg.ctrl.features.03[i,1:115] <- network_features(tom.03.graph,subgraph,subgraph.node.index,subgraph.edge.index)
    neg.ctrl.features.03[i,116:148] <- whole_network_features(tom.03.graph,subgraph.node.index)
    #neg.ctrl.features.03[i,] <- network_features(tom.03.graph,subgraph,subgraph.node.index,subgraph.edge.index)
    rownames(neg.ctrl.features.03)[i] <- as.character(synced.tomato.non.tomato.list[[i]][[1]][[1]]$pathway)
  }
  else{
    gr.feat.a <- network_features(tom.03.graph,subgraph,subgraph.node.index,subgraph.edge.index)
    gr.feat.b <-  whole_network_features(tom.03.graph,subgraph.node.index)
    gr.feat <- cbind(gr.feat.a,gr.feat.b)
    neg.ctrl.features.03  <- rbind(neg.ctrl.features.03 ,gr.feat) 
    #neg.ctrl.features.03 <- rbind(neg.ctrl.features.03, 
    #                              network_features(tom.03.graph,subgraph,subgraph.node.index,subgraph.edge.index))
    rownames(neg.ctrl.features.03)[i] <- as.character(synced.tomato.non.tomato.list[[i]][[1]][[1]]$pathway)
  }
} # end for loop
###############################################################################################################
# for 2004

neg.ctrl.features.04 <- matrix(ncol = length(network.feature.colnames [[1]]))
colnames(neg.ctrl.features.04) <- network.feature.colnames[[1]]
# preparing graph objects
for(i in 1:length(synced.tomato.non.tomato.list)){
  print(paste (synced.tomato.non.tomato.list[[i]][[1]][[1]]$pathway," index: ",i,sep=""))
  
  subgraph.node.index <- match(synced.tomato.non.tomato.list[[i]][[1]][[1]]$overlap,V(tom.04.graph)$names)
  #   nas <- which(is.na(subgraph.node.index))
  #   if(length(nas>0))
  #     subgraph.node.index <- subgraph.node.index[-nas]
  # find appropriate subgraph edge index
  sni <- NULL
  for(j in 1:(length(subgraph.node.index)-1)){
    for(k in (j+1):length(subgraph.node.index)){
      sni <- c(sni,subgraph.node.index[j],subgraph.node.index[k])
    } # end for
  } # end for
  subgraph <- induced.subgraph(tom.04.graph,vids = subgraph.node.index)
  nonzero <- which(get.edge.ids(tom.04.graph,sni)!=0)
  subgraph.edge.index <- get.edge.ids(tom.04.graph,sni)[nonzero]
  if(i==1){
    neg.ctrl.features.04[i,1:115] <- network_features(tom.04.graph,subgraph,subgraph.node.index,subgraph.edge.index)
    neg.ctrl.features.04[i,116:148] <- whole_network_features(tom.04.graph,subgraph.node.index)
   # neg.ctrl.features.04[i,] <- network_features(tom.04.graph,subgraph,subgraph.node.index,subgraph.edge.index)
    rownames(neg.ctrl.features.04)[i] <- as.character(synced.tomato.non.tomato.list[[i]][[1]][[1]]$pathway)
  }
  else{
    gr.feat.a <- network_features(tom.04.graph,subgraph,subgraph.node.index,subgraph.edge.index)
    gr.feat.b <-  whole_network_features(tom.04.graph,subgraph.node.index)
    gr.feat <- cbind(gr.feat.a,gr.feat.b)
    neg.ctrl.features.04  <- rbind(neg.ctrl.features.04 ,gr.feat) 
   # neg.ctrl.features.04 <- rbind(neg.ctrl.features.04, 
   #                                network_features(tom.04.graph,subgraph,subgraph.node.index,subgraph.edge.index))
    rownames(neg.ctrl.features.04)[i] <- as.character(synced.tomato.non.tomato.list[[i]][[1]][[1]]$pathway)
  }
} # end for loop
###############################################################################################################
###############################################################################################################
###############################################################################################################
all.season.negative.network.features <- cbind(neg.ctrl.features.01,neg.ctrl.features.03,neg.ctrl.features.04)
all.season.positive.negative.features <- rbind(all.season.network.features,all.season.negative.network.features)
# change colnames
colnames(all.season.positive.negative.features)[1:148] <- paste(colnames(all.season.positive.negative.features)[1:148],".01",sep="")
colnames(all.season.positive.negative.features)[149:296] <- paste(colnames(all.season.positive.negative.features)[149:296],".03",sep="")
colnames(all.season.positive.negative.features)[297:444] <- paste(colnames(all.season.positive.negative.features)[297:444],".04",sep="")
# network features and pathway features in one dataset
all.features.dataset <- cbind(all.season.positive.negative.features,all.species.pathway.features)
write.csv(all.features.dataset,"IL010304_feature_dataset_060915.csv")
###############################################################################################################
###############################################################################################################
###############################################################################################################
# 07.09.2015
# writing out the percentage of pathway coverage
tom.percentages  <- NULL
for(i in overlap.indices){
  tom.percentages <- rbind(tom.percentages,all.season.pathways[[1]][[i]]$percentage)
}
write.csv(tom.percentages,"tompercentages_070915.csv")
tom.absolute <- NULL
for(i in overlap.indices){
    tom.absolute <- rbind(tom.absolute,length(all.season.pathways[[1]][[i]]$overlap))
}
neg.percentages <- NULL
for(i in 1:length(synced.tomato.non.tomato.list)){
  neg.percentages <- rbind(neg.percentages,synced.tomato.non.tomato.list[[i]][[1]][[1]]$percentage)
}
neg.absolute <- NULL
for(i in 1:length(synced.tomato.non.tomato.list)){
    neg.absolute <- rbind(neg.absolute,length(synced.tomato.non.tomato.list[[i]][[1]][[1]]$overlap))
}
write.csv(neg.percentages,"negativecontrol_percentages_070915.csv")