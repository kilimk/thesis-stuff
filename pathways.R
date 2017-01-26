# 06.09.2015
# changes made. not deleting pathways that are the same - we want this in machine learning; function: find.pathways, lines: 168-187
# 25.06.2015
# Author: David Toubiana
# function that finds all compounds of a pathway witin a ***Cyc

pathways <- function(species){
  
  data <- read.delim(file.choose(),check.names=F)
  pwy <- NULL
  # make list objects of compounds of desired species
  # for all species
  if (tolower(species)=="all"){
    pwy <- vector("list",0)
    for(i in 1:nrow(data)){
      compounds <- data[i,"Compounds of pathway"]
      compounds <- as.character(compounds)
      compounds <- strsplit(compounds," // ")
      if(length(compounds[[1]])>1){
          temp <- compounds[[1]][2:length(compounds[[1]])]
          temp <- unlist(strsplit(temp,"\""))
          temp <- temp[which(temp!="")]
      } # end if
      compounds <- compounds[[1]][1]
      if(length(temp)>0)
        compounds <- c(compounds,temp)
      
      pwy[[i]] <- list(data[i,"Pathways"],compounds,species)
      names(pwy[[i]]) <- c("pathway","compounds","species")
    } # end for
  } # end if
  else{ # specific species
    spe <- grep(tolower(species),tolower(as.character(data[,"Species"])))
    if(length(spe)==0){
      print("Species not found")
      return 
    }
    pwy <- vector("list",0)
    count <- 1  # counter
    for(i in spe){
      compounds <- data[i,"Compounds of pathway"]
      compounds <- as.character(compounds)
      compounds <- strsplit(compounds," // ")
      if(length(compounds[[1]])>1){
        temp <- compounds[[1]][2:length(compounds[[1]])]
        temp <- unlist(strsplit(temp,"\""))
        temp <- temp[which(temp!="")]
      } # end if
      compounds <- compounds[[1]][1]
      if(length(temp)>0)
        compounds <- c(compounds,temp)
     # print(paste("i: ",i))
     # print(paste("count: ",count))
      pwy[[count]] <- list(data[i,"Pathways"],compounds,species)
      names(pwy[[count]]) <- c("pathway","compounds","species")
      count <- count+1
    } # end for
  } # end else
  return (pwy)
} # end function

# function listing all compounds of a given collection of pathways
compounds <- function(pathways){
 
  all.compounds <- NULL
  # pathways can be a vector with different references to pathways
  for(i in 1:length(pathways)){
    pwy <- which(ls(envir=.GlobalEnv)==pathways[i])
    if(length(pwy)<1){
      print("No such pathway variable found!")
      return
    } # end if
    # now browse through the pathway file
    pwy <- get(ls(envir=.GlobalEnv)[pwy])
    for(j in 1:length(pwy)){
      cmp <- match(pwy[[j]]$compounds,all.compounds)
      nas <- which(is.na(cmp))
      if(length(nas)>0) { # found some compound(s) that were not listed yet
        all.compounds <- c(all.compounds,pwy[[j]]$compounds[nas])
      } # end if
    } # end nested for loop
  } # end for loop
  all.compounds <- sort(all.compounds)
  return(all.compounds)
} # end function
################################################################################################################################################################
################################################################################################################################################################


# function that will identify compounds of pathways within the tomato metabolite dataset
# it determines the % of overlapp
# it omits the compounds in the omit_compound variable and pays spececial attention to the
# compounds in the variabel special_compounds
find.pathways <- function(pathway.list,special.compounds,omit.compounds,cnames){
  # initiate return variable
  pathway.overlap <- vector("list",0)
  
  for(i in 1:length(pathway.list)){
    if(i==1){
      actual.compounds  <- tolower(pathway.list[[i]]$compounds) 
      # compounds to omit
      omit <- match(actual.compounds,tolower(as.character(omit.compounds[[1]])))
      omit <- which(!(is.na(omit)))
      # special compound to change names of
      special <- match(tolower(actual.compounds),tolower(as.character(special.compounds)))
      special <- which(!(is.na(special)))
      sugar <- which(tolower(as.character(special.compounds))==actual.compounds[special])
      # identify sugar
      if(length(sugar)>0){
        for(j in 1:length(sugar)){
          if(sugar[j]<4) actual.compounds[special[j]]  <- "galactose" # indicating galactose
          else if (sugar[j]>=4 & sugar[j]<=8) actual.compounds[special[j]]  <- "glucose"
          else actual.compounds[special[j]]  <- "rhamnose"
        } # end for
      } # end if
      # now delete compounds to omit
      actual.compounds <- actual.compounds[-omit]
      overlap <- match(tolower(cnames),tolower(actual.compounds))
      overlap <- which(!is.na(overlap))
      percentage <- NULL
      if(length(overlap)>0){
        percentage <- length(overlap)/length(actual.compounds) 
        overlap <- cnames[overlap]
      } # end if
      else{
        overlap <- "no match"
        percentage <- 0
      } # end else
      pathway.overlap[[i]] <- list(pathway.list[[i]]$pathway, pathway.list[[i]]$compounds, pathway.list[[i]]$species,
                                   actual.compounds, overlap, percentage)
      names(pathway.overlap[[i]]) <- c("pathway","compounds","species","actual.compounds","overlap","percentage")
      #print(i)
    } # end if
    else{
      actual.compounds  <- tolower(pathway.list[[i]]$compounds) 
      # compounds to omit
      omit <- match(actual.compounds,tolower(as.character(omit.compounds[[1]])))
      omit <- which(!(is.na(omit)))
      # special compound to change names of
      special <- match(tolower(actual.compounds),tolower(as.character(special.compounds)))
      special <- which(!(is.na(special)))
      sugar <- which(tolower(as.character(special.compounds))==actual.compounds[special])
      # identify sugar
      if(length(sugar)>0){
        for(j in 1:length(sugar)){
          if(sugar[j]<4) actual.compounds[special[j]]  <- "galactose" # indicating galactose
          else if (sugar[j]>=4 & sugar[j]<=8) actual.compounds[special[j]]  <- "glucose"
          else actual.compounds[special[j]]  <- "rhamnose"
        } # end for
      } # end if
      # now delete compounds to omit
      actual.compounds <- actual.compounds[-omit]
      overlap <- match(tolower(cnames),tolower(actual.compounds))
      overlap <- which(!is.na(overlap))
      percentage <- NULL
      if(length(overlap)>0){
        percentage <- length(overlap)/length(actual.compounds) 
        overlap <- cnames[overlap]
      }
      else{
        overlap <- "no match"
        percentage <- 0
      }
      # compare to other list entries
      exist <- NULL
      # 06.09.2015
      # omitting this part, since we want to have doubles of same entries in machine learning
#       if(overlap[1] != "no match"){ 
#         for(k in 1:(i-1)){
#           e1 <- match(overlap,pathway.overlap[[k]]$overlap)
#           e1  <- which(is.na(e1))
#           e2 <- match(pathway.overlap[[k]]$overlap,overlap)
#           e2  <- which(is.na(e2))
#           if(length(e1)==0 & length(e2)==0 ) { # it means this constellation of compounds already exists
#             exist <- c(exist,k)
#           } # end if
#         } # end for
#         if(length(exist)>0){
#           num <- strsplit(as.character(exist),"")
#           str <- NULL
#           for(l in 1:length(num)){
#             if(l==1) str <- paste(num[[l]],"", sep="")
#             else str <- paste(str,", ", num[[l]], sep="")
#           } # end for
#           overlap <- paste("constellation already exisit in pathways: ", str, sep="")
#         }
#       } # end if
     # print(paste("ok",i))
      pathway.overlap[[i]] <- list(pathway.list[[i]]$pathway, pathway.list[[i]]$compounds, pathway.list[[i]]$species,
                                   actual.compounds, overlap, percentage)
      names(pathway.overlap[[i]]) <- c("pathway","compounds","species","actual.compounds","overlap","percentage")
      #print(i)
    } # end else
  } # end for loop
  return(pathway.overlap)
} # end function
################################################################################################################################################################
################################################################################################################################################################


sync.seasons <- function(pathway.seasons){
  # function that identifies wheter a patwhay was identified in every season
  # uses the list variables $actual.compounds and &overlap, where the length of actual.compounds must be >=2
  # and overlap cannot be "no match" nor "constellation already exisit in pathways: "
  # returns a vector pointing to indices where all seasons have vavlid entries
  
    index  <- NULL
    count <- 1
    #for (i in 1:10){
    for (i in 1:length(pathway.seasons[[1]])){
        for(j in 1:length(pathway.seasons)){
          print(paste(i," : ",j))
          if(length(pathway.seasons[[j]][[i]]$overlap) < 2 | pathway.seasons[[j]][[i]]$overlap == "no match" |
               grepl("constellation already exisit in pathways: ", pathway.seasons[[j]][[i]]$overlap) == TRUE){
            count <- 1
            break
          } # end if
          else{
            if (count<length(pathway.seasons)){
              count <- count+1
              next
            } # end if
            else{
              count <- 1
              index <- c(index,i)
            } # end else
          } # end else
      } # end for
  } # end for
  return(index)
} # end function
################################################################################################################################################################
################################################################################################################################################################
sync.pathways.intern <- function(pathways){
  # function that identifies whether there are multiple pathways within in the same 
  # collection of pathways
  # uses the list variables $actual.compounds and &overlap, where the length of actual.compounds must be >=2
  # and overlap cannot be "no match" nor "constellation already exisit in pathways: "
  # returns a vector pointing to indices where all seasons have vavlid entries
  
  index  <- NULL
  zero.index <- NULL
  #for (i in 1:10){
  for (i in 1:length(pathways)){
    
      #print(paste(i," :"))
      if(length(pathways[[i]]$overlap) < 2 | pathways[[i]]$overlap == "no match" |
           grepl("constellation already exisit in pathways: ", pathways[[i]]$overlap) == TRUE){
        zero.index <- c(zero.index,i) # invalid entries
        next
      } # end if
      else{
        index <- c(index,i) # valid entries
      } # end else
  } # end for
  print(length(zero.index))
  print(length(index))
  
  return(index)
} # end function
################################################################################################################################################################
################################################################################################################################################################
cross.reference.pathways.sync <- function(reference.pathways,reference.index,check.pathway){
  # function that identifies whether there a pathway to be checked also exist in the reference pathway set
  # it receives valid entries as parameters
  # uses the list variables $actual.compounds and &overlap, where the length of actual.compounds must be >=2
  # returns a vector pointing to indices where all indices have vavlid entries 
  
  index  <- NULL
  for (i in 1:length(check.pathway)){
    for (j in reference.index) {
      #print(i)
      #print(j)
      x <- match(check.pathway[[i]][[1]]$overlap, reference.pathways[[1]][[j]]$overlap)
#       print(paste("check: ", check.pathway[[i]][[1]]$overlap)); 
#       print(paste("reference: ",reference.pathways[[1]][[j]]$overlap))
        if (j < max(reference.index)){
          if(length(which(is.na(x)))!=0){
            next
          } #end if
          else
            break
        } # end if
        else
          index <- c(index,i)
      
    } # end for
  } # end for
  
  print(length(index))
  
  return(index)
} # end function







