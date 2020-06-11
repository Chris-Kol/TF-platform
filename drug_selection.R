drug_selection <- function(drugs, distance_path, all, one_by_one,thresh1,thresh2){
  # This function uses the output of the deepsiba predictions and the drug finder function, to determine which cold drugs will be effective
  # The input is a vector of drugs' sig_ids to be used found from the drug_finder function + the path to the distance prediction file
  # If you want cold_drugs that are close to all drugs, all should be 1 and one_by_one should be 0.
  # If you want to see the closest cold_drugs for each drug, one_by_one should be 1 and all should be 0.
  # Thresh1 and Thresh2 are the lower and upper limits of the distance values for the search.
############################################################################################################################################
 
   distance <- read.csv(paste0(distance_path))
  #distance <- read.csv('df_cold_a375.csv')
  indx1 <- which(distance$sig_id.x %in% drugs)
  indx2 <- which(distance$sig_id.y %in% drugs)
  indx <- c(indx1,indx2)
  distance <- distance[indx,]  
  if (all==0 & one_by_one==1){
    library(tidyverse)
    a <- vector("list",1)
  for (i in 1:length(drugs)){
    indx1 <- which(distance$sig_id.x == drugs[i])
    indx2 <- which(distance$sig_id.y == drugs[i])
    indx <- c(indx1,indx2)
    dist <- distance[indx,]
    dist <- dist %>% dplyr::select(sig_id.x, sig_id.y,mu) %>% dplyr:: filter(mu<thresh2 & mu>thresh1)
    dist[which(!(dist$sig_id.x %in% drugs)), c("sig_id.x", "sig_id.y")] <- dist[which(!(dist$sig_id.x %in% drugs)), c("sig_id.y", "sig_id.x")]
    colnames(dist) <- c("train_drug", "test_drug", "value")
    print(paste("dist", drugs[i], sep="__"))
    print(dist)
    a[[i]]<- dist
  }
    return(a)} 
  else if(all==1 & one_by_one==0){
    library(tidyverse)
    all_sigs <- NULL
    for (i in 1:length(drugs)){
      indx1 <- which(distance$sig_id.x == drugs[i])
      indx2 <- which(distance$sig_id.y == drugs[i])
      indx <- c(indx1,indx2)
      dist <- distance[indx,]
      dist <- dist %>% dplyr::select(sig_id.x, sig_id.y,mu) %>% dplyr:: filter(mu<thresh2 & mu>thresh1)
      sigs_x <- unique(as.character(dist$sig_id.x))
      sigs_y <- unique(as.character(dist$sig_id.y))
      sigs <- c(sigs_x,sigs_y)
      sig_ids <- unique(sigs)
      all_sigs <- c(all_sigs,sig_ids)
    }
    all_sigs <- all_sigs[-which(all_sigs %in% drugs)]
    good <- names(table(all_sigs)[which(table(all_sigs) ==length(drugs))])
    print("The drugs that have a close distance to both drugs are:")
    print(good)
    
  return(good)
    }
  }
###example
#drugs<-drugos
#distance_path <- 'df_cold_a375.csv'
#all <-1
#one_by_one <-0
#thresh1<-0
#thresh2<-0.22
#drugas<-drug_selection(drugs, distance_path, all, one_by_one,thresh1,thresh2)
###end of example