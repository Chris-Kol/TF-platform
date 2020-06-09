check_dist_rank <- function(path, thres) {
  
  library(readr)
  dataframe <- readRDS(path)
  cmap_data_drugs = readRDS("C:/Users/Christina/Desktop/BioSysLab/train_genes.rds")       #output of get_cmap_signatures
  landmark = read_tsv("C:/Users/Christina/Desktop/BioSysLab/TF-platform/TF-platform/preprocessing/preprocessing_data/cmap_landmark_genes.txt")
  
  
  ## BUILD MATRIX check (keep distances between drug-TF < thres)
  check <- data.frame(matrix(0, nrow(dataframe), ncol(dataframe)))
  colnames(check)<-colnames(dataframe)
  rownames(check)<-rownames(dataframe)
  vec <- (1:ncol(dataframe))*0
  
  for (i in 1:nrow(dataframe)) {
    vec <- which(dataframe[i,]<thres) 
    if (length(vec)!=0) {
      for (j in 1:length(vec)){
        check[i,vec[j]] <- dataframe[i,vec[j]]
      }
    }
  }
  rows_to_keep <- 0
  for (i in 1:nrow(dataframe)) {
    if (sum(check[i,])>0) {
      rows_to_keep <- union(rows_to_keep, i)
    }
  }
  rows_to_keep <- rows_to_keep [! rows_to_keep %in% 0]
  check <- as.data.frame(check[rows_to_keep,])
  cols_to_keep <- which(colSums(check)!=0)
  temp_names_r <- rownames(check)
  temp_names_c <- colnames(check)[cols_to_keep]
  check <-as.data.frame(check[,cols_to_keep])
  colnames(check)<-temp_names_c
  rownames(check)<-temp_names_r

  
  
  tf_enrichment <- function(cmap_data, landmark){
    #landmark = landmark annotation
    library(dorothea)
    library(tidyverse)
    library(cmapR)
    library(rhdf5)
    library(CARNIVAL)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    library(viper)
    load(file = system.file("BEST_viperRegulon.rdata",package="CARNIVAL"))
    ### calculate tf enrichment
    TF_cmap <- runDoRothEA(cmap_data, regulon=viper_regulon, confidence_level=c('A','B','C')) # Estimating TF activities
    return(TF_cmap)
  }
  
  
  runDoRothEA<-function(df, regulon, confidence_level=c('A','B','C'), write2file = NULL){
    # library(tidyverse)
    library(dplyr)
    library(purrr)
    library(viper)
    library(tibble)
    library(tidyr)
    names(regulon) <- sapply(strsplit(names(regulon), split = ' - '), head, 1)
    filtered_regulon <- regulon %>%
      map_df(.f = function(i) {
        tf_target = i$tfmode %>%
          enframe(name = "target", value="mor") %>%
          mutate(likelihood = i$likelihood)
      },.id = "tf")  %>%
      separate(tf, into=c("tf", "conf"), sep="_") %>%
      filter(conf %in% confidence_level) %>%
      arrange(tf)%>%
      split(.$tf) %>%
      map(function(dat) {
        tf = dat %>% distinct(tf) %>% pull()
        targets = setNames(dat$mor, dat$target)
        likelihood = dat$likelihood
        list(tfmode =targets, likelihood = likelihood)})
    TF_activities = as.data.frame(viper::viper(eset = df, regulon = filtered_regulon, nes = T, method = 'none', minsize = 4, eset.filter = F))
    if(!is.null(write2file)){write.csv2(TF_activities, file = write2file)}
    return(TF_activities)
  }
  
  scores_TF_drugs <- tf_enrichment(cmap_data_drugs, landmark)
  
  scores_TF_drugs <- as.data.frame(scores_TF_drugs[rownames(scores_TF_drugs)%in% colnames(check) ,])
  scores_TF_drugs <- as.data.frame(scores_TF_drugs[,colnames(scores_TF_drugs)%in% rownames(check)])
  scores_TF_drugs <- as.data.frame(t(scores_TF_drugs))
  ES <- (1:nrow(scores_TF_drugs))*0
  ranks <- (1:nrow(scores_TF_drugs))*0
  to_return <- vector("list", ncol(check))
  names(to_return) <- colnames(check)
  
  for (i in 1:ncol(scores_TF_drugs)) {
    
    ind_names <- c(which(check[,i]!=0))
    ES <- as.data.frame(scores_TF_drugs[ind_names,i])
    colnames(ES) <- colnames(scores_TF_drugs)[i]
    rownames(ES) <- rownames(scores_TF_drugs)[ind_names]
    TF <- colnames(scores_TF_drugs)[i]
    
    ranks <- apply(-ES,2,rank)
    ranks <- ranks/(nrow(ranks))*100
    
    dist <- data.frame(1:nrow(ES ))*0
    name <- paste(colnames(ES))
    
    for (j in 1:nrow(ES)) {
      ind <- (which(rownames(check)==rownames(ES)[j]))
      dist[j,1] <- check[ind,i]
      rownames(dist)[j] <- rownames(ES)[j]
    }
    
    colnames(ranks) <- TF
    colnames(dist) <- TF
    
    
    #tmp <- list(scores=t(ES), distances= t(dist), rank_p = t(ranks))
    tmp <- list(scores=(ES), distances= (dist), rank_p = (ranks))
    to_return[[name]] <- tmp
    
  
    a <- hist(ranks[,TF], main = paste(TF, "Ranks %"), freq = T, plot = T)
    
    op <- par(mar=c(4,15,4,2))
    barplot(ranks[,TF], main = paste(TF, "Ranks %"), names.arg = rownames(ES), horiz = TRUE, las = 2, cex.names = 0.7)
    rm(op)
    
    op <- par(mar=c(4,15,4,2))
    barplot(dist[,TF],  main = paste("DISTANCES -", TF), names.arg = rownames(ES), horiz = TRUE, las = 2, cex.names = 0.7)
    rm(op) 
    
    #op <- par(mar=c(4,15,4,2))
    #barplot(ES[,TF], main = paste(TF, "ES"), names.arg = rownames(ES), horiz = TRUE, las = 2, cex.names = 0.7)
    #rm(op)
  
  }
  
  return(to_return)
}

scores_distances_ranks <- check_dist_rank("C:/Users/Christina/Desktop/BioSysLab/best_dist_tf_drugs_a375.rds", thres=0.25)