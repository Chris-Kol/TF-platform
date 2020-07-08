available_tfs_line <- function(tf_kd_enrichment,kd_df,cell,thresh){
  library(tidyverse)
  tfs_kd_ranked <- apply(tf_kd_enrichment,2,rank,"random")
  kd_df_line <- kd_df[which(as.character(kd_df$cell_id) %in% cell),]
  kd_df_line <- kd_df_line %>% filter(quality==1)
  
  tfs_cell <- tfs_kd_ranked[,as.character(kd_df_line$sig_id)]
  tfs_cell <- tfs_cell/nrow(tfs_cell)
  avail <- NULL
  for (i in 1:ncol(tfs_cell)) {
    if (tfs_cell[as.character(kd_df_line$pert_iname[i]),i] <= thresh) {
      avail <- c(avail,as.character(kd_df_line$sig_id[i]))
    }
  }
  return(avail)
}
neighbor_train_selection <- function(tf_sig,tf_name, cell, training_set, go_compounds, go_kds, tf_scores, thresh_go, thresh_tf,workers){
  library(doFuture)
  registerDoFuture()
  plan(multiprocess,workers = workers)
  
  
  distance_scores_all <- function(num_table1,num_table2, threshold_count, names1,names2) {
    library(GeneExpressionSignature)
    library(tidyverse)
    
    ### rank the table
    num_table1 <- as.matrix(num_table1)
    num_table2 <- as.matrix(num_table2)
    table_ranked1 <- apply(X = num_table1, MARGIN = 2, FUN = rank, ties.method ="random")
    table_ranked2 <- apply(X = num_table2, MARGIN = 2, FUN = rank, ties.method ="random")
    final_dist <- matrix(nrow=NCOL(num_table1),ncol=NCOL(num_table2))
    for (i in 1:NCOL(num_table1)) {
      for (j in 1:NCOL(num_table2)) {
        ### create expression set
        merged_table <- cbind(table_ranked1[,i],table_ranked2[,j])
        colnames(merged_table) <- c(names1[i],names2[j])
        pheno <- as.data.frame(colnames(merged_table))
        rownames(pheno) <- colnames(merged_table)
        pheno_new <- new("AnnotatedDataFrame",data=pheno)
        expr_set <- new("ExpressionSet",exprs = merged_table, phenoData=pheno_new)
        ### calculate distances
        distances <- ScoreGSEA(expr_set , threshold_count,"avg")
        final_dist[i,j] <-distances[1,2] #or (2,1) just not the primary diagonal 
      }
      #print((i/ncol(num_table1))*100)
    }
    colnames(final_dist) <- names2
    rownames(final_dist) <- names1
    return(final_dist)
  }
  
  sigs_train <- unique(c(as.character(training_set$sig_id.x),as.character(training_set$sig_id.y)))
  go_train <- go_compounds[,sigs_train]
  go_tf <- go_kds[,as.character(tf)]
  
  thresholds <- c(10,20,30,40,50)
  go_distance <- NULL
  print("starting go_distance")
  go_distance <- foreach(thresh = thresholds) %dopar% {
    
    distance_scores_all(num_table1 = go_train,num_table2 =go_tf, 
                        threshold_count = thresh, names1 = as.character(colnames(go_train)),names2 = as.character(tf))
  }
  go_avg <- (go_distance[[1]]+go_distance[[2]]+go_distance[[3]]+go_distance[[4]]+go_distance[[5]])/5
  go_avg <- go_avg/2
  
  compounds_tf <- tf_scores[,sigs_train]
  compounds_tf <- compounds_tf/nrow(compounds_tf)
  print(paste0("the minimum distance is ",min(go_avg[,1])))
  min_sig <-  rownames(go_avg)[which(go_avg[,1] == min(go_avg[,1]))]
  candidates1 <- rownames(go_avg)[which(go_avg[,1]<=thresh_go)]
  
  candidates_tf <- compounds_tf[as.character(tf_name),as.character(candidates1)]
  sig_results <- names(candidates_tf)[which(candidates_tf<=thresh_tf)]
  
  return(sig_results)
}
inference_selection <- function(predictions,neighbors,test_thresh){
  ind <- unique(c(which(as.character(predictions$sig_id.x) %in% neighbors),which(as.character(predictions$sig_id.y) %in% neighbors)))
  filt <- predictions[ind,]
  filt <- filt %>% filter(mu<=test_thresh)
  
  test_sigs <- unique(c(as.character(filt$sig_id.x),as.character(filt$sig_id.y)))
  test_sigs <- test_sigs[-which(test_sigs %in% neighbors)]
  maj <- NULL
  for (i in 1:length(test_sigs)) {
      ind_test <- unique(c(which(as.character(filt$sig_id.x) %in% test_sigs[i]),
                           which(as.character(filt$sig_id.y) %in% test_sigs[i])))
      filt2 <- filt[ind_test,]
      maj <- c(maj,nrow(filt2)/length(neighbors))
  }
  results <- cbind(test_sigs,maj)

  return(as.data.frame(results))
}
evaluate_inferrence <- function(inferred,maj,tf_sig,tf_name,go_compounds,go_kds,tf_scores,workers,dist_thresh,tf_thresh){
  library(doFuture)
  registerDoFuture()
  plan(multiprocess,workers = workers)
  
  
  distance_scores_all <- function(num_table1,num_table2, threshold_count, names1,names2) {
    library(GeneExpressionSignature)
    library(tidyverse)
    
    ### rank the table
    num_table1 <- as.matrix(num_table1)
    num_table2 <- as.matrix(num_table2)
    table_ranked1 <- apply(X = num_table1, MARGIN = 2, FUN = rank, ties.method ="random")
    table_ranked2 <- apply(X = num_table2, MARGIN = 2, FUN = rank, ties.method ="random")
    final_dist <- matrix(nrow=NCOL(num_table1),ncol=NCOL(num_table2))
    for (i in 1:NCOL(num_table1)) {
      for (j in 1:NCOL(num_table2)) {
        ### create expression set
        merged_table <- cbind(table_ranked1[,i],table_ranked2[,j])
        colnames(merged_table) <- c(names1[i],names2[j])
        pheno <- as.data.frame(colnames(merged_table))
        rownames(pheno) <- colnames(merged_table)
        pheno_new <- new("AnnotatedDataFrame",data=pheno)
        expr_set <- new("ExpressionSet",exprs = merged_table, phenoData=pheno_new)
        ### calculate distances
        distances <- ScoreGSEA(expr_set , threshold_count,"avg")
        final_dist[i,j] <-distances[1,2] #or (2,1) just not the primary diagonal 
      }
      #print((i/ncol(num_table1))*100)
    }
    colnames(final_dist) <- names2
    rownames(final_dist) <- names1
    return(final_dist)
  }
  
  go_test <- go_compounds[,as.character(inferred)]
  go_test <- as.matrix(go_test)
  colnames(go_test) <- as.character(inferred)
  go_tf <- go_kds[,as.character(tf_sig)]
  
  thresholds <- c(10,20,30,40,50)
  go_distance <- NULL
  print("starting go_distance")
  go_distance <- foreach(thresh = thresholds) %dopar% {
    
    distance_scores_all(num_table1 = go_test,num_table2 =go_tf, 
                        threshold_count = thresh, names1 = as.character(colnames(go_test)),names2 = tf)
  }
  go_avg <- (go_distance[[1]]+go_distance[[2]]+go_distance[[3]]+go_distance[[4]]+go_distance[[5]])/5
  go_avg <- go_avg/2
  
  compounds_tf <- tf_scores[,inferred]
  compounds_tf <- as.matrix(compounds_tf)
  colnames(compounds_tf) <- as.character(inferred)
  compounds_tf <- compounds_tf/nrow(compounds_tf)
  compounds_tf <- compounds_tf[tf_name,]
  
  results <- cbind(go_avg,compounds_tf)
  colnames(results) <- c("go_distance","tf_rank")
  results <- as.data.frame(results)
  results$maj <- maj
  dist_acc <- length(which(results$go_distance<=dist_thresh))/nrow(results)
  tf_acc <- length(which(results$tf_rank<=tf_thresh))/nrow(results)
  combined_acc <- length(which(results$go_distance<=dist_thresh & results$tf_rank<=tf_thresh))/nrow(results)
  results$dist_acc <- dist_acc
  results$tf_acc <- tf_acc
  results$combined_acc <- combined_acc
  return(results)
}

tf_kd_enrichment <- readRDS("enrichment_calculations/TF_enrichment/TF_enrichment_results/knockdowns/tf_enrichment_kd.rds")
kd_df <- readRDS("preprocessing/preprocessing_results/knockdowns/tf_knockdowns_consensus.rds")



training_set <- read.csv("data/a375train.csv")
go_compounds <- readRDS("enrichment_calculations/GO_enrichment/GO_enrichment_results/compounds/all_compounds_go.rds")
go_kds <- readRDS("enrichment_calculations/GO_enrichment/GO_enrichment_results/knockdowns/go_nes_results_kd.rds")
tf_scores <- readRDS("enrichment_calculations/TF_enrichment/TF_enrichment_results/compounds/all_tf_ranked.rds")

predictions <- read.csv("data/predictions_test.csv")

a375_tfs <- available_tfs_line(tf_kd_enrichment,kd_df,"A375",0.25)
tf <- a375_tfs[2]
tf_name <- "MYC"
neighbors <- neighbor_train_selection(tf_sig = tf,tf_name = tf_name,cell = "A375",
                                training_set = training_set,go_compounds = go_compounds,go_kds = go_kds,
                                tf_scores = tf_scores,thresh_go = 0.15,thresh_tf = 0.2,workers = 12)


inferred <- inference_selection(predictions,neighbors,0.22)

eval <- evaluate_inferrence(inferred = inferred$test_sigs,maj = inferred$maj,tf_sig = tf,tf_name = tf_name,
                            go_compounds = go_compounds,go_kds = go_kds,tf_scores = tf_scores,workers = 12,
                            dist_thresh = 0.23,tf_thresh = 0.2)




