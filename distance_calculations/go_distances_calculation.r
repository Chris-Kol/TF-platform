library(tidyverse)
library(future.apply)

plan(multiprocess, workers = 12)


vectorized_gsea <- function(x1,x2,thres) {
  library(GeneExpressionSignature)
  x1 <- as.matrix(rank(-x1,ties.method ="random"))
  colnames(x1) <- "q"
  x2 <- apply(X = -x2, MARGIN = 2, FUN = rank, ties.method ="random")
  final_dist <- NULL
  for (j in 1:ncol(x2)) {
    ### create expression set
    merged_table <- cbind(x1,x2[,j])
    colnames(merged_table) <- c(colnames(x1),colnames(x2)[j])
    pheno <- as.data.frame(colnames(merged_table))
    rownames(pheno) <- colnames(merged_table)
    pheno_new <- new("AnnotatedDataFrame",data=pheno)
    expr_set <- new("ExpressionSet",exprs = merged_table, phenoData=pheno_new)
    ### calculate distances
    distances <- ScoreGSEA(expr_set , thres,"avg")
    final_dist[j] <-distances[1,2] #or (2,1) just not the primary diagonal 
  }
  #print("lol")
  return(final_dist)
}

go_cpd_nes <- readRDS("data/GOterms/compounds/all_compounds_go.rds")
go_kd_nes <- readRDS("data/GOterms/kd/go_nes_results_kd.rds")

data <- readRDS("data/screening hits/initial_signatures_with_mycranks.rds")
kd <- readRDS("data/knockdowns/tf_kd_cgs.rds")

thresholds <- c(5,10,15,20,25,35,45,60,70) ### change for each line

go_cpd_line <- go_cpd_nes[,as.character(data$sig_id[which(data$cell_id == "A375")])] ### change for each line
go_kd_line <- go_kd_nes[,as.character(kd$sig_id[which(kd$cell_id == "A375")])] ### change for each line

distances <- NULL

for (i in 1:length(thresholds)) {
  distances[[i]] <- future_apply(go_cpd_line,MARGIN = 2,FUN = vectorized_gsea,go_kd_line,thresholds[i])
  rownames(distances[[i]]) <- colnames(go_kd_line) ### added rownames
  print(i)
}

saveRDS(distances,"distances/goterms/A375_go.rds") ### change for each line
