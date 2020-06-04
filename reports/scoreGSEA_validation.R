library(tidyverse)
library(cancerTiming)

drug_sigs_per_line_dups <- function(cell_line,sig_info,sig_metrics) {
  
  # cell_line character of cell line
  # sig_info dataframe of GSE info
  # sig_metrics dataframe of GSE metrics
  
  library(tidyverse)
  options(warn =- 1)
  
  
  cell <- sig_info %>%
    filter(cell_id == cell_line) %>%
    filter(pert_type == "trt_cp") %>%
    group_by(pert_iname) %>%
    mutate(count = n_distinct(sig_id)) %>%
    ungroup()
  
  print(paste0('the unique drugs for ',cell_line,' are ',length(unique(cell$pert_iname))))
  
  ## add the signature metrics
  
  cell <- left_join(cell,sig_metrics)
  
  ## keep the drugs that we have only 1 signature for this cell line
  
  cell_singles <- cell %>%
    filter(count == 1) %>%
    dplyr::select(-count)
  
  print(paste0('the drugs that have only 1 signature for ',cell_line,' are ',length(unique(cell_singles$pert_iname))))
  
  cell_singles$pert_itime <- factor(cell_singles$pert_itime)
  print("time summary")
  print(summary(cell_singles$pert_itime))
  cell_singles$pert_idose <- factor(cell_singles$pert_idose)
  print("dose summary")
  print(summary(cell_singles$pert_idose))
  
  ## add quality column to single perturbations
  
  cell_singles$quality <- 100
  
  cell_singles <- cell_singles %>%
    mutate(quality = if_else(is_exemplar == 1 & tas > 0.4 & distil_nsample>=2 ,true = 1,false = quality),
           quality = if_else(is_exemplar == 1 & tas > 0.2 & tas<=0.4 & distil_nsample>2 ,true = 2,false = quality),
           quality = if_else(is_exemplar == 1 & tas > 0.2 & tas<=0.4 & distil_nsample <=2 ,true = 3,false = quality),
           quality = if_else(is_exemplar == 1 & tas > 0.1 & tas<=0.2 & distil_nsample>2 ,true = 4,false = quality),
           quality = if_else(is_exemplar == 1 & tas > 0.1 & tas<=0.2 & distil_nsample <= 2 ,true = 5,false = quality),
           quality = if_else(is_exemplar == 1 & tas < 0.1 & distil_nsample > 2 ,true = 6,false = quality),
           quality = if_else(is_exemplar == 1 & tas < 0.1 & distil_nsample <= 2 ,true = 7,false = quality),
           quality = if_else(is_exemplar == 0 ,true = 8,false = quality),
           quality = factor(quality))
  
  print("summary of the quality of drugs with only 1 signature")
  print(summary(cell_singles$quality))
  
  ## keep the multiple signature drugs in cell
  
  cell<- anti_join(cell,cell_singles)
  
  ### add priorities to the multiple signatures
  
  cell$priority <- 100
  cell <- cell %>%
    mutate(priority = if_else(pert_dose == "10.0" & pert_time == 24,true = 1,false = priority),
           priority = if_else(pert_idose == "5 ÂµM" & pert_time == 24,true = 2,false = priority),
           priority = if_else(pert_idose != "5 ÂµM" & pert_dose != "10.0" & pert_time == 24,true = 3,false = priority),
           priority = if_else(pert_dose == "10.0" & pert_time == 6,true = 4,false = priority),
           priority = if_else(pert_idose == "5 ÂµM" & pert_time == 6,true = 5,false = priority),
           priority = if_else(pert_idose != "5 ÂµM" & pert_dose != "10.0" & pert_time == 6,true = 6,false = priority),
           priority = factor(priority))
  
  print("priorities for drugs with multiple signatures")
  print(summary(cell$priority))
  ### add quality to the multiple signatures
  
  cell$quality <- 100
  
  cell <- cell %>%
    mutate(quality = if_else(is_exemplar == 1 & tas > 0.4 & distil_nsample>=2 ,true = 1,false = quality),
           quality = if_else(is_exemplar == 1 & tas > 0.2 & tas<=0.4 & distil_nsample>2 ,true = 2,false = quality),
           quality = if_else(is_exemplar == 1 & tas > 0.2 & tas<=0.4 & distil_nsample <=2 ,true = 3,false = quality),
           quality = if_else(is_exemplar == 1 & tas > 0.1 & tas<=0.2 & distil_nsample>2 ,true = 4,false = quality),
           quality = if_else(is_exemplar == 1 & tas > 0.1 & tas<=0.2 & distil_nsample <= 2 ,true = 5,false = quality),
           quality = if_else(is_exemplar == 1 & tas < 0.1 & distil_nsample > 2 ,true = 6,false = quality),
           quality = if_else(is_exemplar == 1 & tas < 0.1 & distil_nsample <= 2 ,true = 7,false = quality),
           quality = if_else(is_exemplar == 0 ,true = 8,false = quality),
           quality = factor(quality))
  
  
  print("summary of the quality of drugs with multiple signatures")
  print(summary(cell$quality))
  
  
  print(paste0('the drugs that have Multiple signatures for ',cell_line,' are ',length(unique(cell$pert_iname))))
  
  
  
  
  #### clean them based on quality for each drug and then solve the equalities with max tas
  
  
  
  
  cell_cleaned <- cell %>%
    group_by(pert_iname) %>%
    filter(quality == min(as.numeric(quality))) %>%
    ungroup %>%
    dplyr::select(-c(count,priority))
  
  
  
  cell_final <- bind_rows(cell_cleaned,cell_singles)
  
  print("summary of final quality of signatures")
  print(summary(cell_final$quality))
  
  return(cell_cleaned)
}

kegg_path_analysis <- function(sig_ids, cmap_path_to_gctx, landmark_df) {
  ### this function calculates the NES and p.adj of the given signature ids
  ### 125 KEGG pathways are used
  library(tidyverse)
  library(fgsea)
  library(gage)
  library(EGSEAdata)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  
  ### first thing is to load the profiles from the GSE file
  
  profiles <- get_cmap_signatures(cmap_path_to_gctx = cmap_path_to_gctx,sig_ids = sig_ids,landmark = T,landmark_df = landmark_df)
  print("profiles loaded")
  ### load the gene sets
  print("loading gene sets")
  egsea.data(species = "human",returnInfo = TRUE)
  rownames(profiles) <- landmark$`Entrez ID`
  rownames(profiles) <- as.character(rownames(profiles))
  
  print("running fgsea")
  ### run the analysis
  pathway_list <- apply(profiles,MARGIN = 2,fgsea,pathways = kegg.pathways$human$kg.sets,
                        minSize=10,
                        maxSize=500,
                        nperm=1000)
  print("fgsea finished")
  ### get the NES and p.adj
  
  print("preparing output")
  NES <- pathway_list[[1]]$NES
  padj <- pathway_list[[1]]$padj
  
  for (i in 2:length(pathway_list)) {
    
    NES <- cbind(NES,pathway_list[[i]]$NES)
    padj <- cbind(padj,pathway_list[[i]]$padj)
  }
  
  colnames(NES) <- names(pathway_list)
  rownames(NES) <- pathway_list[[1]]$pathway
  colnames(padj) <- names(pathway_list)
  rownames(padj) <- pathway_list[[1]]$pathway
  
  comb <- list(NES,padj)
  
  return(comb)
}

get_cmap_signatures <- function(cmap_path_to_gctx, sig_ids, landmark = TRUE, landmark_df = NULL) {
  
  
  library(tidyverse)
  library(cmapR)
  library(rhdf5)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  
  ds_path <- cmap_path_to_gctx
  if (landmark == TRUE) {
    
    cmap_gctx <- parse.gctx(ds_path,rid = as.character(landmark_df$`Entrez ID`), cid = sig_ids)
    cmap <- cmap_gctx@mat
    
    cmap <- cmap[as.character(landmark_df$`Entrez ID`),]
    
    rownames(cmap) <- landmark_df$Symbol
  }
  
  if (landmark == FALSE) {
    
    cmap_gctx <- parse.gctx(ds_path, cid = sig_ids)
    cmap <- cmap_gctx@mat
    
    entrez <- rownames(cmap)
    anno <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = entrez,
                                  columns = c("SYMBOL", "GENENAME","ENTREZID"),
                                  keytype = "ENTREZID")
    
    anno <- anno %>%
      filter(!is.na(SYMBOL))
    
    cmap <- cmap[anno$ENTREZID,]
    
    rownames(cmap) <- anno$SYMBOL
  }
  
  
  return(cmap)
  
}

distance_scores <- function(num_table, threshold_count, names) {
  library(GeneExpressionSignature)
  library(tidyverse)
  
  ### rank the table
  table_ranked <- apply(X = -num_table, MARGIN = 2, FUN = rank, ties.method = "random")
  
  ### create the phenodata
  pheno2 <- as.data.frame(colnames(num_table))
  rownames(pheno2) <- colnames(num_table)
  pheno_new <- new("AnnotatedDataFrame",data=pheno2)
  ### create expression set
  expr_set <- new("ExpressionSet",exprs = table_ranked, phenoData=pheno_new)
  ### calculate distances
  distances <- ScoreGSEA(expr_set , threshold_count,"avg")
  colnames(distances) <- names
  rownames(distances) <- names
  return(distances)
}



sig <- read.delim("data/cmap_data/GSE92742_Broad_LINCS_sig_info.txt")
ds_path <- "C:/Users/user/Documents/phd/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
landmark <- read_tsv(file = "myc_cmap_pathways/cmap_landmark_genes.txt")
sig_metrics <- read.delim("data/cmap_data/GSE92742_Broad_LINCS_sig_metrics.txt")

test <- drug_sigs_per_line_dups(cell_line = "MCF7",sig_info = sig,sig_metrics = sig_metrics)

test <- test %>% 
  filter(pert_time == "24") %>%
  filter(quality == "1" | quality == "2") %>%
  group_by(pert_iname) %>% 
  mutate(count = n_distinct(sig_id)) %>% ungroup() %>% filter(count > 1)

drugs <- unique(as.character(test$pert_iname))

drug <- drugs[1]

paths <- kegg_path_analysis(sig_ids = test$sig_id[which(test$pert_iname %in% drug)],cmap_path_to_gctx = ds_path,landmark_df = landmark)
profiles <- get_cmap_signatures(cmap_path_to_gctx = ds_path,sig_ids = test$sig_id[which(test$pert_iname %in% drug)],landmark = T,landmark_df = landmark)
nes <- paths[[1]]
pval <- paths[[2]]

thres <- c(3,4,5,6,7,8,9,10,11,12)

sam1_tria <- readRDS("random_non_rep.rds")
rep_mean <- matrix(0,nrow = length(drugs),ncol = 1)
rownames(rep_mean) <- drugs

for (i in 1:length(drugs)) {
  drug <- drugs[i]
  dist_rep <- list(0)
  paths <- kegg_path_analysis(sig_ids = test$sig_id[which(test$pert_iname %in% drug)],cmap_path_to_gctx = ds_path,landmark_df = landmark)
  nes2 <- paths[[1]]
  for (j in 1:length(thres)) {
    dist2 <- distance_scores(num_table = nes2,threshold_count = thres[j],names = colnames(nes2))
    dist_rep[[j]] <- dist2
  }
  sam <- 0
  for (k in 1:length(dist_rep)) {
    sam <- sam + dist_rep[[k]]
  }
  sam <- sam/length(dist_rep)
  rep_mean[i,1] <- mean(sam[upper.tri(x = sam,diag = F)])
  #rep_mean[i,1] <- sam
}

multidensity(list(rep_mean[,1],sam1_tria),main= "distance distribution replicates vs random",xlab = "mean distance GSEA")
legend("topright", 
       legend = c("replicates", "random distances"), 
       col = c('black', 
               'red'),
       pch = 'lines',
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))

t.test(x = rep_mean[,1],y=sam1_tria)
