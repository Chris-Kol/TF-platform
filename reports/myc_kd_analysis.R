library(tidyverse)
library(viridis)
library(reshape2)

ds_path <- "C:/Users/user/Documents/phd/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
landmark <- read_tsv(file = "myc_cmap_pathways/cmap_landmark_genes.txt")

### read the kd consensus signatures

kd <- readRDS("data/knockdowns/tf_kd_cgs.rds")

### keep the myc_kd

myc_kd <- kd %>% filter(pert_iname == "MYC")

### load the profiles

profiles_myc <- get_cmap_signatures(cmap_path_to_gctx = ds_path, sig_ids = myc_kd$sig_id,landmark = T,landmark_df = landmark)

### calculate pathway scores

pathways_myc <- kegg_path_analysis(sig_ids = myc_kd$sig_id,cmap_path_to_gctx = ds_path,landmark_df = landmark)

pathways_myc <- pathways_myc[[1]]

### calculate go terms

gos <- readRDS("data/GOterms/goterm_annotation.rds")
go_myc <- go_path_analysis(sig_ids = myc_kd$sig_id,cmap_path_to_gctx = ds_path,landmark_df = landmark,goterms = gos)

go_myc <- go_myc[[1]]


### calculate distances with scoregsea

genes_dist <- distance_scores(num_table = profiles_myc,threshold_count = 25, names= myc_kd$cell_id)

paths_dist <- distance_scores(num_table = pathways_myc,threshold_count = 10, names = myc_kd$cell_id)

go_dist <- distance_scores(num_table = go_myc, threshold_count = 25, names = myc_kd$cell_id)

### calculate correlations

colnames(profiles_myc) <- myc_kd$cell_id
colnames(pathways_myc) <- myc_kd$cell_id
colnames(go_myc) <- myc_kd$cell_id

genes_cor <- round(cor(profiles_myc,method = "spearman"),2)

pathways_cor <- round(cor(pathways_myc,method = "spearman"),2)

go_cor <- round(cor(go_myc,method = "spearman"),2)

### plot correlations

melted_genes_cor <- melt(genes_cor)

ggplot(data = melted_genes_cor, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+scale_fill_viridis(limits = c(0,1)) + labs(title = "gene-level spearman")

melted_pathways_cor <- melt(pathways_cor)

ggplot(data = melted_pathways_cor, aes(x=Var1, y=Var2, fill= value)) + 
  geom_tile()+scale_fill_viridis(limits = c(0,1)) + labs(title = "KEGG pathway level spearman")

melted_go_cor <- melt(go_cor)

ggplot(data = melted_go_cor, aes(x=Var1, y=Var2, fill= value)) + 
  geom_tile()+scale_fill_viridis(limits = c(0,1)) + labs(title = "GO BP level spearman")


### plot distances

melted_genes_dist <- melt(genes_dist)

ggplot(data = melted_genes_dist, aes(x=Var1, y=Var2, fill= 1-value)) + 
  geom_tile()+scale_fill_viridis(limits = c(0,1)) + labs(title = "Gene Level GSEA scores")

melted_paths_dist <- melt(paths_dist)

ggplot(data = melted_paths_dist, aes(x=Var1, y=Var2, fill= 1-value)) + 
  geom_tile()+scale_fill_viridis(limits = c(0,1)) + labs(title = "KEGG Pathway Level GSEA scores")

melted_go_dist <- melt(go_dist)

ggplot(data = melted_go_dist, aes(x=Var1, y=Var2, fill= 1-abs(value))) + 
  geom_tile()+scale_fill_viridis(limits = c(0,1)) + labs(title = "GO BP Level GSEA scores")

