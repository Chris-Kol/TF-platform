library(tidyverse)
library(reshape2)
library(ggplot2)
library(GeneExpressionSignature)
library(fgsea)
library(gage)
library(EGSEAdata)
library(AnnotationDbi)
library(org.Hs.eg.db)

options(warn =- 1)

ds_path <- "C:/Users/user/Documents/phd/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
landmark <- read_tsv(file = "myc_cmap_pathways/cmap_landmark_genes.txt")
# read cmap signature info dataframe

sig_info <- read.delim(file = "data/cmap_data/GSE92742_Broad_LINCS_sig_info.txt")

# read cmap signature metrics

sig_metrics <- read.delim(file = "data/cmap_data/GSE92742_Broad_LINCS_sig_metrics.txt")

### read the tfs of dorothea

tfs_dorothea <- readRDS("TF_enrichment_landmark/tfs_of_dorothea_for_landmark.rds")

sh <- sig_info %>%
  filter(pert_type == "trt_sh") %>%
  group_by(pert_iname) %>%
  mutate(count = n_distinct(sig_id)) %>%
  ungroup()

sh_cgs <- sig_info %>%
  filter(pert_type == "trt_sh.cgs") %>%
  group_by(pert_iname) %>%
  mutate(count = n_distinct(sig_id)) %>%
  ungroup()

### keep the shrna signatures of the TFS of viper

sh_cgs_filtered <- sh_cgs[which(sh_cgs$pert_iname %in% tfs_dorothea),]

### add sig_metrics to shrna_cgs

sh_cgs_filtered <- left_join(sh_cgs_filtered,sig_metrics)

### add quality to sh_cgs

sh_cgs_filtered$quality <- 100

sh_cgs_filtered <- sh_cgs_filtered %>%
  mutate(quality = if_else(is_exemplar == 1 & tas > 0.4 & distil_nsample>=2 ,true = 1,false = quality),
         quality = if_else(is_exemplar == 1 & tas > 0.4 & distil_nsample>=1 ,true = 1,false = quality),
         quality = if_else(is_exemplar == 1 & tas > 0.2 & tas<=0.4 & distil_nsample>2 ,true = 2,false = quality),
         quality = if_else(is_exemplar == 1 & tas > 0.2 & tas<=0.4 & distil_nsample <=2 ,true = 3,false = quality),
         quality = if_else(is_exemplar == 1 & tas > 0.1 & tas<=0.2 & distil_nsample>2 ,true = 4,false = quality),
         quality = if_else(is_exemplar == 1 & tas > 0.1 & tas<=0.2 & distil_nsample <= 2 ,true = 5,false = quality),
         quality = if_else(is_exemplar == 1 & tas < 0.1 & distil_nsample > 2 ,true = 6,false = quality),
         quality = if_else(is_exemplar == 1 & tas < 0.1 & distil_nsample <= 2 ,true = 7,false = quality),
         quality = if_else(is_exemplar == 0 ,true = 8,false = quality),
         quality = factor(quality))



### filter the duplicate shrna for the same cell line to keep the highest quality

sh_cgs_cleaned <- sh_cgs_filtered %>%
  group_by(pert_iname,cell_id) %>%
  filter(quality == min(as.numeric(quality))) %>%
  filter(tas == max(tas)) %>%
  ungroup #%>%
  #dplyr::select(-c(count,priority))


#saveRDS(sh_cgs_cleaned,"data/knockdowns/tf_kd_cgs.rds")





