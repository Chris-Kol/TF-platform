library(tidyverse)

### read the initial signatures of cmap version 2018

data <- readRDS("data/screening hits/initial_signatures_with_mycranks.rds")


### read the GSE of the perturbations

pert_info <- read.delim(file = "data/cmap_data/GSE92742_Broad_LINCS_pert_info.txt")

pert_info <- pert_info %>%
  filter(pert_type == "trt_cp") %>%
  filter(inchi_key != -666)

### drugs in cmap

cmap_drugs <- data %>% select(pert_id,pert_iname) %>% unique()

### add drug info to drugs in cmap

cmap_drugs <- left_join(cmap_drugs,pert_info)

### add moa and target from the broad repo

broad_repo <- read.delim(file = "data/repo hub/Repurposing_Hub_export (1).txt" ,skip = 0)

broad_repo <- broad_repo %>%
  mutate(InChIKey = str_split(InChIKey,pattern = ",")) %>% unnest(InChIKey) %>%
  mutate(InChIKey = str_trim(InChIKey)) %>%
  filter(InChIKey != "") %>% unique()

broad_repo <- broad_repo  %>%
  mutate(broad_id = substr(x = as.character(Id),start = 1,stop = 13)) %>%
  mutate(broad_id_old = substr(x = as.character(Deprecated.ID),start = 1,stop = 13)) %>%
  mutate(Name = toupper(Name))


### lists of indices for mapping between cmap_drugs + broad_repo

new <- which(unique(cmap_drugs$pert_id) %in% broad_repo$broad_id)
old <- which(unique(cmap_drugs$pert_id) %in% broad_repo$broad_id_old)
inchi <- which(cmap_drugs$inchi_key %in% broad_repo$InChIKey)
smiles <- which(cmap_drugs$canonical_smiles %in% broad_repo$SMILES)
names <- which(cmap_drugs$pert_iname %in% broad_repo$Name)

a <- (unique(union(new,inchi)))
b <- unique(union(a,old))
c <- unique(union(b,smiles))
d <- unique(union(c,names))

### 2245 were mapped to broad repo

df1 <- left_join(cmap_drugs,broad_repo,by = c("pert_id"="broad_id"))
df2 <- left_join(cmap_drugs,broad_repo,by = c("pert_id"="broad_id_old"))
df3 <- left_join(cmap_drugs,broad_repo,by = c("inchi_key"="InChIKey"))
df4 <- left_join(cmap_drugs,broad_repo,by = c("canonical_smiles"="SMILES"))
df5 <- left_join(cmap_drugs,broad_repo,by = c("pert_iname"="Name"))

df <- bind_rows(df1,df2,df3,df4,df5) %>% unique()

df_new <- df %>%
  select(pert_id,Target,MOA,Phase,Disease.Area) %>%
  filter(!(is.na(Target)&is.na(MOA))) %>%
  mutate(nchar_target = nchar(as.character(Target))) %>%
  group_by(pert_id) %>% 
  filter(nchar_target == max(nchar_target)) %>% 
  ungroup() %>% 
  unique()


### add the moa and target to cmap_drugs

cmap_drugs <- left_join(cmap_drugs,df_new)

### save the results

saveRDS(cmap_drugs,"data/screening hits/all_cmapdrugs_moa+target_v1.rds")
