library(tidyverse)

cells <- read.delim("data/cmap_data/GSE92742_Broad_LINCS_cell_info.txt")

cells <- cells %>% select(base_cell_id,primary_site) %>% unique()

#### read original

data <- readRDS("data/screening hits/initial_signatures_with_mycranks.rds")

data <- left_join(data,cells,by=c("cell_id"="base_cell_id"))

data$primary_site[which(is.na(data$primary_site))] <- "large intestine"


### replace data with the new one

saveRDS(data, "data/screening hits/initial_signatures_with_mycranks.rds")



