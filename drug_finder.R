drug_finder <- function(KD, tf_knockdowns, cell_line, sig_ids_train, thresh1, thresh2,thresh3,cmap_path){
  # This function finds drugs similar to knockdowns of transcription factors of a selected cell line
  # The input is a knockdown as character, the paths to the train file of the cell line and the tf_knockdowns file, the cell line name as characters
  # Moreover, the input contains 3 thresholds: The first has to do with the rank of the kds that are good to work with
  # The second has to do with the distances of those kds and the drugs
  # The third has to do with the rank of those drugs 
  # Finally, the cmap_path is the file GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx file from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742 on your working directory.
  # There also needs to be the cmap_landmark_genes.txt file from our Drive: https://drive.google.com/drive/folders/1I6UzU9lK-zISwlqBGD7grDFI9y7Maf1-
    
  
  #the extra functions
  tf_enrichment_ranks <- function(tf_knockdowns, cell_line, sig_ids_train) {
    ###function used
    runDoRothEA<-function(df, regulon, confidence_level=c('A','B','C'), write2file = NULL){
      # library(tidyverse)
      library(dplyr)
      library(purrr)
      library(viper)
      library(tibble)
      library(tidyr)
      names(regulon) <- sapply(strsplit(names(viper_regulon), split = ' - '), head, 1)
      # names(regulon) <- sapply(strsplit(names(viper_regulon), split = ' - '), head, 1)
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
    
    ###go_term_anno to annotate the GO terms with gene identifiers
    library(tidyverse)
    library(topGO)
    library(org.Hs.eg.db)
    library(GO.db)
    ### read landmark genes
    
    landmark_df <- read_tsv(file = "cmap_landmark_genes.txt")
    
    ### go terms
    
    genes <- factor(x = rep(1,978),levels = c(0,1))
    names(genes) <- landmark_df$`Entrez ID`
    
    GOobject <- new("topGOdata",ontology = "BP", allGenes = genes, annot=annFUN.org, mapping="org.Hs.eg.db", 
                    ID = "entrez", nodeSize = 10)
    
    term.genes <- genesInTerm(GOobject, GOobject@graph@nodes)
    ###end of go_term_anno
    
    
    ### go_path_analysis function to exact GSEA for GO terms related to biological processes
    
    #it calculates the NES and p.adj of the given signature ids
    #GOterms in a list form are used
    library(fgsea)
    library(gage)
    library(EGSEAdata)
    library(AnnotationDbi)
    library(GeneExpressionSignature)
    
    ### you need to load the profiles from the GSE file
    #
      sig_ids <- read.csv(paste0(sig_ids_train))
      sigs_x <- unique(as.character(sig_ids$sig_id.x))
      sigs_y <- unique(as.character(sig_ids$sig_id.y))
      sigs <- c(sigs_x,sigs_y)
      sig_ids <- unique(sigs)

      kd_sigs <- readRDS(paste0(tf_knockdowns))
      kd_sigs <- kd_sigs %>% dplyr::filter(quality==1) %>% dplyr::filter(cell_id==paste0(toupper(cell_line)))
      ksig_ids <- kd_sigs %>% dplyr::select(sig_id)
      ksig_ids <- as.character(ksig_ids$sig_id)
      sig_ids <- c(sig_ids,ksig_ids)

    
    
    landmark_df <- read_tsv(file = "cmap_landmark_genes.txt")
    goterms <- term.genes
    
    if (cmap_path != 0){
      profiles_train <- get_cmap_signatures(cmap_path_to_gctx = cmap_path,sig_ids = sig_ids,landmark = T,landmark_df = landmark_df)
    } else {
      profiles_train <- kd_train_genes
    }
    ###end of go_path_analysis  
    
    load(file = system.file("BEST_viperRegulon.rdata",package="CARNIVAL")) # loading the viper regulons
    #landmark = landmark annotation
    #sigs = vector of signature IDs
    #cmap_path = path to lvl 5 mod z gctx
    library(tidyverse)
    library(cmapR)
    library(rhdf5)
    library(CARNIVAL)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    library(viper)
    library(matrixStats)
    
    landmark <- read_tsv(file = "cmap_landmark_genes.txt")
    ### load the profiles expression
    cmap_data <- profiles_train 
    ### calculate tf enrichment
    TF_cmap_train <-runDoRothEA(cmap_data, regulon=viper_regulon, confidence_level=c('A','B','C')) # Estimating TF activities
    TF_cmap_train_rank_rowise <-(rowRanks(as.matrix(TF_cmap_train, rows = NULL, cols = 1:ncol(TF_cmap_train), ties.method = "random")))
    colnames(TF_cmap_train_rank_rowise) <- colnames(TF_cmap_train)
    rownames(TF_cmap_train_rank_rowise) <- rownames(TF_cmap_train)
    TF_cmap_train_rank_rowise <- as.data.frame(TF_cmap_train_rank_rowise)
    return(TF_cmap_train_rank_rowise)
  }
  distance_go_kd_train <- function(tf_knockdowns, cell_line, sig_ids_train,cmap_path) {
    
    ###the functions used
    #1)dinstance_scores_all
    distance_scores_all <- function(num_table1,num_table2, threshold_count, names1,names2) {
      library(GeneExpressionSignature)
      library(tidyverse)
      
      ### rank the table
      table_ranked1 <- apply(X = -num_table1, MARGIN = 2, FUN = rank, ties.method ="random")
      table_ranked2 <- apply(X = -num_table2, MARGIN = 2, FUN = rank, ties.method ="random")
      final_dist <- matrix(nrow=NCOL(num_table1),ncol=NCOL(num_table2))
      for (i in 1:NCOL(num_table1)) {
        for (j in 1:NCOL(num_table2)) {
          ### create expression set
          merged_table <- cbind(table_ranked1[,i],table_ranked2[,j])
          colnames(merged_table) <- c(colnames(table_ranked1)[i],colnames(table_ranked2)[j])
          pheno <- as.data.frame(colnames(merged_table))
          rownames(pheno) <- colnames(merged_table)
          pheno_new <- new("AnnotatedDataFrame",data=pheno)
          expr_set <- new("ExpressionSet",exprs = merged_table, phenoData=pheno_new)
          ### calculate distances
          distances <- ScoreGSEA(expr_set , threshold_count,"avg")
          final_dist[i,j] <-distances[1,2] #or (2,1) just not the primary diagonal 
        }
        print((i/ncol(num_table1))*100)
      }
      colnames(final_dist) <- names1
      rownames(final_dist) <- names2
      return(final_dist)
    }
    
    #2)get_cmap_signatures
    get_cmap_signatures <- function(cmap_path, sig_ids, landmark = TRUE, landmark_df = NULL) {
      library(tidyverse)
      library(cmapR)
      library(rhdf5)
      library(AnnotationDbi)
      library(org.Hs.eg.db)
      
      ds_path <- cmap_path
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
    
    ###go_term_anno to annotate the GO terms with gene identifiers
    library(tidyverse)
    library(topGO)
    library(org.Hs.eg.db)
    library(GO.db)
    ### read landmark genes
    
    landmark_df <- read_tsv(file = "cmap_landmark_genes.txt")
    
    ### go terms
    
    genes <- factor(x = rep(1,978),levels = c(0,1))
    names(genes) <- landmark_df$`Entrez ID`
    
    GOobject <- new("topGOdata",ontology = "BP", allGenes = genes, annot=annFUN.org, mapping="org.Hs.eg.db", 
                    ID = "entrez", nodeSize = 10)
    
    term.genes <- genesInTerm(GOobject, GOobject@graph@nodes)
    ###end of go_term_anno
    
    
    ### go_path_analysis function to exact GSEA for GO terms related to biological processes
    
    #it calculates the NES and p.adj of the given signature ids
    #GOterms in a list form are used
    library(fgsea)
    library(gage)
    library(EGSEAdata)
    library(AnnotationDbi)
    library(GeneExpressionSignature)
    
    ### you need to load the profiles from the GSE file
    #
    
    sig_ids <- read.csv(paste0(sig_ids_train))
    sigs_x <- unique(as.character(sig_ids$sig_id.x))
    sigs_y <- unique(as.character(sig_ids$sig_id.y))
    sigs <- c(sigs_x,sigs_y)
    sig_ids <- unique(sigs)
    
    #cmap_path_to_gctx <- "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
    landmark_df <- read_tsv(file = "cmap_landmark_genes.txt")
    goterms <- term.genes
    kd_sigs <- readRDS(paste0(tf_knockdowns))
    kd_sigs <- kd_sigs %>% filter(quality==1) %>% filter(cell_id==paste0(toupper(cell_line))) %>% dplyr::select(sig_id)
    kd_sigs <- as.vector(kd_sigs$sig_id)
    if (cmap_path != 0){
      profiles_train <- get_cmap_signatures(cmap_path_to_gctx = cmap_path,sig_ids = sig_ids,landmark = T,landmark_df = landmark_df)
    } else {
      profiles_train <- readRDS("train_genes.rds")}
    print("train profiles loaded")
    ### change to entrez ids
    rownames(profiles_train) <- landmark_df$`Entrez ID`
    rownames(profiles_train) <- as.character(rownames(profiles_train))
    
    print("running fgsea")
    
    go_list <- apply(profiles_train,MARGIN = 2,fgsea,pathways = goterms,
                     minSize=10,
                     maxSize=500,
                     nperm=1000)
    
    print("fgsea finished")
    
    ### get the NES and p.adj
    
    print("preparing output")
    NES <- go_list[[1]]$NES
    padj <- go_list[[1]]$padj
    
    for (i in 2:length(go_list)) {
      
      NES <- cbind(NES,go_list[[i]]$NES)
      padj <- cbind(padj,go_list[[i]]$padj)
    }
    
    colnames(NES) <- names(go_list)
    rownames(NES) <- go_list[[1]]$pathway
    colnames(padj) <- names(go_list)
    rownames(padj) <- go_list[[1]]$pathway
    
    comb <- list(NES,padj)
    
    train_go <-comb[[1]]
    #do the same for the kd_sigs that are in the tf_kd_cgs.rds file
    if (cmap_path != 0){
      profiles_train <- get_cmap_signatures(cmap_path_to_gctx = cmap_path,sig_ids = kd_sigs,landmark = T,landmark_df = landmark_df)
    } else {
      profiles_tf <- readRDS("kd_genes.rds")}
    print("tf profiles loaded")
    ### change to entrez ids
    rownames(profiles_tf) <- landmark_df$`Entrez ID`
    rownames(profiles_tf) <- as.character(rownames(profiles_tf))
    
    print("running fgsea")
    
    go_list <- apply(profiles_tf,MARGIN = 2,fgsea,pathways = goterms,
                     minSize=10,
                     maxSize=500,
                     nperm=1000)
    
    print("fgsea finished")
    
    ### get the NES and p.adj
    
    print("preparing output")
    NES <- go_list[[1]]$NES
    padj <- go_list[[1]]$padj
    
    for (i in 2:length(go_list)) {
      
      NES <- cbind(NES,go_list[[i]]$NES)
      padj <- cbind(padj,go_list[[i]]$padj)
    }
    
    colnames(NES) <- names(go_list)
    rownames(NES) <- go_list[[1]]$pathway
    colnames(padj) <- names(go_list)
    rownames(padj) <- go_list[[1]]$pathway
    
    comb <- list(NES,padj)
    
    kd_go <-comb[[1]] 
    ### end of go path analysis
    print("end of go path analysis")
    ###calculation of the go term distances of the tfs and the sig_ids
    #set the 5 thresholds
    library(doFuture)
    thresholds <- c(10,20,30,40,50)
    registerDoFuture()
    plan(multiprocess,workers = 4)
    
    go_distance <- NULL
    print("starting go_distance")
    go_distance <- foreach(thresh = thresholds) %dopar% {
      
      distance_scores_all(num_table1 = train_go,num_table2 =kd_go, threshold_count = thresh, names1 = as.character(colnames(profiles_train)),names2 = as.character(colnames(profiles_tf)))
    }
    
    
    
    final_dist <- (go_distance[1]+go_distance[2]+go_distance[3]+go_distance[4]+go_distance[5])/5
    colnames(final_dist) <- as.character(colnames(go_distance[1]))
    rownames(final_dist) <- as.character(rownames(go_distance[1]))
    final_dist <- as.data.frame(final_dist)
    
    return(final_dist)
  }
  
  library(tidyverse)
  kd_sigs <- readRDS(paste0(tf_knockdowns))
  kd_sigs <- kd_sigs %>% filter(quality==1) %>% filter(cell_id==paste0(toupper(cell_line))) %>% dplyr::select(sig_id)
  sig_ids <- read.csv(paste0(sig_ids_train))
  sigs_x <- unique(as.character(sig_ids$sig_id.x))
  sigs_y <- unique(as.character(sig_ids$sig_id.y))
  sigs <- c(sigs_x,sigs_y)
  sig_ids <- unique(sigs)
  if (cmap_path!=0){kd_genes <-get_cmap_signatures <- function(cmap_path, kd_sigs, landmark = TRUE, landmark_df = NULL)
    kd_genes <-as.data.frame(kd_genes)
  } else {  kd_genes <- readRDS('kd_genes.rds')
  kd_genes <-as.data.frame(kd_genes)}
  if (cmap_path!=0){train_genes <-get_cmap_signatures <- function(cmap_path, sig_ids, landmark = TRUE, landmark_df = NULL)
    train_genes <-as.data.frame(kd_genes)
  } else {train_genes <- readRDS('train_genes.rds')
  train_genes <- as.data.frame(train_genes)}

  indx <- grepl(paste0(toupper(cell_line)), colnames(kd_genes))
  kd_genes <- kd_genes[,indx]
  kd_train_genes <- cbind(kd_genes,train_genes)
  
  
  rank_tf <- tf_enrichment_ranks(tf_knockdowns, cell_line, sig_ids_train)
  dimension <- dim(rank_tf)[2]
  kd_sigs <- readRDS(paste0(tf_knockdowns))
  kd_sigs <- kd_sigs %>% dplyr::filter(quality==1) %>% dplyr::filter(cell_id==paste0(toupper(cell_line))) %>% dplyr::select(sig_id,pert_iname)
  rank_tf_kd <-rank_tf[which(rownames(rank_tf) %in% kd_sigs$pert_iname),]
  if (!(KD %in% rownames(rank_tf_kd))){
    print("Knockdown is not contained in chosen cell line. Please select a different cell line or knockdown.")
    print("Knockdowns in this cell line are:")
    print(rownames(rank_tf_kd))
  }else{
    rank_tf_kd<-rank_tf_kd[,which(colnames(rank_tf_kd) %in% kd_sigs$sig_id)]
    ranks_good <- as.numeric(diag(as.matrix(rank_tf_kd)))
    tfs_good <- as.character(rownames(rank_tf_kd))
    good <- cbind(tfs_good, ranks_good)
    good <- as.data.frame(good)
    good$ranks_good <- as.numeric(good$ranks_good)
    good <- good %>% filter(ranks_good < thresh1*dimension)
    if (!(KD %in% good$tfs_good)){
      print("Knockdown is not working for this cell line. Please select a different cell line or knockdown.")
      print("Knockdowns that are working fine in this cell line are:")
      print(good$tfs_good)
    }else{
      print("The Knockdown chosed works well in this cell line!")
      print("Continuing with distance calculations")
    }
  }
  #####Find Drugs with Low Distance to TFs
  if (cmap_path!=0){distances <- distance_go_kd_train(tf_knockdowns, cell_line, sig_ids_train,cmap_path)
  }else{
  distances <- read.csv('go_dist_kd_train_a375.csv')}
  rownames(distances)<-as.character(distances[,1])
  distances <- distances[,-1]
  help <- integer(ncol(distances))
  for (i in 1:nrow(good)){
    help <- help +grepl(good$tfs_good[i], colnames(distances))
  }
  distances <- distances[,which(help==1)]
  indx <- grepl(paste0(toupper(cell_line)), colnames(distances))
  distances <- distances[indx]
  indx <- grepl(KD, colnames(distances))
  distances <- distances[indx]
  drugs <- rownames(distances)
  drugs <- drugs[(which(distances[,1]<thresh2))]
  if (is_empty(drugs)==TRUE){
    print("No drugs were found. Lowest distance to this kd in this cell line is:")
    print(range(distances)[1])
  }else{
  print("The possible drugs are:")
  print(drugs)
  
  ###Check rank of those drugs
  drugs_ranks <- rank_tf
  drugs_ranks <- drugs_ranks[which(rownames(drugs_ranks)==KD),]
  drugs_ranks <- drugs_ranks[,which(colnames(drugs_ranks) %in% drugs)]
  print("")
  print("The drugs predicted are:")
  a <- c()
  for (i in 1:nrow(drugs_ranks)){
    for (j in 1:ncol(drugs_ranks)){
      if (drugs_ranks[i,j]<thresh3*nrow(distances)){
        print(colnames(drugs_ranks[j]))
        a <- c(a,colnames(drugs_ranks[j]))
      }
    }
  }
  return(a)
  }
}
### example###
#tf_knockdowns <- "tf_kd_cgs.rds"
#sig_ids_train <- "train.csv"
#cell_line <- "a375"
#KD <- c("EGR1")
#cmap_path<-0
#drugos <- drug_finder(KD, tf_knockdowns, cell_line, sig_ids_train, thresh1=0.2,thresh2=0.2,thresh3=0.2,0)
###End of example###
