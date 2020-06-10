val_test_drugs_split <- function(path, n_val, n_test) {
  library(tidyverse)
  library(MLmetrics) 
  
  #load data
  dataframe <- read.csv(path, header=TRUE) #or read.csv(file=val_1, header=TRUE)
  
  breaks <- c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9)
  
  #find which drugs are in the cold set
  drugs_cold <- unique(c((dataframe$rdkit.x[which(dataframe$iscoldx == TRUE)]), (dataframe$rdkit.y[which(dataframe$iscoldy == TRUE)])))
  
  #drugs_cold <- sample(drugs_cold,80)
  
  ind_cold <- c(1:length(drugs_cold))                              
  densities <- matrix(0, length(drugs_cold), (length(breaks)-1))   
  
  all_mses <- matrix(0, length(drugs_cold), length(drugs_cold))     
  
  diag(all_mses) <- 10                                             
  
  for (i in 1:(length(drugs_cold))) {                                  
                                                                                                      
    ind <- c(which(dataframe$rdkit.x == drugs_cold[i]), which(dataframe$rdkit.y == drugs_cold[i]))    
    a = hist(dataframe[ind,]$value,breaks = breaks,freq = T,plot = F)                       
    
    densities[i,] <- a$density                                                                             
  }
  
  i <- 1
  l <-1
  small <- min(n_val, n_test)                                      
  big <- max(n_val, n_test)                                        
  ratio <- floor(big/small)                                        
  left <- big%%small               
  group_s <- (1:small)*0                                           
  group_b <- (1:big)*0                                          
  
  for (k in 1:small) {
    
    for (j in 1:length(drugs_cold)) {
      
      if (all_mses[i,j] != 10) {                                  
        
        all_mses[i,j] <- MSE(densities[i,],densities[j,])          
                                                                   
      }
    }
    
    
    for (z in 1:ratio) {                                          
      to_make_ten <- which.min(all_mses[i,])
      all_mses[,to_make_ten] <- 10
      all_mses[to_make_ten,] <- 10
      group_b[l] <- to_make_ten
      l <- l+1
    }
    group_s[k] <- i
    all_mses[,i] <- 10
    all_mses[i,] <- 10
    
    
    repeat {
      i <- i+1
      if ((! i %in% group_s) & (! i %in% group_b)) {
        break
      }
    }
  }


  if (left != 0) {
    ind_cold <- ind_cold [! ind_cold %in% group_s]
    ind_cold <- ind_cold [! ind_cold %in% group_b]
    group_b <- union(group_b, ind_cold )
  }
  group_b <- group_b [! group_b %in% 0]

  
  ind_s <- c(which(dataframe$rdkit.x == drugs_cold[group_s[1]]), which(dataframe$rdkit.y == drugs_cold[group_s[1]]))
  
  for (i in 2:(length(group_s))) {
    ind_s <- union(ind_s, c(which(dataframe$rdkit.x == drugs_cold[group_s[i]]), which(dataframe$rdkit.y == drugs_cold[group_s[i]])))
  }
  s_hist = hist(dataframe[ind_s,]$value,breaks = breaks,freq = T,plot = T)
  
  #group_b
  ind_b <- c(which(dataframe$rdkit.x == drugs_cold[group_b[1]]), which(dataframe$rdkit.y == drugs_cold[group_b[1]]))
  
  for (i in 2:(length(group_b))) {
    ind_b <- union(ind_b, c(which(dataframe$rdkit.x == drugs_cold[group_b[i]]), which(dataframe$rdkit.y == drugs_cold[group_b[i]])))
  }
  
  b_hist = hist(dataframe[ind_b,]$value,breaks = breaks,freq = T,plot = T)
  
  mse_val_test <- MSE(s_hist$density,b_hist$density)
  
  if (n_val == small) {
    my_list <- list(val_set = group_s, test_set = group_b, cold_list = drugs_cold, mse = mse_val_test)

  } else {
    my_list <- list(val_set = group_b, test_set = group_s, cold_list = drugs_cold, mse = mse_val_test)
  }
  
  return(my_list)
  

}  

both_sets <- val_test_drugs_split(path="C:/Users/Christina/Desktop/BioSysLab/val_1.csv", n_val=40, n_test=40)
both_sets$mse
