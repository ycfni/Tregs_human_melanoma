require(tictoc)
require(Rcpp)
require(KRLS)
require(FNN)
require(rARPACK)
require(GenBinomApps)
require(ggplot2)

`%notin%` <- Negate(`%in%`) 

Rcpp::sourceCpp("/banach1/rq25/DA_project_boris_method/Rcpp/multiply.cpp")

########Step1: Construct affinity matrix
get_affinity_matrix <- function(data, 
                                sigma = NULL,
                                sigma_quantile = 0.0025, 
                                K = 25, 
                                adaptive = T,
                                truncated = F,
                                zero_out_diag = T){
  
  ##### Construct the distance matrix
  dist_mat <- dist(data)
  
  
  ##### Convert the distance matrix to the affinity matrix
  ##### Two options: 1) using fixed bandwidth (sigma); 2) using locally adaptive bandwidth
  
  if (!adaptive){
    if (is.null(sigma)) sigma <- quantile(dist_mat, probs = sigma_quantile)^2
    affinity_mat <- KRLS::gausskernel(data, sigma)
  }
  
  if (adaptive){
    dists <- as.matrix(dist_mat)
    
    sigma_list <- c()
    for (i in 1:nrow(dists)){
      sigma_list[i] <- sort(dists[i,])[K]
    }
    
    affinity_matrix <- matrix(0, nrow=nrow(dists), ncol=ncol(dists))
    if (truncated){
      knn_res <- FNN::get.knn(data, k=K)
      for (i in 1:nrow(affinity_matrix)){
        affinity_matrix[i,knn_res$nn.index[i,]] <- exp(-knn_res$nn.dist[i,]^2/(sigma_list[i]^2))
      }
    }
    if (!truncated){
      for (i in 1:nrow(affinity_matrix)){
        affinity_matrix[i,] <- exp(-dists[i,]^2/(sigma_list[i]^2))
      }
    }
    
    ##### To symmetrize the affinity matrix 
    affinity_mat <- sqrt(affinity_matrix * t(affinity_matrix))
    #affinity_mat <- (affinity_matrix + t(affinity_matrix))/2
    
  }
  
  
  ##### To enhance the performance by blocking self-loop diffusion, as recommended by Boris
  if (zero_out_diag){
    for (i in 1:nrow(affinity_mat)){
      affinity_mat[i,i] = 0
    }
  }
  
  
  affinity_mat
}



########Step2: Partition the Seurat object based on cluster labels
get_partition <- function(data_S,
                          partition.by,
                          reduction = "pca", 
                          dims = 1:30,
                          sigma = NULL,
                          sigma_quantile = NULL,
                          K = 25, 
                          adaptive = T,
                          truncated = F,
                          zero_out_diag = T
){
  
  ##### Function to do doubly normalization
  doubly_normalize <- function(data_k, 
                               maxT = 10000, 
                               tol = 1e-6){
    #message("Running doubly normalization...")
    d <- rep(1, nrow(data_k))
    err = max(abs(d *(data_k %*% d) - 1))
    for (i in 1:maxT){
      d = sqrt(d / (data_k%*%d))
      err = max(abs(d *(data_k %*% d) - 1))
      if (is.nan(err)) stop(paste("Halted - Please choose a larger k to increase the connectivity of graph."))
      if (err < tol) break()
      if (i == maxT) warning("Reached to a maximum of iteration.")
    } 
    data_k <- t(as.vector(d) * data_k) * as.vector(d)
    #message(err)
    data_k <- (data_k + t(data_k))/2
    eigenMapMatMult(data_k, data_k)
  }
  
  
  #####
  ##### Partitioning
  #####
  
  partition <- as.character(data_S@meta.data[, partition.by])
  discarded <- names(table(partition))[which(table(partition) <= K)]
  category <- unique(partition)
  
  
  if (length(discarded) > 0) {
    category <- category[-which(category %in% discarded)]
    message(paste(length(discarded), "partition(s) is(are) discarded due to cell number <=", K))
  }
  
  
  message(paste("The original data will be split into", length(category), "partitions."))
  
  
  cell_barcodes_list <- list()
  for (i in category){
    cell_barcodes_list[[i]] <- colnames(data_S)[which(partition == i)]
  }
  
  
  data_obj_list <- list()
  for (i in category){
    sub_data_S <- data_S[, cell_barcodes_list[[i]]]
    data_obj_list[[i]]<- sub_data_S[[reduction]]@cell.embeddings[, dims]
  }
  
  
  message("Constructing affinity matrices...")
  if (!adaptive) message(paste("Fixed bandwidth with sigma =", sigma, "is used."))
  if (adaptive) message(paste("Using locally adaptive kernel with K =", K))
  affinity_mat_list <- list()
  for (i in category){
    affinity_mat_list[[i]] <- get_affinity_matrix(data_obj_list[[i]],
                                                  sigma = sigma,
                                                  sigma_quantile = sigma_quantile, 
                                                  K = K, 
                                                  adaptive = adaptive,
                                                  truncated = truncated,
                                                  zero_out_diag = zero_out_diag)
  }
  
  #####Iteratively remove isolated outliers
  message("Iteratively remove isolated outliers...")
  for (iter in 1:10){
    N_cells_removed <- 0
    for (i in category){
      outliers_to_remove <- which(apply(affinity_mat_list[[i]], 1, sum) < 1e-12 |
                                    apply(affinity_mat_list[[i]] > 0, 1, sum) <= 2)
      N_cells_removed <- N_cells_removed + length(outliers_to_remove)
      #message(paste(length(outliers_to_remove), "outliers are removed from the partition", i))
      if (length(outliers_to_remove) > 0){
        cell_barcodes_list[[i]] <- cell_barcodes_list[[i]][-outliers_to_remove]
        data_obj_list[[i]] <- data_obj_list[[i]][-outliers_to_remove,]
        affinity_mat_list[[i]] <- affinity_mat_list[[i]][-outliers_to_remove, -outliers_to_remove]
      }
    }
    if (N_cells_removed == 0) break
  }
  
  
  
  message("Perform doubly normalization...")
  tic("Doubly Normalization")
  DN_affinity_mat_list <- list()
  for (i in category){
    #message(paste("Processing partition -", i))
    DN_affinity_mat_list[[i]] <- doubly_normalize(affinity_mat_list[[i]])
    DN_affinity_mat_list[[i]] <- DN_affinity_mat_list[[i]] %*% DN_affinity_mat_list[[i]]
  }
  toc()
  
  
  ##### Eigendecomposition
  message("Perform eigendecomposition...")
  tic("Eigendecomposition")
  E_list <- list()
  for (i in category){
    #message(paste("Processing partition -", i))
    #options(warn=0)
    E_list[[i]] <- rARPACK::eigs_sym(DN_affinity_mat_list[[i]], k = nrow(DN_affinity_mat_list[[i]]))
    #options(warn=1)
    if (E_list[[i]]$values[2] == 1) stop("Partition-", i, "consists of disconnected components and needs to be further partitioned.")
    if (E_list[[i]]$values[2] >= 0.999) print(paste0("WARNING: Please increase K in partition-", i, " to improve computation efficiency."))
  }
  toc()
  
  
  ##### Store data_obj_list, DN_affinity_mat_list, E_list, cell_barcodes_list
  partition_object_list <- list()
  partition_object_list[["data_obj_list"]] <- data_obj_list
  partition_object_list[["DN_affinity_mat_list"]] <- DN_affinity_mat_list
  partition_object_list[["E_list"]] <- E_list
  partition_object_list[["cell_barcodes_list"]] <- cell_barcodes_list
  
  
  partition_object_list
}



get_scan_statistics <- function(data_S,
                                partition_object_list, 
                                q = 0.01, 
                                condition,
                                COND_1){
  
  find_eps <- function(E_list, q){
    message("Find the optimal epsilon...")
    
    L = length(E_list)
    lambda2_vec = rep(0, L)
    n_vec = rep(0, L)
    for (l in 1:L){
      lambda2_vec[l] = E_list[[l]]$values[2]
      n_vec[l] = length(E_list[[l]]$values)
    }
    gridRes = 1e-5
    varEps = seq(from = gridRes, to = 1, by = gridRes)
    
    n_vec_mat <- matrix(n_vec, ncol=L, nrow=length(varEps), byrow=TRUE)
    lambda2_vec_mat <- matrix(lambda2_vec, ncol=L, nrow=length(varEps), byrow=TRUE)
    
    T_mat = t(ceiling(log(n_vec_mat/varEps)/log((1/lambda2_vec_mat))))
    H_mat = t(ceiling(log(1+n_vec_mat)/log(1+(1/n_vec_mat)*varEps^2)))
    
    riskBound = varEps + sqrt(0.5*log(t(n_vec/q) %*% pmin(T_mat,H_mat)))
    idx = which.min(riskBound)
    varEps_best = varEps[idx]
    min(riskBound)
    
    varEps_best
  }
  
  
  get_time_series <- function(P_vec, scale){
    findclosest <- function(P_vec, val){
      ind <- which.min(abs(P_vec-val))
      ifelse(P_vec[ind]<=val, ind-1, ind)
    }
    
    t_vec <- c(length(P_vec))
    t_tmp <- length(P_vec)
    t_flag <- length(P_vec)-1
    for (i in 1:length(P_vec)){
      thr <- P_vec[t_tmp]*scale
      t_tmp <- findclosest(P_vec[1:t_flag], thr)
      t_tmp <- min(t_tmp, t_flag)
      if(t_tmp <= 1){
        t_vec <- c(t_vec,1)
        break
      }
      if(t_tmp > 1) {
        t_vec <- c(t_vec,t_tmp)
        t_flag <- t_tmp-1
      }
    }
    t_vec
  }
  
  
  partition_object_list[["data_obj_list"]] -> data_obj_list
  partition_object_list[["E_list"]] -> E_list
  partition_object_list[["cell_barcodes_list"]] -> cell_barcodes_list
  
  
  ##### Find the optimal epsilon
  epsilon <- find_eps(E_list, q = q)
  message(sprintf("The epsilon is %s", as.character(epsilon)))
  
  ##### Estimate confidence interval for binomial proportion
  p <- GenBinomApps::clopper.pearson.ci(k = table(data_S@meta.data[, condition])[COND_1], 
                                        n = ncol(data_S), 
                                        alpha = q, 
                                        CI = "two.sided")$Upper.limit
  message(sprintf("The upper bound of p is %s", as.character(p)))
  
  
  
  message("Running local two-sample tests...")
  
  statistic_obj_list <- list()
  M <- c()
  
  for (ind in names(data_obj_list)){
    label <- as.character(data_S[,cell_barcodes_list[[ind]]]@meta.data[, condition])
    dss_y <- t(E_list[[ind]]$vectors) %*% (as.integer(label == COND_1) - p)
    max_T <- ceiling(log(length(E_list[[ind]]$values)/epsilon)/log(1/E_list[[ind]]$values[2]))
    
    
    ##### Construct intermediate matrices and values
    lambda_mat <- matrix(E_list[[ind]]$values, ncol = 1)
    lambda_mat <- cbind(lambda_mat, lambda_mat^2)
    for (i in 2:(as.integer(log2(max_T))+1)){
      lambda_mat <- cbind(lambda_mat, E_list[[ind]]$values^(2^(i-1)) * lambda_mat)
    }
    lambda_mat <- lambda_mat[,1:max_T]
    
    P_t <- sqrt(eigenMapMatMult(E_list[[ind]]$vectors^2, lambda_mat^2))
    P_t2 <- P_t^2
    
    
    ##### Find time steps
    scale <- 1+epsilon^2/length(label)
    time_series <- list()
    for (i in 1:nrow(P_t2)){
      time_series[[i]] <- get_time_series(P_t2[i,], scale)
    }
    
    
    ##### Calculate statistics
    P_ty <- eigenMapMatMult(E_list[[ind]]$vectors, lambda_mat * dss_y[,1])
    statistic <- P_ty/P_t
    
    
    statistic_obj <- list()
    statistic_obj$statistic <- statistic
    statistic_obj$P_t <- P_t
    statistic_obj$time_series <- time_series
    
    
    statistic_obj_list[[ind]] <- statistic_obj
    M <- c(M, length(unlist(time_series)))
  }
  
  
  ##### Calculate the threshold/cutoff
  thr = sqrt(0.5*log(sum(M)/(q/2)))
  
  
  ##### Summarize
  scan_statistic <- c()
  barcodeIdx <- c()
  timeIdx <- c()
  partitionIdx <- c()
  scan_statistic_max_per_point <- c()
  
  
  for (ind in names(data_obj_list)){
    statistic <- statistic_obj_list[[ind]]$statistic
    P_t2 <- statistic_obj_list[[ind]]$P_t^2
    time_series <- statistic_obj_list[[ind]]$time_series
    statistic_time_series <- list()
    
    for (i in 1:nrow(P_t2)){
      statistic_time_series[[i]] <- statistic[i,rev(time_series[[i]])]
      idx <- which(statistic_time_series[[i]] > thr)
      scan_statistic <- c(scan_statistic, statistic_time_series[[i]][idx])
      barcodeIdx <- c(barcodeIdx, rep(i, length(idx)))
      timeIdx <- c(timeIdx,rev(time_series[[i]])[idx])
      partitionIdx <- c(partitionIdx, rep(ind, length(idx)))
      scan_statistic_max_per_point <- c(scan_statistic_max_per_point, max(statistic[i,rev(time_series[[i]])]))
    }
  }
  
  RWSS_obj <- list()
  RWSS_obj[["scan_statistic"]] <- scan_statistic
  RWSS_obj[["barcodeIdx"]] <- barcodeIdx
  RWSS_obj[["timeIdx"]] <- timeIdx
  RWSS_obj[["partitionIdx"]] <- partitionIdx
  RWSS_obj[["scan_statistic_max_per_point"]] <- scan_statistic_max_per_point
  RWSS_obj[["all_intermediate_metrics"]] <- statistic_obj_list
  
  if (length(scan_statistic) == 0){
    message("No significant DA cell is found.")
    RWSS_obj[["p_values"]] <- c()
    RWSS_obj[["Ws_lower_bound"]] <- NA
    RWSS_obj[["originalBarcodeIdx"]] <- c()
    RWSS_obj[["partitionBarcodeIdx"]] <- c()
    RWSS_obj[["W_norm"]] <- c()
  }
  
  else {
    ##### Compute local p-values and lower bounds for <w,s>
    scanRes.p_values = sum(M)*exp(-2*(scan_statistic)^2)
    scanRes.W_norm = rep(0,length(scan_statistic))
    scanRes.originalBarcodeIdx = rep(0,length(scan_statistic))
    scanRes.partitionBarcodeIdx = rep(0,length(scan_statistic))
    for (i in 1:length(scan_statistic)){
      if (length(statistic_obj_list[[partitionIdx[i]]]$P_t[barcodeIdx[i] ,timeIdx[i]]) == 0) message(c(dim(statistic_obj_list[[partitionIdx[i]]]$P_t), barcodeIdx[i] ,timeIdx[i]))
      scanRes.W_norm[i] = statistic_obj_list[[partitionIdx[i]]]$P_t[barcodeIdx[i] ,timeIdx[i]]
      scanRes.originalBarcodeIdx[i] = cell_barcodes_list[[partitionIdx[i]]][barcodeIdx[i]]
      scanRes.partitionBarcodeIdx[i] = barcodeIdx[i]
    }
    scanRes.Ws_lower_bound = (scan_statistic - thr)*scanRes.W_norm
    
    
    RWSS_obj[["p_values"]] <- scanRes.p_values
    RWSS_obj[["Ws_lower_bound"]] <- scanRes.Ws_lower_bound
    RWSS_obj[["originalBarcodeIdx"]] <- scanRes.originalBarcodeIdx
    RWSS_obj[["partitionBarcodeIdx"]] <- scanRes.partitionBarcodeIdx
    RWSS_obj[["W_norm"]] <- scanRes.W_norm
  }
  
  
  RWSS_obj
}





perform_RWSS_test <- function(data_S,
                              partition_object_list, 
                              q = 0.01, 
                              condition){
  
  RWSS_res_list <- list()
  conditions <- as.character(unique(data_S@meta.data[, condition]))
  
  
  message(paste("Now look at", conditions[1], ">", conditions[2]))
  RWSS_res <- get_scan_statistics(data_S,
                                  partition_object_list,
                                  q, 
                                  condition,
                                  COND_1 = conditions[1])
  RWSS_res_list[[conditions[1]]] <- RWSS_res
  
  message(paste("Now look at", conditions[2], ">", conditions[1]))
  RWSS_res <- get_scan_statistics(data_S,
                                  partition_object_list,
                                  q, 
                                  condition,
                                  COND_1 = conditions[2])
  RWSS_res_list[[conditions[2]]] <- RWSS_res
  
  
  RWSS_res_list
}




get_DA_score <- function(data_S, RWSS_res_list){
  conditions <- names(RWSS_res_list)
  
  summary_mat <- as.data.frame(matrix(0, nrow = ncol(data_S), ncol = 3))
  rownames(summary_mat) <- colnames(data_S)
  colnames(summary_mat) <- c(paste0("DA_score_", conditions), "DA_score")
  
  for (j in 1:2){
    cd <- conditions[j]
    RWSS_res_list[[cd]]$Ws_lower_bound
    RWSS_res_list[[cd]]$originalBarcodeIdx
    DA_cells <- unique(RWSS_res_list[[cd]]$originalBarcodeIdx)
    for (i in DA_cells){
      summary_mat[i, j] <- max(RWSS_res_list[[cd]]$Ws_lower_bound[which(RWSS_res_list[[cd]]$originalBarcodeIdx == i)])
    }
  }
  
  summary_mat[, 3] <- summary_mat[, 1] - summary_mat[, 2]
  summary_mat[, 3][which(summary_mat[, 1] != 0 & summary_mat[, 2] != 0 )] <- NA
  
  data_S <- AddMetaData(data_S, summary_mat, col.name = colnames(summary_mat))
  
  data_S
}




get_DA_cell <- function(data_S,
                        RWSS_res_list,
                        cutoff = 0){
  
  conditions <- names(RWSS_res_list)
  
  summary_mat <- as.data.frame(matrix("Other", nrow = ncol(data_S), ncol = 3))
  rownames(summary_mat) <- colnames(data_S)
  colnames(summary_mat) <- c(paste0("DA_", conditions), "DA")
  
  res <- RWSS_res_list[[conditions[1]]]
  summary_mat[unique(res$originalBarcodeIdx[which(res$Ws_lower_bound > cutoff)]), 1] <- paste0(conditions[1], "-enriched")
  res <- RWSS_res_list[[conditions[2]]]
  summary_mat[unique(res$originalBarcodeIdx[which(res$Ws_lower_bound > cutoff)]), 2] <- paste0(conditions[2], "-enriched")
  
  summary_mat[, "DA"] <- summary_mat[, 1]
  summary_mat[which(summary_mat[, 2] != "Other"), "DA"] <- paste0(conditions[2], "-enriched")
  summary_mat[which(summary_mat[, 1] != "Other" & summary_mat[, 2] != "Other"), "DA"] <- "Ambiguous"
  
  data_S <- AddMetaData(data_S, summary_mat, col.name = colnames(summary_mat))
  data_S$DA <- factor(data_S$DA, levels = c(paste0(conditions[1], "-enriched"),
                                            paste0(conditions[2], "-enriched"),
                                            "Other",
                                            "Ambiguous"))
  
  
  data_S
}


get_DA_cluster <- function(data_S, 
                           resolution = 0.1,
                           partition.by = "cluster",
                           assay = "RNA"){
  
  data_S_all <- data_S
  data_S <- subset(data_S_all, DA %notin% c("Other",
                                            "Ambiguous"))
  
  data_S$tmp <- paste0(data_S$DA, "|", data_S@meta.data[, partition.by])
  table(data_S$tmp)
  data_S$da_cluster <- NA
  for (i in unique(data_S$tmp)){
    sub_data_S <- subset(data_S, tmp == i)
    sub_data_S <- FindNeighbors(sub_data_S, dims = 1:30, verbose = F)
    sub_data_S <- FindClusters(sub_data_S, resolution = resolution, verbose = F)
    data_S$da_cluster[colnames(sub_data_S)] <- paste0(i, "-", sub_data_S@meta.data[, paste0(assay, "_snn_res.", as.character(resolution))])
  }
  data_S$da_cluster_indexed <- as.integer(as.factor(data_S$da_cluster))
  
  data_S_all$DA_cluster <- NA
  data_S_all$DA_cluster[colnames(data_S)] <- data_S$da_cluster_indexed
  
  data_S_all
  
}
