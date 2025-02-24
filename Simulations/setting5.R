#!/usr/bin/env Rscript
#SBATCH --mem-per-cpu=2g  --time=4:00:00 --mail-type=ALL --mail-user=yisha.yao@yale.edu


load("/gpfs/gibbs/project/zhang_heping/yy634/cluster_gene_association/real_data/chr22/chr22_design.RData")
rm(n_X, n_X_group_sizes)  
alpha <- 500;  opt_lam <- 2800; lam_0 <- 3800;  lam_thred <- 1000;  p <- 200;  del <- 0.05*n/n_q  # tuning lam_0, lam_thred, alpha=(20,100,200,500)
 
library(plyr)
library(MASS)
library(Matrix)
library(plus)
library(mclust)
library(grplasso)
rand_vect <- function(N, M, sd = 1, pos.only = TRUE) {
  vec <- rnorm(N, M/N, sd)
  if (abs(sum(vec)) < 0.01) vec <- vec + 1
  vec <- round(vec / sum(vec) * M)
  deviation <- M - sum(vec)
  for (. in seq_len(abs(deviation))) {
    vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
  }
  if (pos.only) while (any(vec < 0)) {
    negs <- vec < 0
    pos  <- vec > 0
    vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
    vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
  }
  vec
}

Beta <- matrix(0, n_pc_d, p)
Beta_cluster_num <- rep(0, n_q)
Beta_labels <- matrix(0, n_q, p)
true_means <- matrix(0, n_pc_d, p)
set.seed(58)   
for (i in 1:n_q){
  r <- rbinom(1,1,del) 
  if (r ==0){} else{
    x <- sample(2:6, 1)
    Beta_cluster_num[i] <- x  #the number of clusters in the ith row-block of B
    gaussian_mix_size <- sample(c(ceiling(p/2), ceiling(2*p/3)), 1)
    cluster_sizes <- rand_vect(x, gaussian_mix_size) #the number of elements in each cluster
    label_temp <- as.matrix(rep(0:x, c((p-gaussian_mix_size),cluster_sizes)))  
    Beta_temp <- matrix(0, n_PC_group_sizes[i], (p-gaussian_mix_size))
    v <- runif(x, 0.5, 1.5)
    means <- matrix(0, n_PC_group_sizes[i], x)
    for (k in 1:x){
      means[,k] <- runif(n_PC_group_sizes[i], -3, 3)
    }
    for (j in 1:x){
      multigaussian <- mvrnorm(n = cluster_sizes[j], means[,j], v[j]*diag(n_PC_group_sizes[i]))
      Beta_temp <- cbind(Beta_temp, t(multigaussian))
    }
    index <- sample(nrow(label_temp))
    Beta_labels[i,] <- label_temp[index]
    Beta[n_pc_low[i]:n_pc_upp[i],] <- Beta_temp[,index]
    for (l in 1:p){
      if (Beta_labels[i,l] !=0){true_means[n_pc_low[i]:n_pc_upp[i],l] <- means[,Beta_labels[i,l]]}
    }
  }
}
rm(gaussian_mix_size, r, i, j, k, l, x, cluster_sizes, label_temp, Beta_temp, means, multigaussian, index)
E <- mvrnorm(n, mu=rep(0, p), Sigma=diag(p))
Y <- n_s_PC_design %*% Beta + E

Betas <- list()  #Betas[[t]] <- matrix(0, n_pc_d, p)
Blabels <- list()  #Blabels[[t]] <- matrix(rep(1:p, n_q), ncol=p, byrow = T)
means <- list()  #means[[t]] <- matrix(0, n_pc_d, p)
classifications <- list()  #classifications[[t]] <- matrix(0, n_q, p)
cluster_nums <- list()  #cluster_nums[[t]] <- rep(0, n_q)
best_models <- list()  #best_models[[t]] <- c()
covariances <- list()  #covariances[[t]] is also a list.
#to obtain the starting point
Beta_0 <- matrix(0, n_pc_d, p)
Beta_0_labels <- matrix(rep(1:p, n_q), ncol=p, byrow = T)
for (j in 1:p){
  fit <- grplasso(n_s_PC_design, y=Y[,j], index=n_PC_grouping, lambda=lam_0, model=LinReg(), penscale=sqrt, 
                  control=grpl.control(update.hess="lambda", trace=0), center=FALSE)
  grp_est <- fit$coefficients[,1]
  if (all(grp_est==0)){}else {
    indicator <- c(1:n_q)
    for (k in 1:n_q){
      if (all(grp_est[n_pc_low[k]:n_pc_upp[k]]==0)){indicator[k] <- 0}
    }  
    indicator <- indicator[which(indicator!=0)]
    grp_list <- list()
    for (l in indicator){
      grp_list[[l]] <- c(n_pc_low[l]:n_pc_upp[l])
    } 
    grp_ind <- unlist(grp_list)
    mcp <- plus(n_s_PC_design[,grp_ind], Y[,j], method = "mc", normalize = TRUE, intercept = FALSE)
    Beta_0[grp_ind,j] <- mcp$beta[dim(mcp$beta)[1],]
  }
  for (i in 1:n_q){
    if (all(Beta_0[n_pc_low[i]:n_pc_upp[i],j]==0)){Beta_0_labels[i,j] <- 0}
  }
}
Betas[[1]] <- Beta_0
Blabels[[1]] <- Beta_0_labels 
rm(i,j,k,l, fit, mcp, grp_list, grp_ind, grp_est, indicator, Beta_0, Beta_0_labels)

Beta_temp <- Betas[[1]]
Blabels_temp <- Blabels[[1]]
means_temp <- matrix(0, n_pc_d, p)
covs_temp <- list()
best_models_temp <- c()
classifications_temp <- matrix(0, n_q, p)
cluster_nums_temp <- rep(0, n_q)
for (i in 1:n_q){
  if (length(which(Blabels_temp[i,]!=0))<5){ 
    classifications_temp[i,] <- rep(0, p)
    best_models_temp <- c("nocluster")
    means_temp[n_pc_low[i]:n_pc_upp[i],] <- 0
    covs_temp[[i]] <- 0
    cluster_nums_temp[i] <- 0
  } else{
    temp <- Beta_temp[(n_pc_low[i]:n_pc_upp[i]),Blabels_temp[i,]]
    mod <- Mclust(t(temp), G=1:10, control=emControl(eps=0.001))
    index1 <- which(Blabels_temp[i,] ==0)
    part1 <- as.data.frame(cbind(index1, 0))
    colnames(part1) <- c("index", "label")
    index2 <- which(Blabels_temp[i,] !=0)
    part2 <- as.data.frame(cbind(index2, mod$classification))
    colnames(part2) <- c("index", "label")
    est <- rbind(part1, part2)
    est <- est[order(est$index),]
    classifications_temp[i,] <- est$label
    best_models_temp[i] <- mod$modelName
    for (j in 1:p){
      if (classifications_temp[i,j]==0){means_temp[n_pc_low[i]:n_pc_upp[i],j] <- 0}else {
        if (n_PC_group_sizes[i]==1){
          means_temp[n_pc_low[i]:n_pc_upp[i],j] <- mod$parameters$mean[classifications_temp[i,j]]
        }else {
          means_temp[n_pc_low[i]:n_pc_upp[i],j] <- mod$parameters$mean[,classifications_temp[i,j]]
        }
      }
    }
    covs_temp[[i]] <- mod$parameters$variance$sigma
    cluster_nums_temp[i] <- mod$G
  }
}
means[[1]] <- means_temp
covariances[[1]] <- covs_temp
classifications[[1]] <- classifications_temp
cluster_nums[[1]] <- cluster_nums_temp
best_models[[1]] <- best_models_temp
rm(i, temp, mod, index1, part1, index2, part2, est, Beta_temp, Blabels_temp) 
#  means_temp, classifications_temp, cluster_nums_temp, best_models_temp, covs_temp


t <- 1
repeat{
  #update Beta
  alpha <- alpha*0.8
  opt_lam <- max(opt_lam*0.9, lam_thred) #sometimes opt_lam decreasing would lower the TNR
  Beta_temp <- matrix(0, n_pc_d, p)
  Blabels_temp <- matrix(rep(1:p, n_q), ncol=p, byrow = T)
  for (j in 1:p){
    #update the on-support part
    n_g_ind <- which(classifications_temp[,j]!=0)
    z_g_ind <- which(classifications_temp[,j]==0)  
    nonzero_list <- list()
    cov_list <- list()
    for (i in n_g_ind){ 
      nonzero_list[[i]] <- c(n_pc_low[i]:n_pc_upp[i])
      covs <- covs_temp[[i]]
      class_ind <- classifications_temp[i,j]  
      if (n_PC_group_sizes[i] == 1){
        if (length(covs)==1){cov_list[[i]] <- 1/covs} else{cov_list[[i]] <- 1/covs[class_ind]}
      } else{cov_list[[i]] <- ginv(covs[,,class_ind])}
    }
    nonzero_ind <- unlist(nonzero_list)
    Xn <- n_s_PC_design[, nonzero_ind]
    muj <- means_temp[nonzero_ind,j]
    cov_list[sapply(cov_list, is.null)] <- NULL
    covj <- as.matrix(bdiag(cov_list))
    
    zero_list <- list()
    for (k in z_g_ind){
      zero_list[[k]] <- c(n_pc_low[k]:n_pc_upp[k]) 
    }
    zero_ind <- unlist(zero_list)
    Xz <- n_s_PC_design[, zero_ind]
    rm(i,k)
    leftmatrix <- (t(Xn) %*% Xn)/n + alpha*covj/n_pc_d
    rightmatrix <- t(Xn) %*% (Y[,j]- n_s_PC_design %*% Beta_temp[,j])/n + alpha*covj%*%muj/n_pc_d
    Beta_temp[nonzero_ind, j] <- ginv(leftmatrix) %*% rightmatrix
    
    #update the off-support part
    indicator <- c()
    res <- as.vector(Y[,j] - n_s_PC_design %*% Beta_temp[, j])
    fit1 <- grplasso(Xz, y=res, index=n_PC_grouping[zero_ind], lambda=opt_lam, model=LinReg(), penscale=sqrt, 
                     control=grpl.control(update.hess="lambda", trace=0), center=FALSE)
    grp_est1 <- fit1$coefficients[,1]
    if (all(grp_est1==0)){}else {
      indicator <- z_g_ind
      count <- 0
      for (k in z_g_ind){
        l <- count + 1
        u <- count + n_PC_group_sizes[k]
        count <- count + n_PC_group_sizes[k]
        if (all(grp_est1[l:u]==0)){indicator[which(indicator==k)] <- 0}
      }  
      indicator <- indicator[which(indicator!=0)] #newly entered groups, that is previously zero groups but now nonzero
    }
    
    #update the support
    new_n_g_ind <- sort(c(n_g_ind, indicator))
    if (length(new_n_g_ind) >0){
      new_nonzero_list <- list()
      for (k in new_n_g_ind){
        new_nonzero_list[[k]] <- c(n_pc_low[k]:n_pc_upp[k]) 
      }
      new_nonzero_ind <- unlist(new_nonzero_list)
      Xnn <- n_s_PC_design[, new_nonzero_ind]
      fit2 <- grplasso(Xnn, y=Y[,j], index=n_PC_grouping[new_nonzero_ind], lambda=opt_lam, model=LinReg(), penscale=sqrt, 
                       control=grpl.control(update.hess="lambda", trace=0), center=FALSE)
      grp_est2 <- fit2$coefficients[,1]
      indicatornew <- new_n_g_ind
      count <- 0
      for (k in new_n_g_ind){
        l <- count + 1
        u <- count + n_PC_group_sizes[k]
        count <- count + n_PC_group_sizes[k]
        if (all(grp_est2[l:u]==0)){indicatornew[which(indicatornew==k)] <- 0}
      }  
      indicatornew <- indicatornew[which(indicatornew !=0)] 
      Beta_temp[, j] <- rep(0, n_pc_d)
      grp_list <- list()
      for (i in indicatornew){
        grp_list[[i]] <- c(n_pc_low[i]:n_pc_upp[i])
      }
      grp_ind <- unlist(grp_list)
      mcp <- plus(n_s_PC_design[,grp_ind], Y[,j], method = "mc", normalize = TRUE, intercept = FALSE)
      Beta_temp[grp_ind, j] <- mcp$beta[dim(mcp$beta)[1],]
    }
    
    for (i in 1:n_q){
      if (all(Beta_temp[n_pc_low[i]:n_pc_upp[i],j]==0)){Blabels_temp[i,j] <- 0}
    }
  }
  t <- t+1
  Betas[[t]] <- Beta_temp
  Blabels[[t]] <- Blabels_temp
  rm(n_g_ind, nonzero_list, cov_list, covs, 
     class_ind, nonzero_ind, Xn, muj, covj, z_g_ind, zero_list, zero_ind, Xz, 
     leftmatrix, rightmatrix, res, fit1, fit2, mcp, grp_list, grp_ind, indicator, indicatornew, i, j,
     count, grp_est1, grp_est2)
  #clustering on Beta_temp
  means_temp <- matrix(0, n_pc_d, p)
  covs_temp <- list()
  best_models_temp <- c()
  classifications_temp <- matrix(0, n_q, p)
  cluster_nums_temp <- rep(0, n_q)
  for (i in 1:n_q){
    if (length(which(Blabels_temp[i,]!=0))<5){
      classifications_temp[i,] <- rep(0, p)
      best_models_temp <- c("nocluster")
      means_temp[n_pc_low[i]:n_pc_upp[i],] <- 0
      covs_temp[[i]] <- 0
      cluster_nums_temp[i] <- 0
    } else{
      temp <- Beta_temp[(n_pc_low[i]:n_pc_upp[i]),Blabels_temp[i,]]
      mod <- Mclust(t(temp), G=1:10, control=emControl(eps=0.001))
      #if do not specify control=emControl(eps=0.001), the resulting cov might be singular, and then resort to ginv(con), not accurate in updating Beta
      index1 <- which(Blabels_temp[i,] ==0)
      part1 <- as.data.frame(cbind(index1, 0))
      colnames(part1) <- c("index", "label")
      index2 <- which(Blabels_temp[i,] !=0)
      part2 <- as.data.frame(cbind(index2, mod$classification))
      colnames(part2) <- c("index", "label")
      est <- rbind(part1, part2)
      est <- est[order(est$index),]
      classifications_temp[i,] <- est$label
      best_models_temp[i] <- mod$modelName
      for (j in 1:p){
        if (classifications_temp[i,j]==0){means_temp[n_pc_low[i]:n_pc_upp[i],j] <- 0}else {
          if (n_PC_group_sizes[i]==1){
            means_temp[n_pc_low[i]:n_pc_upp[i],j] <- mod$parameters$mean[classifications_temp[i,j]]
          }else {
            means_temp[n_pc_low[i]:n_pc_upp[i],j] <- mod$parameters$mean[,classifications_temp[i,j]]
          }
        }
      }
      covs_temp[[i]] <- mod$parameters$variance$sigma
      cluster_nums_temp[i] <- mod$G
    }
  }
  means[[t]] <- means_temp
  covariances[[t]] <- covs_temp
  classifications[[t]] <- classifications_temp
  cluster_nums[[t]] <- cluster_nums_temp
  best_models[[t]] <- best_models_temp
  rm(i, temp, mod, index1, part1, index2, part2, est) 
  if(sum((Beta_temp-Betas[[t-1]])^2)/(n_pc_d*p) < 0.0002 & t>5){break}
}

# Beta_labels records the true classifications
save(n_pc_d, p, n_q, n_PC_group_sizes, Beta, Beta_labels, Beta_cluster_num, true_means,
     Betas, Blabels, classifications, means, covariances, cluster_nums, best_models,
     file = "/gpfs/gibbs/project/zhang_heping/yy634/cluster_gene_association/submit2/simulation2/setting5_alpha500_lamthred1000.RData")
